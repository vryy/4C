/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_artery_coupling_linebased.cpp

 \brief base algorithm for line-based (non-conforming) coupling between
        poromultiphase_scatra-framework and flow in artery networks
        including scalar transport

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#include "poromultiphase_scatra_artery_coupling_linebased.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_porofluidmultiphase_ele/porofluidmultiphase_ele_parameter.H"
#include "../linalg/linalg_multiply.H"
#include "poromultiphase_scatra_artery_coupling_pair.H"
#include "poromultiphase_scatra_artery_coupling_defines.H"


/*----------------------------------------------------------------------*
 | constructor                                         kremheller 05/18 |
 *----------------------------------------------------------------------*/
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::PoroMultiPhaseScaTraArtCouplLineBased(
    Teuchos::RCP<DRT::Discretization>  arterydis,
    Teuchos::RCP<DRT::Discretization>  contdis,
    const Teuchos::ParameterList&      meshtyingparams,
    const std::string&                 condname,
    const std::string&                 artcoupleddofname,
    const std::string&                 contcoupleddofname
    ):
    PoroMultiPhaseScaTraArtCouplBase(arterydis, contdis, meshtyingparams, condname, artcoupleddofname, contcoupleddofname),
    meshtyingparams_(meshtyingparams),
    porofluidmanagersset_(false),
    coupling_method_(DRT::INPUT::IntegralValue<INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod>(meshtyingparams,"ARTERY_COUPLING_METHOD")),
    pp_(meshtyingparams.get<double>("PENALTY")),
    timersearch_(Comm())
{

  // user info
  if(myrank_ == 0)
  {
    std::cout << "<                                                  >" << std::endl;
    PrintOutCouplingMethod();
    std::cout << "<                                                  >" << std::endl;
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "\n";
  }

  return;
}

/*----------------------------------------------------------------------*
 | init the strategy                                 kremheller 07/18   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::Init()
{
  // we do not have a moving mesh
  if(DRT::Problem::Instance()->ProblemType() == prb_porofluidmultiphase)
    evaluate_in_ref_config_ = true;

  // fill the vectors
  FillFunctionAndScaleVectors();

  // find possible pairs
  BruteForceSearch();

  // create the pairs
  CreateMeshTyingPairs();

  // pre-evaluate the pairs
  PreEvaluateMeshTyingPairs();

  // fill length of artery elements that is not influenced if the underlying
  // 2D/3D mesh moves (basically protruding artery elements or segments)
  FillUnaffectedArteryLength();

  // create the GID to segment vector
  CreateGIDToSegmentVector();

  // print out summary of pairs
  OutputSummary();

  // initialize phinp for continuous dis
  phinp_cont_ = Teuchos::rcp(new Epetra_Vector(*contdis_->DofRowMap(), true));
  // initialize phinp for artery dis
  phinp_art_ = Teuchos::rcp(new Epetra_Vector(*arterydis_->DofRowMap(), true));

  // initialize phinp for continuous dis (colmap)
  phinp_cont_colmap_ = Teuchos::rcp(new Epetra_Vector(*contdis_->DofColMap(), true));
  // initialize phinp for artery dis (colmap)
  phinp_art_colmap_ = Teuchos::rcp(new Epetra_Vector(*arterydis_->DofColMap(), true));

  // initialize phinp for continuous dis
  zeros_cont_ = Teuchos::rcp(new Epetra_Vector(*contdis_->DofRowMap(), true));
  // initialize phinp for artery dis
  zeros_art_ = Teuchos::rcp(new Epetra_Vector(*arterydis_->DofRowMap(), true));

  // -------------------------------------------------------------------
  // create empty D and M matrices (27 adjacent nodes as 'good' guess)
  // -------------------------------------------------------------------
  D_ = Teuchos::rcp(new LINALG::SparseMatrix(*(arterydis_->DofRowMap()),27,false,true));
  M_ = Teuchos::rcp(new LINALG::SparseMatrix(*(arterydis_->DofRowMap()),27,false,true));

  // full map of continous and artery dofs
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(Teuchos::rcp(new Epetra_Map(*contdis_->DofRowMap())));
  maps.push_back(Teuchos::rcp(new Epetra_Map(*arterydis_->DofRowMap())));

  fullmap_ = LINALG::MultiMapExtractor::MergeMaps(maps);
  /// dof row map of coupled problem splitted in (field) blocks
  globalex_ = Teuchos::rcp(new LINALG::MultiMapExtractor());
  globalex_->Setup(*fullmap_,maps);

  // check global map extractor
  globalex_->CheckForValidMapExtractor();

  return;
}

/*----------------------------------------------------------------------*
 | setup the linear system of equations                kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::SetupSystem(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat,
    Teuchos::RCP<Epetra_Vector>                 rhs,
    Teuchos::RCP<LINALG::SparseMatrix>          sysmat_cont,
    Teuchos::RCP<LINALG::SparseMatrix>          sysmat_art,
    Teuchos::RCP<const Epetra_Vector>           rhs_cont,
    Teuchos::RCP<const Epetra_Vector>           rhs_art,
    Teuchos::RCP<const LINALG::MapExtractor>    dbcmap_cont,
    Teuchos::RCP<const LINALG::MapExtractor>    dbcmap_art
    )
{

  if(!porofluidmanagersset_)
  {
    for(unsigned i = 0; i < mesht_elepairs_.size(); i++)
      mesht_elepairs_[i]->SetupFluidManagersAndMaterials(contdis_->Name());
    porofluidmanagersset_ = true;
  }

  // reset
  rhs->PutScalar(0.0);
  sysmat->Matrix(0,1).Zero();
  sysmat->Matrix(1,0).Zero();

  // evaluate and assemble the pairs
  // note: OD terms, i.e., (0,1) and (1,0) directly assembled into sysmat
  //       main-diag terms, i.e., (0,0) and (1,1) added into _cont and _art
  //       rhs is added into global rhs
  EvaluateMeshTyingPairs(sysmat, rhs, sysmat_cont, sysmat_art);

  // add normal part to rhs
  rhs->Update(1.0, *globalex_->InsertVector(rhs_cont,0), 1.0);
  rhs->Update(1.0, *globalex_->InsertVector(rhs_art,1), 1.0);

  // apply DBCs
  // 1) on vector
  LINALG::ApplyDirichlettoSystem(rhs, zeros_cont_, *(dbcmap_cont->CondMap()));
  LINALG::ApplyDirichlettoSystem(rhs, zeros_art_,  *(dbcmap_art->CondMap()));
  // 2) on main-diag-matrices
  sysmat_cont->ApplyDirichlet(*(dbcmap_cont->CondMap()),true);
  sysmat_art->ApplyDirichlet(*(dbcmap_art->CondMap()),true);
  // 3) on OD-matrices
  sysmat->Matrix(0,1).Complete(sysmat_art->RangeMap(), sysmat_cont->RangeMap());
  sysmat->Matrix(1,0).Complete(sysmat_cont->RangeMap(), sysmat_art->RangeMap());
  sysmat->Matrix(0,1).ApplyDirichlet(*(dbcmap_cont->CondMap()),false);
  sysmat->Matrix(1,0).ApplyDirichlet(*(dbcmap_art->CondMap()),false);

  // get also the main-diag terms into the global sysmat
  sysmat->Assign(0,0,LINALG::View,*sysmat_cont);
  sysmat->Assign(1,1,LINALG::View,*sysmat_art);

  bool matlab = false;
  if (matlab)
  {
    //sparse_matrix
    std::string filename = "../o/mymatrix.dat";
    std::string filename_vc = "../o/myvec.dat";
    LINALG::PrintBlockMatrixInMatlabFormat(filename, *(sysmat));
    LINALG::PrintVectorInMatlabFormat(filename_vc, *rhs, true);
    dserror("exit");
  }

  return;
}

/*----------------------------------------------------------------------*
 | brute force search for near elements                kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::BruteForceSearch()
{
  if(myrank_ == 0)
    std::cout << "Starting with brute force search for coupling ... ";
  // reset timer
  timersearch_.ResetStartTime();
  // *********** time measurement ***********
  double dtcpu = timersearch_.WallTime();
  // *********** time measurement ***********



  // compute the search radius
  const double searchradius = ComputeSearchRadius();

  // scheme is elecolmap -- nodecolmap
  for (int i=0;i<arterydis_->ElementColMap()->NumMyElements();++i)
  {
    const int artelegid = arterydis_->ElementColMap()->GID(i);
    DRT::Element* artele = arterydis_->gElement(artelegid);

    DRT::Node** artnodes = artele->Nodes();

    // loop over all nodes of this element
    for(int jnode = 0; jnode < artele->NumNode(); jnode++)
    {

      DRT::Node* artnode = artnodes[jnode];

      static LINALG::Matrix<3,1> artpos;
      artpos(0) = artnode->X()[0];
      artpos(1) = artnode->X()[1];
      artpos(2) = artnode->X()[2];

      // loop over all column nodes of 2D/3D discretization
      for (int j=0;j<contdis_->NodeColMap()->NumMyElements();++j)
      {
        const int contgid = contdis_->NodeColMap()->GID(j);
        DRT::Node* contnode = contdis_->gNode(contgid);

        static LINALG::Matrix<3,1> contpos;
        contpos(0) = contnode->X()[0];
        contpos(1) = contnode->X()[1];
        contpos(2) = contnode->X()[2];

        static LINALG::Matrix<3,1> dist;
        dist.Update(1.0, contpos, -1.0, artpos, 0.0);

        if(dist.Norm2() < searchradius)
        {
          DRT::Element** contneighboureles = contnode->Elements();

          for(int l = 0; l < contnode->NumElement(); l++)
          {
            DRT::Element* thiscontele = contneighboureles[l];
            const int conteleid = thiscontele->Id();

            nearbyelepairs_[artelegid].insert(conteleid);
          }
        }
      } // col-nodes 2D/3D
    } // col-nodes of artery
  } // col-eles of artery

//  Debug output
//  std::map<int, std::set<int>>::iterator it;
//
//  for ( it = nearbyelepairs_.begin(); it != nearbyelepairs_.end(); it++ )
//  {
//    std::cout << "on proc " << myrank_ << " artery ele " << it->first  // string (key)
//               << ": ";
//    for(auto f : it->second) {
//      std::cout << f << ", ";
//    }
//    std::cout << "\n";
//  }



  // *********** time measurement ***********
  double mydtsearch = timersearch_.WallTime() - dtcpu;
  double maxdtsearch = 0.0;
  Comm().MaxAll(&mydtsearch,&maxdtsearch,1);
  // *********** time measurement ***********
  if(myrank_ == 0)
    std::cout << "Completed in " << maxdtsearch << "s" << std::endl;

  return;
}

/*----------------------------------------------------------------------*
 | create the pairs                                    kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::CreateMeshTyingPairs()
{
  // loop over pairs found by search
  std::map<int, std::set<int> >::const_iterator nearbyeleiter;
  for ( nearbyeleiter = nearbyelepairs_.begin(); nearbyeleiter != nearbyelepairs_.end(); ++nearbyeleiter )
  {
    const int artelegid = nearbyeleiter->first;
    std::vector< DRT::Element const *> ele_ptrs(2);
    ele_ptrs[0] = arterydis_->gElement(artelegid);

    std::set<int>::const_iterator secondeleiter;
    for ( secondeleiter = nearbyeleiter->second.begin(); secondeleiter != nearbyeleiter->second.end(); ++secondeleiter )
    {
      const int contelegid = *secondeleiter;
      ele_ptrs[1] = contdis_->gElement(contelegid);

      // construct, init and setup coupling pairs
      Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase> newpair =
          POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::CreateNewArteryCouplingPair(ele_ptrs);
      newpair->Init(ele_ptrs, meshtyingparams_, coupleddofs_cont_, coupleddofs_art_, scale_vec_, funct_vec_);

      // add to list of current contact pairs
      mesht_elepairs_.push_back(newpair);
    }
  }

  // output
  int total_numactive_pairs = 0;
  int numactive_pairs = static_cast<int>(mesht_elepairs_.size());
  Comm().SumAll(&numactive_pairs, &total_numactive_pairs, 1);

  if(myrank_ == 0)
    std::cout << "\nFound "<< total_numactive_pairs << " Artery-to-PoroMultiphaseScatra coupling pairs (segments)" << std::endl;

}

/*----------------------------------------------------------------------*
 | pre-evaluate the pairs and sort out duplicates      kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::PreEvaluateMeshTyingPairs()
{

  // pre-evaluate
  for(unsigned i = 0; i < mesht_elepairs_.size(); i++)
    mesht_elepairs_[i]->PreEvaluate();

  // delete the inactive and duplicate pairs
  std::vector<Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase> > active_msht_elepairs;
  for(unsigned i = 0; i < mesht_elepairs_.size(); i++)
  {
    if(mesht_elepairs_[i]->IsActive() && !IsDuplicateSegment(active_msht_elepairs, mesht_elepairs_[i]))
      active_msht_elepairs.push_back(mesht_elepairs_[i]);
  }

  mesht_elepairs_ = active_msht_elepairs;

  // output
  int total_numactive_pairs = 0;
  int numactive_pairs = static_cast<int>(mesht_elepairs_.size());
  Comm().SumAll(&numactive_pairs, &total_numactive_pairs, 1);
  if(myrank_ == 0)
    std::cout << "Only "<< total_numactive_pairs << " Artery-to-PoroMultiphaseScatra coupling pairs (segments) are active" << std::endl;

  // evaluate kappa
  if(coupling_method_ == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp and num_coupled_dofs_ > 0)
    EvaluateKappa();

  return;
}

/*----------------------------------------------------------------------*
 | evaluate kappa                                      kremheller 07/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::EvaluateKappa()
{

  LINALG::SerialDenseVector eleKappa;
  kappaInv_ = Teuchos::rcp(new Epetra_Vector(*arterydis_->DofRowMap(), true));
  for(unsigned i = 0; i < mesht_elepairs_.size(); i++)
  {
    // evaluate element kappa
    mesht_elepairs_[i]->EvaluateKappa(&eleKappa);
    const int artelegid = mesht_elepairs_[i]->Ele1GID();

    const DRT::Element* ele1 = arterydis_->gElement(artelegid);
    // get element location vector and ownerships
    std::vector<int> lmrow1;
    std::vector<int> lmrowowner1;
    std::vector<int> lmstride;

    // assemble
    ele1->LocationVector(*arterydis_,lmrow1,lmrowowner1,lmstride);
    kappaInv_->SumIntoGlobalValues(eleKappa.Length(),eleKappa.Values(),&lmrow1[0]);
  }
  // invert (pay attention to protruding elements)
  for (int i=0;i<arterydis_->DofRowMap()->NumMyElements();++i)
  {
    const int artdofgid = arterydis_->DofRowMap()->GID(i);
    const double kappaVal = (*kappaInv_)[kappaInv_->Map().LID(artdofgid)];
    if(fabs(kappaVal) > KAPPAINVTOL)
      kappaInv_->ReplaceGlobalValue(artdofgid,0,1.0/kappaVal);
    else
      kappaInv_->ReplaceGlobalValue(artdofgid,0,0.0);
  }

  return;
}

/*------------------------------------------------------------------------*
 | fill the unaffected lenght and initialize curr length kremheller 05/18 |
 *------------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::FillUnaffectedArteryLength()
{

  // the unaffected lenght is the length of 1D elements not changed through deformation,
  // basically if these elements protrude.
  // for each element this length is computed as: ele_length - sum_segments seg_length
  // if the above quantity is bigger than zero, a 1D element protrudes

  // initialize the unaffected and current lengths
  unaffected_length_artery_ = Teuchos::rcp(new Epetra_Vector(*arterydis_->ElementRowMap(), true));
  current_length_artery_    = Teuchos::rcp(new Epetra_Vector(*arterydis_->ElementRowMap(), true));

  // initialize the unaffected length
  for (int i=0;i<arterydis_->ElementRowMap()->NumMyElements();++i)
  {
    const int artelegid = arterydis_->ElementRowMap()->GID(i);
    DRT::Element* artele = arterydis_->gElement(artelegid);

    DRT::Node** artnodes = artele->Nodes();

    DRT::Node* artnode0 = artnodes[0];

    static LINALG::Matrix<3,1> artpos0;
    artpos0(0) = artnode0->X()[0];
    artpos0(1) = artnode0->X()[1];
    artpos0(2) = artnode0->X()[2];

    DRT::Node* artnode1 = artnodes[1];

    static LINALG::Matrix<3,1> artpos1;
    artpos1(0) = artnode1->X()[0];
    artpos1(1) = artnode1->X()[1];
    artpos1(2) = artnode1->X()[2];

    static LINALG::Matrix<3,1> dist;
    dist.Update(1.0, artpos0, -1.0, artpos1, 0.0);
    const double init_length = dist.Norm2();

    unaffected_length_artery_->SumIntoGlobalValues(1,&init_length,&artelegid);
  }

  // subtract the segment lengths only if we evaluate in current configuration
  if(!evaluate_in_ref_config_)
  {
    for(unsigned i = 0; i < mesht_elepairs_.size(); i++)
    {
      double init_segment_length = mesht_elepairs_[i]->ApplyMeshMovement(Teuchos::null, contdis_);
      init_segment_length *= -1.0;
      const int artelegid = mesht_elepairs_[i]->Ele1GID();
      if((arterydis_->gElement(artelegid))->Owner() == myrank_)
        unaffected_length_artery_->SumIntoGlobalValues(1,&(init_segment_length),&artelegid);
    }
  }
  else
    current_length_artery_->Update(1.0,*unaffected_length_artery_,0.0);

  return;
}

/*----------------------------------------------------------------------*
 | create the GID to segment vector                    kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::CreateGIDToSegmentVector()
{

  // fill the GID-to-segment vector
  for(unsigned i = 0; i < mesht_elepairs_.size(); i++)
  {
    const int artelegid = mesht_elepairs_[i]->Ele1GID();
    const double etaA = mesht_elepairs_[i]->EtaA();
    const double etaB = mesht_elepairs_[i]->EtaB();

    gid_to_segment_[artelegid].push_back(etaA);
    gid_to_segment_[artelegid].push_back(etaB);
  }

  // sort and take care of special cases
  for (int i=0;i<arterydis_->ElementColMap()->NumMyElements();++i)
  {
    const int artelegid = arterydis_->ElementColMap()->GID(i);
    if(gid_to_segment_[artelegid].size() > 0) // check if element projects
    {
      // sort
      std::sort(gid_to_segment_[artelegid].begin(), gid_to_segment_[artelegid].end());
      const int end = gid_to_segment_[artelegid].size();
      const double valueAtEnd = gid_to_segment_[artelegid][end-1];
      // end of element lies outside domain
      if(fabs(valueAtEnd - 1.0) > ETAABTOL)
      {
        gid_to_segment_[artelegid].push_back(valueAtEnd);
        gid_to_segment_[artelegid].push_back(1.0);
      }
      const double valueAtBegin = gid_to_segment_[artelegid][0];
      // beginning of element lies outside domain
      if(fabs(valueAtBegin + 1.0) > ETAABTOL)
      {
        gid_to_segment_[artelegid].insert(gid_to_segment_[artelegid].begin(),valueAtBegin);
        gid_to_segment_[artelegid].insert(gid_to_segment_[artelegid].begin(),-1.0);
      }
    }
    // this element does not project
    else
    {
      gid_to_segment_[artelegid].push_back(-1.0);
      gid_to_segment_[artelegid].push_back( 1.0);
    }
    // create the gid to seglength vector
    gid_to_seglength_[artelegid] = std::vector<double>((int)(gid_to_segment_[artelegid].size()/2), -1.0);
  }

  // set segment ID on coupling pairs and fill GID-to-segment length vector
  for (int iele=0;iele<arterydis_->ElementColMap()->NumMyElements();++iele)
  {
    const int artelegid = arterydis_->ElementColMap()->GID(iele);
    std::vector<double> segmentboundaries = gid_to_segment_[artelegid];
    DRT::Element* thisele = arterydis_->gElement(artelegid);
    // TODO: this will not work for higher order artery eles
    const double initlength = GetMaxNodalDistance(thisele, arterydis_);
    for(unsigned int iseg = 0; iseg < segmentboundaries.size()/2; iseg++)
    {
      int id = -1;
      // get EtaA and etaB and calculate initial length
      const double etaA = segmentboundaries[iseg*2];
      const double etaB = segmentboundaries[iseg*2+1];
      gid_to_seglength_[artelegid][iseg] = initlength*(etaB - etaA)/2.0;

      // return also id -> index in mesht_elepairs_ of this segment
      if(IsIdenticalSegment(mesht_elepairs_,
          artelegid,
          etaA,
          etaB,
          id)
          )
        mesht_elepairs_[id]->SetSegmentID(iseg);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate the pairs                                  kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::EvaluateMeshTyingPairs(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat,
    Teuchos::RCP<Epetra_Vector>                 rhs,
    Teuchos::RCP<LINALG::SparseMatrix>          sysmat_cont,
    Teuchos::RCP<LINALG::SparseMatrix>          sysmat_art
    )
{

  // reset
  if(coupling_method_ == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp)
  {
    D_->Zero();
    M_->Zero();
  }

  // resulting discrete element force vectors of the two interacting elements
  std::vector< LINALG::SerialDenseVector > eleforce(2);

  // linearizations
  std::vector< std::vector< LINALG::SerialDenseMatrix > > elestiff ( 2,
      std::vector<LINALG::SerialDenseMatrix>(2) );

  // element mortar coupling matrices
  LINALG::SerialDenseMatrix D_ele;
  LINALG::SerialDenseMatrix M_ele;

  // export from row to col
  if(contdis_->Name() == "porofluid")
    contdis_->SetState("phinp_fluid",phinp_cont_colmap_);
  else if(contdis_->Name() == "scatra")
    contdis_->SetState("phinp",phinp_cont_colmap_);
  else
    dserror("should not happen");

  // evaluate all pairs
  for(unsigned i = 0; i < mesht_elepairs_.size(); i++)
  {
    // get element quantities
    std::vector<double> contelephinp = GetElePhinp(contdis_, mesht_elepairs_[i]->Ele2GID(), phinp_cont_colmap_);
    std::vector<double> artelephinp  = GetElePhinp(arterydis_, mesht_elepairs_[i]->Ele1GID(), phinp_art_colmap_);

    // set on pairs
    mesht_elepairs_[i]->ResetState(contelephinp, artelephinp, contdis_, arterydis_);

    // evaluate
    mesht_elepairs_[i]->Evaluate(
            &eleforce[0],
            &eleforce[1],
            &elestiff[0][0],
            &elestiff[0][1],
            &elestiff[1][0],
            &elestiff[1][1],
            &D_ele,
            &M_ele,
            gid_to_seglength_[mesht_elepairs_[i]->Ele1GID()]);

    // assemble
    AssembleEleForceStiffIntoSystemVectorMatrix(mesht_elepairs_[i]->Ele1GID(),
                                                mesht_elepairs_[i]->Ele2GID(),
                                                eleforce,
                                                elestiff,
                                                sysmat,
                                                rhs,
                                                sysmat_cont,
                                                sysmat_art);

    // in case of MP, assemble D and M
    if(coupling_method_ == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp and num_coupled_dofs_ > 0)
      AssembleDM(mesht_elepairs_[i]->Ele1GID(),
                 mesht_elepairs_[i]->Ele2GID(),
                 D_ele,
                 M_ele);

  }

  // assemble D and M contributions into global force and stiffness
  if(coupling_method_ == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp and num_coupled_dofs_ > 0)
    SumDMIntoGlobalForceStiff(sysmat,
                              rhs,
                              sysmat_cont,
                              sysmat_art);

  return;
}

/*----------------------------------------------------------------------*
 | assemble into global force and stiffness            kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::AssembleEleForceStiffIntoSystemVectorMatrix(
    const int& ele1gid,
    const int& ele2gid,
    std::vector< LINALG::SerialDenseVector > const& elevec,
    std::vector< std::vector< LINALG::SerialDenseMatrix > > const& elemat,
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat,
    Teuchos::RCP<Epetra_Vector>                 rhs,
    Teuchos::RCP<LINALG::SparseMatrix>          sysmat_cont,
    Teuchos::RCP<LINALG::SparseMatrix>          sysmat_art)
{

  const DRT::Element* ele1 = arterydis_->gElement(ele1gid);
  const DRT::Element* ele2 = contdis_->gElement(ele2gid);

  // get element location vector and ownerships
  std::vector<int> lmrow1;
  std::vector<int> lmrow2;
  std::vector<int> lmrowowner1;
  std::vector<int> lmrowowner2;
  std::vector<int> lmstride;

  ele1->LocationVector(*arterydis_,lmrow1,lmrowowner1,lmstride);
  ele2->LocationVector(*contdis_,lmrow2,lmrowowner2,lmstride);

  // assemble rhs
  rhs->SumIntoGlobalValues(elevec[0].Length(),elevec[0].Values(),&lmrow1[0]);
  rhs->SumIntoGlobalValues(elevec[1].Length(),elevec[1].Values(),&lmrow2[0]);

  // assemble matrices
  sysmat_art->Assemble(0, elemat[0][0], lmrow1, lmrowowner1, lmrow1);
  sysmat->Matrix(1,0).Assemble(0, elemat[0][1], lmrow1, lmrowowner1, lmrow2);
  sysmat->Matrix(0,1).Assemble(0, elemat[1][0], lmrow2, lmrowowner2, lmrow1);
  sysmat_cont->Assemble(0, elemat[1][1], lmrow2, lmrowowner2, lmrow2);

  return;

}

/*----------------------------------------------------------------------*
 | assemble D and M into global D and M                kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::AssembleDM(
    const int& ele1gid,
    const int& ele2gid,
    const LINALG::SerialDenseMatrix& D_ele,
    const LINALG::SerialDenseMatrix& M_ele)
{

  const DRT::Element* ele1 = arterydis_->gElement(ele1gid);
  const DRT::Element* ele2 = contdis_->gElement(ele2gid);

  // get element location vector and ownerships
  std::vector<int> lmrow1;
  std::vector<int> lmrow2;
  std::vector<int> lmrowowner1;
  std::vector<int> lmrowowner2;
  std::vector<int> lmstride;

  ele1->LocationVector(*arterydis_,lmrow1,lmrowowner1,lmstride);
  ele2->LocationVector(*contdis_,lmrow2,lmrowowner2,lmstride);

  D_->Assemble(0, D_ele, lmrow1, lmrowowner1, lmrow1);
  M_->Assemble(0, M_ele, lmrow1, lmrowowner1, lmrow2);

  return;
}

/*----------------------------------------------------------------------*
 | sum global D and M into global force and stiff      kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::SumDMIntoGlobalForceStiff(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat,
    Teuchos::RCP<Epetra_Vector>                 rhs,
    Teuchos::RCP<LINALG::SparseMatrix>          sysmat_cont,
    Teuchos::RCP<LINALG::SparseMatrix>          sysmat_art)
{

  bool matlab = false;
  if (matlab)
  {
    //sparse_matrix
    std::string filename_D = "../o/D.dat";
    std::string filename_M = "../o/M.dat";
    LINALG::PrintMatrixInMatlabFormat(filename_D, *(D_->EpetraMatrix()));
    LINALG::PrintMatrixInMatlabFormat(filename_M, *(M_->EpetraMatrix()));
    dserror("exit");
  }

  // complete
  D_->Complete();
  M_->Complete(*contdis_->DofRowMap(), *arterydis_->DofRowMap());

  // get kappa matrix
  Teuchos::RCP<LINALG::SparseMatrix> kappaInvMat = Teuchos::rcp(new LINALG::SparseMatrix(*kappaInv_));
  kappaInvMat->Complete();

  // kappa^{-1}*M
  Teuchos::RCP<LINALG::SparseMatrix> km = LINALG::MLMultiply(*kappaInvMat,false,*M_,false,false,false,true);
  // kappa^{-1}*D
  Teuchos::RCP<LINALG::SparseMatrix> kd = LINALG::MLMultiply(*kappaInvMat,false,*D_,false,false,false,true);

  // D^T*kappa^{-1}*D
  Teuchos::RCP<LINALG::SparseMatrix> dtkd = LINALG::MLMultiply(*D_, true, *kd, false, false, false, true);
  // D^T*kappa^{-1}*M
  Teuchos::RCP<LINALG::SparseMatrix> dtkm = LINALG::MLMultiply(*D_,true,*km,false,false,false,true);
  // M^T*kappa^{-1}*M
  Teuchos::RCP<LINALG::SparseMatrix> mtkm = LINALG::MLMultiply(*M_,true,*km,false,false,false,true);

  // add matrices
  sysmat_art->Add(*dtkd, false, pp_, 1.0);
  sysmat->Matrix(1,0).Add(*dtkm, false, -pp_, 1.0);
  sysmat->Matrix(0,1).Add(*dtkm, true, -pp_, 1.0);
  sysmat_cont->Add(*mtkm, false, pp_, 1.0);

  // add vector
  Teuchos::RCP<Epetra_Vector> art_contribution = Teuchos::rcp(new Epetra_Vector(*arterydis_->DofRowMap()));
  Teuchos::RCP<Epetra_Vector> cont_contribution = Teuchos::rcp(new Epetra_Vector(*contdis_->DofRowMap()));

  // Note: negative since rhs
  // pp*D^T*kappa^{-1}*D*phi_np^art
  dtkd->Multiply(false, *phinp_art_, *art_contribution);
  rhs->Update(-pp_, *globalex_->InsertVector(art_contribution,1), 1.0);

  // -pp*D^T*kappa^{-1}*M*phi_np^cont
  dtkm->Multiply(false, *phinp_cont_, *art_contribution);
  rhs->Update(pp_, *globalex_->InsertVector(art_contribution,1), 1.0);

  // pp*M^T*kappa^{-1}*M*phi_np^cont
  mtkm->Multiply(false, *phinp_cont_, *cont_contribution);
  rhs->Update(-pp_, *globalex_->InsertVector(cont_contribution,0), 1.0);

  // -pp*M^T*kappa^{-1}*D*phi_np^art = -pp*(D^T*kappa^{-1}*M)^T*phi_np^art
  dtkm->Multiply(true, *phinp_art_, *cont_contribution);
  rhs->Update(pp_, *globalex_->InsertVector(cont_contribution,0), 1.0);

  return;
}

/*----------------------------------------------------------------------*
 | get element primary variable                        kremheller 05/18 |
 *----------------------------------------------------------------------*/
std::vector<double> POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::GetElePhinp(
    Teuchos::RCP<DRT::Discretization>  dis,
    const int elegid,
    Teuchos::RCP<const Epetra_Vector> phinp
    )
{
  std::vector<int> lm, lmowner, lmstride;
  std::vector<double> elephinp;

  (dis->gElement(elegid))->LocationVector(*dis,lm,lmowner,lmstride);

  DRT::UTILS::ExtractMyValues(*phinp,elephinp,lm);

  return elephinp;

}

/*----------------------------------------------------------------------*
 | check for duplicate segment                         kremheller 05/18 |
 *----------------------------------------------------------------------*/
bool POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::IsDuplicateSegment(
    const std::vector<Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase> >& msht_elepairs,
    const Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>                possible_duplicate
    )
{

  // we have to sort out duplicate segments, these might occur if the artery element
  // lies exactly between two different 2D/3D-elements

  const double eta_a = possible_duplicate->EtaA();
  const double eta_b = possible_duplicate->EtaB();
  const int ele1gid = possible_duplicate->Ele1GID();
  int elepairID = -1;

  return IsIdenticalSegment(msht_elepairs, ele1gid, eta_a, eta_b, elepairID);
}

/*----------------------------------------------------------------------*
 | check for identical segment                         kremheller 05/18 |
 *----------------------------------------------------------------------*/
bool POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::IsIdenticalSegment(
    const std::vector<Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase> >& msht_elepairs,
    const int&    ele1gid,
    const double& etaA,
    const double& etaB,
    int&          elepairID
    )
{

  for(unsigned i = 0; i < msht_elepairs.size(); i++)
  {
    // first check if ele1-Gid is identical
    if(ele1gid == msht_elepairs[i]->Ele1GID())
      // check if integration segment is the same
      if(fabs(etaA-msht_elepairs[i]->EtaA()) < ETAABTOL && fabs(etaB-msht_elepairs[i]->EtaB()) < ETAABTOL)
      {
        if(PROJOUTPUT)
          std::cout << "found duplicate integration segment" << std::endl;
        elepairID = i;
        return true;
      }
  }

  elepairID = -1;
  return false;
}

/*----------------------------------------------------------------------*
 | factory method to create single pairs               kremheller 05/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::CreateNewArteryCouplingPair(
    std::vector< DRT::Element const *> const& ele_ptrs)
{
  const DRT::Element::DiscretizationType distypeart = ele_ptrs[0]->Shape();
  switch (distypeart)
  {
    case DRT::Element::line2:
    {
      const DRT::Element::DiscretizationType distypecont = ele_ptrs[1]->Shape();
      switch (distypecont)
      {
        case DRT::Element::quad4:
          return Teuchos::rcp (new POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<DRT::Element::line2,DRT::Element::quad4>());
        case DRT::Element::hex8:
          return Teuchos::rcp (new POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<DRT::Element::line2,DRT::Element::hex8>());
        default:
          dserror("only quad4 and hex8 elements supported for continuous elements so far");
      }
      break;
    }
    default:
      dserror("only line 2 elements supported for artery elements so far");
  }

  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 | compute search radius                               kremheller 05/18 |
 *----------------------------------------------------------------------*/
double POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::ComputeSearchRadius()
{

  double maxdist = 0.0;

  // loop over all elements in row map
  for (int i=0;i<arterydis_->ElementRowMap()->NumMyElements();i++)
  {
    // get pointer onto element
    int gid = arterydis_->ElementRowMap()->GID(i);
    DRT::Element* thisele = arterydis_->gElement(gid);

    maxdist = std::max(maxdist, GetMaxNodalDistance(thisele, arterydis_));
  }

  // loop over all elements in row map
  for (int i=0;i<contdis_->ElementRowMap()->NumMyElements();i++)
  {
    // get pointer onto element
    int gid = contdis_->ElementRowMap()->GID(i);
    DRT::Element* thisele = contdis_->gElement(gid);

    maxdist = std::max(maxdist, GetMaxNodalDistance(thisele, contdis_));
  }

  double globalmaxdist = 0.0;
  Comm().MaxAll(&maxdist,&globalmaxdist,1);

  // safety factor
  return BRUTEFORCESAFETY/2.0*globalmaxdist;

}

/*----------------------------------------------------------------------*
 | get maximum nodal distance                          kremheller 05/18 |
 *----------------------------------------------------------------------*/
double POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::GetMaxNodalDistance(
    DRT::Element*                     ele,
    Teuchos::RCP<DRT::Discretization> dis
    )
{

  double maxdist = 0.0;

  // get first node and its position
  int node0_gid = ele->NodeIds()[0];
  DRT::Node* node0 = dis->gNode(node0_gid);

  static LINALG::Matrix<3,1> pos0;
  pos0(0) = node0->X()[0];
  pos0(1) = node0->X()[1];
  pos0(2) = node0->X()[2];

  // loop over second node to numnode to compare distances with first node
  for(int inode = 1; inode < ele->NumNode(); inode++)
  {
    int node1_gid = ele->NodeIds()[inode];
    DRT::Node* node1 = dis->gNode(node1_gid);

    static LINALG::Matrix<3,1> pos1;
    pos1(0) = node1->X()[0];
    pos1(1) = node1->X()[1];
    pos1(2) = node1->X()[2];

    static LINALG::Matrix<3,1> dist;
    dist.Update(1.0, pos0, -1.0, pos1, 0.0);

    maxdist = std::max(maxdist, dist.Norm2());
  }

  return maxdist;
}

/*----------------------------------------------------------------------*
 | setup a global vector                               kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::SetupVector(
Teuchos::RCP<Epetra_Vector>                 vec,
Teuchos::RCP<const Epetra_Vector>           vec_cont,
Teuchos::RCP<const Epetra_Vector>           vec_art
)
{
  // zero out
  vec->PutScalar(0.0);
  // set up global vector
  globalex_->InsertVector(*vec_cont,0,*vec);
  globalex_->InsertVector(*vec_art,1,*vec);


  return;
}

/*----------------------------------------------------------------------*
 | extract single field vectors                        kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::ExtractSingleFieldVectors(
  Teuchos::RCP<const Epetra_Vector>       globalvec,
  Teuchos::RCP<const Epetra_Vector>&      vec_cont,
  Teuchos::RCP<const Epetra_Vector>&      vec_art
  )
{

  // process first field (continuous)
  vec_cont = globalex_->ExtractVector(globalvec, 0);
  // process second field (artery)
  vec_art = globalex_->ExtractVector(globalvec, 1);

  return;
}

/*----------------------------------------------------------------------*
 | check if initial fields on coupled DOFs match       kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::CheckInitialFields(
    Teuchos::RCP<const Epetra_Vector>      vec_cont,
    Teuchos::RCP<const Epetra_Vector>      vec_art
    )
{

  // not performed here since penalty approach will force solution to be
  // equal anyway

  return;
}

/*----------------------------------------------------------------------*
 | artery dof row map                                  kremheller 05/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::ArteryDofRowMap() const
{

  return globalex_->Map(1);
}

/*----------------------------------------------------------------------*
 | artery dof row map                                  kremheller 05/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::DofRowMap() const
{

  return fullmap_;
}

/*----------------------------------------------------------------------*
 | set solution vectors of single fields               kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::SetSolutionVectors(
    Teuchos::RCP<const Epetra_Vector>            phinp_cont,
    Teuchos::RCP<const Epetra_Vector>            phinp_art
    )
{

  phinp_cont_ = phinp_cont;
  phinp_art_ = phinp_art;
  // export from row to col map
  LINALG::Export( *phinp_cont, *phinp_cont_colmap_ );
  LINALG::Export( *phinp_art, *phinp_art_colmap_  );

  return;
}

/*----------------------------------------------------------------------*
 | apply mesh movement                                 kremheller 06/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::ApplyMeshMovement(
    Teuchos::RCP<const Epetra_Vector> disp
    )
{

  // only if we evalute in current configuration
  if(!evaluate_in_ref_config_)
  {
    if(disp == Teuchos::null)
      dserror("Received null pointer on displacement");

    // export from row to col map
    Teuchos::RCP<Epetra_Vector> dispnp_col = Teuchos::rcp(new Epetra_Vector(*contdis_->DofColMap(1), true));
    LINALG::Export( *disp, *dispnp_col );
    current_length_artery_->PutScalar(0.0);
    // apply movement on pairs and fill gid-to-seglength and current_length_artery_
    for(unsigned i = 0; i < mesht_elepairs_.size(); i++)
    {
      const double newsegmentlength = mesht_elepairs_[i]->ApplyMeshMovement(dispnp_col, contdis_);
      const int artelegid = mesht_elepairs_[i]->Ele1GID();
      const int segid = mesht_elepairs_[i]->GetSegmentID();

      gid_to_seglength_[artelegid][segid] = newsegmentlength;

      if((arterydis_->gElement(artelegid))->Owner() == myrank_)
        current_length_artery_->SumIntoGlobalValues(1,&(newsegmentlength),&artelegid);
    }

    // update with unaffected length
    current_length_artery_->Update(1.0, *unaffected_length_artery_, 1.0);
  }

  // set state on artery dis
  arterydis_->SetState(1,"curr_ele_length", current_length_artery_);

  return;
}

/*----------------------------------------------------------------------*
 | print out method                                    kremheller 06/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::PrintOutCouplingMethod() const
{

  std::string name;
  if(coupling_method_ == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp)
    name = "Mortar Penalty";
  else if(coupling_method_ == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::gpts)
    name = "Gauss-Point-To-Segment";
  else
    dserror("unknown coupling method");

  std::cout << "<   Coupling-Method : " << std::left << std::setw(22) << name << "       >" << std::endl;
  std::cout << "<   Penalty         : " << std::left << std::setw(6) << pp_ << "                       >" << std::endl;
  if(evaluate_in_ref_config_)
    std::cout << "<   Moving arteries : No                           >" << std::endl;
  else
    std::cout << "<   Moving arteries : Yes                          >" << std::endl;


  return;

}

/*----------------------------------------------------------------------*
 | print out summary                                   kremheller 06/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::OutputSummary() const
{
  if(myrank_ == 0)
  {
    std::cout << "\nSummary of coupling pairs (segments):" << std::endl;
    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
  }
  Comm().Barrier();
  for(unsigned i = 0; i < mesht_elepairs_.size(); i++)
  {
    std::cout << "Proc " << std::right << std::setw(2) << myrank_ << ": Artery-ele " << std::right << std::setw(5) << mesht_elepairs_[i]->Ele1GID() << ":   [" << std::left << std::setw(11) << mesht_elepairs_[i]->EtaA()
       << "," << std::right << std::setw(11) << mesht_elepairs_[i]->EtaB() << "] <---> continuous-ele " << std::right << std::setw(7) << mesht_elepairs_[i]->Ele2GID() << std::endl;
  }
  Comm().Barrier();
  if(myrank_ == 0)
    std::cout << "\n";

  return;
}

/*----------------------------------------------------------------------*
 | fill vectors as read from input                     kremheller 06/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::FillFunctionAndScaleVectors()
{
  scale_vec_.resize(2);
  funct_vec_.resize(2);

  // get the actual coupled DOFs  ----------------------------------------------------
  // 1) 1D artery discretization
  int    word1;
  std::istringstream scale_art_stream(Teuchos::getNumericStringParameter(meshtyingparams_,"SCALEREAC_ART"));
  while (scale_art_stream >> word1)
    scale_vec_[0].push_back((int)(word1));

  std::istringstream funct_art_stream(Teuchos::getNumericStringParameter(meshtyingparams_,"REACFUNCT_ART"));
  while (funct_art_stream >> word1)
    funct_vec_[0].push_back((int)(word1-1));

  // 2) 2D, 3D continuous field discretization
  std::istringstream scale_cont_stream(Teuchos::getNumericStringParameter(meshtyingparams_,"SCALEREAC_CONT"));
  while (scale_cont_stream >> word1)
    scale_vec_[1].push_back((int)(word1));

  std::istringstream funct_cont_stream(Teuchos::getNumericStringParameter(meshtyingparams_,"REACFUNCT_CONT"));
  while (funct_cont_stream >> word1)
    funct_vec_[1].push_back((int)(word1-1));

  return;
}
