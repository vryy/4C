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
#include <Epetra_FEVector.h>

#include "../drt_porofluidmultiphase/porofluidmultiphase_utils.H"
#include "../drt_porofluidmultiphase_ele/porofluidmultiphase_ele_parameter.H"
#include "../drt_scatra_ele/scatra_ele_parameter_timint.H"
#include "../linalg/linalg_multiply.H"
#include "poromultiphase_scatra_artery_coupling_pair.H"
#include "poromultiphase_scatra_artery_coupling_defines.H"


/*----------------------------------------------------------------------*
 | constructor                                         kremheller 05/18 |
 *----------------------------------------------------------------------*/
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::PoroMultiPhaseScaTraArtCouplLineBased(
    Teuchos::RCP<DRT::Discretization> arterydis, Teuchos::RCP<DRT::Discretization> contdis,
    const Teuchos::ParameterList& couplingparams, const std::string& condname,
    const std::string& artcoupleddofname, const std::string& contcoupleddofname)
    : PoroMultiPhaseScaTraArtCouplBase(
          arterydis, contdis, couplingparams, condname, artcoupleddofname, contcoupleddofname),
      couplingparams_(couplingparams),
      porofluidmanagersset_(false),
      issetup_(false),
      porofluidprob_(false),
      timefacrhs_art_(0.0),
      timefacrhs_cont_(0.0),
      coupling_method_(
          DRT::INPUT::IntegralValue<INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod>(
              couplingparams, "ARTERY_COUPLING_METHOD")),
      pp_(couplingparams_.get<double>("PENALTY")),
      timersearch_(Comm())
{
  // user info
  if (myrank_ == 0)
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
  if (DRT::Problem::Instance()->ProblemType() == prb_porofluidmultiphase)
  {
    evaluate_in_ref_config_ = true;
    porofluidprob_ = true;
  }

  // fill the vectors
  FillFunctionAndScaleVectors();

  // initialize phinp for continuous dis
  phinp_cont_ = Teuchos::rcp(new Epetra_Vector(*contdis_->DofRowMap(), true));
  // initialize phin for continuous dis
  phin_cont_ = Teuchos::rcp(new Epetra_Vector(*contdis_->DofRowMap(), true));
  // initialize phinp for artery dis
  phinp_art_ = Teuchos::rcp(new Epetra_Vector(*arterydis_->DofRowMap(), true));

  // initialize phinp for continuous dis
  zeros_cont_ = Teuchos::rcp(new Epetra_Vector(*contdis_->DofRowMap(), true));
  // initialize phinp for artery dis
  zeros_art_ = Teuchos::rcp(new Epetra_Vector(*arterydis_->DofRowMap(), true));

  // -------------------------------------------------------------------
  // create empty D and M matrices (27 adjacent nodes as 'good' guess)
  // -------------------------------------------------------------------
  D_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *(arterydis_->DofRowMap()), 27, false, true, LINALG::SparseMatrix::FE_MATRIX));
  M_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *(arterydis_->DofRowMap()), 27, false, true, LINALG::SparseMatrix::FE_MATRIX));
  kappaInv_ = Teuchos::rcp(new Epetra_FEVector(*arterydis_->DofRowMap(), true));

  // full map of continous and artery dofs
  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  maps.push_back(Teuchos::rcp(new Epetra_Map(*contdis_->DofRowMap())));
  maps.push_back(Teuchos::rcp(new Epetra_Map(*arterydis_->DofRowMap())));

  fullmap_ = LINALG::MultiMapExtractor::MergeMaps(maps);
  /// dof row map of coupled problem splitted in (field) blocks
  globalex_ = Teuchos::rcp(new LINALG::MultiMapExtractor());
  globalex_->Setup(*fullmap_, maps);

  FEmat_ = Teuchos::rcp(
      new LINALG::SparseMatrix(*fullmap_, 81, true, true, LINALG::SparseMatrix::FE_MATRIX));

  FErhs_ = Teuchos::rcp(new Epetra_FEVector(*fullmap_));

  // check global map extractor
  globalex_->CheckForValidMapExtractor();

  return;
}

/*----------------------------------------------------------------------*
 | setup the strategy                                kremheller 03/19   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::Setup()
{
  // create the pairs
  CreateCouplingPairs();

  // pre-evaluate the pairs
  PreEvaluateCouplingPairs();

  // create the GID to segment vector
  CreateGIDToSegmentVector();

  // fill length of artery elements that is not influenced if the underlying
  // 2D/3D mesh moves (basically protruding artery elements or segments)
  FillUnaffectedArteryLength();

  // print out summary of pairs
  OutputSummary();

  issetup_ = true;

  return;
}

/*----------------------------------------------------------------------*
 | setup the linear system of equations                kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::SetupSystem(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs,
    Teuchos::RCP<LINALG::SparseMatrix> sysmat_cont, Teuchos::RCP<LINALG::SparseMatrix> sysmat_art,
    Teuchos::RCP<const Epetra_Vector> rhs_cont, Teuchos::RCP<const Epetra_Vector> rhs_art,
    Teuchos::RCP<const LINALG::MapExtractor> dbcmap_cont,
    Teuchos::RCP<const LINALG::MapExtractor> dbcmap_art)
{
  if (!issetup_) dserror("Setup() has not been called");

  if (!porofluidmanagersset_)
  {
    // set the right-hand side time factors (we assume constant time step size here)
    SetTimeFacRhs();
    for (unsigned i = 0; i < coupl_elepairs_.size(); i++)
      coupl_elepairs_[i]->SetupFluidManagersAndMaterials(
          contdis_->Name(), timefacrhs_art_, timefacrhs_cont_);
    porofluidmanagersset_ = true;
  }

  // evaluate and assemble the pairs
  // note: OD terms, i.e., (0,1) and (1,0) directly assembled into sysmat
  //       main-diag terms, i.e., (0,0) and (1,1) added into _cont and _art
  //       rhs is added into global rhs
  EvaluateCouplingPairs(sysmat, rhs, sysmat_cont, sysmat_art);

  // add normal part to rhs
  rhs->Update(1.0, *globalex_->InsertVector(rhs_cont, 0), 1.0);
  rhs->Update(1.0, *globalex_->InsertVector(rhs_art, 1), 1.0);

  // apply DBCs
  // 1) on vector
  LINALG::ApplyDirichlettoSystem(rhs, zeros_cont_, *(dbcmap_cont->CondMap()));
  LINALG::ApplyDirichlettoSystem(rhs, zeros_art_, *(dbcmap_art->CondMap()));
  // 2) on main-diag-matrices
  sysmat_cont->ApplyDirichlet(*(dbcmap_cont->CondMap()), true);
  sysmat_art->ApplyDirichlet(*(dbcmap_art->CondMap()), true);
  // 3) on OD-matrices
  sysmat->Matrix(0, 1).Complete(sysmat_art->RangeMap(), sysmat_cont->RangeMap());
  sysmat->Matrix(1, 0).Complete(sysmat_cont->RangeMap(), sysmat_art->RangeMap());
  sysmat->Matrix(0, 1).ApplyDirichlet(*(dbcmap_cont->CondMap()), false);
  sysmat->Matrix(1, 0).ApplyDirichlet(*(dbcmap_art->CondMap()), false);

  // get also the main-diag terms into the global sysmat
  sysmat->Assign(0, 0, LINALG::View, *sysmat_cont);
  sysmat->Assign(1, 1, LINALG::View, *sysmat_art);

  bool matlab = false;
  if (matlab)
  {
    // sparse_matrix
    std::string filename = "../o/mymatrix.dat";
    std::string filename_vc = "../o/myvec.dat";
    LINALG::PrintBlockMatrixInMatlabFormat(filename, *(sysmat));
    LINALG::PrintVectorInMatlabFormat(filename_vc, *rhs, true);
    dserror("exit");
  }

  return;
}

/*----------------------------------------------------------------------*
 | create the pairs                                    kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::CreateCouplingPairs()
{
  // loop over pairs found by search
  std::map<int, std::set<int>>::const_iterator nearbyeleiter;
  for (nearbyeleiter = nearbyelepairs_.begin(); nearbyeleiter != nearbyelepairs_.end();
       ++nearbyeleiter)
  {
    const int artelegid = nearbyeleiter->first;
    std::vector<DRT::Element const*> ele_ptrs(2);
    ele_ptrs[0] = arterydis_->gElement(artelegid);

    std::set<int>::const_iterator secondeleiter;
    for (secondeleiter = nearbyeleiter->second.begin();
         secondeleiter != nearbyeleiter->second.end(); ++secondeleiter)
    {
      const int contelegid = *secondeleiter;
      ele_ptrs[1] = contdis_->gElement(contelegid);

      // construct, init and setup coupling pairs
      Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase> newpair =
          POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::CreateNewArteryCouplingPair(
              ele_ptrs);
      newpair->Init(
          ele_ptrs, couplingparams_, coupleddofs_cont_, coupleddofs_art_, scale_vec_, funct_vec_);

      // add to list of current contact pairs
      coupl_elepairs_.push_back(newpair);
    }
  }

  // output
  int total_numactive_pairs = 0;
  int numactive_pairs = static_cast<int>(coupl_elepairs_.size());
  Comm().SumAll(&numactive_pairs, &total_numactive_pairs, 1);

  if (myrank_ == 0)
    std::cout << "\nFound " << total_numactive_pairs
              << " Artery-to-PoroMultiphaseScatra coupling pairs (segments)" << std::endl;
}

/*----------------------------------------------------------------------*
 | pre-evaluate the pairs and sort out duplicates      kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::PreEvaluateCouplingPairs()
{
  // pre-evaluate
  for (unsigned i = 0; i < coupl_elepairs_.size(); i++) coupl_elepairs_[i]->PreEvaluate();

  // delete the inactive and duplicate pairs
  std::vector<Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>>
      active_coupl_elepairs;
  for (unsigned i = 0; i < coupl_elepairs_.size(); i++)
  {
    const int contelegid = coupl_elepairs_[i]->Ele2GID();
    DRT::Element* contele = contdis_->gElement(contelegid);

    if (coupl_elepairs_[i]->IsActive() &&
        !IsDuplicateSegment(active_coupl_elepairs, coupl_elepairs_[i]) &&
        contele->Owner() == myrank_)
      active_coupl_elepairs.push_back(coupl_elepairs_[i]);
  }

  coupl_elepairs_ = active_coupl_elepairs;

  // output
  int total_numactive_pairs = 0;
  int numactive_pairs = static_cast<int>(coupl_elepairs_.size());
  Comm().SumAll(&numactive_pairs, &total_numactive_pairs, 1);
  if (myrank_ == 0)
    std::cout << "Only " << total_numactive_pairs
              << " Artery-to-PoroMultiphaseScatra coupling pairs (segments) are active"
              << std::endl;

  return;
}

/*------------------------------------------------------------------------*
 | fill the unaffected length and initialize curr length kremheller 05/18 |
 *------------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::FillUnaffectedArteryLength()
{
  // no need to do this
  if (porofluidprob_)
  {
    for (int i = 0; i < arterydis_->ElementColMap()->NumMyElements(); ++i)
    {
      const int artelegid = arterydis_->ElementColMap()->GID(i);
      DRT::Element* artele = arterydis_->gElement(artelegid);

      // TODO: this will not work for higher order artery eles
      const double initlength = POROFLUIDMULTIPHASE::UTILS::GetMaxNodalDistance(artele, arterydis_);
      const int numseg = (int)(gid_to_segment_[artelegid].size() / 2);
      gid_to_seglength_[artelegid].resize(numseg);
      for (int iseg = 0; iseg < numseg; iseg++)
      {
        const double etaA = gid_to_segment_[artelegid][2 * iseg];
        const double etaB = gid_to_segment_[artelegid][2 * iseg + 1];
        gid_to_seglength_[artelegid][iseg] = initlength * (etaB - etaA) / 2.0;

        int id = -1;
        // return also id -> index in coupl_elepairs_ of this segment
        // and set iseg as segment id of the coupling pairs
        if (IsIdenticalSegment(coupl_elepairs_, artelegid, etaA, etaB, id))
          coupl_elepairs_[id]->SetSegmentID(iseg);
      }
    }

    return;
  }

  // the unaffected length is the length of 1D elements not changed through deformation,
  // basically if these elements protrude.
  // for each element this length is computed as: ele_length - sum_segments seg_length
  // if the above quantity is bigger than zero, a 1D element protrudes

  // initialize the unaffected and current lengths
  unaffected_seg_lengths_artery_ =
      Teuchos::rcp(new Epetra_FEVector(*arterydis_->DofRowMap(1), true));
  current_seg_lengths_artery_ = Teuchos::rcp(new Epetra_FEVector(*arterydis_->DofRowMap(1)));

  // set segment ID on coupling pairs and fill the unaffected artery length
  for (int iele = 0; iele < arterydis_->ElementColMap()->NumMyElements(); ++iele)
  {
    const int artelegid = arterydis_->ElementColMap()->GID(iele);
    DRT::Element* thisele = arterydis_->gElement(artelegid);

    // TODO: this will not work for higher order artery eles
    const double initlength = POROFLUIDMULTIPHASE::UTILS::GetMaxNodalDistance(thisele, arterydis_);

    std::vector<double> segmentboundaries = gid_to_segment_[artelegid];
    for (unsigned int iseg = 0; iseg < segmentboundaries.size() / 2; iseg++)
    {
      int id = -1;
      // get EtaA and etaB and calculate initial length
      const double etaA = segmentboundaries[iseg * 2];
      const double etaB = segmentboundaries[iseg * 2 + 1];
      const double seglength = initlength * (etaB - etaA) / 2.0;

      // since we use an FE-vector
      if (thisele->Owner() == myrank_)
      {
        // build the location array
        std::vector<int> seglengthdofs = arterydis_->Dof(1, thisele);
        unaffected_seg_lengths_artery_->SumIntoGlobalValues(1, &seglengthdofs[iseg], &seglength);
      }

      // return also id -> index in coupl_elepairs_ of this segment
      // and set iseg as segment id of the coupling pairs
      if (IsIdenticalSegment(coupl_elepairs_, artelegid, etaA, etaB, id))
        coupl_elepairs_[id]->SetSegmentID(iseg);
    }
  }

  if (unaffected_seg_lengths_artery_->GlobalAssemble(Add, false) != 0)
    dserror("GlobalAssemble of unaffected_seg_lengths_artery_ failed");

  // subtract the segment lengths only if we evaluate in current configuration
  if (!evaluate_in_ref_config_)
  {
    for (unsigned i = 0; i < coupl_elepairs_.size(); i++)
    {
      // get the initial lengths
      double init_segment_length = coupl_elepairs_[i]->ApplyMeshMovement(true, contdis_);
      init_segment_length *= -1.0;

      const int artelegid = coupl_elepairs_[i]->Ele1GID();
      DRT::Element* artele = arterydis_->gElement(artelegid);

      std::vector<int> seglengthdofs = arterydis_->Dof(1, artele);
      const int segid = coupl_elepairs_[i]->GetSegmentID();

      unaffected_seg_lengths_artery_->SumIntoGlobalValues(
          1, &seglengthdofs[segid], &(init_segment_length));
    }
    if (unaffected_seg_lengths_artery_->GlobalAssemble(Add, false) != 0)
      dserror("GlobalAssemble of unaffected_seg_lengths_artery_ failed");
  }
  // the current length is simply the unaffected length
  else
    current_seg_lengths_artery_->Update(1.0, *unaffected_seg_lengths_artery_, 0.0);

  return;
}

/*----------------------------------------------------------------------*
 | create the GID to segment vector                    kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::CreateGIDToSegmentVector()
{
  // fill the GID-to-segment vector
  for (unsigned i = 0; i < coupl_elepairs_.size(); i++)
  {
    const int artelegid = coupl_elepairs_[i]->Ele1GID();
    const int contelegid = coupl_elepairs_[i]->Ele2GID();

    const DRT::Element* contele = contdis_->gElement(contelegid);

    const double etaA = coupl_elepairs_[i]->EtaA();
    const double etaB = coupl_elepairs_[i]->EtaB();

    if (contele->Owner() == myrank_)
    {
      gid_to_segment_[artelegid].push_back(etaA);
      gid_to_segment_[artelegid].push_back(etaB);
    }
    else
      dserror(
          "Something went wrong here, pair in coupling ele pairs where continuous-discretization "
          "element is not owned by this proc.");
  }

  // communicate it to all procs.
  std::vector<int> allproc(Comm().NumProc());
  for (int i = 0; i < Comm().NumProc(); ++i) allproc[i] = i;
  LINALG::Gather<double>(
      gid_to_segment_, gid_to_segment_, (int)allproc.size(), &allproc[0], Comm());

  // sort and take care of special cases
  for (int i = 0; i < arterydis_->ElementColMap()->NumMyElements(); ++i)
  {
    const int artelegid = arterydis_->ElementColMap()->GID(i);
    if (gid_to_segment_[artelegid].size() > 0)  // check if element projects
    {
      // sort
      std::sort(gid_to_segment_[artelegid].begin(), gid_to_segment_[artelegid].end());
      const int end = gid_to_segment_[artelegid].size();
      const double valueAtEnd = gid_to_segment_[artelegid][end - 1];
      // end of element lies outside domain
      if (fabs(valueAtEnd - 1.0) > ETAABTOL)
      {
        gid_to_segment_[artelegid].push_back(valueAtEnd);
        gid_to_segment_[artelegid].push_back(1.0);
      }
      const double valueAtBegin = gid_to_segment_[artelegid][0];
      // beginning of element lies outside domain
      if (fabs(valueAtBegin + 1.0) > ETAABTOL)
      {
        gid_to_segment_[artelegid].insert(gid_to_segment_[artelegid].begin(), valueAtBegin);
        gid_to_segment_[artelegid].insert(gid_to_segment_[artelegid].begin(), -1.0);
      }
    }
    // this element does not project
    else
    {
      gid_to_segment_[artelegid].push_back(-1.0);
      gid_to_segment_[artelegid].push_back(1.0);
    }
  }

  // safety checks
  for (int i = 0; i < arterydis_->ElementColMap()->NumMyElements(); ++i)
  {
    // 1) check if artery element has more than MAXNUMSEGPERELE segments
    const int artelegid = arterydis_->ElementColMap()->GID(i);
    if (gid_to_segment_[artelegid].size() > 2 * MAXNUMSEGPERELE)
      dserror(
          "Artery element %i has %i segments, which is more than the maximum allowed number of %i "
          "segments per artery element, increase MAXNUMSEGPERELE",
          artelegid, (int)(gid_to_segment_[artelegid].size() / 2), MAXNUMSEGPERELE);
    // 2) check if segment has been overlooked
    for (int iseg = 0; iseg < (int)(gid_to_segment_[artelegid].size() / 2) - 1; iseg++)
    {
      if (fabs(gid_to_segment_[artelegid][2 * iseg + 1] -
               gid_to_segment_[artelegid][2 * iseg + 2]) > SMALLESTSEGMENT)
      {
        std::cout << "Problem with segments of artery-element " << artelegid << ":" << std::endl;
        for (int jseg = 0; jseg < (int)(gid_to_segment_[artelegid].size() / 2); jseg++)
        {
          std::cout << "[" << gid_to_segment_[artelegid][2 * jseg] << ", "
                    << gid_to_segment_[artelegid][2 * jseg + 1] << "]" << std::endl;
        }
        dserror(
            "artery element %i has probably not found all possible segments, increase "
            "NUMPROJCHECKS",
            artelegid);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate the pairs                                  kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::EvaluateCouplingPairs(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs,
    Teuchos::RCP<LINALG::SparseMatrix> sysmat_cont, Teuchos::RCP<LINALG::SparseMatrix> sysmat_art)
{
  // reset
  if (coupling_method_ == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp)
  {
    D_->Zero();
    M_->Zero();
    kappaInv_->PutScalar(0.0);
  }

  FEmat_->Zero();
  FErhs_->PutScalar(0.0);

  // resulting discrete element force vectors of the two interacting elements
  std::vector<LINALG::SerialDenseVector> eleforce(2);

  // linearizations
  std::vector<std::vector<LINALG::SerialDenseMatrix>> elestiff(
      2, std::vector<LINALG::SerialDenseMatrix>(2));

  // element mortar coupling matrices
  LINALG::SerialDenseMatrix D_ele;
  LINALG::SerialDenseMatrix M_ele;
  LINALG::SerialDenseVector Kappa_ele;

  // set states
  if (contdis_->Name() == "porofluid")
  {
    contdis_->SetState("phinp_fluid", phinp_cont_);
    contdis_->SetState("phin_fluid", phin_cont_);
    arterydis_->SetState("one_d_artery_pressure", phinp_art_);
    if (not evaluate_in_ref_config_ && not contdis_->HasState(1, "velocity field"))
      dserror("evaluation in current configuration wanted but solid phase velocity not available!");
  }
  else if (contdis_->Name() == "scatra")
  {
    contdis_->SetState("phinp", phinp_cont_);
    arterydis_->SetState("one_d_artery_phinp", phinp_art_);
  }
  else
    dserror(
        "Only porofluid and scatra-discretizations are supported for linebased-coupling so far");

  // evaluate all pairs
  for (unsigned i = 0; i < coupl_elepairs_.size(); i++)
  {
    // reset state on pairs
    coupl_elepairs_[i]->ResetState(contdis_, arterydis_);

    // get the segment lengths
    const std::vector<double> seglengths = GetEleSegmentLengths(coupl_elepairs_[i]->Ele1GID());

    // evaluate
    coupl_elepairs_[i]->Evaluate(&eleforce[0], &eleforce[1], &elestiff[0][0], &elestiff[0][1],
        &elestiff[1][0], &elestiff[1][1], &D_ele, &M_ele, &Kappa_ele, seglengths);

    // assemble
    FEAssembleEleForceStiffIntoSystemVectorMatrix(coupl_elepairs_[i]->Ele1GID(),
        coupl_elepairs_[i]->Ele2GID(), eleforce, elestiff, sysmat, rhs, sysmat_cont, sysmat_art);

    // in case of MP, assemble D, M and Kappa
    if (coupling_method_ == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp and
        num_coupled_dofs_ > 0)
      FEAssembleDMKappa(
          coupl_elepairs_[i]->Ele1GID(), coupl_elepairs_[i]->Ele2GID(), D_ele, M_ele, Kappa_ele);
  }

  if (FErhs_->GlobalAssemble(Add, false) != 0) dserror("GlobalAssemble of right hand side failed");
  rhs->Update(1.0, *FErhs_, 0.0);

  FEmat_->Complete();
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockartery =
      FEmat_->Split<LINALG::DefaultBlockMatrixStrategy>(*globalex_, *globalex_);

  blockartery->Complete();
  sysmat->Matrix(1, 0).Add(blockartery->Matrix(1, 0), false, 1.0, 0.0);
  sysmat->Matrix(0, 1).Add(blockartery->Matrix(0, 1), false, 1.0, 0.0);
  sysmat_cont->Add(blockartery->Matrix(0, 0), false, 1.0, 1.0);
  sysmat_art->Add(blockartery->Matrix(1, 1), false, 1.0, 1.0);

  // assemble D and M contributions into global force and stiffness
  if (coupling_method_ == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp and
      num_coupled_dofs_ > 0)
    SumDMIntoGlobalForceStiff(sysmat, rhs, sysmat_cont, sysmat_art);

  return;
}

/*----------------------------------------------------------------------*
 | FE-assemble into global force and stiffness         kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::
    FEAssembleEleForceStiffIntoSystemVectorMatrix(const int& ele1gid, const int& ele2gid,
        std::vector<LINALG::SerialDenseVector> const& elevec,
        std::vector<std::vector<LINALG::SerialDenseMatrix>> const& elemat,
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs,
        Teuchos::RCP<LINALG::SparseMatrix> sysmat_cont,
        Teuchos::RCP<LINALG::SparseMatrix> sysmat_art)
{
  const DRT::Element* ele1 = arterydis_->gElement(ele1gid);
  const DRT::Element* ele2 = contdis_->gElement(ele2gid);

  // get element location vector and ownerships
  std::vector<int> lmrow1;
  std::vector<int> lmrow2;
  std::vector<int> lmrowowner1;
  std::vector<int> lmrowowner2;
  std::vector<int> lmstride;

  ele1->LocationVector(*arterydis_, lmrow1, lmrowowner1, lmstride);
  ele2->LocationVector(*contdis_, lmrow2, lmrowowner2, lmstride);

  FEmat_->FEAssemble(elemat[0][0], lmrow1, lmrow1);
  FEmat_->FEAssemble(elemat[0][1], lmrow1, lmrow2);
  FEmat_->FEAssemble(elemat[1][0], lmrow2, lmrow1);
  FEmat_->FEAssemble(elemat[1][1], lmrow2, lmrow2);

  FErhs_->SumIntoGlobalValues(elevec[0].Length(), &lmrow1[0], elevec[0].Values());
  FErhs_->SumIntoGlobalValues(elevec[1].Length(), &lmrow2[0], elevec[1].Values());

  return;
}

/*----------------------------------------------------------------------*
 | assemble D, M and kappa into global D, M and kappa  kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::FEAssembleDMKappa(
    const int& ele1gid, const int& ele2gid, const LINALG::SerialDenseMatrix& D_ele,
    const LINALG::SerialDenseMatrix& M_ele, const LINALG::SerialDenseVector& Kappa_ele)
{
  const DRT::Element* ele1 = arterydis_->gElement(ele1gid);
  const DRT::Element* ele2 = contdis_->gElement(ele2gid);

  // get element location vector and ownerships
  std::vector<int> lmrow1;
  std::vector<int> lmrow2;
  std::vector<int> lmrowowner1;
  std::vector<int> lmrowowner2;
  std::vector<int> lmstride;

  ele1->LocationVector(*arterydis_, lmrow1, lmrowowner1, lmstride);
  ele2->LocationVector(*contdis_, lmrow2, lmrowowner2, lmstride);

  D_->FEAssemble(D_ele, lmrow1, lmrow1);
  M_->FEAssemble(M_ele, lmrow1, lmrow2);
  kappaInv_->SumIntoGlobalValues(Kappa_ele.Length(), &lmrow1[0], Kappa_ele.Values());

  return;
}

/*----------------------------------------------------------------------*
 | sum global D and M into global force and stiff      kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::SumDMIntoGlobalForceStiff(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs,
    Teuchos::RCP<LINALG::SparseMatrix> sysmat_cont, Teuchos::RCP<LINALG::SparseMatrix> sysmat_art)
{
  bool matlab = false;
  if (matlab)
  {
    // sparse_matrix
    std::string filename_D = "../o/D.dat";
    std::string filename_M = "../o/M.dat";
    LINALG::PrintMatrixInMatlabFormat(filename_D, *(D_->EpetraMatrix()));
    LINALG::PrintMatrixInMatlabFormat(filename_M, *(M_->EpetraMatrix()));
    dserror("exit");
  }

  // invert
  InvertKappa();

  // complete
  D_->Complete();
  M_->Complete(*contdis_->DofRowMap(), *arterydis_->DofRowMap());

  // get kappa matrix
  Teuchos::RCP<LINALG::SparseMatrix> kappaInvMat =
      Teuchos::rcp(new LINALG::SparseMatrix(*new Epetra_Vector(Copy, *kappaInv_, 0)));
  kappaInvMat->Complete();

  // kappa^{-1}*M
  Teuchos::RCP<LINALG::SparseMatrix> km =
      LINALG::MLMultiply(*kappaInvMat, false, *M_, false, false, false, true);
  // kappa^{-1}*D
  Teuchos::RCP<LINALG::SparseMatrix> kd =
      LINALG::MLMultiply(*kappaInvMat, false, *D_, false, false, false, true);

  // D^T*kappa^{-1}*D
  Teuchos::RCP<LINALG::SparseMatrix> dtkd =
      LINALG::MLMultiply(*D_, true, *kd, false, false, false, true);
  // D^T*kappa^{-1}*M
  Teuchos::RCP<LINALG::SparseMatrix> dtkm =
      LINALG::MLMultiply(*D_, true, *km, false, false, false, true);
  // M^T*kappa^{-1}*M
  Teuchos::RCP<LINALG::SparseMatrix> mtkm =
      LINALG::MLMultiply(*M_, true, *km, false, false, false, true);

  // add matrices
  sysmat_art->Add(*dtkd, false, pp_ * timefacrhs_art_, 1.0);
  sysmat->Matrix(1, 0).Add(*dtkm, false, -pp_ * timefacrhs_art_, 1.0);
  sysmat->Matrix(0, 1).Add(*dtkm, true, -pp_ * timefacrhs_cont_, 1.0);
  sysmat_cont->Add(*mtkm, false, pp_ * timefacrhs_cont_, 1.0);

  // add vector
  Teuchos::RCP<Epetra_Vector> art_contribution =
      Teuchos::rcp(new Epetra_Vector(*arterydis_->DofRowMap()));
  Teuchos::RCP<Epetra_Vector> cont_contribution =
      Teuchos::rcp(new Epetra_Vector(*contdis_->DofRowMap()));

  // Note: all terms are negative since rhs
  // pp*D^T*kappa^{-1}*D*phi_np^art
  dtkd->Multiply(false, *phinp_art_, *art_contribution);
  rhs->Update(-pp_ * timefacrhs_art_, *globalex_->InsertVector(art_contribution, 1), 1.0);

  // -pp*D^T*kappa^{-1}*M*phi_np^cont
  dtkm->Multiply(false, *phinp_cont_, *art_contribution);
  rhs->Update(pp_ * timefacrhs_art_, *globalex_->InsertVector(art_contribution, 1), 1.0);

  // pp*M^T*kappa^{-1}*M*phi_np^cont
  mtkm->Multiply(false, *phinp_cont_, *cont_contribution);
  rhs->Update(-pp_ * timefacrhs_cont_, *globalex_->InsertVector(cont_contribution, 0), 1.0);

  // -pp*M^T*kappa^{-1}*D*phi_np^art = -pp*(D^T*kappa^{-1}*M)^T*phi_np^art
  dtkm->Multiply(true, *phinp_art_, *cont_contribution);
  rhs->Update(pp_ * timefacrhs_cont_, *globalex_->InsertVector(cont_contribution, 0), 1.0);

  return;
}

/*----------------------------------------------------------------------*
 | invert kappa vector                                 kremheller 08/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::InvertKappa()
{
  // global assemble
  if (kappaInv_->GlobalAssemble(Add, false) != 0) dserror("GlobalAssemble of kappaInv_ failed");

  // invert (pay attention to protruding elements)
  for (int i = 0; i < arterydis_->DofRowMap()->NumMyElements(); ++i)
  {
    const int artdofgid = arterydis_->DofRowMap()->GID(i);
    const double kappaVal = (*kappaInv_)[0][kappaInv_->Map().LID(artdofgid)];
    if (fabs(kappaVal) > KAPPAINVTOL)
      kappaInv_->ReplaceGlobalValue(artdofgid, 0, 1.0 / kappaVal);
    else
      kappaInv_->ReplaceGlobalValue(artdofgid, 0, 0.0);
  }

  return;
}

/*----------------------------------------------------------------------*
 | get element segment lengths                         kremheller 09/18 |
 *----------------------------------------------------------------------*/
std::vector<double>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::GetEleSegmentLengths(
    const int artelegid)
{
  if (porofluidprob_) return gid_to_seglength_[artelegid];

  // safety checks
  if (!arterydis_->HasState(1, "curr_seg_lengths")) dserror("cannot get state curr_seg_lengths");

  // build the location array
  DRT::Element* artele = arterydis_->gElement(artelegid);
  std::vector<int> seglengthdofs = arterydis_->Dof(1, artele);

  Teuchos::RCP<const Epetra_Vector> curr_seg_lengths = arterydis_->GetState(1, "curr_seg_lengths");

  std::vector<double> seglengths(MAXNUMSEGPERELE);
  DRT::UTILS::ExtractMyValues(*curr_seg_lengths, seglengths, seglengthdofs);

  return seglengths;
}
/*----------------------------------------------------------------------*
 | check for duplicate segment                         kremheller 05/18 |
 *----------------------------------------------------------------------*/
bool POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::IsDuplicateSegment(
    const std::vector<Teuchos::RCP<
        POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>>& coupl_elepairs,
    const Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>
        possible_duplicate)
{
  // we have to sort out duplicate segments, these might occur if the artery element
  // lies exactly between two different 2D/3D-elements

  const double eta_a = possible_duplicate->EtaA();
  const double eta_b = possible_duplicate->EtaB();
  const int ele1gid = possible_duplicate->Ele1GID();
  int elepairID = -1;

  return IsIdenticalSegment(coupl_elepairs, ele1gid, eta_a, eta_b, elepairID);
}

/*----------------------------------------------------------------------*
 | check for identical segment                         kremheller 05/18 |
 *----------------------------------------------------------------------*/
bool POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::IsIdenticalSegment(
    const std::vector<Teuchos::RCP<
        POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>>& coupl_elepairs,
    const int& ele1gid, const double& etaA, const double& etaB, int& elepairID)
{
  for (unsigned i = 0; i < coupl_elepairs.size(); i++)
  {
    // first check if ele1-Gid is identical
    if (ele1gid == coupl_elepairs[i]->Ele1GID())
      // check if integration segment is the same
      if (fabs(etaA - coupl_elepairs[i]->EtaA()) < ETAABTOL &&
          fabs(etaB - coupl_elepairs[i]->EtaB()) < ETAABTOL)
      {
        if (PROJOUTPUT) std::cout << "found duplicate integration segment" << std::endl;
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
    std::vector<DRT::Element const*> const& ele_ptrs)
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
          return Teuchos::rcp(
              new POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<DRT::Element::line2,
                  DRT::Element::quad4>());
        case DRT::Element::hex8:
          return Teuchos::rcp(
              new POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<DRT::Element::line2,
                  DRT::Element::hex8>());
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
 | setup a global vector                               kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::SetupVector(
    Teuchos::RCP<Epetra_Vector> vec, Teuchos::RCP<const Epetra_Vector> vec_cont,
    Teuchos::RCP<const Epetra_Vector> vec_art)
{
  // zero out
  vec->PutScalar(0.0);
  // set up global vector
  globalex_->InsertVector(*vec_cont, 0, *vec);
  globalex_->InsertVector(*vec_art, 1, *vec);


  return;
}

/*----------------------------------------------------------------------*
 | extract single field vectors                        kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::ExtractSingleFieldVectors(
    Teuchos::RCP<const Epetra_Vector> globalvec, Teuchos::RCP<const Epetra_Vector>& vec_cont,
    Teuchos::RCP<const Epetra_Vector>& vec_art)
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
    Teuchos::RCP<const Epetra_Vector> vec_cont, Teuchos::RCP<const Epetra_Vector> vec_art)
{
  // not performed here since penalty approach will force solution to be
  // equal anyway

  return;
}

/*----------------------------------------------------------------------*
 | artery dof row map                                  kremheller 05/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::ArteryDofRowMap() const
{
  return globalex_->Map(1);
}

/*----------------------------------------------------------------------*
 | artery dof row map                                  kremheller 05/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::DofRowMap() const
{
  return fullmap_;
}

/*----------------------------------------------------------------------*
 | set solution vectors of single fields               kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::SetSolutionVectors(
    Teuchos::RCP<const Epetra_Vector> phinp_cont, Teuchos::RCP<const Epetra_Vector> phin_cont,
    Teuchos::RCP<const Epetra_Vector> phinp_art)
{
  phinp_cont_ = phinp_cont;
  if (phin_cont != Teuchos::null) phin_cont_ = phin_cont;
  phinp_art_ = phinp_art;

  return;
}

/*----------------------------------------------------------------------*
 | apply mesh movement                                 kremheller 06/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::ApplyMeshMovement()
{
  // no need to do this
  if (porofluidprob_) return;

  // only if we evalute in current configuration
  if (!evaluate_in_ref_config_)
  {
    // safety
    if (!contdis_->HasState(1, "dispnp")) dserror("cannot get displacement state");

    // update with unaffected length
    current_seg_lengths_artery_->Update(1.0, *unaffected_seg_lengths_artery_, 0.0);

    // apply movement on pairs and fill gid-to-seglength and current_seg_lengths_artery_
    for (unsigned i = 0; i < coupl_elepairs_.size(); i++)
    {
      const double newsegmentlength = coupl_elepairs_[i]->ApplyMeshMovement(false, contdis_);
      const int artelegid = coupl_elepairs_[i]->Ele1GID();
      const int segid = coupl_elepairs_[i]->GetSegmentID();

      DRT::Element* artele = arterydis_->gElement(artelegid);
      // build the location array
      std::vector<int> seglengthdofs = arterydis_->Dof(1, artele);

      current_seg_lengths_artery_->SumIntoGlobalValues(
          1, &seglengthdofs[segid], &(newsegmentlength));
    }

    if (current_seg_lengths_artery_->GlobalAssemble(Add, false) != 0)
      dserror("GlobalAssemble of current_seg_lengths_artery_ failed");
  }

  // set state on artery dis
  arterydis_->SetState(1, "curr_seg_lengths",
      Teuchos::rcp(new Epetra_Vector(Copy, *current_seg_lengths_artery_, 0)));

  return;
}

/*----------------------------------------------------------------------*
 | print out method                                    kremheller 06/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::PrintOutCouplingMethod() const
{
  std::string name;
  if (coupling_method_ == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp)
    name = "Mortar Penalty";
  else if (coupling_method_ == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::gpts)
    name = "Gauss-Point-To-Segment";
  else
    dserror("unknown coupling method");

  std::cout << "<   Coupling-Method : " << std::left << std::setw(22) << name << "       >"
            << std::endl;
  std::cout << "<   Penalty         : " << std::left << std::setw(6) << pp_
            << "                       >" << std::endl;
  if (evaluate_in_ref_config_)
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
  if (myrank_ == 0)
  {
    std::cout << "\nSummary of coupling pairs (segments):" << std::endl;
    std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << std::endl;
  }
  Comm().Barrier();
  for (unsigned i = 0; i < coupl_elepairs_.size(); i++)
  {
    std::cout << "Proc " << std::right << std::setw(2) << myrank_ << ": Artery-ele " << std::right
              << std::setw(5) << coupl_elepairs_[i]->Ele1GID() << ":   [" << std::left
              << std::setw(11) << coupl_elepairs_[i]->EtaA() << "," << std::right << std::setw(11)
              << coupl_elepairs_[i]->EtaB() << "] <---> continuous-ele " << std::right
              << std::setw(7) << coupl_elepairs_[i]->Ele2GID() << std::endl;
  }
  Comm().Barrier();
  if (myrank_ == 0) std::cout << "\n";

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
  int word1;
  std::istringstream scale_art_stream(
      Teuchos::getNumericStringParameter(couplingparams_, "SCALEREAC_ART"));
  while (scale_art_stream >> word1) scale_vec_[0].push_back((int)(word1));

  std::istringstream funct_art_stream(
      Teuchos::getNumericStringParameter(couplingparams_, "REACFUNCT_ART"));
  while (funct_art_stream >> word1) funct_vec_[0].push_back((int)(word1 - 1));

  // 2) 2D, 3D continuous field discretization
  std::istringstream scale_cont_stream(
      Teuchos::getNumericStringParameter(couplingparams_, "SCALEREAC_CONT"));
  while (scale_cont_stream >> word1) scale_vec_[1].push_back((int)(word1));

  std::istringstream funct_cont_stream(
      Teuchos::getNumericStringParameter(couplingparams_, "REACFUNCT_CONT"));
  while (funct_cont_stream >> word1) funct_vec_[1].push_back((int)(word1 - 1));

  return;
}

/*----------------------------------------------------------------------*
 | set factor for right hand side                      kremheller 03/19 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::SetTimeFacRhs()
{
  // set the right hand side factor
  if (contdis_->Name() == "porofluid")
  {
    DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter* eleparams =
        DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter::Instance("porofluid");
    // artery
    timefacrhs_art_ = 1.0;
    // continuous
    timefacrhs_cont_ = eleparams->TimeFacRhs();
  }
  else if (contdis_->Name() == "scatra")
  {
    DRT::ELEMENTS::ScaTraEleParameterTimInt* eleparams =
        DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra");
    // artery
    timefacrhs_art_ = eleparams->TimeFacRhs();
    // continuous
    timefacrhs_cont_ = eleparams->TimeFacRhs();
  }
  else
    dserror(
        "Only porofluid and scatra-discretizations are supported for linebased-coupling so far");

  return;
}

/*-------------------------------------------------------------------------*
 | set element pairs that are close                       kremheller 03/19 |
 *------------------------------------------------------------------------ */
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::SetNearbyElePairs(
    const std::map<int, std::set<int>>* nearbyelepairs)
{
  nearbyelepairs_ = *nearbyelepairs;
  return;
}
