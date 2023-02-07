/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for line-based (non-conforming) coupling between
        poromultiphase_scatra-framework and flow in artery networks
        including scalar transport

   \level 3

 *----------------------------------------------------------------------*/

#include "poromultiphase_scatra_artery_coupling_linebased.H"
#include "lib_globalproblem.H"
#include <Epetra_FEVector.h>
#include <Epetra_IntVector.h>

#include "porofluidmultiphase_utils.H"
#include "linalg_utils_densematrix_communication.H"
#include "poromultiphase_scatra_artery_coupling_pair.H"
#include "poromultiphase_scatra_artery_coupling_defines.H"
#include "mat_cnst_1d_art.H"

#include "linalg_utils_sparse_algebra_manipulation.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::PoroMultiPhaseScaTraArtCouplLineBased(
    Teuchos::RCP<DRT::Discretization> arterydis, Teuchos::RCP<DRT::Discretization> contdis,
    const Teuchos::ParameterList& couplingparams, const std::string& condname,
    const std::string& artcoupleddofname, const std::string& contcoupleddofname)
    : PoroMultiPhaseScaTraArtCouplNonConforming(
          arterydis, contdis, couplingparams, condname, artcoupleddofname, contcoupleddofname),
      maxnumsegperartele_(DRT::Problem::Instance()
                              ->PoroFluidMultiPhaseDynamicParams()
                              .sublist("ARTERY COUPLING")
                              .get<int>("MAXNUMSEGPERARTELE"))
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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::Setup()
{
  // call base class
  POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::Setup();

  // pre-evaluate the pairs
  PreEvaluateCouplingPairs();

  // create the GID to segment vector
  CreateGIDToSegmentVector();

  // fill length of artery elements that is not influenced if the underlying
  // 2D/3D mesh moves (basically protruding artery elements or segments)
  FillUnaffectedArteryLength();

  // fill unaffected integrated diam (basically protruding artery elements or segments)
  if (contdis_->Name() == "porofluid" && has_varying_diam_) FillUnaffectedIntegratedDiam();

  // calculate blood vessel volume fraction (only porofluid needs to do this)
  if (contdis_->Name() == "porofluid" &&
      (DRT::INPUT::IntegralValue<int>(couplingparams_, "OUTPUT_BLOODVESSELVOLFRAC")))
    CalculateBloodVesselVolumeFraction();

  // print out summary of pairs
  if (contdis_->Name() == "porofluid" &&
      (DRT::INPUT::IntegralValue<int>(couplingparams_, "PRINT_OUT_SUMMARY_PAIRS")))
    OutputSummary();

  issetup_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::SetupSystem(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs,
    Teuchos::RCP<LINALG::SparseMatrix> sysmat_cont, Teuchos::RCP<LINALG::SparseMatrix> sysmat_art,
    Teuchos::RCP<const Epetra_Vector> rhs_cont, Teuchos::RCP<const Epetra_Vector> rhs_art,
    Teuchos::RCP<const LINALG::MapExtractor> dbcmap_cont,
    Teuchos::RCP<const LINALG::MapExtractor> dbcmap_art)
{
  // copy vector
  Teuchos::RCP<Epetra_Vector> rhs_art_with_collapsed = Teuchos::rcp(new Epetra_Vector(*rhs_art));
  Teuchos::RCP<Epetra_Map> dbcmap_art_with_collapsed =
      GetAdditionalDBCForCollapsedEles(dbcmap_art, rhs_art_with_collapsed);

  // call base class
  POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::SetupSystem(sysmat, rhs,
      sysmat_cont, sysmat_art, rhs_cont, rhs_art_with_collapsed, dbcmap_cont, dbcmap_art->CondMap(),
      dbcmap_art_with_collapsed);
}

Teuchos::RCP<Epetra_Map>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::GetAdditionalDBCForCollapsedEles(
    Teuchos::RCP<const LINALG::MapExtractor> dbcmap_art,
    Teuchos::RCP<Epetra_Vector> rhs_art_with_collapsed)
{
  // Zero flux is automatically assumed for nodes which border a collapsed element
  // since the respective collapsed element is not evaluated
  // nodes which only border collapsed elements are not evaluated at all, hence, leading to zero
  // rows in global stiffness matrix and to singularity of this matrix
  // here we identify these nodes and set a zero dirichlet boundary condition on them
  // Note that this procedure is equivalent to taking collapsed elements out of the simulation
  // entirely

  int artelematerial = contdis_->Name() == "scatra" ? 1 : 0;
  std::vector<int> mydirichdofs(0);

  const int numrownodes = arterydis_->NumMyRowNodes();
  const Epetra_Map* dofrowmap = arterydis_->DofRowMap();

  for (int inode = 0; inode < numrownodes; ++inode)
  {
    DRT::Node* actnode = arterydis_->lRowNode(inode);
    DRT::Element** eles = actnode->Elements();
    bool all_eles_collapsed = true;
    for (int iele = 0; iele < (actnode->NumElement()); iele++)
    {
      DRT::Element* actele = eles[iele];
      const auto& arterymat =
          Teuchos::rcp_dynamic_cast<const MAT::Cnst_1d_art>(actele->Material(artelematerial));
      if (not arterymat->IsCollapsed())
      {
        all_eles_collapsed = false;
        break;
      }
    }

    // all elements of this node are collapsed
    if (all_eles_collapsed)
    {
      // 1) insert all dofs of this node into dirichlet dof vector
      std::vector<int> dofs = arterydis_->Dof(0, actnode);
      mydirichdofs.insert(mydirichdofs.end(), dofs.begin(), dofs.end());
      // 2) insert the negative value of all dofs of this node into the rhs, with the employed
      // incremental form this will force the value to zero
      for (const auto& mydof : dofs)
        rhs_art_with_collapsed->ReplaceGlobalValue(mydof, 0, -(*phinp_art_)[dofrowmap->LID(mydof)]);
    }
  }

  // build map
  int nummydirichvals = mydirichdofs.size();
  Teuchos::RCP<Epetra_Map> dirichmap =
      Teuchos::rcp(new Epetra_Map(-1, nummydirichvals, mydirichdofs.data(), 0, arterydis_->Comm()));

  // build vector of maps
  std::vector<Teuchos::RCP<const Epetra_Map>> condmaps;
  condmaps.push_back(dirichmap);
  condmaps.push_back(dbcmap_art->CondMap());

  // combined map
  Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);

  return condmerged;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::PreEvaluateCouplingPairs()
{
  // pre-evaluate
  for (unsigned i = 0; i < coupl_elepairs_.size(); i++)
    coupl_elepairs_[i]->PreEvaluate(Teuchos::null);

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

  // the following case takes care of the special occurence where the 1D element lies exactly in
  // between two 2D/3D-elements which are owned by different processors

  // fill the GID-to-segment vector
  std::map<int, std::vector<double>> gid_to_seglength;
  FillGIDToSegmentVector(active_coupl_elepairs, gid_to_seglength);

  // dummy map to collect duplicates in form [ele2gid, eta_a, eta_b, ... ];
  std::map<int, std::vector<double>> duplicates;

  // loop over all artery elements
  for (int i = 0; i < arterydis_->ElementColMap()->NumMyElements(); ++i)
  {
    const int artelegid = arterydis_->ElementColMap()->GID(i);
    if (gid_to_seglength[artelegid].size() > 0)  // check if element projects
    {
      // compare all segment with each other if it might be identical
      for (int iseg = 0; iseg < (int)(gid_to_seglength[artelegid].size() / 2); iseg++)
      {
        const double eta_a = gid_to_seglength[artelegid][2 * iseg];
        const double eta_b = gid_to_seglength[artelegid][2 * iseg + 1];
        for (int jseg = iseg + 1; jseg < (int)(gid_to_seglength[artelegid].size() / 2); jseg++)
        {
          const double eta_a_jseg = gid_to_seglength[artelegid][2 * jseg];
          const double eta_b_jseg = gid_to_seglength[artelegid][2 * jseg + 1];
          // identical segment found
          if (fabs(eta_a - eta_a_jseg) < XIETATOL && fabs(eta_b - eta_b_jseg) < XIETATOL)
          {
            // we need this to get the ele2gid
            int id = -1;
            if (IsIdenticalSegment(active_coupl_elepairs, artelegid, eta_a, eta_b, id))
            {
              const int ele2gid = active_coupl_elepairs[id]->Ele2GID();
              duplicates[artelegid].push_back((double)(ele2gid));
              duplicates[artelegid].push_back(eta_a);
              duplicates[artelegid].push_back(eta_b);
            }
          }
        }
      }
    }
  }

  // communicate the dummy map to all procs.
  std::vector<int> allproc(Comm().NumProc());
  for (int i = 0; i < Comm().NumProc(); ++i) allproc[i] = i;
  LINALG::Gather<double>(duplicates, duplicates, (int)allproc.size(), allproc.data(), Comm());

  // loop over duplicates and delete one duplicate (the one where the 2D/3D element has the larger
  // id)
  std::map<int, std::vector<double>>::iterator it;
  for (it = duplicates.begin(); it != duplicates.end(); it++)
  {
    const int artelegid = it->first;
    std::vector<double> myduplicates = it->second;
    // should always be a multiple of six because we should always find exactly two/four, etc.
    // duplicates
    if (myduplicates.size() % 6 != 0)
      dserror("duplicate vector has size %i, should be multiple of six", myduplicates.size());
    // compare the possible duplicates
    for (int idupl = 0; idupl < (int)(myduplicates.size() / 3); idupl++)
    {
      const double eta_a = myduplicates[3 * idupl + 1];
      const double eta_b = myduplicates[3 * idupl + 2];
      for (int jdupl = idupl + 1; jdupl < (int)(myduplicates.size() / 3); jdupl++)
      {
        const double eta_a_jdupl = myduplicates[3 * jdupl + 1];
        const double eta_b_jdupl = myduplicates[3 * jdupl + 2];
        // duplicate found
        if (fabs(eta_a - eta_a_jdupl) < XIETATOL && fabs(eta_b - eta_b_jdupl) < XIETATOL)
        {
          const int ele_i = (int)myduplicates[3 * idupl];
          const int ele_j = (int)myduplicates[3 * jdupl];
          const int ele_to_be_erased = std::max(ele_i, ele_j);
          int id = -1;
          // delete the duplicate with the larger ele2gid
          if (IsIdenticalSegment(active_coupl_elepairs, artelegid, eta_a, eta_b, id))
          {
            if (active_coupl_elepairs[id]->Ele2GID() == ele_to_be_erased)
            {
              active_coupl_elepairs.erase(active_coupl_elepairs.begin() + id);
            }
          }
        }
      }
    }
  }

  // overwrite the coupling pairs
  coupl_elepairs_ = active_coupl_elepairs;

  // output
  int total_numactive_pairs = 0;
  int numactive_pairs = static_cast<int>(coupl_elepairs_.size());
  Comm().SumAll(&numactive_pairs, &total_numactive_pairs, 1);
  if (myrank_ == 0)
  {
    std::cout << "Only " << total_numactive_pairs
              << " Artery-to-PoroMultiphaseScatra coupling pairs (segments) are active"
              << std::endl;
  }
}

/*------------------------------------------------------------------------*
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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::FillUnaffectedIntegratedDiam()
{
  Teuchos::RCP<Epetra_FEVector> unaffected_diams_artery_row =
      Teuchos::rcp(new Epetra_FEVector(*arterydis_->ElementRowMap(), true));

  for (int i = 0; i < arterydis_->ElementRowMap()->NumMyElements(); ++i)
  {
    const int artelegid = arterydis_->ElementRowMap()->GID(i);
    DRT::Element* artele = arterydis_->gElement(artelegid);

    // TODO: this will not work for higher order artery eles
    const double initlength = POROFLUIDMULTIPHASE::UTILS::GetMaxNodalDistance(artele, arterydis_);

    // first add all contributions int unaffected_diams_artery_row-vector
    Teuchos::RCP<MAT::Cnst_1d_art> arterymat =
        Teuchos::rcp_dynamic_cast<MAT::Cnst_1d_art>(artele->Material());
    if (arterymat == Teuchos::null) dserror("cast to artery material failed");
    const double length_diam = initlength * arterymat->Diam();
    unaffected_diams_artery_row->SumIntoGlobalValues(1, &artelegid, &length_diam);
  }
  // then subtract the coupling pairs to detect protruding parts
  for (unsigned i = 0; i < coupl_elepairs_.size(); i++)
  {
    // get the initial lengths
    double init_segment_length = coupl_elepairs_[i]->ApplyMeshMovement(true, contdis_);
    init_segment_length *= -1.0;

    const int artelegid = coupl_elepairs_[i]->Ele1GID();
    DRT::Element* artele = arterydis_->gElement(artelegid);

    Teuchos::RCP<MAT::Cnst_1d_art> arterymat =
        Teuchos::rcp_dynamic_cast<MAT::Cnst_1d_art>(artele->Material());
    if (arterymat == Teuchos::null) dserror("cast to artery material failed");
    const double length_diam = init_segment_length * arterymat->Diam();
    unaffected_diams_artery_row->SumIntoGlobalValues(1, &artelegid, &length_diam);
  }

  // global assembly and export
  if (unaffected_diams_artery_row->GlobalAssemble(Add, false) != 0)
    dserror("GlobalAssemble of unaffected_seg_lengths_artery_ failed");
  LINALG::Export(*unaffected_diams_artery_row, *unaffected_integrated_diams_artery_col_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::
    CalculateBloodVesselVolumeFraction()
{
  bloodvesselvolfrac_ = Teuchos::rcp(new Epetra_Vector(*contdis_->ElementRowMap(), true));

  double totalvolblood = 0.0;
  // evaluate all pairs
  for (unsigned i = 0; i < coupl_elepairs_.size(); i++)
  {
    const int artelegid = coupl_elepairs_[i]->Ele1GID();
    const int contelegid = coupl_elepairs_[i]->Ele2GID();

    DRT::Element* artele = arterydis_->gElement(artelegid);

    Teuchos::RCP<MAT::Cnst_1d_art> arterymat =
        Teuchos::rcp_dynamic_cast<MAT::Cnst_1d_art>(artele->Material());
    if (arterymat == Teuchos::null) dserror("cast to artery material failed");

    // TODO: this will not work for higher order artery eles
    const double etaA = coupl_elepairs_[i]->EtaA();
    const double etaB = coupl_elepairs_[i]->EtaB();
    const double length = POROFLUIDMULTIPHASE::UTILS::GetMaxNodalDistance(artele, arterydis_);

    const double vol_cont = coupl_elepairs_[i]->CalculateVol2D3D();
    const double vol_art =
        (etaB - etaA) / 2.0 * length * arterymat->Diam() * arterymat->Diam() * M_PI / 4.0;

    totalvolblood += vol_art;

    const double volfrac = vol_art / vol_cont;

    // note: this works since we 2D/3D continuous element of each pair is always owned by this proc.
    int err = bloodvesselvolfrac_->SumIntoGlobalValues(1, &volfrac, &contelegid);
    if (err) dserror("SumIntoGlobalValues failed!");
  }

  // user output
  double vol_sumall = 0.0;
  Comm().SumAll(&totalvolblood, &vol_sumall, 1);
  if (myrank_ == 0)
  {
    std::cout << "\n<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "<    Calculating blood vessel volume fraction      >" << std::endl;
    std::cout << "<    total volume blood:    " << std::setw(5) << vol_sumall
              << "                 >" << std::endl;
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::SetVaryingDiamFlag()
{
  PoroMultiPhaseScaTraArtCouplNonConforming::SetVaryingDiamFlag();

  // set up the required vectors
  if (has_varying_diam_)
  {
    integrated_diams_artery_row_ =
        Teuchos::rcp(new Epetra_FEVector(*arterydis_->ElementRowMap(), true));
    unaffected_integrated_diams_artery_col_ =
        Teuchos::rcp(new Epetra_Vector(*arterydis_->ElementColMap(), true));
    integrated_diams_artery_col_ =
        Teuchos::rcp(new Epetra_Vector(*arterydis_->ElementColMap(), true));
    ele_diams_artery_col_ = Teuchos::rcp(new Epetra_Vector(*arterydis_->ElementColMap(), true));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::CreateGIDToSegmentVector()
{
  // fill the GID-to-segment vector
  FillGIDToSegmentVector(coupl_elepairs_, gid_to_segment_);

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
      if (fabs(valueAtEnd - 1.0) > XIETATOL)
      {
        gid_to_segment_[artelegid].push_back(valueAtEnd);
        gid_to_segment_[artelegid].push_back(1.0);
      }
      const double valueAtBegin = gid_to_segment_[artelegid][0];
      // beginning of element lies outside domain
      if (fabs(valueAtBegin + 1.0) > XIETATOL)
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
    // 1) check if artery element has more than MAXNUMSEGPERARTELE segments
    const int artelegid = arterydis_->ElementColMap()->GID(i);
    if ((int)gid_to_segment_[artelegid].size() > 2 * maxnumsegperartele_)
    {
      dserror(
          "Artery element %i has %i segments, which is more than the maximum allowed number of %i "
          "segments per artery element, increase MAXNUMSEGPERARTELE",
          artelegid, (int)(gid_to_segment_[artelegid].size() / 2), maxnumsegperartele_);
    }
    // 2) check if segment has been overlooked
    for (int iseg = 0; iseg < (int)(gid_to_segment_[artelegid].size() / 2) - 1; iseg++)
    {
      if (fabs(gid_to_segment_[artelegid][2 * iseg + 1] -
               gid_to_segment_[artelegid][2 * iseg + 2]) > XIETATOL)
      {
        std::cout << "Problem with segments of artery-element " << artelegid << ":" << std::endl;
        for (int jseg = 0; jseg < (int)(gid_to_segment_[artelegid].size() / 2); jseg++)
        {
          std::cout << "[" << gid_to_segment_[artelegid][2 * jseg] << ", "
                    << gid_to_segment_[artelegid][2 * jseg + 1] << "]" << std::endl;
        }
        dserror("artery element %i has probably not found all possible segments", artelegid);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::FillGIDToSegmentVector(
    const std::vector<Teuchos::RCP<
        POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>>& coupl_elepairs,
    std::map<int, std::vector<double>>& gid_to_seglength)
{
  // fill the GID-to-segment vector
  for (unsigned i = 0; i < coupl_elepairs.size(); i++)
  {
    const int artelegid = coupl_elepairs[i]->Ele1GID();
    const int contelegid = coupl_elepairs[i]->Ele2GID();

    const DRT::Element* contele = contdis_->gElement(contelegid);

    const double etaA = coupl_elepairs[i]->EtaA();
    const double etaB = coupl_elepairs[i]->EtaB();

    if (contele->Owner() == myrank_)
    {
      gid_to_seglength[artelegid].push_back(etaA);
      gid_to_seglength[artelegid].push_back(etaB);
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
      gid_to_seglength, gid_to_seglength, (int)allproc.size(), allproc.data(), Comm());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::
    FEAssembleEleForceStiffIntoSystemVectorMatrix(const int& ele1gid, const int& ele2gid,
        const double& integrated_diam, std::vector<LINALG::SerialDenseVector> const& elevec,
        std::vector<std::vector<LINALG::SerialDenseMatrix>> const& elemat,
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs)
{
  // call base class
  POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::
      FEAssembleEleForceStiffIntoSystemVectorMatrix(
          ele1gid, ele2gid, integrated_diam, elevec, elemat, sysmat, rhs);

  // also assemble the diameter if necessary
  if (contdis_->Name() == "porofluid" && has_varying_diam_)
    integrated_diams_artery_row_->SumIntoGlobalValues(1, &ele1gid, &(integrated_diam));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::SetArteryDiamInMaterial()
{
  // assemble
  if (integrated_diams_artery_row_->GlobalAssemble(Add, false) != 0)
    dserror("GlobalAssemble of integrated_integrated_diams_artery_row_ failed");

  // export to column format
  LINALG::Export(*integrated_diams_artery_row_, *integrated_diams_artery_col_);

  // fill the vector collecting the element diameter
  FillArteryEleDiamCol();

  // find the free-hanging elements which will be deleted
  std::vector<int> eles_to_be_deleted;
  if (delete_free_hanging_eles_) FindFreeHanging1DElements(eles_to_be_deleted);

  // set the diameter in material
  for (int i = 0; i < arterydis_->NumMyColElements(); ++i)
  {
    // pointer to current element
    DRT::Element* actele = arterydis_->lColElement(i);
    const int elegid = actele->Id();

    double diam = (*ele_diams_artery_col_)[i];

    // set to zero for free-hanging elements
    if (delete_free_hanging_eles_)
    {
      if (std::find(eles_to_be_deleted.begin(), eles_to_be_deleted.end(), elegid) !=
          eles_to_be_deleted.end())
        diam = 0.0;
    }

    // get the artery-material
    Teuchos::RCP<MAT::Cnst_1d_art> arterymat =
        Teuchos::rcp_dynamic_cast<MAT::Cnst_1d_art>(actele->Material());
    if (arterymat == Teuchos::null) dserror("cast to artery material failed");

    // set to zero if collapsed
    if (diam < arterymat->CollapseThreshold())
    {
      // Collapse happens for first time --> inform user
      if (arterymat->Diam() >= arterymat->CollapseThreshold() && actele->Owner() == myrank_)
        std::cout << ">>>>>> Artery element " << actele->Id() << " just collapsed <<<<<<"
                  << std::endl;
      arterymat->SetDiam(0.0);
    }
    else  // otherwise set to calculated diameter
      arterymat->SetDiam(diam);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::ResetIntegratedDiamToZero()
{
  integrated_diams_artery_row_->PutScalar(0.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::FillArteryEleDiamCol()
{
  // reset
  ele_diams_artery_col_->PutScalar(0.0);
  // set the diameter in the vector
  for (int i = 0; i < arterydis_->NumMyColElements(); ++i)
  {
    // pointer to current element
    DRT::Element* actele = arterydis_->lColElement(i);
    const int elegid = actele->Id();

    const std::vector<double> seglengths = GetEleSegmentLengths(elegid);
    const double curr_ele_length = std::accumulate(seglengths.begin(), seglengths.end(), 0.0);
    // diam = int(diam)/length_element
    // also add the unaffected diameter --> diameter of artery elements which protrude
    const double diam =
        ((*integrated_diams_artery_col_)[i] + (*unaffected_integrated_diams_artery_col_)[i]) /
        curr_ele_length;

    ele_diams_artery_col_->ReplaceMyValue(i, 0, diam);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::FindFreeHanging1DElements(
    std::vector<int>& eles_to_be_deleted)
{
  // user info
  if (myrank_ == 0)
  {
    std::cout << "\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<"
                 "<<<<<<<<<<<<<<<<<<<<<<<"
              << std::endl;
    std::cout << ">>>>>>                               Find free-hanging 1D elements               "
                 "               <<<<<<"
              << std::endl;
  }
  // get fully-overlapping discretization
  Teuchos::RCP<DRT::Discretization> artconncompdis =
      POROFLUIDMULTIPHASE::UTILS::CreateFullyOverlappingArteryDiscretization(
          arterydis_, "conn_comp_dis", true);

  // vector to mark visited nodes
  Teuchos::RCP<Epetra_IntVector> visited =
      Teuchos::rcp(new Epetra_IntVector(*artconncompdis->NodeColMap(), true));

  // get fully-overlapping diams vector
  Teuchos::RCP<Epetra_Vector> ele_diams_artery_full_overlap =
      Teuchos::rcp(new Epetra_Vector(*artconncompdis->ElementColMap(), true));
  Teuchos::RCP<Epetra_Vector> ele_diams_artery_row =
      Teuchos::rcp(new Epetra_Vector(*arterydis_->ElementRowMap(), true));
  LINALG::Export(*ele_diams_artery_col_, *ele_diams_artery_row);
  LINALG::Export(*ele_diams_artery_row, *ele_diams_artery_full_overlap);

  // vector of connected components of 1D graph
  std::vector<std::vector<int>> connected_components;
  int num_conn_components = 0;
  int num_conn_components_wo_single_nodes = 0;

  // loop over fully-overlapping discretization
  for (int i = 0; i < artconncompdis->NumMyColNodes(); ++i)
  {
    // pointer to current node
    DRT::Node* actnode = artconncompdis->lColNode(i);
    // if not visited start a new connected component
    if ((*visited)[actnode->LID()] == 0)
    {
      connected_components.push_back(std::vector<int>());
      // recursive call to depth-first search
      DepthFirstSearchUtil(actnode, visited, artconncompdis, ele_diams_artery_full_overlap,
          connected_components[num_conn_components]);
      // single nodes are not of interest as they are detected (and taken out of simulation) anyways
      if (connected_components[num_conn_components].size() > 1)
        num_conn_components_wo_single_nodes++;

      num_conn_components++;
    }
  }

  // user info
  if (myrank_ == 0 && num_conn_components_wo_single_nodes > 1)
  {
    std::cout << "found " << num_conn_components_wo_single_nodes << " connected components"
              << std::endl;
  }

  // loop over all connected components
  for (unsigned int i = 0; i < connected_components.size(); ++i)
  {
    int conn_comp_size = connected_components[i].size();
    // single nodes are not of interest as they are detected anyways
    if (conn_comp_size > 1)
    {
      // user info
      if (myrank_ == 0)
        std::cout << "connected_component with ID " << i << " of size: " << conn_comp_size
                  << std::endl;

      // check if any of the nodes of this connected component has a Dirichlet BC
      DRT::Condition* dirich = nullptr;
      for (int j = 0; j < conn_comp_size; ++j)
      {
        DRT::Node* mynode = artconncompdis->gNode((connected_components[i])[j]);
        dirich = mynode->GetCondition("Dirichlet");
        if (dirich != nullptr)
        {
          if (myrank_ == 0)
            std::cout << "   ---> has at least one Dirichlet boundary condition" << std::endl;
          break;
        }
      }

      // if no node of this connected component has a DBC or if it is smaller than the
      // user-specified threshold, all its elements are taken out
      if (dirich == nullptr or conn_comp_size < (int)(delete_free_hanging_eles_threshold_ *
                                                      artconncompdis->NumGlobalNodes()))
      {
        // get the elements which have to be deleted
        for (int j = 0; j < conn_comp_size; ++j)
        {
          DRT::Node* mynode = artconncompdis->gNode((connected_components[i])[j]);
          DRT::Element** myeles = mynode->Elements();
          for (int i_element = 0; i_element < mynode->NumElement(); i_element++)
            eles_to_be_deleted.push_back(myeles[i_element]->Id());
        }
        // user info
        if (myrank_ == 0)
        {
          if (dirich == nullptr)
            std::cout
                << "   ---> has no Dirichlet boundary condition --> its elements will be taken out"
                << std::endl;
          if (dirich != nullptr and conn_comp_size < (int)(delete_free_hanging_eles_threshold_ *
                                                           artconncompdis->NumGlobalNodes()))
            std::cout << "   ---> smaller than threshold size of "
                      << (int)(delete_free_hanging_eles_threshold_ *
                               artconncompdis->NumGlobalNodes())
                      << " --> its elements will be taken out" << std::endl;
        }
      }
    }
  }

  // user info
  if (myrank_ == 0)
  {
    std::cout << "\n>>>>>>                           End of Find free-hanging 1D elements          "
                 "                 <<<<<<"
              << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<"
                 "<<<<<<<<<<<<<<<<<<<<<<<\n"
              << std::endl;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::DepthFirstSearchUtil(
    DRT::Node* actnode, Teuchos::RCP<Epetra_IntVector> visited,
    Teuchos::RCP<DRT::Discretization> artconncompdis,
    Teuchos::RCP<const Epetra_Vector> ele_diams_artery_full_overlap,
    std::vector<int>& this_connected_comp)
{
  // mark this node visited and add it to this connected component
  const int lid = visited->Map().LID(actnode->Id());
  (*visited)[lid] = 1;
  this_connected_comp.push_back(actnode->Id());

  // check all adjacent elements (edges)
  DRT::Element** Elements = actnode->Elements();
  for (int i_element = 0; i_element < actnode->NumElement(); i_element++)
  {
    DRT::Node** Nodes = Elements[i_element]->Nodes();

    // get diameter
    const double diam = (*ele_diams_artery_full_overlap)[Elements[i_element]->LID()];

    // get the artery-material
    Teuchos::RCP<MAT::Cnst_1d_art> arterymat =
        Teuchos::rcp_dynamic_cast<MAT::Cnst_1d_art>(Elements[i_element]->Material());
    if (arterymat == Teuchos::null) dserror("cast to artery material failed");

    // if the element is not collapsed it is connected to this node and we continue with the
    // depth-first search with all nodes of this element
    if (diam >= arterymat->CollapseThreshold())
    {
      for (int i_node = 0; i_node < Elements[i_element]->NumNode(); i_node++)
      {
        if ((*visited)[Nodes[i_node]->LID()] == 0)
          DepthFirstSearchUtil(Nodes[i_node], visited, artconncompdis,
              ele_diams_artery_full_overlap, this_connected_comp);
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::
    EvaluateAdditionalLinearizationofIntegratedDiam()
{
  // linearizations
  std::vector<LINALG::SerialDenseMatrix> elestiff(2);

  // evaluate all pairs
  for (unsigned i = 0; i < coupl_elepairs_.size(); i++)
  {
    // only needed if varying diameter is set for this pair
    if (coupl_elepairs_[i]->DiamFunctionActive())
    {
      // evaluate
      coupl_elepairs_[i]->EvaluateAdditionalLinearizationofIntegratedDiam(
          &(elestiff[0]), &(elestiff[1]));

      // and FE-Assemble
      const int ele1gid = coupl_elepairs_[i]->Ele1GID();
      const int ele2gid = coupl_elepairs_[i]->Ele2GID();
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

      FEmat_->FEAssemble(elestiff[0], lmrow1, lmrow1);
      FEmat_->FEAssemble(elestiff[1], lmrow1, lmrow2);
    }
  }
}

/*----------------------------------------------------------------------*
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

  std::vector<double> seglengths(maxnumsegperartele_);
  DRT::UTILS::ExtractMyValues(*curr_seg_lengths, seglengths, seglengthdofs);

  return seglengths;
}

/*----------------------------------------------------------------------*
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
      if (fabs(etaA - coupl_elepairs[i]->EtaA()) < XIETATOL &&
          fabs(etaB - coupl_elepairs[i]->EtaB()) < XIETATOL)
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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::PrintOutCouplingMethod() const
{
  std::cout << "<   Line-based formulation                         >" << std::endl;
  PoroMultiPhaseScaTraArtCouplNonConforming::PrintOutCouplingMethod();
}

/*----------------------------------------------------------------------*
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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplLineBased::BloodVesselVolumeFraction()
{
  return bloodvesselvolfrac_;
}
