/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for non-conforming coupling between poromultiphase_scatra-
        framework and flow in artery networks including scalar transport

   \level 3

 *----------------------------------------------------------------------*/

#include "baci_poromultiphase_scatra_artery_coupling_nonconforming.hpp"

#include "baci_global_data.hpp"
#include "baci_lib_discret.hpp"
#include "baci_linalg_multiply.hpp"
#include "baci_linalg_serialdensevector.hpp"
#include "baci_linalg_utils_sparse_algebra_assemble.hpp"
#include "baci_linalg_utils_sparse_algebra_print.hpp"
#include "baci_mat_cnst_1d_art.hpp"
#include "baci_porofluidmultiphase_ele_parameter.hpp"
#include "baci_porofluidmultiphase_utils.hpp"
#include "baci_poromultiphase_scatra_artery_coupling_defines.hpp"
#include "baci_poromultiphase_scatra_artery_coupling_pair.hpp"
#include "baci_rebalance_binning_based.hpp"
#include "baci_scatra_ele_parameter_timint.hpp"

#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::
    PoroMultiPhaseScaTraArtCouplNonConforming(Teuchos::RCP<DRT::Discretization> arterydis,
        Teuchos::RCP<DRT::Discretization> contdis, const Teuchos::ParameterList& couplingparams,
        const std::string& condname, const std::string& artcoupleddofname,
        const std::string& contcoupleddofname)
    : PoroMultiPhaseScaTraArtCouplBase(
          arterydis, contdis, couplingparams, condname, artcoupleddofname, contcoupleddofname),
      couplingparams_(couplingparams),
      condname_(condname),
      porofluidmanagersset_(false),
      issetup_(false),
      porofluidprob_(false),
      has_varying_diam_(false),
      delete_free_hanging_eles_(CORE::UTILS::IntegralValue<int>(
          GLOBAL::Problem::Instance()->PoroFluidMultiPhaseDynamicParams().sublist(
              "ARTERY COUPLING"),
          "DELETE_FREE_HANGING_ELES")),
      delete_free_hanging_eles_threshold_(GLOBAL::Problem::Instance()
                                              ->PoroFluidMultiPhaseDynamicParams()
                                              .sublist("ARTERY COUPLING")
                                              .get<double>("DELETE_SMALL_FREE_HANGING_COMPS")),
      coupling_method_(
          CORE::UTILS::IntegralValue<INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod>(
              couplingparams, "ARTERY_COUPLING_METHOD")),
      timefacrhs_art_(0.0),
      timefacrhs_cont_(0.0),
      pp_(couplingparams_.get<double>("PENALTY"))
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::Init()
{
  // we do not have a moving mesh
  if (GLOBAL::Problem::Instance()->GetProblemType() == GLOBAL::ProblemType::porofluidmultiphase)
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
  d_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
      *(arterydis_->DofRowMap()), 27, false, true, CORE::LINALG::SparseMatrix::FE_MATRIX));
  m_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
      *(arterydis_->DofRowMap()), 27, false, true, CORE::LINALG::SparseMatrix::FE_MATRIX));
  kappa_inv_ = Teuchos::rcp(new Epetra_FEVector(*arterydis_->DofRowMap(), true));

  // full map of continous and artery dofs
  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  maps.push_back(Teuchos::rcp(new Epetra_Map(*contdis_->DofRowMap())));
  maps.push_back(Teuchos::rcp(new Epetra_Map(*arterydis_->DofRowMap())));

  fullmap_ = CORE::LINALG::MultiMapExtractor::MergeMaps(maps);
  /// dof row map of coupled problem splitted in (field) blocks
  globalex_ = Teuchos::rcp(new CORE::LINALG::MultiMapExtractor());
  globalex_->Setup(*fullmap_, maps);

  FEmat_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
      *fullmap_, 81, true, true, CORE::LINALG::SparseMatrix::FE_MATRIX));

  fe_rhs_ = Teuchos::rcp(new Epetra_FEVector(*fullmap_));

  // check global map extractor
  globalex_->CheckForValidMapExtractor();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::Setup()
{
  // get the coupling method
  auto arterycoupl =
      CORE::UTILS::IntegralValue<INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod>(
          GLOBAL::Problem::Instance()->PoroFluidMultiPhaseDynamicParams().sublist(
              "ARTERY COUPLING"),
          "ARTERY_COUPLING_METHOD");

  // create the pairs
  if (arterycoupl == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::ntp)
  {
    GetCouplingIdsfromInput();
    if (couplingnodes_ntp_.size() == 0)
      FOUR_C_THROW("No 1D Coupling Node Ids found for NTP Coupling");
    CreateCouplingPairsNTP();
  }
  else
  {
    CreateCouplingPairsLineSurfBased();
  }


  // check if varying diameter is used
  if (contdis_->Name() == "porofluid") SetVaryingDiamFlag();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::GetCouplingIdsfromInput()
{
  // get 1D coupling IDs from Input
  std::vector<DRT::Condition*> artCoupcond;

  arterydis_->GetCondition(condname_, artCoupcond);

  couplingnodes_ntp_.resize(artCoupcond.size());

  for (unsigned iter = 0; iter < artCoupcond.size(); ++iter)
  {
    const std::vector<int>* ArteryNodeIds = (artCoupcond[iter])->GetNodes();
    for (auto couplingids : *ArteryNodeIds)
    {
      couplingnodes_ntp_[iter] = couplingids;
      if (myrank_ == 0)
      {
        std::cout << "Artery Coupling Node Id " << iter + 1
                  << " from Input = " << couplingnodes_ntp_[iter] << "\n";
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::Evaluate(
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs)
{
  if (!issetup_) FOUR_C_THROW("Setup() has not been called");

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
  EvaluateCouplingPairs(sysmat, rhs);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::SetupSystem(
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs,
    Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_cont,
    Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_art, Teuchos::RCP<const Epetra_Vector> rhs_cont,
    Teuchos::RCP<const Epetra_Vector> rhs_art,
    Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmap_cont,
    Teuchos::RCP<const Epetra_Map> dbcmap_art,
    Teuchos::RCP<const Epetra_Map> dbcmap_art_with_collapsed)
{
  // add normal part to rhs
  rhs->Update(1.0, *globalex_->InsertVector(rhs_cont, 0), 1.0);
  rhs->Update(1.0, *globalex_->InsertVector(rhs_art, 1), 1.0);

  // apply DBCs
  // 1) on vector
  CORE::LINALG::ApplyDirichletToSystem(*rhs, *zeros_cont_, *(dbcmap_cont->CondMap()));
  CORE::LINALG::ApplyDirichletToSystem(*rhs, *zeros_art_, *(dbcmap_art));
  // 2) on OD-matrices
  sysmat->Matrix(0, 1).Complete(sysmat_art->RangeMap(), sysmat_cont->RangeMap());
  sysmat->Matrix(1, 0).Complete(sysmat_cont->RangeMap(), sysmat_art->RangeMap());
  sysmat->Matrix(0, 1).ApplyDirichlet(*(dbcmap_cont->CondMap()), false);
  sysmat->Matrix(1, 0).ApplyDirichlet(*(dbcmap_art_with_collapsed), false);

  // 3) get also the main-diag terms into the global sysmat
  sysmat->Matrix(0, 0).Add(*sysmat_cont, false, 1.0, 1.0);
  sysmat->Matrix(1, 1).Add(*sysmat_art, false, 1.0, 1.0);
  sysmat->Matrix(0, 0).Complete();
  sysmat->Matrix(1, 1).Complete();
  // and apply DBC
  sysmat->Matrix(0, 0).ApplyDirichlet(*(dbcmap_cont->CondMap()), true);
  sysmat->Matrix(1, 1).ApplyDirichlet(*(dbcmap_art_with_collapsed), true);
  // Assign view to 3D system matrix (such that it now includes also contributions from coupling)
  // this is important! Monolithic algorithms use this matrix
  sysmat_cont->Assign(CORE::LINALG::View, sysmat->Matrix(0, 0));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::
    CreateCouplingPairsLineSurfBased()
{
  const Teuchos::ParameterList& fluidcouplingparams =
      GLOBAL::Problem::Instance()->PoroFluidMultiPhaseDynamicParams().sublist("ARTERY COUPLING");
  // loop over pairs found by search
  std::map<int, std::set<int>>::const_iterator nearbyeleiter;
  int numactive_pairs = 0;
  for (nearbyeleiter = nearbyelepairs_.begin(); nearbyeleiter != nearbyelepairs_.end();
       ++nearbyeleiter)
    numactive_pairs += nearbyeleiter->second.size();

  coupl_elepairs_.resize(numactive_pairs);

  int mypair = 0;
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
      if (ele_ptrs[1]->Owner() == myrank_)
      {
        // construct, init and setup coupling pairs
        Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase> newpair =
            POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::
                CreateNewArteryCouplingPair(ele_ptrs);
        newpair->Init(ele_ptrs, couplingparams_, fluidcouplingparams, coupleddofs_cont_,
            coupleddofs_art_, scale_vec_, funct_vec_, condname_,
            couplingparams_.get<double>("PENALTY"));

        // add to list of current contact pairs
        coupl_elepairs_[mypair] = newpair;
        mypair++;
      }
    }
  }
  coupl_elepairs_.resize(mypair);

  // output
  int total_numactive_pairs = 0;
  numactive_pairs = static_cast<int>(coupl_elepairs_.size());
  Comm().SumAll(&numactive_pairs, &total_numactive_pairs, 1);


  if (myrank_ == 0)
  {
    std::cout << "\nFound " << total_numactive_pairs
              << " Artery-to-PoroMultiphaseScatra coupling pairs (segments)" << std::endl;
  }

  // not needed any more
  nearbyelepairs_.clear();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::CreateCouplingPairsNTP()
{
  const Teuchos::ParameterList& fluidcouplingparams =
      GLOBAL::Problem::Instance()->PoroFluidMultiPhaseDynamicParams().sublist("ARTERY COUPLING");

  int numactive_pairs = std::accumulate(nearbyelepairs_.begin(), nearbyelepairs_.end(), 0,
      [](int a, auto b) { return a + (static_cast<int>(b.second.size())); });

  coupl_elepairs_.resize(numactive_pairs);

  // loop over pairs found by search
  int mypair = 0;
  for (const auto& nearbyeleiter : nearbyelepairs_)
  {
    // create vector of active coupling pairs
    std::vector<DRT::Element const*> ele_ptrs(2);
    // assign artery element
    ele_ptrs[0] = arterydis_->gElement(nearbyeleiter.first);

    // get nodes of artery element
    const DRT::Node* const* artnodes = ele_ptrs[0]->Nodes();

    // loop over nodes of artery element
    for (int i = 0; i < ele_ptrs[0]->NumNode(); i++)
    {
      // loop over prescribed couplings nodes from input
      for (unsigned int j = 0; j < couplingnodes_ntp_.size(); j++)
      {
        // check if artery node is prescribed coupling node
        if (artnodes[i]->Id() == couplingnodes_ntp_[j])
        {
          // get coupling type (ARTERY or AIRWAY ?)
          std::vector<DRT::Condition*> coupcond;
          arterydis_->GetCondition(condname_, coupcond);
          std::string coupling_element_type_ = *(coupcond[j])->Get<std::string>("coupling_type");

          // recompute coupling dofs
          RecomputeCoupledDOFsForNTP(coupcond, j);

          // get penalty parameter
          const auto penalty = *coupcond[j]->Get<double>("PENALTY");

          // get eta (parameter coordinate of corresponding node)
          const int eta_ntp = (i == 0) ? -1 : 1;

          // loop over assigned 2D/3D elements
          for (const auto continuouseleiter : nearbyeleiter.second)
          {
            // assign 2D/3D element
            ele_ptrs[1] = contdis_->gElement(continuouseleiter);

            // only those pairs, where the 3D element is owned by this proc actually evaluated by
            // this proc
            if (ele_ptrs[1]->Owner() == myrank_)
            {
              // construct, init and setup coupling pairs
              Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>
                  newpair = POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::
                      CreateNewArteryCouplingPair(ele_ptrs);
              newpair->Init(ele_ptrs, couplingparams_, fluidcouplingparams, coupleddofs_cont_,
                  coupleddofs_art_, scale_vec_, funct_vec_, condname_, penalty,
                  coupling_element_type_, eta_ntp);
              // add to list of current contact pairs
              coupl_elepairs_[mypair] = newpair;
              mypair++;
            }
          }
        }
      }
    }
  }
  coupl_elepairs_.resize(mypair);

  // output
  int total_numactive_pairs = 0;
  numactive_pairs = static_cast<int>(coupl_elepairs_.size());
  Comm().SumAll(&numactive_pairs, &total_numactive_pairs, 1);


  if (myrank_ == 0)
  {
    std::cout << "\nFound " << total_numactive_pairs
              << " Artery-to-PoroMultiphaseScatra coupling pairs (segments)" << std::endl;
  }

  // not needed any more
  nearbyelepairs_.clear();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::SetVaryingDiamFlag()
{
  int has_varying_diam = 0;
  // check all column elements if one of them uses the diameter law by function
  for (int i = 0; i < arterydis_->NumMyColElements(); ++i)
  {
    // pointer to current element
    DRT::Element* actele = arterydis_->lColElement(i);

    // get the artery-material
    Teuchos::RCP<MAT::Cnst1dArt> arterymat =
        Teuchos::rcp_dynamic_cast<MAT::Cnst1dArt>(actele->Material());
    if (arterymat == Teuchos::null) FOUR_C_THROW("cast to artery material failed");

    if (arterymat->DiameterLaw() == MAT::PAR::ArteryDiameterLaw::diameterlaw_by_function)
    {
      has_varying_diam = 1;
      break;
    }
  }

  // sum over all procs.
  int sum_has_varying_diam = 0;
  Comm().SumAll(&has_varying_diam, &sum_has_varying_diam, 1);
  // if one has a varying diameter set the flag to true
  if (sum_has_varying_diam > 0) has_varying_diam_ = true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::EvaluateCouplingPairs(
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs)
{
  // reset
  if (coupling_method_ == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp)
  {
    d_->Zero();
    m_->Zero();
    kappa_inv_->PutScalar(0.0);
  }

  FEmat_->Zero();
  fe_rhs_->PutScalar(0.0);

  // resulting discrete element force vectors of the two interacting elements
  std::vector<CORE::LINALG::SerialDenseVector> eleforce(2);

  // linearizations
  std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> elestiff(
      2, std::vector<CORE::LINALG::SerialDenseMatrix>(2));

  // element mortar coupling matrices
  CORE::LINALG::SerialDenseMatrix D_ele;
  CORE::LINALG::SerialDenseMatrix M_ele;
  CORE::LINALG::SerialDenseVector Kappa_ele;

  // set states
  if (contdis_->Name() == "porofluid")
  {
    contdis_->SetState("phinp_fluid", phinp_cont_);
    contdis_->SetState("phin_fluid", phin_cont_);
    arterydis_->SetState("one_d_artery_pressure", phinp_art_);
    if (not evaluate_in_ref_config_ && not contdis_->HasState(1, "velocity field"))
      FOUR_C_THROW(
          "evaluation in current configuration wanted but solid phase velocity not available!");
    if (has_varying_diam_) ResetIntegratedDiamToZero();
  }
  else if (contdis_->Name() == "scatra")
  {
    contdis_->SetState("phinp", phinp_cont_);
    arterydis_->SetState("one_d_artery_phinp", phinp_art_);
  }
  else
    FOUR_C_THROW(
        "Only porofluid and scatra-discretizations are supported for linebased-coupling so far");

  // evaluate all pairs
  for (unsigned i = 0; i < coupl_elepairs_.size(); i++)
  {
    // reset state on pairs
    coupl_elepairs_[i]->ResetState(contdis_, arterydis_);

    // get the segment lengths
    const std::vector<double> seglengths = GetEleSegmentLengths(coupl_elepairs_[i]->Ele1GID());

    // evaluate
    const double integrated_diam = coupl_elepairs_[i]->Evaluate(&(eleforce[0]), &(eleforce[1]),
        &(elestiff[0][0]), &(elestiff[0][1]), &(elestiff[1][0]), &(elestiff[1][1]), &D_ele, &M_ele,
        &Kappa_ele, seglengths);

    // assemble
    FEAssembleEleForceStiffIntoSystemVectorMatrix(coupl_elepairs_[i]->Ele1GID(),
        coupl_elepairs_[i]->Ele2GID(), integrated_diam, eleforce, elestiff, sysmat, rhs);

    // in case of MP, assemble D, M and Kappa
    if (coupling_method_ == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp and
        num_coupled_dofs_ > 0)
      FEAssembleDMKappa(
          coupl_elepairs_[i]->Ele1GID(), coupl_elepairs_[i]->Ele2GID(), D_ele, M_ele, Kappa_ele);
  }

  // set artery diameter in material to be able to evalute the 1D elements with varying diameter
  // and evaluate additional linearization of (integrated) element diameters
  if (contdis_->Name() == "porofluid" && has_varying_diam_)
  {
    SetArteryDiamInMaterial();
    EvaluateAdditionalLinearizationofIntegratedDiam();
  }

  if (fe_rhs_->GlobalAssemble(Add, false) != 0)
    FOUR_C_THROW("GlobalAssemble of right hand side failed");
  rhs->Update(1.0, *fe_rhs_, 0.0);

  FEmat_->Complete();
  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> blockartery =
      FEmat_->Split<CORE::LINALG::DefaultBlockMatrixStrategy>(*globalex_, *globalex_);

  blockartery->Complete();
  sysmat->Matrix(1, 0).Add(blockartery->Matrix(1, 0), false, 1.0, 0.0);
  sysmat->Matrix(0, 1).Add(blockartery->Matrix(0, 1), false, 1.0, 0.0);
  sysmat->Matrix(0, 0).Add(blockartery->Matrix(0, 0), false, 1.0, 0.0);
  sysmat->Matrix(1, 1).Add(blockartery->Matrix(1, 1), false, 1.0, 0.0);

  // assemble D and M contributions into global force and stiffness
  if (coupling_method_ == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp and
      num_coupled_dofs_ > 0)
    SumDMIntoGlobalForceStiff(sysmat, rhs);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::
    FEAssembleEleForceStiffIntoSystemVectorMatrix(const int& ele1gid, const int& ele2gid,
        const double& integrated_diam, std::vector<CORE::LINALG::SerialDenseVector> const& elevec,
        std::vector<std::vector<CORE::LINALG::SerialDenseMatrix>> const& elemat,
        Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs)
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

  fe_rhs_->SumIntoGlobalValues(elevec[0].length(), lmrow1.data(), elevec[0].values());
  fe_rhs_->SumIntoGlobalValues(elevec[1].length(), lmrow2.data(), elevec[1].values());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::FEAssembleDMKappa(
    const int& ele1gid, const int& ele2gid, const CORE::LINALG::SerialDenseMatrix& D_ele,
    const CORE::LINALG::SerialDenseMatrix& M_ele, const CORE::LINALG::SerialDenseVector& Kappa_ele)
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

  d_->FEAssemble(D_ele, lmrow1, lmrow1);
  m_->FEAssemble(M_ele, lmrow1, lmrow2);
  kappa_inv_->SumIntoGlobalValues(Kappa_ele.length(), lmrow1.data(), Kappa_ele.values());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::SumDMIntoGlobalForceStiff(
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs)
{
  // invert
  InvertKappa();

  // complete
  d_->Complete();
  m_->Complete(*contdis_->DofRowMap(), *arterydis_->DofRowMap());

  // get kappa matrix
  Teuchos::RCP<CORE::LINALG::SparseMatrix> kappaInvMat =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*new Epetra_Vector(Copy, *kappa_inv_, 0)));
  kappaInvMat->Complete();

  // kappa^{-1}*M
  Teuchos::RCP<CORE::LINALG::SparseMatrix> km =
      CORE::LINALG::MLMultiply(*kappaInvMat, false, *m_, false, false, false, true);
  // kappa^{-1}*D
  Teuchos::RCP<CORE::LINALG::SparseMatrix> kd =
      CORE::LINALG::MLMultiply(*kappaInvMat, false, *d_, false, false, false, true);

  // D^T*kappa^{-1}*D
  Teuchos::RCP<CORE::LINALG::SparseMatrix> dtkd =
      CORE::LINALG::MLMultiply(*d_, true, *kd, false, false, false, true);
  // D^T*kappa^{-1}*M
  Teuchos::RCP<CORE::LINALG::SparseMatrix> dtkm =
      CORE::LINALG::MLMultiply(*d_, true, *km, false, false, false, true);
  // M^T*kappa^{-1}*M
  Teuchos::RCP<CORE::LINALG::SparseMatrix> mtkm =
      CORE::LINALG::MLMultiply(*m_, true, *km, false, false, false, true);

  // add matrices
  sysmat->Matrix(0, 0).Add(*mtkm, false, pp_ * timefacrhs_cont_, 1.0);
  sysmat->Matrix(1, 1).Add(*dtkd, false, pp_ * timefacrhs_art_, 1.0);
  sysmat->Matrix(1, 0).Add(*dtkm, false, -pp_ * timefacrhs_art_, 1.0);
  sysmat->Matrix(0, 1).Add(*dtkm, true, -pp_ * timefacrhs_cont_, 1.0);

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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::InvertKappa()
{
  // global assemble
  if (kappa_inv_->GlobalAssemble(Add, false) != 0)
    FOUR_C_THROW("GlobalAssemble of kappaInv_ failed");

  // invert (pay attention to protruding elements)
  for (int i = 0; i < arterydis_->DofRowMap()->NumMyElements(); ++i)
  {
    const int artdofgid = arterydis_->DofRowMap()->GID(i);
    const double kappaVal = (*kappa_inv_)[0][kappa_inv_->Map().LID(artdofgid)];
    if (fabs(kappaVal) > KAPPAINVTOL)
      kappa_inv_->ReplaceGlobalValue(artdofgid, 0, 1.0 / kappaVal);
    else
      kappa_inv_->ReplaceGlobalValue(artdofgid, 0, 0.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::CreateNewArteryCouplingPair(
    std::vector<DRT::Element const*> const& ele_ptrs)
{
  const CORE::FE::CellType distypeart = ele_ptrs[0]->Shape();
  switch (distypeart)
  {
    case CORE::FE::CellType::line2:
    {
      const CORE::FE::CellType distypecont = ele_ptrs[1]->Shape();
      switch (distypecont)
      {
        case CORE::FE::CellType::quad4:
        {
          switch (GLOBAL::Problem::Instance()->NDim())
          {
            case 1:
              return Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<
                  CORE::FE::CellType::line2, CORE::FE::CellType::quad4, 1>());
            case 2:
              return Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<
                  CORE::FE::CellType::line2, CORE::FE::CellType::quad4, 2>());
            case 3:
              return Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<
                  CORE::FE::CellType::line2, CORE::FE::CellType::quad4, 3>());
            default:
              FOUR_C_THROW("Unsupported dimension %d.", GLOBAL::Problem::Instance()->NDim());
          }
        }
        case CORE::FE::CellType::hex8:
        {
          switch (GLOBAL::Problem::Instance()->NDim())
          {
            case 1:
              return Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<
                  CORE::FE::CellType::line2, CORE::FE::CellType::hex8, 1>());
            case 2:
              return Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<
                  CORE::FE::CellType::line2, CORE::FE::CellType::hex8, 2>());
            case 3:
              return Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<
                  CORE::FE::CellType::line2, CORE::FE::CellType::hex8, 3>());
            default:
              FOUR_C_THROW("Unsupported dimension %d.", GLOBAL::Problem::Instance()->NDim());
          }
        }
        case CORE::FE::CellType::tet4:
        {
          switch (GLOBAL::Problem::Instance()->NDim())
          {
            case 1:
              return Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<
                  CORE::FE::CellType::line2, CORE::FE::CellType::tet4, 1>());
            case 2:
              return Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<
                  CORE::FE::CellType::line2, CORE::FE::CellType::tet4, 2>());
            case 3:
              return Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<
                  CORE::FE::CellType::line2, CORE::FE::CellType::tet4, 3>());
            default:
              FOUR_C_THROW("Unsupported dimension %d.", GLOBAL::Problem::Instance()->NDim());
          }
        }
        case CORE::FE::CellType::tet10:
        {
          switch (GLOBAL::Problem::Instance()->NDim())
          {
            case 1:
              return Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<
                  CORE::FE::CellType::line2, CORE::FE::CellType::tet10, 1>());
            case 2:
              return Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<
                  CORE::FE::CellType::line2, CORE::FE::CellType::tet10, 2>());
            case 3:
              return Teuchos::rcp(new POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<
                  CORE::FE::CellType::line2, CORE::FE::CellType::tet10, 3>());
            default:
              FOUR_C_THROW("Unsupported dimension %d.", GLOBAL::Problem::Instance()->NDim());
          }
        }
        default:
          FOUR_C_THROW(
              "only quad4, hex8, tet4 and tet10 elements supported for continuous elements so far");
      }
    }
    default:
      FOUR_C_THROW("only line 2 elements supported for artery elements so far");
  }

  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::SetupVector(
    Teuchos::RCP<Epetra_Vector> vec, Teuchos::RCP<const Epetra_Vector> vec_cont,
    Teuchos::RCP<const Epetra_Vector> vec_art)
{
  // zero out
  vec->PutScalar(0.0);
  // set up global vector
  globalex_->InsertVector(*vec_cont, 0, *vec);
  globalex_->InsertVector(*vec_art, 1, *vec);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::ExtractSingleFieldVectors(
    Teuchos::RCP<const Epetra_Vector> globalvec, Teuchos::RCP<const Epetra_Vector>& vec_cont,
    Teuchos::RCP<const Epetra_Vector>& vec_art)
{
  // process first field (continuous)
  vec_cont = globalex_->ExtractVector(globalvec, 0);
  // process second field (artery)
  vec_art = globalex_->ExtractVector(globalvec, 1);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::ArteryDofRowMap() const
{
  return globalex_->Map(1);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::DofRowMap() const
{
  return fullmap_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::SetSolutionVectors(
    Teuchos::RCP<const Epetra_Vector> phinp_cont, Teuchos::RCP<const Epetra_Vector> phin_cont,
    Teuchos::RCP<const Epetra_Vector> phinp_art)
{
  phinp_cont_ = phinp_cont;
  if (phin_cont != Teuchos::null) phin_cont_ = phin_cont;
  phinp_art_ = phinp_art;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::BloodVesselVolumeFraction()
{
  FOUR_C_THROW("Not implemented in base class");
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::PrintOutCouplingMethod() const
{
  std::string name;
  if (coupling_method_ == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp)
    name = "Mortar Penalty";
  else if (coupling_method_ == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::gpts)
    name = "Gauss-Point-To-Segment";
  else
    FOUR_C_THROW("unknown coupling method");

  std::cout << "<   Coupling-Method : " << std::left << std::setw(22) << name << "       >"
            << std::endl;
  std::cout << "<   Penalty         : " << std::left << std::setw(6) << pp_
            << "                       >" << std::endl;
  if (evaluate_in_ref_config_)
    std::cout << "<   Moving arteries : No                           >" << std::endl;
  else
    std::cout << "<   Moving arteries : Yes                          >" << std::endl;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::FillFunctionAndScaleVectors()
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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::SetTimeFacRhs()
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
  {
    FOUR_C_THROW(
        "Only porofluid and scatra-discretizations are supported for non-conforming coupling so "
        "far");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::SetNearbyElePairs(
    const std::map<int, std::set<int>>* nearbyelepairs)
{
  nearbyelepairs_ = *nearbyelepairs;
}

FOUR_C_NAMESPACE_CLOSE
