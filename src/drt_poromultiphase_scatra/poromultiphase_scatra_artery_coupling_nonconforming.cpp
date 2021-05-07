/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for non-conforming coupling between poromultiphase_scatra-
        framework and flow in artery networks including scalar transport

   \level 3

 *----------------------------------------------------------------------*/

#include "poromultiphase_scatra_artery_coupling_nonconforming.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_serialdensevector.H"
#include <Epetra_FEVector.h>
#include <Epetra_IntVector.h>

#include "../drt_porofluidmultiphase/porofluidmultiphase_utils.H"
#include "../drt_porofluidmultiphase_ele/porofluidmultiphase_ele_parameter.H"
#include "../drt_scatra_ele/scatra_ele_parameter_timint.H"
#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_densematrix_communication.H"
#include "../linalg/linalg_utils_sparse_algebra_print.H"
#include "poromultiphase_scatra_artery_coupling_pair.H"
#include "poromultiphase_scatra_artery_coupling_defines.H"
#include "../drt_mat/cnst_1d_art.H"

#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

/*----------------------------------------------------------------------*
 | constructor                                         kremheller 05/18 |
 *----------------------------------------------------------------------*/
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::
    PoroMultiPhaseScaTraArtCouplNonConforming(Teuchos::RCP<DRT::Discretization> arterydis,
        Teuchos::RCP<DRT::Discretization> contdis, const Teuchos::ParameterList& couplingparams,
        const std::string& condname, const std::string& artcoupleddofname,
        const std::string& contcoupleddofname)
    : PoroMultiPhaseScaTraArtCouplBase(
          arterydis, contdis, couplingparams, condname, artcoupleddofname, contcoupleddofname),
      couplingparams_(couplingparams),
      porofluidmanagersset_(false),
      issetup_(false),
      porofluidprob_(false),
      has_varying_diam_(false),
      delete_free_hanging_eles_(DRT::INPUT::IntegralValue<int>(
          DRT::Problem::Instance()->PoroFluidMultiPhaseDynamicParams().sublist("ARTERY COUPLING"),
          "DELETE_FREE_HANGING_ELES")),
      timefacrhs_art_(0.0),
      timefacrhs_cont_(0.0),
      delete_free_hanging_eles_threshold_(DRT::Problem::Instance()
                                              ->PoroFluidMultiPhaseDynamicParams()
                                              .sublist("ARTERY COUPLING")
                                              .get<double>("DELETE_SMALL_FREE_HANGING_COMPS")),
      coupling_method_(
          DRT::INPUT::IntegralValue<INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod>(
              couplingparams, "ARTERY_COUPLING_METHOD")),
      pp_(couplingparams_.get<double>("PENALTY"))
{
}

/*----------------------------------------------------------------------*
 | init the strategy                                 kremheller 07/18   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::Init()
{
  // we do not have a moving mesh
  if (DRT::Problem::Instance()->GetProblemType() == prb_porofluidmultiphase)
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
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::Setup()
{
  // create the pairs
  CreateCouplingPairs();

  // check if varying diameter is used
  if (contdis_->Name() == "porofluid") SetVaryingDiamFlag();
}

/*----------------------------------------------------------------------*
 | evaluate the coupling                               kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::Evaluate(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs)
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
  EvaluateCouplingPairs(sysmat, rhs);
  return;
}

/*----------------------------------------------------------------------*
 | setup the linear system of equations                kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::SetupSystem(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs,
    Teuchos::RCP<LINALG::SparseMatrix> sysmat_cont, Teuchos::RCP<LINALG::SparseMatrix> sysmat_art,
    Teuchos::RCP<const Epetra_Vector> rhs_cont, Teuchos::RCP<const Epetra_Vector> rhs_art,
    Teuchos::RCP<const LINALG::MapExtractor> dbcmap_cont, Teuchos::RCP<const Epetra_Map> dbcmap_art,
    Teuchos::RCP<const Epetra_Map> dbcmap_art_with_collapsed)
{
  // add normal part to rhs
  rhs->Update(1.0, *globalex_->InsertVector(rhs_cont, 0), 1.0);
  rhs->Update(1.0, *globalex_->InsertVector(rhs_art, 1), 1.0);

  // apply DBCs
  // 1) on vector
  LINALG::ApplyDirichlettoSystem(rhs, zeros_cont_, *(dbcmap_cont->CondMap()));
  LINALG::ApplyDirichlettoSystem(rhs, zeros_art_, *(dbcmap_art));
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
  sysmat_cont->Assign(LINALG::View, sysmat->Matrix(0, 0));

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
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::CreateCouplingPairs()
{
  const Teuchos::ParameterList& fluidcouplingparams =
      DRT::Problem::Instance()->PoroFluidMultiPhaseDynamicParams().sublist("ARTERY COUPLING");
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
            coupleddofs_art_, scale_vec_, funct_vec_);

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
    std::cout << "\nFound " << total_numactive_pairs
              << " Artery-to-PoroMultiphaseScatra coupling pairs (segments)" << std::endl;

  // not needed any more
  nearbyelepairs_.clear();
}

/*------------------------------------------------------------------------*
 | set flag if varying diameter has to be calculated     kremheller 04/21 |
 *------------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::SetVaryingDiamFlag()
{
  int has_varying_diam = 0;
  // check all column elements if one of them uses the diameter law by function
  for (int i = 0; i < arterydis_->NumMyColElements(); ++i)
  {
    // pointer to current element
    DRT::Element* actele = arterydis_->lColElement(i);

    // get the artery-material
    Teuchos::RCP<MAT::Cnst_1d_art> arterymat =
        Teuchos::rcp_dynamic_cast<MAT::Cnst_1d_art>(actele->Material());
    if (arterymat == Teuchos::null) dserror("cast to artery material failed");

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
 | evaluate the pairs                                  kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::EvaluateCouplingPairs(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs)
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
    if (has_varying_diam_) ResetIntegratedDiamToZero();
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
    const double integrated_diam =
        coupl_elepairs_[i]->Evaluate(&eleforce[0], &eleforce[1], &elestiff[0][0], &elestiff[0][1],
            &elestiff[1][0], &elestiff[1][1], &D_ele, &M_ele, &Kappa_ele, seglengths);

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

  if (FErhs_->GlobalAssemble(Add, false) != 0) dserror("GlobalAssemble of right hand side failed");
  rhs->Update(1.0, *FErhs_, 0.0);

  FEmat_->Complete();
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockartery =
      FEmat_->Split<LINALG::DefaultBlockMatrixStrategy>(*globalex_, *globalex_);

  blockartery->Complete();
  sysmat->Matrix(1, 0).Add(blockartery->Matrix(1, 0), false, 1.0, 0.0);
  sysmat->Matrix(0, 1).Add(blockartery->Matrix(0, 1), false, 1.0, 0.0);
  sysmat->Matrix(0, 0).Add(blockartery->Matrix(0, 0), false, 1.0, 0.0);
  sysmat->Matrix(1, 1).Add(blockartery->Matrix(1, 1), false, 1.0, 0.0);

  // assemble D and M contributions into global force and stiffness
  if (coupling_method_ == INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp and
      num_coupled_dofs_ > 0)
    SumDMIntoGlobalForceStiff(sysmat, rhs);

  return;
}

/*----------------------------------------------------------------------*
 | FE-assemble into global force and stiffness         kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::
    FEAssembleEleForceStiffIntoSystemVectorMatrix(const int& ele1gid, const int& ele2gid,
        const double& integrated_diam, std::vector<LINALG::SerialDenseVector> const& elevec,
        std::vector<std::vector<LINALG::SerialDenseMatrix>> const& elemat,
        Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs)
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
}

/*----------------------------------------------------------------------*
 | assemble D, M and kappa into global D, M and kappa  kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::FEAssembleDMKappa(
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
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::SumDMIntoGlobalForceStiff(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat, Teuchos::RCP<Epetra_Vector> rhs)
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

  return;
}

/*----------------------------------------------------------------------*
 | invert kappa vector                                 kremheller 08/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::InvertKappa()
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
 | factory method to create single pairs               kremheller 05/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPairBase>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::CreateNewArteryCouplingPair(
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
        case DRT::Element::tet4:
          return Teuchos::rcp(
              new POROMULTIPHASESCATRA::PoroMultiPhaseScatraArteryCouplingPair<DRT::Element::line2,
                  DRT::Element::tet4>());
        default:
          dserror("only quad4, hex8 and tet4 elements supported for continuous elements so far");
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
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::SetupVector(
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
 | artery dof row map                                  kremheller 05/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::ArteryDofRowMap() const
{
  return globalex_->Map(1);
}

/*----------------------------------------------------------------------*
 | artery dof row map                                  kremheller 05/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::DofRowMap() const
{
  return fullmap_;
}

/*----------------------------------------------------------------------*
 | set solution vectors of single fields               kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::SetSolutionVectors(
    Teuchos::RCP<const Epetra_Vector> phinp_cont, Teuchos::RCP<const Epetra_Vector> phin_cont,
    Teuchos::RCP<const Epetra_Vector> phinp_art)
{
  phinp_cont_ = phinp_cont;
  if (phin_cont != Teuchos::null) phin_cont_ = phin_cont;
  phinp_art_ = phinp_art;

  return;
}

/*----------------------------------------------------------------------*
 | access to blood vessel volume fraction              kremheller 08/19 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector>
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::BloodVesselVolumeFraction()
{
  dserror("Not implemented in base class");
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 | print out method                                    kremheller 06/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::PrintOutCouplingMethod() const
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
}

/*----------------------------------------------------------------------*
 | fill vectors as read from input                     kremheller 06/18 |
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
 | set factor for right hand side                      kremheller 03/19 |
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
    dserror(
        "Only porofluid and scatra-discretizations are supported for non-conforming coupling so "
        "far");
}

/*-------------------------------------------------------------------------*
 | set element pairs that are close                       kremheller 03/19 |
 *------------------------------------------------------------------------ */
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraArtCouplNonConforming::SetNearbyElePairs(
    const std::map<int, std::set<int>>* nearbyelepairs)
{
  nearbyelepairs_ = *nearbyelepairs;
}
