/*----------------------------------------------------------------------*/
/*! \file
\brief mortar coupling terms of ehl

\level 3

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  headers                                                             |
 *----------------------------------------------------------------------*/
#include "4C_adapter_coupling_ehl_mortar.hpp"

#include "4C_contact_friction_node.hpp"
#include "4C_contact_interface.hpp"
#include "4C_contact_lagrange_strategy_tsi.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

FOUR_C_NAMESPACE_OPEN

ADAPTER::CouplingEhlMortar::CouplingEhlMortar(int spatial_dimension,
    Teuchos::ParameterList mortar_coupling_params, Teuchos::ParameterList contact_dynamic_params,
    CORE::FE::ShapeFunctionType shape_function_type)
    : CouplingNonLinMortar(
          spatial_dimension, mortar_coupling_params, contact_dynamic_params, shape_function_type),
      contact_regularization_(CORE::UTILS::IntegralValue<int>(
          GLOBAL::Problem::Instance()->ContactDynamicParams(), "REGULARIZED_NORMAL_CONTACT")),
      regularization_thickness_(GLOBAL::Problem::Instance()->ContactDynamicParams().get<double>(
          "REGULARIZATION_THICKNESS")),
      regularization_compliance_(GLOBAL::Problem::Instance()->ContactDynamicParams().get<double>(
          "REGULARIZATION_STIFFNESS"))
{
  if (Teuchos::getIntegralValue<INPAR::MORTAR::ParallelRedist>(
          GLOBAL::Problem::Instance()->MortarCouplingParams().sublist("PARALLEL REDISTRIBUTION"),
          "PARALLEL_REDIST") != INPAR::MORTAR::ParallelRedist::redist_none)
    FOUR_C_THROW(
        "EHL does not support parallel redistribution. Set \"PARALLEL_REDIST none\" in section "
        "\"MORTAR COUPLING\"");

  if (contact_regularization_)
    if (regularization_compliance_ <= 0. || regularization_thickness_ <= 0.)
      FOUR_C_THROW("need positive REGULARIZATION_THICKNESS and REGULARIZATION_STIFFNESS");
  if (contact_regularization_) regularization_compliance_ = 1. / regularization_compliance_;
  if (CORE::UTILS::IntegralValue<int>(GLOBAL::Problem::Instance()->ContactDynamicParams(),
          "REGULARIZED_NORMAL_CONTACT") == true &&
      CORE::UTILS::IntegralValue<bool>(
          GLOBAL::Problem::Instance()->ElastoHydroDynamicParams(), "DRY_CONTACT_MODEL") == false)
    FOUR_C_THROW("for dry contact model you need REGULARIZED_NORMAL_CONTACT and DRY_CONTACT_MODEL");
}

/*----------------------------------------------------------------------*
 |  read mortar condition                                               |
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingEhlMortar::ReadMortarCondition(Teuchos::RCP<DRT::Discretization> masterdis,
    Teuchos::RCP<DRT::Discretization> slavedis, std::vector<int> coupleddof,
    const std::string& couplingcond, Teuchos::ParameterList& input,
    std::map<int, DRT::Node*>& mastergnodes, std::map<int, DRT::Node*>& slavegnodes,
    std::map<int, Teuchos::RCP<DRT::Element>>& masterelements,
    std::map<int, Teuchos::RCP<DRT::Element>>& slaveelements)
{
  ADAPTER::CouplingNonLinMortar::ReadMortarCondition(masterdis, slavedis, coupleddof, couplingcond,
      input, mastergnodes, slavegnodes, masterelements, slaveelements);

  input.set<int>("PROBTYPE", INPAR::CONTACT::ehl);
}

void ADAPTER::CouplingEhlMortar::Setup(Teuchos::RCP<DRT::Discretization> masterdis,
    Teuchos::RCP<DRT::Discretization> slavedis, std::vector<int> coupleddof,
    const std::string& couplingcond)
{
  ADAPTER::CouplingNonLinMortar::Setup(masterdis, slavedis, coupleddof, couplingcond);
  z_ = Teuchos::rcp(new Epetra_Vector(*interface_->SlaveRowDofs(), true));
  fscn_ = Teuchos::rcp(new Epetra_Vector(*interface_->SlaveRowDofs(), true));

  INPAR::CONTACT::FrictionType ftype = CORE::UTILS::IntegralValue<INPAR::CONTACT::FrictionType>(
      GLOBAL::Problem::Instance()->ContactDynamicParams(), "FRICTION");

  std::vector<DRT::Condition*> ehl_conditions(0);
  masterdis->GetCondition(couplingcond, ehl_conditions);
  std::vector<const std::string*> sides((int)ehl_conditions.size());
  double fr_coeff = -1.;
  for (int i = 0; i < (int)ehl_conditions.size(); ++i)
  {
    [[maybe_unused]] const int group1id = ehl_conditions[i]->Get<int>("Interface ID");
    const auto fr = ehl_conditions[i]->Get<double>("FrCoeffOrBound");
    if (fr != ehl_conditions[0]->Get<double>("FrCoeffOrBound"))
      FOUR_C_THROW("inconsistency in friction coefficients");
    fr_coeff = fr;
  }

  switch (ftype)
  {
    case INPAR::CONTACT::friction_tresca:
      FOUR_C_THROW("no tresca friction supported");
      break;
    case INPAR::CONTACT::friction_none:
      break;
    case INPAR::CONTACT::friction_coulomb:
      interface_->InterfaceParams().set<double>("FRCOEFF", fr_coeff);
      interface_->InterfaceParams().set<double>("FRBOUND", -1.);
      break;
    default:
      FOUR_C_THROW("don't know what to do with this friction type");
      break;
  }
}


/*----------------------------------------------------------------------*
 |  perform interface integration and assembly                          |
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingEhlMortar::Integrate(Teuchos::RCP<const Epetra_Vector> disp, const double dt)
{
  // safety check
  CheckSetup();

  // return if this state has already been evaluated
  if (AlreadyEvaluated(disp)) return;

  // set current displ state
  interface_->SetState(MORTAR::state_new_displacement, *disp);

  // init internal data
  interface_->Initialize();
  interface_->SetElementAreas();
  // call interface evaluate (d,m,gap...)
  interface_->Evaluate();

  // some first assemblies, that don't require any additional states
  D_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(*slavedofrowmap_, 81, false, false));
  M_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(*slavedofrowmap_, 81, false, false));
  interface_->AssembleDM(*D_, *M_);
  D_->Complete();
  M_->Complete(*masterdofrowmap_, *slavedofrowmap_);
  N_->Complete(*smdofrowmap_, *slavedofrowmap_);
  AssembleRealGap();
  AssembleRealGapDeriv();
  AssembleNormals();
  AssembleNormalsDeriv();
  AssembleSurfGrad();
  AssembleInterfaceVelocities(dt);

  // save that state as the last evaluated one
  evaluated_state_ = Teuchos::rcp(new Epetra_Vector(*disp));

  // all done
  return;
}

/*----------------------------------------------------------------------*
 |  perform interface integration and assembly                          |
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingEhlMortar::CondenseContact(
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> sysmat,
    Teuchos::RCP<Epetra_Vector>& combined_RHS, Teuchos::RCP<const Epetra_Vector> disp,
    const double dt)
{
  const double alphaf_ = 0.;  // statics!
  const INPAR::CONTACT::ConstraintDirection& constr_direction_ =
      CORE::UTILS::IntegralValue<INPAR::CONTACT::ConstraintDirection>(
          Interface()->InterfaceParams(), "CONSTRAINT_DIRECTIONS");

  // return if this state has already been evaluated
  if (not AlreadyEvaluated(disp)) Integrate(disp, dt);

  // get the relative movement for frictional contact
  EvaluateRelMov();

  // update active set
  as_converged_ = interface_->UpdateActiveSetSemiSmooth();
  interface_->BuildActiveSet();

  // assemble the constraint lines for the active contact nodes
  Teuchos::RCP<CORE::LINALG::SparseMatrix> dcsdd = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
      *interface_->ActiveDofs(), 100, true, false, CORE::LINALG::SparseMatrix::FE_MATRIX));
  Teuchos::RCP<CORE::LINALG::SparseMatrix> dcsdLMc = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
      *interface_->ActiveDofs(), 100, true, false, CORE::LINALG::SparseMatrix::FE_MATRIX));
  Teuchos::RCP<Epetra_Vector> fcsa = CORE::LINALG::CreateVector(*interface_->ActiveDofs(), true);
  Teuchos::RCP<Epetra_Vector> g_all;
  if (constr_direction_ == INPAR::CONTACT::constr_xyz)
    g_all = CORE::LINALG::CreateVector(*interface_->SlaveRowDofs(), true);
  else
    g_all = CORE::LINALG::CreateVector(*interface_->SlaveRowNodes(), true);

  Teuchos::RCP<CORE::LINALG::SparseMatrix> dmatrix =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*interface_->SlaveRowDofs(), 10));
  Teuchos::RCP<CORE::LINALG::SparseMatrix> mmatrix =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*interface_->SlaveRowDofs(), 100));
  interface_->AssembleDM(*dmatrix, *mmatrix);
  dmatrix->Complete();
  mmatrix->Complete(*masterdofrowmap_, *slavedofrowmap_);

  // setup some linearizations
  CORE::LINALG::SparseMatrix linDcontactLM(
      *interface_->SlaveRowDofs(), 100, true, false, CORE::LINALG::SparseMatrix::FE_MATRIX);
  CORE::LINALG::SparseMatrix linMcontactLM(
      *interface_->MasterRowDofs(), 100, true, false, CORE::LINALG::SparseMatrix::FE_MATRIX);
  interface_->AssembleLinDM(linDcontactLM, linMcontactLM);

  // D and M matrix for the active nodes
  Teuchos::RCP<CORE::LINALG::SparseMatrix> dInv =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*interface_->SlaveRowDofs(), 100, true, false));

  // linearized normal contact
  interface_->AssembleS(*dcsdd);
  interface_->AssembleG(*g_all);

  if (contact_regularization_)
  {
    interface_->AssembleNormalContactRegularization(*dcsdd, *dcsdLMc, *fcsa);

    // linearized tangential contact (friction)
    if (interface_->IsFriction())
    {
      Teuchos::RCP<Epetra_Vector> rcsa_fr =
          CORE::LINALG::CreateVector(*interface_->ActiveDofs(), true);
      interface_->AssembleLinSlipNormalRegularization(*dcsdLMc, *dcsdd, *rcsa_fr);
      interface_->AssembleLinStick(*dcsdLMc, *dcsdd, *rcsa_fr);
      rcsa_fr->Scale(-1.);
      CONTACT::UTILS::AddVector(*rcsa_fr, *fcsa);
    }
    else
    {
      Teuchos::RCP<Epetra_Vector> rcsa_fr =
          CORE::LINALG::CreateVector(*interface_->ActiveDofs(), true);
      interface_->AssembleTN(dcsdLMc, Teuchos::null);
      interface_->AssembleTNderiv(dcsdd, Teuchos::null);
      interface_->AssembleTangrhs(*rcsa_fr);
      rcsa_fr->Scale(-1.);
      CONTACT::UTILS::AddVector(*rcsa_fr, *fcsa);
    }
  }
  else
    FOUR_C_THROW("stop");

  // complete all those linearizations
  //                             colmap        rowmap
  linDcontactLM.Complete(*SMdofMap(), *interface_->SlaveRowDofs());
  linMcontactLM.Complete(*SMdofMap(), *interface_->MasterRowDofs());

  // normal contact
  Teuchos::RCP<Epetra_Vector> gact;
  if (constr_direction_ == INPAR::CONTACT::constr_xyz)
  {
    gact = CORE::LINALG::CreateVector(*interface_->ActiveDofs(), true);
    if (gact->GlobalLength()) CORE::LINALG::Export(*g_all, *gact);
  }
  else
  {
    gact = CORE::LINALG::CreateVector(*interface_->ActiveNodes(), true);
    if (gact->GlobalLength())
    {
      CORE::LINALG::Export(*g_all, *gact);
      if (gact->ReplaceMap(*interface_->ActiveNDofs())) FOUR_C_THROW("replaceMap went wrong");
    }
  }
  CONTACT::UTILS::AddVector(*gact, *fcsa);
  fcsa->Norm2(&contact_rhs_norm_);

  // complete all the new matrix blocks
  // Note: since the contact interace assemled them, they are all based
  //       on displacement row and col maps. Hence, some still need to be transformed
  dcsdd->Complete(*SMdofMap(), *interface_->ActiveDofs());
  dcsdLMc->Complete(*interface_->ActiveDofs(), *interface_->ActiveDofs());

  // get the seperate blocks of the 2x2 TSI block system
  // View mode!!! Since we actually want to add things there
  Teuchos::RCP<CORE::LINALG::SparseMatrix> kss =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(sysmat->Matrix(0, 0), CORE::LINALG::Copy));
  Teuchos::RCP<CORE::LINALG::SparseMatrix> kst =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(sysmat->Matrix(0, 1), CORE::LINALG::Copy));
  Teuchos::RCP<CORE::LINALG::SparseMatrix> kts =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(sysmat->Matrix(1, 0), CORE::LINALG::Copy));
  Teuchos::RCP<CORE::LINALG::SparseMatrix> ktt =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(sysmat->Matrix(1, 1), CORE::LINALG::Copy));

  // get some maps
  Teuchos::RCP<Epetra_Map> gdisp_DofRowMap = Teuchos::rcp(new Epetra_Map(kss->RowMap()));
  Teuchos::RCP<Epetra_Map> gpres_DofRowMap = Teuchos::rcp(new Epetra_Map(ktt->RowMap()));
  Teuchos::RCP<Epetra_Map> gmdof = Teuchos::rcp(new Epetra_Map(*interface_->MasterRowDofs()));
  Teuchos::RCP<Epetra_Map> active_dofs = Teuchos::rcp(new Epetra_Map(*interface_->ActiveDofs()));

  // split rhs
  Teuchos::RCP<Epetra_Vector> rs = Teuchos::rcp(new Epetra_Vector(kss->RowMap(), true));
  Teuchos::RCP<Epetra_Vector> rt = Teuchos::rcp(new Epetra_Vector(ktt->RowMap(), true));
  CORE::LINALG::Export(*combined_RHS, *rs);
  CORE::LINALG::Export(*combined_RHS, *rt);

  // we don't want the rhs but the residual
  rs->Scale(-1.);
  rt->Scale(-1.);

  // add last time step contact forces to rhs
  if (fscn_ != Teuchos::null)  // in the first time step, we don't have any history of the
                               // contact force, after that, fscn_ should be initialized propperly
  {
    Epetra_Vector tmp(kss->RowMap());
    CORE::LINALG::Export(*fscn_, tmp);
    if (rs->Update(alphaf_, tmp, 1.) != 0)  // fscn already scaled with alphaf_ in update
      FOUR_C_THROW("update went wrong");
  }


  // map containing the inactive and non-contact structural dofs
  Teuchos::RCP<Epetra_Map> str_gni_dofs =
      CORE::LINALG::SplitMap(*CORE::LINALG::SplitMap(kss->RowMap(), *interface_->MasterRowDofs()),
          *interface_->ActiveDofs());

  // add to kss
  kss->UnComplete();
  kss->Add(linDcontactLM, false, 1. - alphaf_, 1.);
  kss->Add(linMcontactLM, false, 1. - alphaf_, 1.);

  // complete the matrix blocks again, now that we have added
  // the additional displacement linearizations
  kss->Complete();

  // now we have added the additional linearizations.
  // if there are no active nodes, we can leave now
  if (interface_->ActiveNodes()->NumGlobalElements() == 0)
  {
    sysmat->Reset();
    sysmat->Assign(0, 0, CORE::LINALG::Copy, *kss);
    sysmat->Assign(0, 1, CORE::LINALG::Copy, *kst);
    sysmat->Assign(1, 0, CORE::LINALG::Copy, *kts);
    sysmat->Assign(1, 1, CORE::LINALG::Copy, *ktt);
    return;
  }

  // split matrix blocks in 3 rows: Active, Master and (Inactive+others)
  Teuchos::RCP<CORE::LINALG::SparseMatrix> kss_ni, kss_m, kss_a, kst_ni, kst_m, kst_a, kts_ni,
      kts_m, kts_a, ktt_ni, ktt_m, ktt_a, dummy1, dummy2, dummy3;

  // temporary matrix
  Teuchos::RCP<CORE::LINALG::SparseMatrix> tmp;
  Teuchos::RCP<Epetra_Vector> tmpv;

  // an empty dummy map
  Teuchos::RCP<Epetra_Map> dummy_map1, dummy_map2;

  // ****************************************************
  // split kss block*************************************
  // ****************************************************
  // split first row
  CORE::LINALG::SplitMatrix2x2(
      kss, str_gni_dofs, dummy_map1, gdisp_DofRowMap, dummy_map2, kss_ni, dummy1, tmp, dummy2);

  // this shoud be a split in rows, so that two blocks should have zero columns
  if (dummy1->DomainMap().NumGlobalElements() != 0 || dummy2->DomainMap().NumGlobalElements() != 0)
    FOUR_C_THROW("this split should only split rows, no columns expected for this matrix blocks");

  // reset
  dummy1 = Teuchos::null;
  dummy2 = Teuchos::null;
  dummy_map1 = Teuchos::null;
  dummy_map2 = Teuchos::null;

  // split the remaining two rows
  CORE::LINALG::SplitMatrix2x2(
      tmp, gmdof, dummy_map1, gdisp_DofRowMap, dummy_map2, kss_m, dummy1, kss_a, dummy2);

  // this shoud be a split in rows, so that two blocks should have zero columns
  if (dummy1->DomainMap().NumGlobalElements() != 0 || dummy2->DomainMap().NumGlobalElements() != 0)
    FOUR_C_THROW("this split should only split rows, no columns expected for this matrix blocks");

  // reset
  dummy1 = Teuchos::null;
  dummy2 = Teuchos::null;
  dummy_map1 = Teuchos::null;
  dummy_map2 = Teuchos::null;
  tmp = Teuchos::null;
  // ****************************************************
  // split kss block*************************************
  // ****************************************************

  // ****************************************************
  // split kst block*************************************
  // ****************************************************
  // split first row
  CORE::LINALG::SplitMatrix2x2(
      kst, str_gni_dofs, dummy_map1, gpres_DofRowMap, dummy_map2, kst_ni, dummy1, tmp, dummy2);

  // this shoud be a split in rows, so that two blocks should have zero columns
  if (dummy1->DomainMap().NumGlobalElements() != 0 || dummy2->DomainMap().NumGlobalElements() != 0)
    FOUR_C_THROW("this split should only split rows, no columns expected for this matrix blocks");

  // reset
  dummy1 = Teuchos::null;
  dummy2 = Teuchos::null;
  dummy_map1 = Teuchos::null;
  dummy_map2 = Teuchos::null;

  // split the remaining two rows
  CORE::LINALG::SplitMatrix2x2(
      tmp, gmdof, dummy_map1, gpres_DofRowMap, dummy_map2, kst_m, dummy1, kst_a, dummy2);

  // this shoud be a split in rows, so that two blocks should have zero columns
  if (dummy1->DomainMap().NumGlobalElements() != 0 || dummy2->DomainMap().NumGlobalElements() != 0)
    FOUR_C_THROW("this split should only split rows, no columns expected for this matrix blocks");

  // reset
  dummy1 = Teuchos::null;
  dummy2 = Teuchos::null;
  dummy_map1 = Teuchos::null;
  dummy_map2 = Teuchos::null;
  tmp = Teuchos::null;
  // ****************************************************
  // split kst block*************************************
  // ****************************************************

  // ****************************************************
  // split rhs vectors***********************************
  // ****************************************************
  // split structural rhs
  Epetra_Vector rsni(*str_gni_dofs);
  CORE::LINALG::Export(*rs, rsni);
  Epetra_Vector rsm(*interface_->MasterRowDofs());
  CORE::LINALG::Export(*rs, rsm);
  Teuchos::RCP<Epetra_Vector> rsa = Teuchos::rcp(new Epetra_Vector(*interface_->ActiveDofs()));
  CORE::LINALG::Export(*rs, *rsa);
  // ****************************************************
  // split rhs vectors***********************************
  // ****************************************************

  // D and M matrix for the active nodes
  Teuchos::RCP<CORE::LINALG::SparseMatrix> dInvA =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*interface_->ActiveDofs(), 100, true, false));
  Teuchos::RCP<CORE::LINALG::SparseMatrix> mA =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*interface_->ActiveDofs(), 100, true, false));

  dummy_map1 = dummy_map2 = Teuchos::null;
  dummy1 = dummy2 = dummy3 = Teuchos::null;
  CORE::LINALG::SplitMatrix2x2(
      dmatrix, active_dofs, dummy_map1, active_dofs, dummy_map2, dInvA, dummy1, dummy2, dummy3);
  dInvA->Complete(*interface_->ActiveDofs(), *interface_->ActiveDofs());
  // invert D-matrix
  Epetra_Vector dDiag(*interface_->ActiveDofs());
  dInvA->ExtractDiagonalCopy(dDiag);
  if (dDiag.Reciprocal(dDiag)) FOUR_C_THROW("inversion of diagonal D matrix failed");
  dInvA->ReplaceDiagonalValues(dDiag);

  dummy_map1 = dummy_map2 = Teuchos::null;
  dummy1 = dummy2 = dummy3 = Teuchos::null;
  CORE::LINALG::SplitMatrix2x2(
      mmatrix, active_dofs, dummy_map1, gmdof, dummy_map2, mA, dummy1, dummy2, dummy3);
  mA->Complete(*interface_->MasterRowDofs(), *interface_->ActiveDofs());

  // get dinv * M
  Teuchos::RCP<CORE::LINALG::SparseMatrix> dInvMa =
      CORE::LINALG::MLMultiply(*dInvA, false, *mA, false, false, false, true);

  // we need to add another term, since AssembleLinStick/Slip assumes that we solve
  // for the Lagrange multiplier increments. However, we solve for the LM directly.
  // We can do that, since the system is linear in the LMs.
  tmpv = Teuchos::rcp(new Epetra_Vector(*interface_->ActiveDofs()));
  Teuchos::RCP<Epetra_Vector> tmpv2 = Teuchos::rcp(new Epetra_Vector(*interface_->ActiveDofs()));
  CORE::LINALG::Export(*z_, *tmpv2);
  dcsdLMc->Multiply(false, *tmpv2, *tmpv);
  tmpv->Scale(-1.);
  CONTACT::UTILS::AddVector(*tmpv, *fcsa);
  tmpv = Teuchos::null;
  tmpv2 = Teuchos::null;


  // save some matrix blocks for recovery
  dinvA_ = dInvA;
  kss_a_ = kss_a;
  kst_a_ = kst_a;
  rs_a_ = rsa;
  // apply contact symmetry conditions
  if (sdirichtoggle_.is_null()) FOUR_C_THROW("you didn't call StoreDirichletStatus");
  if (constr_direction_ == INPAR::CONTACT::constr_xyz)
  {
    double haveDBC = 0;
    sdirichtoggle_->Norm1(&haveDBC);
    if (haveDBC > 0.)
    {
      Teuchos::RCP<Epetra_Vector> diag =
          CORE::LINALG::CreateVector(*interface_->ActiveDofs(), true);
      dInvA->ExtractDiagonalCopy(*diag);
      Teuchos::RCP<Epetra_Vector> lmDBC =
          CORE::LINALG::CreateVector(*interface_->ActiveDofs(), true);
      CORE::LINALG::Export(*sdirichtoggle_, *lmDBC);
      Teuchos::RCP<Epetra_Vector> tmp = CORE::LINALG::CreateVector(*interface_->ActiveDofs(), true);
      tmp->Multiply(1., *diag, *lmDBC, 0.);
      diag->Update(-1., *tmp, 1.);
      dInvA->ReplaceDiagonalValues(*diag);
      dInvMa = CORE::LINALG::MLMultiply(*dInvA, false, *mA, false, false, false, true);
    }
  }

  // reset the tangent stiffness
  // (for the condensation we have constructed copies above)
  sysmat->UnComplete();

  // need diagonal block kss with explicitdirichtlet_=true
  // to be able to apply dirichlet values for contact symmetry condition
  CORE::LINALG::SparseMatrix tmpkss(
      *gdisp_DofRowMap, 100, false, false, CORE::LINALG::SparseMatrix::FE_MATRIX);
  sysmat->Assign(0, 0, CORE::LINALG::Copy, tmpkss);

  // get references to the blocks (just for convenience)
  CORE::LINALG::SparseMatrix& kss_new = sysmat->Matrix(0, 0);
  CORE::LINALG::SparseMatrix& kst_new = sysmat->Matrix(0, 1);
  kss_new.Reset();
  kst_new.Reset();
  // reynolds equation blocks remain untouched

  // reset rhs
  combined_RHS->PutScalar(0.);
  CONTACT::UTILS::AddVector(*rt, *combined_RHS);

  // **********************************************************************
  // **********************************************************************
  // BUILD CONDENSED SYSTEM
  // **********************************************************************
  // **********************************************************************

  // (1) add the blocks, we do nothing with (i.e. (Inactive+others))
  kss_new.Add(*kss_ni, false, 1., 1.);
  kst_new.Add(*kst_ni, false, 1., 1.);
  CONTACT::UTILS::AddVector(rsni, *combined_RHS);

  // (2) add the 'uncondensed' blocks (i.e. everything w/o a D^-1
  // (2)a actual stiffness blocks of the master-rows
  kss_new.Add(*kss_m, false, 1., 1.);
  kst_new.Add(*kst_m, false, 1., 1.);
  CONTACT::UTILS::AddVector(rsm, *combined_RHS);

  // (2)b active constraints in the active slave rows
  kss_new.Add(*dcsdd, false, 1., 1.);
  CONTACT::UTILS::AddVector(*fcsa, *combined_RHS);

  // (3) condensed parts
  // second row
  kss_new.Add(
      *CORE::LINALG::MLMultiply(*dInvMa, true, *kss_a, false, false, false, true), false, 1., 1.);
  kst_new.Add(
      *CORE::LINALG::MLMultiply(*dInvMa, true, *kst_a, false, false, false, true), false, 1., 1.);
  tmpv = Teuchos::rcp(new Epetra_Vector(*interface_->MasterRowDofs()));
  if (dInvMa->Multiply(true, *rsa, *tmpv)) FOUR_C_THROW("multiply failed");
  CONTACT::UTILS::AddVector(*tmpv, *combined_RHS);
  tmpv = Teuchos::null;

  // third row
  Teuchos::RCP<CORE::LINALG::SparseMatrix> wDinv =
      CORE::LINALG::MLMultiply(*dcsdLMc, false, *dInvA, true, false, false, true);
  kss_new.Add(*CORE::LINALG::MLMultiply(*wDinv, false, *kss_a, false, false, false, true), false,
      -1. / (1. - alphaf_), 1.);
  kst_new.Add(*CORE::LINALG::MLMultiply(*wDinv, false, *kst_a, false, false, false, true), false,
      -1. / (1. - alphaf_), 1.);
  tmpv = Teuchos::rcp(new Epetra_Vector(*interface_->ActiveDofs()));
  wDinv->Multiply(false, *rsa, *tmpv);
  tmpv->Scale(-1. / (1. - alphaf_));
  CONTACT::UTILS::AddVector(*tmpv, *combined_RHS);
  tmpv = Teuchos::null;
  wDinv = Teuchos::null;

  // and were done with the system matrix
  sysmat->Complete();

  // we need to return the rhs, not the residual
  combined_RHS->Scale(-1.);

  return;
}


void ADAPTER::CouplingEhlMortar::EvaluateRelMov()
{
  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    DRT::Node* node = Interface()->Discret().lRowNode(i);
    if (!node) FOUR_C_THROW("node not found");
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);
    if (!cnode) FOUR_C_THROW("not a contact node");

    cnode->FriData().GetDerivJump().resize(3);
    // write it to nodes
    for (int dim = 0; dim < interface_->Dim(); dim++)
    {
      cnode->FriData().jump()[dim] = cnode->EhlData().GetWeightedRelTangVel()(dim);
      for (auto p = cnode->EhlData().GetWeightedRelTangVelDeriv().begin();
           p != cnode->EhlData().GetWeightedRelTangVelDeriv().end(); ++p)
        cnode->FriData().GetDerivJump()[dim][p->first] = p->second(dim);
    }
  }
}

void ADAPTER::CouplingEhlMortar::RecoverCoupled(
    Teuchos::RCP<Epetra_Vector> sinc, Teuchos::RCP<Epetra_Vector> tinc)
{
  const double alphaf_ = 0.;  // statics!

  Teuchos::RCP<Epetra_Vector> z_old = Teuchos::null;
  if (z_ != Teuchos::null) z_old = Teuchos::rcp(new Epetra_Vector(*z_));

  // recover contact LM
  if (interface_->ActiveNodes()->NumGlobalElements() > 0)
  {
    // do we have everything we need?
    if (rs_a_ == Teuchos::null || kss_a_ == Teuchos::null || kst_a_ == Teuchos::null ||
        dinvA_ == Teuchos::null)
      FOUR_C_THROW("some data for LM recovery is missing");

    Epetra_Vector lmc_a_new(*interface_->ActiveDofs(), false);
    Epetra_Vector tmp(*interface_->ActiveDofs(), false);
    lmc_a_new.Update(1., *rs_a_, 0.);
    kss_a_->Multiply(false, *sinc, tmp);
    lmc_a_new.Update(1., tmp, 1.);
    kst_a_->Multiply(false, *tinc, tmp);
    lmc_a_new.Update(1., tmp, 1.);
    dinvA_->Multiply(false, lmc_a_new, tmp);
    tmp.Scale(-1. / (1. - alphaf_));
    z_ = Teuchos::rcp(new Epetra_Vector(*interface_->SlaveRowDofs()));

    CORE::LINALG::Export(tmp, *z_);
  }

  else
    z_ = Teuchos::rcp(new Epetra_Vector(*interface_->SlaveRowDofs()));

  if (z_old != Teuchos::null)
  {
    z_old->Update(-1., *z_, 1.);
    z_old->Norm2(&contact_LM_incr_norm_);
  }

  // store updated LM into nodes
  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(interface_->Discret().lRowNode(i));
    for (int dof = 0; dof < interface_->Dim(); ++dof)
      cnode->MoData().lm()[dof] = z_->operator[](z_->Map().LID(cnode->Dofs()[dof]));
  }

  return;
}



/*----------------------------------------------------------------------*
 |  Store dirichlet B.C. status into CNode                    popp 06/09|
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingEhlMortar::StoreDirichletStatus(
    Teuchos::RCP<const CORE::LINALG::MapExtractor> dbcmaps)
{
  // loop over all slave row nodes on the current interface
  for (int j = 0; j < interface_->SlaveRowNodes()->NumMyElements(); ++j)
  {
    int gid = interface_->SlaveRowNodes()->GID(j);
    DRT::Node* node = interface_->Discret().gNode(gid);
    if (!node) FOUR_C_THROW("ERROR: Cannot find node with gid %", gid);
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);

    // check if this node's dofs are in dbcmap
    for (int k = 0; k < cnode->NumDof(); ++k)
    {
      int currdof = cnode->Dofs()[k];
      int lid = (dbcmaps->CondMap())->LID(currdof);

      // store dbc status if found
      if (lid >= 0 && cnode->DbcDofs()[k] == false) cnode->SetDbc() = true;

      // check compatibility of contact symmetry condition and displacement dirichlet conditions
      if (lid < 0 && cnode->DbcDofs()[k] == true)
      {
        std::cout << "node " << cnode->Id() << " at: " << cnode->X()[0] << " " << cnode->X()[1]
                  << " " << cnode->X()[2] << std::endl;
        std::cout << "dbcdofs: " << cnode->DbcDofs()[0] << cnode->DbcDofs()[1]
                  << cnode->DbcDofs()[2] << std::endl;
        FOUR_C_THROW(
            "Inconsistency in structure Dirichlet conditions and Mortar symmetry conditions");
      }
    }
  }
  // create old style dirichtoggle vector (supposed to go away)
  sdirichtoggle_ = Teuchos::rcp(new Epetra_Vector(*interface_->SlaveRowDofs(), true));
  Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*(dbcmaps->CondMap())));
  temp->PutScalar(1.0);
  CORE::LINALG::Export(*temp, *sdirichtoggle_);

  return;
}



bool ADAPTER::CouplingEhlMortar::AlreadyEvaluated(Teuchos::RCP<const Epetra_Vector> disp)
{
  if (evaluated_state_.is_null()) return false;
  Teuchos::RCP<Epetra_Vector> diff = Teuchos::rcp(new Epetra_Vector(*disp));
  if (diff->Update(-1., *evaluated_state_, 1.)) FOUR_C_THROW("update failed");
  double inf_diff = -1.;
  if (diff->NormInf(&inf_diff)) FOUR_C_THROW("NormInf failed");
  if (inf_diff < 1.e-13) return true;

  return false;
}

Teuchos::RCP<CORE::LINALG::SparseMatrix> ADAPTER::CouplingEhlMortar::AssembleEHLLinD(
    const Teuchos::RCP<Epetra_Vector> x  // slave dof vector
)
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> DLinEHL = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
      *slavedofrowmap_, 81, true, false, CORE::LINALG::SparseMatrix::FE_MATRIX));
  DLinEHL->Zero();
  DLinEHL->UnComplete();

  interface_->AssembleCoupLinD(*DLinEHL, x);

  DLinEHL->Complete(*smdofrowmap_, *slavedofrowmap_);

  return DLinEHL;
}

Teuchos::RCP<CORE::LINALG::SparseMatrix> ADAPTER::CouplingEhlMortar::AssembleEHLLinM(
    const Teuchos::RCP<Epetra_Vector> x  // slave dof vector
)
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> MLinEHL = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
      *masterdofrowmap_, 81, true, false, CORE::LINALG::SparseMatrix::FE_MATRIX));
  MLinEHL->Zero();
  MLinEHL->UnComplete();

  interface_->AssembleCoupLinM(*MLinEHL, x);

  MLinEHL->Complete(*smdofrowmap_, *masterdofrowmap_);

  return MLinEHL;
}

void ADAPTER::CouplingEhlMortar::AssembleNormals()
{
  normals_ = Teuchos::rcp(new Epetra_Vector(*SlaveDofMap(), true));

  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) FOUR_C_THROW("node not found");
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);
    if (!cnode) FOUR_C_THROW("not a contact node");

    for (int d = 0; d < interface_->Dim(); ++d)
      normals_->ReplaceGlobalValue(cnode->Dofs()[d], 0, cnode->MoData().n()[d]);
  }
}


void ADAPTER::CouplingEhlMortar::AssembleNormalsDeriv()
{
  Nderiv_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(*slavedofrowmap_, 81, false, false));
  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) FOUR_C_THROW("node not found");
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);
    if (!cnode) FOUR_C_THROW("not a contact node");

    for (int d = 0; d < Interface()->Dim(); ++d)
      for (auto p = cnode->Data().GetDerivN()[d].begin(); p != cnode->Data().GetDerivN()[d].end();
           ++p)
        Nderiv_->Assemble(p->second, cnode->Dofs()[d], p->first);
  }
  Nderiv_->Complete();
}

void ADAPTER::CouplingEhlMortar::AssembleRealGap()
{
  nodal_gap_ = Teuchos::rcp(new Epetra_Vector(*slavenoderowmap_, true));

  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) FOUR_C_THROW("node not found");
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);
    if (!cnode) FOUR_C_THROW("not a contact node");
    double real_gap = cnode->Data().Getg();
    switch (cnode->MoData().GetD().size())
    {
      case 0:
        break;
      case 1:
        if (cnode->MoData().GetD().begin()->first != cnode->Id())
          FOUR_C_THROW("something is wrong. Here should by my own Id");
        real_gap /= cnode->MoData().GetD().at(cnode->Id());
        break;
      default:
        FOUR_C_THROW(
            "GetD should be of size 0 (unprojectable) or 1 (projectable). Are you not using "
            "duals?");
    }
    nodal_gap_->ReplaceGlobalValue(cnode->Id(), 0, real_gap);
  }

  static const double offset =
      GLOBAL::Problem::Instance()->LubricationDynamicParams().get<double>("GAP_OFFSET");
  for (int i = 0; i < nodal_gap_->Map().NumMyElements(); ++i) nodal_gap_->operator[](i) += offset;
}

void ADAPTER::CouplingEhlMortar::AssembleRealGapDeriv()
{
  deriv_nodal_gap_ =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*slavedofrowmap_, 81, false, false));

  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) FOUR_C_THROW("node not found");
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);
    if (!cnode) FOUR_C_THROW("not a contact node");

    if (cnode->Data().GetDerivD().size() != cnode->MoData().GetD().size())
      FOUR_C_THROW("size inconsistency");

    const double w_gap = cnode->Data().Getg();
    double d = -1.;
    switch (cnode->Data().GetDerivD().size())
    {
      case 0:
        break;
      case 1:
        if (cnode->Data().GetDerivD().begin()->first != cnode->Id())
          FOUR_C_THROW("something is wrong. Here should by my own Id");
        d = cnode->MoData().GetD().at(cnode->Id());
        break;
      default:
        FOUR_C_THROW(
            "GetDerivD should be of size 0 (unprojectable) or 1 (projectable). Are you not using "
            "duals?");
    }

    if (cnode->Data().GetDerivD().size())
      for (auto p = cnode->Data().GetDerivD().at(cnode->Id()).begin();
           p != cnode->Data().GetDerivD().at(cnode->Id()).end(); ++p)
      {
        const double val = -w_gap / (d * d) * p->second;
        for (int d = 0; d < interface_->Dim(); ++d)
          deriv_nodal_gap_->Assemble(val, cnode->Dofs()[d], p->first);
      }

    if (d == -1 && cnode->Data().GetDerivG().size() != 0) FOUR_C_THROW("inconsistency");

    if (cnode->Data().GetDerivG().size())
      for (auto p = cnode->Data().GetDerivG().begin(); p != cnode->Data().GetDerivG().end(); ++p)
      {
        const double val = p->second / d;
        for (int d = 0; d < interface_->Dim(); ++d)
          deriv_nodal_gap_->Assemble(val, cnode->Dofs()[d], p->first);
      }
  }
  deriv_nodal_gap_->Complete(*smdofrowmap_, *slavedofrowmap_);
}

void ADAPTER::CouplingEhlMortar::AssembleInterfaceVelocities(const double dt)
{
  relTangVel_ = Teuchos::rcp(new Epetra_Vector(*slavedofrowmap_));
  avTangVel_ = Teuchos::rcp(new Epetra_Vector(*slavedofrowmap_));
  relTangVel_deriv_ =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*slavedofrowmap_, 81, false, false));
  avTangVel_deriv_ =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*slavedofrowmap_, 81, false, false));

  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) FOUR_C_THROW("node not found");
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);
    if (!cnode) FOUR_C_THROW("not a contact node");


    double d_val = 0.;
    switch (cnode->MoData().GetD().size())
    {
      case 0:
        break;
      case 1:
        if (cnode->MoData().GetD().begin()->first != cnode->Id())
          FOUR_C_THROW("something is wrong. Here should by my own Id");
        d_val = cnode->MoData().GetD().at(cnode->Id());
        break;
      default:
        FOUR_C_THROW(
            "GetD should be of size 0 (unprojectable) or 1 (projectable). Are you not using "
            "duals?");
    }

    if (d_val == 0.) continue;

    for (int d = 0; d < Interface()->Dim(); ++d)
    {
      relTangVel_->ReplaceGlobalValue(
          cnode->Dofs()[d], 0, cnode->EhlData().GetWeightedRelTangVel()(d) / d_val);
      avTangVel_->ReplaceGlobalValue(
          cnode->Dofs()[d], 0, cnode->EhlData().GetWeightedAvTangVel()(d) / d_val);
    }

    for (auto p = cnode->Data().GetDerivD().at(cnode->Id()).begin();
         p != cnode->Data().GetDerivD().at(cnode->Id()).end(); ++p)
    {
      const int col = p->first;
      for (int d = 0; d < Interface()->Dim(); ++d)
      {
        const int row = cnode->Dofs()[d];
        const double rel_val =
            -cnode->EhlData().GetWeightedRelTangVel()(d) / (d_val * d_val) * p->second;
        const double av_val =
            -cnode->EhlData().GetWeightedAvTangVel()(d) / (d_val * d_val) * p->second;
        relTangVel_deriv_->Assemble(rel_val, row, col);
        avTangVel_deriv_->Assemble(av_val, row, col);
      }
    }
    for (auto p = cnode->EhlData().GetWeightedAvTangVelDeriv().begin();
         p != cnode->EhlData().GetWeightedAvTangVelDeriv().end(); ++p)
    {
      const int col = p->first;
      for (int d = 0; d < Interface()->Dim(); ++d)
      {
        const int row = cnode->Dofs()[d];
        const double val = p->second(d) / d_val;
        avTangVel_deriv_->Assemble(val, row, col);
      }
    }
    for (auto p = cnode->EhlData().GetWeightedRelTangVelDeriv().begin();
         p != cnode->EhlData().GetWeightedRelTangVelDeriv().end(); ++p)
    {
      const int col = p->first;
      for (int d = 0; d < Interface()->Dim(); ++d)
      {
        const int row = cnode->Dofs()[d];
        const double val = p->second(d) / d_val;
        relTangVel_deriv_->Assemble(val, row, col);
      }
    }
  }

  relTangVel_->Scale(1. / dt);
  avTangVel_->Scale(1. / dt);
  relTangVel_deriv_->Complete(*smdofrowmap_, *slavedofrowmap_);
  avTangVel_deriv_->Complete(*smdofrowmap_, *slavedofrowmap_);
  relTangVel_deriv_->Scale(1. / dt);
  avTangVel_deriv_->Scale(1. / dt);
}

void ADAPTER::CouplingEhlMortar::AssembleSurfGrad()
{
  SurfGrad_ = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
      *slavedofrowmap_, 81, false, false, CORE::LINALG::SparseMatrix::FE_MATRIX));

  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) FOUR_C_THROW("ERROR: Cannot find node");
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);
    if (!cnode) FOUR_C_THROW("this is not a contact node");

    double dval = 1.;
    switch (cnode->MoData().GetD().size())
    {
      case 0:
        dval = 1.e32;
        break;  // large number so no tangential gradient
      case 1:
        if (cnode->MoData().GetD().begin()->first != cnode->Id())
          FOUR_C_THROW("something is wrong. Here should by my own Id");
        dval = cnode->MoData().GetD().at(cnode->Id());
        break;
      default:
        FOUR_C_THROW(
            "GetD should be of size 0 (unprojectable) or 1 (projectable). Are you not using "
            "duals?");
    }

    for (auto p = cnode->EhlData().GetSurfGrad().begin(); p != cnode->EhlData().GetSurfGrad().end();
         ++p)
      for (int d = 0; d < Interface()->Dim(); ++d)
        SurfGrad_->Assemble(p->second(d) / dval, cnode->Dofs()[d], p->first);
  }

  SurfGrad_->Complete();
}

Teuchos::RCP<CORE::LINALG::SparseMatrix> ADAPTER::CouplingEhlMortar::AssembleSurfGradDeriv(
    const Teuchos::RCP<const Epetra_Vector> x)
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> SurfGradDeriv =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(
          *slavedofrowmap_, 81, false, false, CORE::LINALG::SparseMatrix::FE_MATRIX));

  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) FOUR_C_THROW("ERROR: Cannot find node");
    CONTACT::Node* cnode = dynamic_cast<CONTACT::Node*>(node);
    if (!cnode) FOUR_C_THROW("this is not a contact node");

    double dval = 1.;
    switch (cnode->MoData().GetD().size())
    {
      case 0:
        dval = 1.e32;
        break;  // large number so no tangential gradient
      case 1:
        if (cnode->MoData().GetD().begin()->first != cnode->Id())
          FOUR_C_THROW("something is wrong. Here should by my own Id");
        dval = cnode->MoData().GetD().at(cnode->Id());
        break;
      default:
        FOUR_C_THROW(
            "GetD should be of size 0 (unprojectable) or 1 (projectable). Are you not using "
            "duals?");
    }

    for (auto p = cnode->EhlData().GetSurfGradDeriv().begin();
         p != cnode->EhlData().GetSurfGradDeriv().end(); ++p)
    {
      const int col = p->first;
      for (auto q = p->second.begin(); q != p->second.end(); ++q)
      {
        const int lid = x->Map().LID(q->first);
        if (lid < 0) FOUR_C_THROW("not my gid");
        const double x_val = x->operator[](lid);
        for (int d = 0; d < Interface()->Dim(); ++d)
        {
          const double val = x_val * q->second(d) / dval;
          SurfGradDeriv->Assemble(val, cnode->Dofs()[d], col);
        }
      }
    }

    if (cnode->Data().GetDerivD().size())
      for (auto p = cnode->Data().GetDerivD().at(cnode->Id()).begin();
           p != cnode->Data().GetDerivD().at(cnode->Id()).end(); ++p)
      {
        const int col = p->first;

        for (auto q = cnode->EhlData().GetSurfGrad().begin();
             q != cnode->EhlData().GetSurfGrad().end(); ++q)
          for (int d = 0; d < Interface()->Dim(); ++d)
          {
            const int row = cnode->Dofs()[d];
            const int x_gid = q->first;
            const int x_lid = x->Map().LID(x_gid);
            if (x_lid < 0) FOUR_C_THROW("not my gid");
            double x_val = x->operator[](x_lid);
            const double val = -x_val * q->second(d) / (dval * dval) * p->second;
            SurfGradDeriv->Assemble(val, row, col);
          }
      }
  }
  SurfGradDeriv->Complete(*smdofrowmap_, *slavedofrowmap_);
  return SurfGradDeriv;
}


void ADAPTER::CouplingEhlMortar::CreateForceVec(
    Teuchos::RCP<Epetra_Vector>& n, Teuchos::RCP<Epetra_Vector>& t)
{
  n = Teuchos::rcp(new Epetra_Vector(*interface_->SlaveRowDofs()));
  t = Teuchos::rcp(new Epetra_Vector(*interface_->SlaveRowDofs()));
  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(interface_->Discret().lRowNode(i));
    if (!cnode) FOUR_C_THROW("cast failed");
    const CORE::LINALG::Matrix<3, 1> lm(cnode->MoData().lm(), true);
    const CORE::LINALG::Matrix<3, 1> nor(cnode->MoData().n(), true);
    CORE::LINALG::Matrix<3, 3> nn;
    nn.MultiplyNT(nor, nor);
    CORE::LINALG::Matrix<3, 1> lmn;
    lmn.Multiply(nn, lm);
    CORE::LINALG::Matrix<3, 1> lmt(lm);
    lmt.Update(-1., lmn, 1.);
    for (int d = 0; d < 3; ++d)
    {
      n->operator[](n->Map().LID(cnode->Dofs()[d])) = lmn(d);
      t->operator[](t->Map().LID(cnode->Dofs()[d])) = lmt(d);
    }
  }
}

void ADAPTER::CouplingEhlMortar::CreateActiveSlipToggle(Teuchos::RCP<Epetra_Vector>* active,
    Teuchos::RCP<Epetra_Vector>* slip, Teuchos::RCP<Epetra_Vector>* active_old)
{
  *active = Teuchos::rcp(new Epetra_Vector(*interface_->SlaveRowNodes()));
  *slip = Teuchos::rcp(new Epetra_Vector(*interface_->SlaveRowNodes()));
  if (active_old != nullptr)
    *active_old = Teuchos::rcp(new Epetra_Vector(*interface_->SlaveRowNodes()));
  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(interface_->Discret().lRowNode(i));
    if (!cnode) FOUR_C_THROW("cast failed");
    if (cnode->Active())
      (*active)->operator[](i) = 1.;
    else
      (*active)->operator[](i) = 0.;
    if (cnode->FriData().Slip())
      (*slip)->operator[](i) = 1.;
    else
      (*slip)->operator[](i) = 0.;

    if (active_old != nullptr)
    {
      if (cnode->Data().ActiveOld())
        (*active_old)->operator[](i) = 1.;
      else
        (*active_old)->operator[](i) = 0.;
    }
  }
}

void ADAPTER::CouplingEhlMortar::WriteRestart(IO::DiscretizationWriter& output)
{
  if (!contact_regularization_) return;

  output.WriteVector("last_contact_force", fscn_);
  output.WriteVector("contact_lm", z_);

  Teuchos::RCP<Epetra_Vector> active_toggle, active_old_toggle, slip_toggle;
  CreateActiveSlipToggle(&active_toggle, &slip_toggle, &active_old_toggle);

  output.WriteVector("active_toggle", active_toggle);
  output.WriteVector("active_old_toggle", active_old_toggle);
  output.WriteVector("slip_toggle", slip_toggle);
}

void ADAPTER::CouplingEhlMortar::ReadRestart(IO::DiscretizationReader& reader)
{
  if (!contact_regularization_) return;

  reader.ReadVector(fscn_, "last_contact_force");
  reader.ReadVector(z_, "contact_lm");

  Teuchos::RCP<Epetra_Vector> active_toggle =
      Teuchos::rcp(new Epetra_Vector(*interface_->SlaveRowNodes()));
  Teuchos::RCP<Epetra_Vector> active_old_toggle =
      Teuchos::rcp(new Epetra_Vector(*interface_->SlaveRowNodes()));
  Teuchos::RCP<Epetra_Vector> slip_toggle =
      Teuchos::rcp(new Epetra_Vector(*interface_->SlaveRowNodes()));
  reader.ReadVector(active_toggle, "active_toggle");
  reader.ReadVector(active_old_toggle, "active_old_toggle");
  reader.ReadVector(slip_toggle, "slip_toggle");

  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(interface_->Discret().lRowNode(i));
    if (!cnode) FOUR_C_THROW("cast failed");
    cnode->Active() = active_toggle->operator[](i);
    cnode->FriData().Slip() = slip_toggle->operator[](i);
    cnode->Data().ActiveOld() = active_old_toggle->operator[](i);
    for (int d = 0; d < interface_->Dim(); ++d)
      cnode->MoData().lm()[d] = z_->operator[](z_->Map().LID(cnode->Dofs()[d]));
  }
}

int ADAPTER::CouplingEhlMortar::ActiveContact()
{
  return interface_->ActiveNodes()->NumGlobalElements();
}

int ADAPTER::CouplingEhlMortar::SlipContact()
{
  return interface_->SlipNodes()->NumGlobalElements();
}

FOUR_C_NAMESPACE_CLOSE
