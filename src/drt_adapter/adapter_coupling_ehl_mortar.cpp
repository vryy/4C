/*!----------------------------------------------------------------------
\file adapter_coupling_ehl_mortar.cpp
\brief mortar coupling terms of ehl

\level 3
\maintainer Alexander Seitz

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  headers                                                             |
 *----------------------------------------------------------------------*/
#include "adapter_coupling_ehl_mortar.H"

#include "../drt_contact/contact_interface.H"
#include "../drt_contact/friction_node.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_contact/contact_tsi_lagrange_strategy.H"

#include "../linalg/linalg_multiply.H"
#include "../drt_io/io.H"

ADAPTER::CouplingEhlMortar::CouplingEhlMortar()
    : CouplingNonLinMortar(),
      contact_regularization_(DRT::INPUT::IntegralValue<int>(
          DRT::Problem::Instance()->ContactDynamicParams(), "REGULARIZED_NORMAL_CONTACT")),
      regularization_thickness_(
          DRT::Problem::Instance()->ContactDynamicParams().get<double>("REGULARIZATION_THICKNESS")),
      regularization_compliance_(
          DRT::Problem::Instance()->ContactDynamicParams().get<double>("REGULARIZATION_STIFFNESS"))
{
  if (DRT::INPUT::IntegralValue<INPAR::MORTAR::ParRedist>(
          DRT::Problem::Instance()->MortarCouplingParams(), "PARALLEL_REDIST") !=
      INPAR::MORTAR::parredist_none)
    dserror(
        "EHL does not support parallel redistribution. Set \"PARALLEL_REDIST none\" in section "
        "\"MORTAR COUPLING\"");

  if (contact_regularization_)
    if (regularization_compliance_ <= 0. || regularization_thickness_ <= 0.)
      dserror("need positive REGULARIZATION_THICKNESS and REGULARIZATION_STIFFNESS");
  if (contact_regularization_) regularization_compliance_ = 1. / regularization_compliance_;
  if (DRT::INPUT::IntegralValue<int>(
          DRT::Problem::Instance()->ContactDynamicParams(), "REGULARIZED_NORMAL_CONTACT") == true &&
      DRT::INPUT::IntegralValue<bool>(
          DRT::Problem::Instance()->ElastoHydroDynamicParams(), "DRY_CONTACT_MODEL") == false)
    dserror("for dry contact model you need REGULARIZED_NORMAL_CONTACT and DRY_CONTACT_MODEL");
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

  INPAR::CONTACT::FrictionType ftype = DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(
      DRT::Problem::Instance()->ContactDynamicParams(), "FRICTION");

  std::vector<DRT::Condition*> ehl_conditions(0);
  masterdis->GetCondition(couplingcond, ehl_conditions);
  std::vector<const std::string*> sides((int)ehl_conditions.size());
  double fr_coeff = -1.;
  for (int i = 0; i < (int)ehl_conditions.size(); ++i)
  {
    const std::vector<int>* group1v = ehl_conditions[i]->Get<std::vector<int>>("Interface ID");
    if (!group1v) dserror("no interface Id found");
    if (group1v->operator[](0) != 1) dserror("only one interface id expected");
    const double fr = ehl_conditions[i]->GetDouble("FrCoeffOrBound");
    if (fr != ehl_conditions[0]->GetDouble("FrCoeffOrBound"))
      dserror("inconsistency in friction coefficients");
    fr_coeff = fr;
  }

  switch (ftype)
  {
    case INPAR::CONTACT::friction_tresca:
      dserror("no tresca friction supported");
      break;
    case INPAR::CONTACT::friction_none:
      break;
    case INPAR::CONTACT::friction_coulomb:
      interface_->IParams().set<double>("FRCOEFF", fr_coeff);
      interface_->IParams().set<double>("FRBOUND", -1.);
      break;
    default:
      dserror("don't know what to do with this friction type");
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
  D_ = Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap_, 81, false, false));
  M_ = Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap_, 81, false, false));
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
void ADAPTER::CouplingEhlMortar::CondenseContact(Teuchos::RCP<LINALG::BlockSparseMatrixBase> sysmat,
    Teuchos::RCP<Epetra_Vector>& combined_RHS, Teuchos::RCP<const Epetra_Vector> disp,
    const double dt)
{
  const double alphaf_ = 0.;  // statics!
  const INPAR::CONTACT::ConstraintDirection& constr_direction_ =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::ConstraintDirection>(
          Interface()->IParams(), "CONSTRAINT_DIRECTIONS");

  // return if this state has already been evaluated
  if (not AlreadyEvaluated(disp)) Integrate(disp, dt);

  // get the relative movement for frictional contact
  EvaluateRelMov();

  // update active set
  as_converged_ = interface_->UpdateActiveSetSemiSmooth();
  interface_->BuildActiveSet();

  // assemble the constraint lines for the active contact nodes
  Teuchos::RCP<LINALG::SparseMatrix> dcsdd = Teuchos::rcp(new LINALG::SparseMatrix(
      *interface_->ActiveDofs(), 100, true, false, LINALG::SparseMatrix::FE_MATRIX));
  Teuchos::RCP<LINALG::SparseMatrix> dcsdLMc = Teuchos::rcp(new LINALG::SparseMatrix(
      *interface_->ActiveDofs(), 100, true, false, LINALG::SparseMatrix::FE_MATRIX));
  Teuchos::RCP<Epetra_Vector> fcsa = LINALG::CreateVector(*interface_->ActiveDofs(), true);
  Teuchos::RCP<Epetra_Vector> g_all;
  if (constr_direction_ == INPAR::CONTACT::constr_xyz)
    g_all = LINALG::CreateVector(*interface_->SlaveRowDofs(), true);
  else
    g_all = LINALG::CreateVector(*interface_->SlaveRowNodes(), true);

  Teuchos::RCP<LINALG::SparseMatrix> dmatrix =
      Teuchos::rcp(new LINALG::SparseMatrix(*interface_->SlaveRowDofs(), 10));
  Teuchos::RCP<LINALG::SparseMatrix> mmatrix =
      Teuchos::rcp(new LINALG::SparseMatrix(*interface_->SlaveRowDofs(), 100));
  interface_->AssembleDM(*dmatrix, *mmatrix);
  dmatrix->Complete();
  mmatrix->Complete(*masterdofrowmap_, *slavedofrowmap_);

  // setup some linearizations
  LINALG::SparseMatrix linDcontactLM(
      *interface_->SlaveRowDofs(), 100, true, false, LINALG::SparseMatrix::FE_MATRIX);
  LINALG::SparseMatrix linMcontactLM(
      *interface_->MasterRowDofs(), 100, true, false, LINALG::SparseMatrix::FE_MATRIX);
  interface_->AssembleLinDM(linDcontactLM, linMcontactLM);

  // D and M matrix for the active nodes
  Teuchos::RCP<LINALG::SparseMatrix> dInv =
      Teuchos::rcp(new LINALG::SparseMatrix(*interface_->SlaveRowDofs(), 100, true, false));

  // linearized normal contact
  interface_->AssembleS(*dcsdd);
  interface_->AssembleG(*g_all);

  if (contact_regularization_)
  {
    interface_->AssembleNormalContactRegularization(*dcsdd, *dcsdLMc, *fcsa);

    // linearized tangential contact (friction)
    if (interface_->IsFriction())
    {
      Teuchos::RCP<Epetra_Vector> rcsa_fr = LINALG::CreateVector(*interface_->ActiveDofs(), true);
      interface_->AssembleLinSlipNormalRegularization(*dcsdLMc, *dcsdd, *rcsa_fr);
      interface_->AssembleLinStick(*dcsdLMc, *dcsdd, *rcsa_fr);
      rcsa_fr->Scale(-1.);
      CONTACT::UTILS::AddVector(*rcsa_fr, *fcsa);
    }
    else
    {
      Teuchos::RCP<Epetra_Vector> rcsa_fr = LINALG::CreateVector(*interface_->ActiveDofs(), true);
      interface_->AssembleTN(dcsdLMc, Teuchos::null);
      interface_->AssembleTNderiv(dcsdd, Teuchos::null);
      interface_->AssembleTangrhs(*rcsa_fr);
      rcsa_fr->Scale(-1.);
      CONTACT::UTILS::AddVector(*rcsa_fr, *fcsa);
    }
  }
  else
    dserror("stop");

  // complete all those linearizations
  //                             colmap        rowmap
  linDcontactLM.Complete(*SMdofMap(), *interface_->SlaveRowDofs());
  linMcontactLM.Complete(*SMdofMap(), *interface_->MasterRowDofs());

  // normal contact
  Teuchos::RCP<Epetra_Vector> gact;
  if (constr_direction_ == INPAR::CONTACT::constr_xyz)
  {
    gact = LINALG::CreateVector(*interface_->ActiveDofs(), true);
    if (gact->GlobalLength()) LINALG::Export(*g_all, *gact);
  }
  else
  {
    gact = LINALG::CreateVector(*interface_->ActiveNodes(), true);
    if (gact->GlobalLength())
    {
      LINALG::Export(*g_all, *gact);
      if (gact->ReplaceMap(*interface_->ActiveNDofs())) dserror("replaceMap went wrong");
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
  Teuchos::RCP<LINALG::SparseMatrix> kss =
      Teuchos::rcp(new LINALG::SparseMatrix(sysmat->Matrix(0, 0), LINALG::Copy));
  Teuchos::RCP<LINALG::SparseMatrix> kst =
      Teuchos::rcp(new LINALG::SparseMatrix(sysmat->Matrix(0, 1), LINALG::Copy));
  Teuchos::RCP<LINALG::SparseMatrix> kts =
      Teuchos::rcp(new LINALG::SparseMatrix(sysmat->Matrix(1, 0), LINALG::Copy));
  Teuchos::RCP<LINALG::SparseMatrix> ktt =
      Teuchos::rcp(new LINALG::SparseMatrix(sysmat->Matrix(1, 1), LINALG::Copy));

  // get some maps
  Teuchos::RCP<Epetra_Map> gdisp_DofRowMap = Teuchos::rcp(new Epetra_Map(kss->RowMap()));
  Teuchos::RCP<Epetra_Map> gpres_DofRowMap = Teuchos::rcp(new Epetra_Map(ktt->RowMap()));
  Teuchos::RCP<Epetra_Map> gmdof = Teuchos::rcp(new Epetra_Map(*interface_->MasterRowDofs()));
  Teuchos::RCP<Epetra_Map> active_dofs = Teuchos::rcp(new Epetra_Map(*interface_->ActiveDofs()));

  // split rhs
  Teuchos::RCP<Epetra_Vector> rs = Teuchos::rcp(new Epetra_Vector(kss->RowMap(), true));
  Teuchos::RCP<Epetra_Vector> rt = Teuchos::rcp(new Epetra_Vector(ktt->RowMap(), true));
  LINALG::Export(*combined_RHS, *rs);
  LINALG::Export(*combined_RHS, *rt);

  // we don't want the rhs but the residual
  rs->Scale(-1.);
  rt->Scale(-1.);

  // add last time step contact forces to rhs
  if (fscn_ != Teuchos::null)  // in the first time step, we don't have any history of the
                               // contact force, after that, fscn_ should be initialized propperly
  {
    Epetra_Vector tmp(kss->RowMap());
    LINALG::Export(*fscn_, tmp);
    if (rs->Update(alphaf_, tmp, 1.) != 0)  // fscn already scaled with alphaf_ in update
      dserror("update went wrong");
  }


  // map containing the inactive and non-contact structural dofs
  Teuchos::RCP<Epetra_Map> str_gni_dofs = LINALG::SplitMap(
      *LINALG::SplitMap(kss->RowMap(), *interface_->MasterRowDofs()), *interface_->ActiveDofs());

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
    sysmat->Assign(0, 0, LINALG::Copy, *kss);
    sysmat->Assign(0, 1, LINALG::Copy, *kst);
    sysmat->Assign(1, 0, LINALG::Copy, *kts);
    sysmat->Assign(1, 1, LINALG::Copy, *ktt);
    return;
  }

  // split matrix blocks in 3 rows: Active, Master and (Inactive+others)
  Teuchos::RCP<LINALG::SparseMatrix> kss_ni, kss_m, kss_a, kst_ni, kst_m, kst_a, kts_ni, kts_m,
      kts_a, ktt_ni, ktt_m, ktt_a, dummy1, dummy2, dummy3;

  // temporary matrix
  Teuchos::RCP<LINALG::SparseMatrix> tmp;
  Teuchos::RCP<Epetra_Vector> tmpv;

  // an empty dummy map
  Teuchos::RCP<Epetra_Map> dummy_map1, dummy_map2;

  // ****************************************************
  // split kss block*************************************
  // ****************************************************
  // split first row
  LINALG::SplitMatrix2x2(
      kss, str_gni_dofs, dummy_map1, gdisp_DofRowMap, dummy_map2, kss_ni, dummy1, tmp, dummy2);

  // this shoud be a split in rows, so that two blocks should have zero columns
  if (dummy1->DomainMap().NumGlobalElements() != 0 || dummy2->DomainMap().NumGlobalElements() != 0)
    dserror("this split should only split rows, no columns expected for this matrix blocks");

  // reset
  dummy1 = Teuchos::null;
  dummy2 = Teuchos::null;
  dummy_map1 = Teuchos::null;
  dummy_map2 = Teuchos::null;

  // split the remaining two rows
  LINALG::SplitMatrix2x2(
      tmp, gmdof, dummy_map1, gdisp_DofRowMap, dummy_map2, kss_m, dummy1, kss_a, dummy2);

  // this shoud be a split in rows, so that two blocks should have zero columns
  if (dummy1->DomainMap().NumGlobalElements() != 0 || dummy2->DomainMap().NumGlobalElements() != 0)
    dserror("this split should only split rows, no columns expected for this matrix blocks");

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
  LINALG::SplitMatrix2x2(
      kst, str_gni_dofs, dummy_map1, gpres_DofRowMap, dummy_map2, kst_ni, dummy1, tmp, dummy2);

  // this shoud be a split in rows, so that two blocks should have zero columns
  if (dummy1->DomainMap().NumGlobalElements() != 0 || dummy2->DomainMap().NumGlobalElements() != 0)
    dserror("this split should only split rows, no columns expected for this matrix blocks");

  // reset
  dummy1 = Teuchos::null;
  dummy2 = Teuchos::null;
  dummy_map1 = Teuchos::null;
  dummy_map2 = Teuchos::null;

  // split the remaining two rows
  LINALG::SplitMatrix2x2(
      tmp, gmdof, dummy_map1, gpres_DofRowMap, dummy_map2, kst_m, dummy1, kst_a, dummy2);

  // this shoud be a split in rows, so that two blocks should have zero columns
  if (dummy1->DomainMap().NumGlobalElements() != 0 || dummy2->DomainMap().NumGlobalElements() != 0)
    dserror("this split should only split rows, no columns expected for this matrix blocks");

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
  LINALG::Export(*rs, rsni);
  Epetra_Vector rsm(*interface_->MasterRowDofs());
  LINALG::Export(*rs, rsm);
  Teuchos::RCP<Epetra_Vector> rsa = Teuchos::rcp(new Epetra_Vector(*interface_->ActiveDofs()));
  LINALG::Export(*rs, *rsa);
  // ****************************************************
  // split rhs vectors***********************************
  // ****************************************************

  // D and M matrix for the active nodes
  Teuchos::RCP<LINALG::SparseMatrix> dInvA =
      Teuchos::rcp(new LINALG::SparseMatrix(*interface_->ActiveDofs(), 100, true, false));
  Teuchos::RCP<LINALG::SparseMatrix> mA =
      Teuchos::rcp(new LINALG::SparseMatrix(*interface_->ActiveDofs(), 100, true, false));

  dummy_map1 = dummy_map2 = Teuchos::null;
  dummy1 = dummy2 = dummy3 = Teuchos::null;
  LINALG::SplitMatrix2x2(
      dmatrix, active_dofs, dummy_map1, active_dofs, dummy_map2, dInvA, dummy1, dummy2, dummy3);
  dInvA->Complete(*interface_->ActiveDofs(), *interface_->ActiveDofs());
  // invert D-matrix
  Epetra_Vector dDiag(*interface_->ActiveDofs());
  dInvA->ExtractDiagonalCopy(dDiag);
  if (dDiag.Reciprocal(dDiag)) dserror("inversion of diagonal D matrix failed");
  dInvA->ReplaceDiagonalValues(dDiag);

  dummy_map1 = dummy_map2 = Teuchos::null;
  dummy1 = dummy2 = dummy3 = Teuchos::null;
  LINALG::SplitMatrix2x2(
      mmatrix, active_dofs, dummy_map1, gmdof, dummy_map2, mA, dummy1, dummy2, dummy3);
  mA->Complete(*interface_->MasterRowDofs(), *interface_->ActiveDofs());

  // get dinv * M
  Teuchos::RCP<LINALG::SparseMatrix> dInvMa =
      LINALG::MLMultiply(*dInvA, false, *mA, false, false, false, true);

  // we need to add another term, since AssembleLinStick/Slip assumes that we solve
  // for the Lagrange multiplier increments. However, we solve for the LM directly.
  // We can do that, since the system is linear in the LMs.
  tmpv = Teuchos::rcp(new Epetra_Vector(*interface_->ActiveDofs()));
  Teuchos::RCP<Epetra_Vector> tmpv2 = Teuchos::rcp(new Epetra_Vector(*interface_->ActiveDofs()));
  LINALG::Export(*z_, *tmpv2);
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
  if (sdirichtoggle_.is_null()) dserror("you didn't call StoreDirichletStatus");
  if (constr_direction_ == INPAR::CONTACT::constr_xyz)
  {
    double haveDBC = 0;
    sdirichtoggle_->Norm1(&haveDBC);
    if (haveDBC > 0.)
    {
      Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(*interface_->ActiveDofs(), true);
      dInvA->ExtractDiagonalCopy(*diag);
      Teuchos::RCP<Epetra_Vector> lmDBC = LINALG::CreateVector(*interface_->ActiveDofs(), true);
      LINALG::Export(*sdirichtoggle_, *lmDBC);
      Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*interface_->ActiveDofs(), true);
      tmp->Multiply(1., *diag, *lmDBC, 0.);
      diag->Update(-1., *tmp, 1.);
      dInvA->ReplaceDiagonalValues(*diag);
      dInvMa = LINALG::MLMultiply(*dInvA, false, *mA, false, false, false, true);
    }
  }

  // reset the tangent stiffness
  // (for the condensation we have constructed copies above)
  sysmat->UnComplete();

  // need diagonal block kss with explicitdirichtlet_=true
  // to be able to apply dirichlet values for contact symmetry condition
  LINALG::SparseMatrix tmpkss(*gdisp_DofRowMap, 100, false, false, LINALG::SparseMatrix::FE_MATRIX);
  sysmat->Assign(0, 0, LINALG::Copy, tmpkss);

  // get references to the blocks (just for convenience)
  LINALG::SparseMatrix& kss_new = sysmat->Matrix(0, 0);
  LINALG::SparseMatrix& kst_new = sysmat->Matrix(0, 1);
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
  kss_new.Add(*LINALG::MLMultiply(*dInvMa, true, *kss_a, false, false, false, true), false, 1., 1.);
  kst_new.Add(*LINALG::MLMultiply(*dInvMa, true, *kst_a, false, false, false, true), false, 1., 1.);
  tmpv = Teuchos::rcp(new Epetra_Vector(*interface_->MasterRowDofs()));
  if (dInvMa->Multiply(true, *rsa, *tmpv)) dserror("multiply failed");
  CONTACT::UTILS::AddVector(*tmpv, *combined_RHS);
  tmpv = Teuchos::null;

  // third row
  Teuchos::RCP<LINALG::SparseMatrix> wDinv =
      LINALG::MLMultiply(*dcsdLMc, false, *dInvA, true, false, false, true);
  kss_new.Add(*LINALG::MLMultiply(*wDinv, false, *kss_a, false, false, false, true), false,
      -1. / (1. - alphaf_), 1.);
  kst_new.Add(*LINALG::MLMultiply(*wDinv, false, *kst_a, false, false, false, true), false,
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
    if (!node) dserror("node not found");
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(node);
    if (!cnode) dserror("not a contact node");

    cnode->FriData().GetDerivJump().resize(3);
    // write it to nodes
    for (int dim = 0; dim < interface_->Dim(); dim++)
    {
      cnode->FriData().jump()[dim] = cnode->CoEhlData().GetWeightedRelTangVel()(dim);
      for (auto p = cnode->CoEhlData().GetWeightedRelTangVelDeriv().begin();
           p != cnode->CoEhlData().GetWeightedRelTangVelDeriv().end(); ++p)
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
      dserror("some data for LM recovery is missing");

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

    LINALG::Export(tmp, *z_);
  }

  else
    z_ = Teuchos::rcp(new Epetra_Vector(*interface_->SlaveRowDofs()));

  if (z_old != Teuchos::null)
  {
    z_old->Update(-1., *z_, 1.);
    z_old->Norm2(&contact_LM_incr_norm_);
  }

  // store updated LM into nodes
  for (int i = 0; i < interface_->SlaveColNodes()->NumMyElements(); ++i)
  {
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(interface_->Discret().lColNode(i));
    for (int dof = 0; dof < interface_->Dim(); ++dof)
      cnode->MoData().lm()[dof] = z_->operator[](z_->Map().LID(cnode->Dofs()[dof]));
  }

  return;
}



/*----------------------------------------------------------------------*
 |  Store dirichlet B.C. status into CNode                    popp 06/09|
 *----------------------------------------------------------------------*/
void ADAPTER::CouplingEhlMortar::StoreDirichletStatus(
    Teuchos::RCP<const LINALG::MapExtractor> dbcmaps)
{
  // loop over all slave row nodes on the current interface
  for (int j = 0; j < interface_->SlaveRowNodes()->NumMyElements(); ++j)
  {
    int gid = interface_->SlaveRowNodes()->GID(j);
    DRT::Node* node = interface_->Discret().gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);

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
        dserror("Inconsistency in structure Dirichlet conditions and Mortar symmetry conditions");
      }
    }
  }
  // create old style dirichtoggle vector (supposed to go away)
  sdirichtoggle_ = Teuchos::rcp(new Epetra_Vector(*interface_->SlaveRowDofs(), true));
  Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*(dbcmaps->CondMap())));
  temp->PutScalar(1.0);
  LINALG::Export(*temp, *sdirichtoggle_);

  return;
}



bool ADAPTER::CouplingEhlMortar::AlreadyEvaluated(Teuchos::RCP<const Epetra_Vector> disp)
{
  if (evaluated_state_.is_null()) return false;
  Teuchos::RCP<Epetra_Vector> diff = Teuchos::rcp(new Epetra_Vector(*disp));
  if (diff->Update(-1., *evaluated_state_, 1.)) dserror("update failed");
  double inf_diff = -1.;
  if (diff->NormInf(&inf_diff)) dserror("NormInf failed");
  if (inf_diff < 1.e-13) return true;

  return false;
}

Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::CouplingEhlMortar::AssembleEHLLinD(
    const Teuchos::RCP<Epetra_Vector> x  // slave dof vector
)
{
  Teuchos::RCP<LINALG::SparseMatrix> DLinEHL = Teuchos::rcp(
      new LINALG::SparseMatrix(*slavedofrowmap_, 81, true, false, LINALG::SparseMatrix::FE_MATRIX));
  DLinEHL->Zero();
  DLinEHL->UnComplete();

  interface_->AssembleCoupLinD(*DLinEHL, x);

  DLinEHL->Complete(*smdofrowmap_, *slavedofrowmap_);

  return DLinEHL;
}

Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::CouplingEhlMortar::AssembleEHLLinM(
    const Teuchos::RCP<Epetra_Vector> x  // slave dof vector
)
{
  Teuchos::RCP<LINALG::SparseMatrix> MLinEHL = Teuchos::rcp(new LINALG::SparseMatrix(
      *masterdofrowmap_, 81, true, false, LINALG::SparseMatrix::FE_MATRIX));
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
    if (!node) dserror("node not found");
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    if (!cnode) dserror("not a contact node");

    for (int d = 0; d < interface_->Dim(); ++d)
      normals_->ReplaceGlobalValue(cnode->Dofs()[d], 0, cnode->MoData().n()[d]);
  }
}


void ADAPTER::CouplingEhlMortar::AssembleNormalsDeriv()
{
  Nderiv_ = Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap_, 81, false, false));
  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) dserror("node not found");
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    if (!cnode) dserror("not a contact node");

    for (int d = 0; d < Interface()->Dim(); ++d)
      for (auto p = cnode->CoData().GetDerivN()[d].begin();
           p != cnode->CoData().GetDerivN()[d].end(); ++p)
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
    if (!node) dserror("node not found");
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    if (!cnode) dserror("not a contact node");
    double real_gap = cnode->CoData().Getg();
    switch (cnode->MoData().GetD().size())
    {
      case 0:
        break;
      case 1:
        if (cnode->MoData().GetD().begin()->first != cnode->Id())
          dserror("something is wrong. Here should by my own Id");
        real_gap /= cnode->MoData().GetD().at(cnode->Id());
        break;
      default:
        dserror(
            "GetD should be of size 0 (unprojectable) or 1 (projectable). Are you not using "
            "duals?");
    }
    nodal_gap_->ReplaceGlobalValue(cnode->Id(), 0, real_gap);
  }

  static const double offset =
      DRT::Problem::Instance()->LubricationDynamicParams().get<double>("GAP_OFFSET");
  for (int i = 0; i < nodal_gap_->Map().NumMyElements(); ++i) nodal_gap_->operator[](i) += offset;
}

void ADAPTER::CouplingEhlMortar::AssembleRealGapDeriv()
{
  deriv_nodal_gap_ = Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap_, 81, false, false));

  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) dserror("node not found");
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    if (!cnode) dserror("not a contact node");

    if (cnode->CoData().GetDerivD().size() != cnode->MoData().GetD().size())
      dserror("size inconsistency");

    const double w_gap = cnode->CoData().Getg();
    double d = -1.;
    switch (cnode->CoData().GetDerivD().size())
    {
      case 0:
        break;
      case 1:
        if (cnode->CoData().GetDerivD().begin()->first != cnode->Id())
          dserror("something is wrong. Here should by my own Id");
        d = cnode->MoData().GetD().at(cnode->Id());
        break;
      default:
        dserror(
            "GetDerivD should be of size 0 (unprojectable) or 1 (projectable). Are you not using "
            "duals?");
    }

    if (cnode->CoData().GetDerivD().size())
      for (auto p = cnode->CoData().GetDerivD().at(cnode->Id()).begin();
           p != cnode->CoData().GetDerivD().at(cnode->Id()).end(); ++p)
      {
        const double val = -w_gap / (d * d) * p->second;
        for (int d = 0; d < interface_->Dim(); ++d)
          deriv_nodal_gap_->Assemble(val, cnode->Dofs()[d], p->first);
      }

    if (d == -1 && cnode->CoData().GetDerivG().size() != 0) dserror("inconsistency");

    if (cnode->CoData().GetDerivG().size())
      for (auto p = cnode->CoData().GetDerivG().begin(); p != cnode->CoData().GetDerivG().end();
           ++p)
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
  relTangVel_deriv_ = Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap_, 81, false, false));
  avTangVel_deriv_ = Teuchos::rcp(new LINALG::SparseMatrix(*slavedofrowmap_, 81, false, false));

  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) dserror("node not found");
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    if (!cnode) dserror("not a contact node");


    double d_val = 0.;
    switch (cnode->MoData().GetD().size())
    {
      case 0:
        break;
      case 1:
        if (cnode->MoData().GetD().begin()->first != cnode->Id())
          dserror("something is wrong. Here should by my own Id");
        d_val = cnode->MoData().GetD().at(cnode->Id());
        break;
      default:
        dserror(
            "GetD should be of size 0 (unprojectable) or 1 (projectable). Are you not using "
            "duals?");
    }

    if (d_val == 0.) continue;

    for (int d = 0; d < Interface()->Dim(); ++d)
    {
      relTangVel_->ReplaceGlobalValue(
          cnode->Dofs()[d], 0, cnode->CoEhlData().GetWeightedRelTangVel()(d) / d_val);
      avTangVel_->ReplaceGlobalValue(
          cnode->Dofs()[d], 0, cnode->CoEhlData().GetWeightedAvTangVel()(d) / d_val);
    }

    for (auto p = cnode->CoData().GetDerivD().at(cnode->Id()).begin();
         p != cnode->CoData().GetDerivD().at(cnode->Id()).end(); ++p)
    {
      const int col = p->first;
      for (int d = 0; d < Interface()->Dim(); ++d)
      {
        const int row = cnode->Dofs()[d];
        const double rel_val =
            -cnode->CoEhlData().GetWeightedRelTangVel()(d) / (d_val * d_val) * p->second;
        const double av_val =
            -cnode->CoEhlData().GetWeightedAvTangVel()(d) / (d_val * d_val) * p->second;
        relTangVel_deriv_->Assemble(rel_val, row, col);
        avTangVel_deriv_->Assemble(av_val, row, col);
      }
    }
    for (auto p = cnode->CoEhlData().GetWeightedAvTangVelDeriv().begin();
         p != cnode->CoEhlData().GetWeightedAvTangVelDeriv().end(); ++p)
    {
      const int col = p->first;
      for (int d = 0; d < Interface()->Dim(); ++d)
      {
        const int row = cnode->Dofs()[d];
        const double val = p->second(d) / d_val;
        avTangVel_deriv_->Assemble(val, row, col);
      }
    }
    for (auto p = cnode->CoEhlData().GetWeightedRelTangVelDeriv().begin();
         p != cnode->CoEhlData().GetWeightedRelTangVelDeriv().end(); ++p)
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
  SurfGrad_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *slavedofrowmap_, 81, false, false, LINALG::SparseMatrix::FE_MATRIX));

  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) dserror("ERROR: Cannot find node");
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    if (!cnode) dserror("this is not a contact node");

    double dval = 1.;
    switch (cnode->MoData().GetD().size())
    {
      case 0:
        dval = 1.e32;
        break;  // large number so no tangential gradient
      case 1:
        if (cnode->MoData().GetD().begin()->first != cnode->Id())
          dserror("something is wrong. Here should by my own Id");
        dval = cnode->MoData().GetD().at(cnode->Id());
        break;
      default:
        dserror(
            "GetD should be of size 0 (unprojectable) or 1 (projectable). Are you not using "
            "duals?");
    }

    for (auto p = cnode->CoEhlData().GetSurfGrad().begin();
         p != cnode->CoEhlData().GetSurfGrad().end(); ++p)
      for (int d = 0; d < Interface()->Dim(); ++d)
        SurfGrad_->Assemble(p->second(d) / dval, cnode->Dofs()[d], p->first);
  }

  SurfGrad_->Complete();
}

Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::CouplingEhlMortar::AssembleSurfGradDeriv(
    const Teuchos::RCP<const Epetra_Vector> x)
{
  Teuchos::RCP<LINALG::SparseMatrix> SurfGradDeriv = Teuchos::rcp(new LINALG::SparseMatrix(
      *slavedofrowmap_, 81, false, false, LINALG::SparseMatrix::FE_MATRIX));

  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    DRT::Node* node = Interface()->Discret().gNode(interface_->SlaveRowNodes()->GID(i));
    if (!node) dserror("ERROR: Cannot find node");
    CONTACT::CoNode* cnode = dynamic_cast<CONTACT::CoNode*>(node);
    if (!cnode) dserror("this is not a contact node");

    double dval = 1.;
    switch (cnode->MoData().GetD().size())
    {
      case 0:
        dval = 1.e32;
        break;  // large number so no tangential gradient
      case 1:
        if (cnode->MoData().GetD().begin()->first != cnode->Id())
          dserror("something is wrong. Here should by my own Id");
        dval = cnode->MoData().GetD().at(cnode->Id());
        break;
      default:
        dserror(
            "GetD should be of size 0 (unprojectable) or 1 (projectable). Are you not using "
            "duals?");
    }

    for (auto p = cnode->CoEhlData().GetSurfGradDeriv().begin();
         p != cnode->CoEhlData().GetSurfGradDeriv().end(); ++p)
    {
      const int col = p->first;
      for (auto q = p->second.begin(); q != p->second.end(); ++q)
      {
        const int lid = x->Map().LID(q->first);
        if (lid < 0) dserror("not my gid");
        const double x_val = x->operator[](lid);
        for (int d = 0; d < Interface()->Dim(); ++d)
        {
          const double val = x_val * q->second(d) / dval;
          SurfGradDeriv->Assemble(val, cnode->Dofs()[d], col);
        }
      }
    }

    if (cnode->CoData().GetDerivD().size())
      for (auto p = cnode->CoData().GetDerivD().at(cnode->Id()).begin();
           p != cnode->CoData().GetDerivD().at(cnode->Id()).end(); ++p)
      {
        const int col = p->first;

        for (auto q = cnode->CoEhlData().GetSurfGrad().begin();
             q != cnode->CoEhlData().GetSurfGrad().end(); ++q)
          for (int d = 0; d < Interface()->Dim(); ++d)
          {
            const int row = cnode->Dofs()[d];
            const int x_gid = q->first;
            const int x_lid = x->Map().LID(x_gid);
            if (x_lid < 0) dserror("not my gid");
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
    if (!cnode) dserror("cast failed");
    const LINALG::Matrix<3, 1> lm(cnode->MoData().lm(), true);
    const LINALG::Matrix<3, 1> nor(cnode->MoData().n(), true);
    LINALG::Matrix<3, 3> nn;
    nn.MultiplyNT(nor, nor);
    LINALG::Matrix<3, 1> lmn;
    lmn.Multiply(nn, lm);
    LINALG::Matrix<3, 1> lmt(lm);
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
  if (active_old != NULL)
    *active_old = Teuchos::rcp(new Epetra_Vector(*interface_->SlaveRowNodes()));
  for (int i = 0; i < interface_->SlaveRowNodes()->NumMyElements(); ++i)
  {
    CONTACT::FriNode* cnode = dynamic_cast<CONTACT::FriNode*>(interface_->Discret().lRowNode(i));
    if (!cnode) dserror("cast failed");
    if (cnode->Active())
      (*active)->operator[](i) = 1.;
    else
      (*active)->operator[](i) = 0.;
    if (cnode->FriData().Slip())
      (*slip)->operator[](i) = 1.;
    else
      (*slip)->operator[](i) = 0.;

    if (active_old != NULL)
    {
      if (cnode->CoData().ActiveOld())
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
    if (!cnode) dserror("cast failed");
    cnode->Active() = active_toggle->operator[](i);
    cnode->FriData().Slip() = slip_toggle->operator[](i);
    cnode->CoData().ActiveOld() = active_old_toggle->operator[](i);
    for (int d = 0; d < interface_->Dim(); ++d)
      cnode->MoData().lm()[d] = z_->operator[](z_->Map().LID(cnode->Dofs()[d]));
  }
}
