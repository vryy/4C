/*----------------------------------------------------------------------------*/
/** \file
\brief NOX interface implementation for the xcontact level-set elliptical
       reinitialization algorithm


\level 3
*/
/*----------------------------------------------------------------------------*/

#include "xcontact_levelset_reinit_elliptic.H"
#include "xcontact_levelset_algorithm.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_sparsematrix.H"

#include "../drt_lib/drt_utils_discret.H"

#include "../solver_nonlin_nox/nox_nln_aux.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool XCONTACT::LEVELSET::REINIT::Elliptic::computeF(
    const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  dsassert(fillFlag == NOX::Epetra::Interface::Required::Residual, "Unsupported fill flag!");

  PreEvaluate(x);

  // set parameters
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", SCATRA::calc_rhs);
  SetElementParameters(eleparams);
  Algorithm().AddTimeIntegrationSpecificVectors();

  // time measurement: element calls
  TEUCHOS_FUNC_TIME_MONITOR("XCONTACT::LEVELSET::REINIT::Elliptic::computeF");

  // reset right-hand-side
  F.PutScalar(0.0);

  DRT::UTILS::Evaluate(*Algorithm().Discretization(), eleparams, Teuchos::null,
      Teuchos::rcpFromRef(F), &Algorithm().ActiveColEleMap());

  Algorithm().Discretization()->ClearState();

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool XCONTACT::LEVELSET::REINIT::Elliptic::computeFandJacobian(
    const Epetra_Vector& x, Epetra_Vector& F, Epetra_Operator& jac)
{
  PreEvaluate(x);

  // set parameters
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", SCATRA::calc_mat_and_rhs);
  SetElementParameters(eleparams);
  Algorithm().AddTimeIntegrationSpecificVectors();

  // time measurement: element calls
  TEUCHOS_FUNC_TIME_MONITOR("XCONTACT::LEVELSET::REINIT::Elliptic::computeFandJacobian");

  // reset jacobian and right-hand-side
  LINALG::SparseOperator* jac_linalg = dynamic_cast<LINALG::SparseOperator*>(&jac);
  dsassert(jac_linalg != NULL, "Dynamic cast failed");

  jac_linalg->Zero();
  F.PutScalar(0.0);

  DRT::UTILS::Evaluate(*Algorithm().Discretization(), eleparams, Teuchos::rcp(jac_linalg, false),
      Teuchos::rcpFromRef(F), &Algorithm().ActiveColEleMap());

  if (not jac_linalg->Filled()) jac_linalg->Complete();

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool XCONTACT::LEVELSET::REINIT::Elliptic::computeJacobian(
    const Epetra_Vector& x, Epetra_Operator& jac)
{
  PreEvaluate(x);

  // set parameters
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", SCATRA::calc_mat);
  SetElementParameters(eleparams);
  Algorithm().AddTimeIntegrationSpecificVectors();

  // time measurement: element calls
  TEUCHOS_FUNC_TIME_MONITOR("XCONTACT::LEVELSET::REINIT::Elliptic::computeJacobian");

  // reset jacobian and right-hand-side
  LINALG::SparseOperator* jac_linalg = dynamic_cast<LINALG::SparseOperator*>(&jac);
  dsassert(jac_linalg != NULL, "Dynamic cast failed");

  jac_linalg->Zero();

  DRT::UTILS::Evaluate(*Algorithm().Discretization(), eleparams, Teuchos::rcp(jac_linalg, false),
      Teuchos::null, &Algorithm().ActiveColEleMap());

  if (not jac_linalg->Filled()) jac_linalg->Complete();

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::REINIT::Elliptic::PreEvaluate(const Epetra_Vector& x)
{
  if (Algorithm().UseProjection())
  {
    // calculate the l2-projected gradient
    Algorithm().NodalGradBasedValue() =
        Algorithm().ComputeNodalGradientL2Projection(x, INPAR::SCATRA::l2_proj_system_dual);
  }

  // set the state in the current active time integration scheme
  Algorithm().Discretization()->ClearState();
  Algorithm().SetState(x);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::REINIT::Elliptic::SetElementParameters(Teuchos::ParameterList& eleparams)
{
  // perform a l2-projection of the level-set gradient, if desired
  if (Algorithm().UseProjection())
  {
    Algorithm().Discretization()->AddMultiVectorToParameterList(
        eleparams, "gradphi", Algorithm().NodalGradBasedValue());
    Algorithm().Discretization()->SetState(
        "l2_proj_system_mat_diag", Algorithm().L2ProjSysMatDiagonalPtr());
  }

  // activate the reinitialization routines
  eleparams.set<bool>("solve reinit eq", true);

  // set the nodal velocities
  eleparams.set<int>("ndsvel", Algorithm().nds_vel_);

  // add interface integration cells
  eleparams.set<Teuchos::RCP<const std::map<int, GEO::BoundaryIntCellPtrs>>>(
      "boundary cells", Algorithm().GetZeroIsoLinePtr());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double XCONTACT::LEVELSET::REINIT::Elliptic::GetPrimaryRHSNorms(const Epetra_Vector& F,
    const NOX::NLN::StatusTest::QuantityType& checkquantity,
    const NOX::Abstract::Vector::NormType& type, const bool& isscaled) const
{
  double rhsnorm = -1.0;

  if (checkquantity != NOX::NLN::StatusTest::quantity_levelset_reinit)
    dserror(
        "How did you come here? We support only one quantity type at this "
        "point! ( given quantity = % s )",
        NOX::NLN::StatusTest::QuantityType2String(checkquantity).c_str());

  // wrap the given rhs vector as const NOX::Epetra::Vector
  Epetra_Vector& rhs = const_cast<Epetra_Vector&>(F);
  Teuchos::RCP<Epetra_Vector> rhs_ptr = Teuchos::rcpFromRef(rhs);

  const NOX::Epetra::Vector rhs_nox = NOX::Epetra::Vector(rhs_ptr, NOX::Epetra::Vector::CreateView);

  rhsnorm = rhs_nox.norm(type);
  // do the scaling if desired
  if (isscaled) rhsnorm /= static_cast<double>(rhs_nox.length());

  return rhsnorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double XCONTACT::LEVELSET::REINIT::Elliptic::GetPrimarySolutionUpdateRMS(const Epetra_Vector& xnew,
    const Epetra_Vector& xold, const double& atol, const double& rtol,
    const NOX::NLN::StatusTest::QuantityType& checkquantity,
    const bool& disable_implicit_weighting) const
{
  double rms = -1.0;

  if (checkquantity != NOX::NLN::StatusTest::quantity_levelset_reinit)
    dserror(
        "How did you come here? We support only one quantity type at this "
        "point! ( given quantity = % s )",
        NOX::NLN::StatusTest::QuantityType2String(checkquantity).c_str());

  Epetra_Vector incr = Epetra_Vector(xold);

  incr.Update(1.0, xnew, -1.0);
  rms = NOX::NLN::AUX::RootMeanSquareNorm(
      atol, rtol, Teuchos::rcpFromRef(xnew), Teuchos::rcpFromRef(incr), disable_implicit_weighting);

  return rms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double XCONTACT::LEVELSET::REINIT::Elliptic::GetPrimarySolutionUpdateNorms(
    const Epetra_Vector& xnew, const Epetra_Vector& xold,
    const NOX::NLN::StatusTest::QuantityType& checkquantity,
    const NOX::Abstract::Vector::NormType& type, const bool& isscaled) const
{
  double updatenorm = -1.0;

  if (checkquantity != NOX::NLN::StatusTest::quantity_levelset_reinit)
    dserror(
        "How did you come here? We support only one quantity type at this "
        "point! ( given quantity = % s )",
        NOX::NLN::StatusTest::QuantityType2String(checkquantity).c_str());

  Epetra_Vector incr = Epetra_Vector(xold);

  incr.Update(1.0, xnew, -1.0);
  const NOX::Epetra::Vector incr_nox =
      NOX::Epetra::Vector(Teuchos::rcpFromRef(incr), NOX::Epetra::Vector::CreateView);

  updatenorm = incr_nox.norm(type);
  // do the scaling if desired
  if (isscaled) updatenorm /= static_cast<double>(incr_nox.length());

  return updatenorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double XCONTACT::LEVELSET::REINIT::Elliptic::GetPreviousPrimarySolutionNorms(
    const Epetra_Vector& xold, const NOX::NLN::StatusTest::QuantityType& checkquantity,
    const NOX::Abstract::Vector::NormType& type, const bool& isscaled) const
{
  double xoldnorm = -1.0;

  if (checkquantity != NOX::NLN::StatusTest::quantity_levelset_reinit)
    dserror(
        "How did you come here? We support only one quantity type at this "
        "point! ( given quantity = % s )",
        NOX::NLN::StatusTest::QuantityType2String(checkquantity).c_str());

  // export the displacement solution if necessary
  Epetra_Vector& xold_mutable = const_cast<Epetra_Vector&>(xold);
  const NOX::Epetra::Vector xold_nox =
      NOX::Epetra::Vector(Teuchos::rcpFromRef(xold_mutable), NOX::Epetra::Vector::CreateView);

  xoldnorm = xold_nox.norm(type);
  // do the scaling if desired
  if (isscaled) xoldnorm /= static_cast<double>(xold_nox.length());

  return xoldnorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double XCONTACT::LEVELSET::REINIT::Elliptic::GetModelValue(const Epetra_Vector& x,
    const Epetra_Vector& F, enum NOX::NLN::MeritFunction::MeritFctName mrt_func_type) const
{
  dserror("Currently unsupported!");
  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double XCONTACT::LEVELSET::REINIT::Elliptic::CalcRefNormForce()
{
  dserror("Currently unsupported!");
  exit(EXIT_FAILURE);
}
