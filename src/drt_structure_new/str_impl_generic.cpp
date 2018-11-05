/*-----------------------------------------------------------*/
/*!
\file str_impl_generic.cpp

\brief Generic class for all implicit time integrators.

\maintainer Michael Hiermeier

\date Mar 11, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_impl_generic.H"
#include "str_timint_implicit.H"
#include "str_model_evaluator.H"
#include "str_model_evaluator_data.H"

#include "../solver_nonlin_nox/nox_nln_group.H"
#include "../solver_nonlin_nox/nox_nln_aux.H"
#include "../solver_nonlin_nox/nox_nln_solver_linesearchbased.H"
#include "../solver_nonlin_nox/nox_nln_group_prepostoperator.H"

#include "../drt_io/io_pstream.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::IMPLICIT::Generic::Generic() : ispredictor_state_(false)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Generic::Setup()
{
  CheckInit();
  // call base class first
  STR::Integrator::Setup();
  // ---------------------------------------------------------------------------
  // set the new pre/post operator for the nox nln group in the parameter list
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList& p_grp_opt = SDyn().GetMutableNoxParams().sublist("Group Options");

  // create the new generic pre/post operator
  Teuchos::RCP<NOX::NLN::Abstract::PrePostOperator> prepost_generic_ptr =
      Teuchos::rcp(new NOX::NLN::PrePostOp::IMPLICIT::Generic(*this));

  // Get the current map. If there is no map, return a new empty one. (reference)
  NOX::NLN::GROUP::PrePostOperator::Map& prepostgroup_map =
      NOX::NLN::GROUP::PrePostOp::GetMutableMap(p_grp_opt);

  // insert/replace the old pointer in the map
  prepostgroup_map[NOX::NLN::GROUP::prepost_impl_generic] = prepost_generic_ptr;

  // ---------------------------------------------------------------------------
  // set the new pre/post operator for the nox nln solver in the parameter list
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList& p_sol_opt = SDyn().GetMutableNoxParams().sublist("Solver Options");

  NOX::NLN::AUX::AddToPrePostOpVector(p_sol_opt, prepost_generic_ptr);

  // No issetup_ = true, since the Setup() functions of the derived classes
  // have to be called and finished first!

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Generic::SetIsPredictorState(const bool& ispredictor_state)
{
  ispredictor_state_ = ispredictor_state;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const bool& STR::IMPLICIT::Generic::IsPredictorState() const { return ispredictor_state_; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::ParameterList& STR::IMPLICIT::Generic::GetNoxParams()
{
  return SDyn().GetMutableNoxParams();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::Generic::GetDefaultStepLength() const
{
  const Teuchos::ParameterList& p_nox = TimInt().GetDataSDyn().GetNoxParams();
  const std::string& nln_solver = p_nox.get<std::string>("Nonlinear Solver");
  // The pseudo transient implementation holds also a line search object!
  if (nln_solver == "Line Search Based" or nln_solver == "Pseudo Transient")
  {
    const Teuchos::ParameterList& p_ls = p_nox.sublist("Line Search");
    const std::string& method = p_ls.get<std::string>("Method");
    const Teuchos::ParameterList& p_method = p_ls.sublist(method);
    if (p_method.isParameter("Default Step")) return p_method.get<double>("Default Step");
  }
  // default: return a step length of 1.0
  return 1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Generic::ResetEvalParams()
{
  // set the time step dependent parameters for the element evaluation
  EvalData().SetTotalTime(GlobalState().GetTimeNp());
  EvalData().SetDeltaTime((*GlobalState().GetDeltaTime())[0]);
  EvalData().SetIsTolerateError(true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Generic::PrintJacobianInMatlabFormat(const NOX::NLN::Group& curr_grp) const
{
  const STR::TIMINT::Implicit& timint_impl = dynamic_cast<const STR::TIMINT::Implicit&>(TimInt());

  timint_impl.PrintJacobianInMatlabFormat(curr_grp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::Generic::ApplyCorrectionSystem(const enum NOX::NLN::CorrectionType type,
    const std::vector<INPAR::STR::ModelType>& constraint_models, const Epetra_Vector& x,
    Epetra_Vector& f, LINALG::SparseOperator& jac)
{
  CheckInitSetup();

  ResetEvalParams();

  EvalData().SetCorrectionType(type);

  bool ok = false;
  switch (type)
  {
    case NOX::NLN::CorrectionType::soc_full:
    {
      // Do a standard full step.
      /* Note that there is a difference, since we tagged this evaluation by
       * setting it to a non-default step. */
      ok = ApplyForceStiff(x, f, jac);
      break;
    }
    case NOX::NLN::CorrectionType::soc_cheap:
    {
      ok = ModelEval().ApplyCheapSOCRhs(type, constraint_models, x, f, 1.0);
      if (not jac.Filled()) dserror("The jacobian is supposed to be filled at this point!");

      break;
    }
    default:
    {
      dserror(
          "No action defined for the given second order correction type: "
          "\"%s\"",
          NOX::NLN::CorrectionType2String(type).c_str());
      exit(EXIT_FAILURE);
    }
  }

  if (not ok) return false;

  if (not jac.Filled()) jac.Complete();

  return ok;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Generic::ConditionNumber(const NOX::NLN::Group& grp) const
{
  const STR::TIMINT::Implicit& timint_impl = dynamic_cast<const STR::TIMINT::Implicit&>(TimInt());

  timint_impl.ComputeConditionNumber(grp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::PrePostOp::IMPLICIT::Generic::runPreComputeX(const NOX::NLN::Group& input_grp,
    const Epetra_Vector& dir, const double& step, const NOX::NLN::Group& curr_grp)
{
  // set the evaluation parameters
  const Epetra_Vector& xold =
      dynamic_cast<const NOX::Epetra::Vector&>(input_grp.getX()).getEpetraVector();
  Epetra_Vector& dir_mutable = const_cast<Epetra_Vector&>(dir);

  const bool isdefaultstep = (step == default_step_);
  impl_.ModelEval().RunPreComputeX(xold, dir_mutable, step, curr_grp, isdefaultstep);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::PrePostOp::IMPLICIT::Generic::runPostComputeX(const NOX::NLN::Group& input_grp,
    const Epetra_Vector& dir, const double& step, const NOX::NLN::Group& curr_grp)
{
  // set the evaluation parameters
  const Epetra_Vector& xold =
      dynamic_cast<const NOX::Epetra::Vector&>(input_grp.getX()).getEpetraVector();
  const Epetra_Vector& xnew =
      dynamic_cast<const NOX::Epetra::Vector&>(curr_grp.getX()).getEpetraVector();

  bool isdefaultstep = (step == default_step_);
  impl_.ModelEval().RunPostComputeX(xold, dir, step, xnew, isdefaultstep);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::PrePostOp::IMPLICIT::Generic::runPostIterate(const NOX::Solver::Generic& solver)
{
  double step = 0.0;
  const bool isdefaultstep = getStep(step, solver);
  const int num_corrs = getNumberOfModifiedNewtonCorrections(solver);

  impl_.ModelEval().RunPostIterate(solver, step, isdefaultstep, num_corrs);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::PrePostOp::IMPLICIT::Generic::runPreSolve(const NOX::Solver::Generic& solver)
{
  double step = 0.0;
  const bool isdefaultstep = getStep(step, solver);

  impl_.ModelEval().RunPreSolve(solver, step, isdefaultstep);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::PrePostOp::IMPLICIT::Generic::runPostApplyJacobianInverse(
    const NOX::Abstract::Vector& rhs, NOX::Abstract::Vector& result,
    const NOX::Abstract::Vector& xold, const NOX::NLN::Group& grp)
{
  impl_.ModelEval().RunPostApplyJacobianInverse(
      convert2EpetraVector(rhs), convert2EpetraVector(result), convert2EpetraVector(xold), grp);

  impl_.PrintJacobianInMatlabFormat(grp);
  impl_.ConditionNumber(grp);

  // reset any possible set correction type at this point
  const STR::MODELEVALUATOR::Data& eval_data = impl_.EvalData();
  const_cast<STR::MODELEVALUATOR::Data&>(eval_data).SetCorrectionType(
      NOX::NLN::CorrectionType::vague);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_Vector& NOX::NLN::PrePostOp::IMPLICIT::Generic::convert2EpetraVector(
    NOX::Abstract::Vector& vec) const
{
  NOX::Epetra::Vector* epetra_vec = dynamic_cast<NOX::Epetra::Vector*>(&vec);
  if (epetra_vec == NULL) dserror("The given NOX::Abstract::Vector is no NOX::Epetra::Vector!");

  return epetra_vec->getEpetraVector();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Vector& NOX::NLN::PrePostOp::IMPLICIT::Generic::convert2EpetraVector(
    const NOX::Abstract::Vector& vec) const
{
  const NOX::Epetra::Vector* epetra_vec = dynamic_cast<const NOX::Epetra::Vector*>(&vec);
  if (epetra_vec == NULL) dserror("The given NOX::Abstract::Vector is no NOX::Epetra::Vector!");

  return epetra_vec->getEpetraVector();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::PrePostOp::IMPLICIT::Generic::getStep(
    double& step, const NOX::Solver::Generic& solver) const
{
  // try to cast the given solver object
  const NOX::NLN::Solver::LineSearchBased* ls_solver =
      dynamic_cast<const NOX::NLN::Solver::LineSearchBased*>(&solver);

  bool isdefaultstep = false;

  if (not ls_solver)
  {
    step = default_step_;
    isdefaultstep = true;
  }
  else
  {
    step = ls_solver->getStepSize();
    isdefaultstep = (step == default_step_);
  }

  return isdefaultstep;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int NOX::NLN::PrePostOp::IMPLICIT::Generic::getNumberOfModifiedNewtonCorrections(
    const NOX::Solver::Generic& solver) const
{
  int number_of_corr = 0;
  const Teuchos::ParameterList& pmod =
      solver.getList().sublist("Direction").sublist("Newton").sublist("Modified");

  if (pmod.isParameter("Number of Corrections"))
    number_of_corr = pmod.get<int>("Number of Corrections");

  return number_of_corr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Generic::ComputeJacobianContributionsFromElementLevelForPTC(
    Teuchos::RCP<LINALG::SparseMatrix>& scalingMatrixOpPtr)
{
  ModelEval().ComputeJacobianContributionsFromElementLevelForPTC(scalingMatrixOpPtr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::Generic::RemoveCondensedContributionsFromRhs(Epetra_Vector& rhs) const
{
  ModelEval().RemoveCondensedContributionsFromRhs(rhs);
}
