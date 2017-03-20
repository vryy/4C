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
#include "str_timint_base.H"
#include "str_model_evaluator.H"
#include "str_model_evaluator_data.H"

#include "../solver_nonlin_nox/nox_nln_group.H"
#include "../solver_nonlin_nox/nox_nln_aux.H"
#include "../solver_nonlin_nox/nox_nln_solver_linesearchbased.H"
#include "../solver_nonlin_nox/nox_nln_group_prepostoperator.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::IMPLICIT::Generic::Generic()
    : ispredictor_state_(false)
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
  Teuchos::ParameterList& p_grp_opt =
      SDyn().GetMutableNoxParams().sublist("Group Options");

  // create the new generic pre/post operator
  Teuchos::RCP<NOX::NLN::Abstract::PrePostOperator> prepost_group_ptr =
      Teuchos::rcp(new NOX::NLN::GROUP::PrePostOp::IMPLICIT::Generic(
          *this ) );

  // Get the current map. If there is no map, return a new empty one. (reference)
  NOX::NLN::GROUP::PrePostOperator::Map& prepostgroup_map =
      NOX::NLN::GROUP::PrePostOp::GetMutableMap(p_grp_opt);

  // insert/replace the old pointer in the map
  prepostgroup_map[NOX::NLN::GROUP::prepost_impl_generic] = prepost_group_ptr;

  // ---------------------------------------------------------------------------
  // set the new pre/post operator for the nox nln solver in the parameter list
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList& p_sol_opt =
      SDyn().GetMutableNoxParams().sublist("Solver Options");

  Teuchos::RCP<NOX::Abstract::PrePostOperator> prepost_solver_ptr =
      Teuchos::rcp( new NOX::NLN::Solver::PrePostOp::IMPLICIT::Generic(
          *this ) );

  NOX::NLN::AUX::AddToPrePostOpVector( p_sol_opt, prepost_solver_ptr );

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
const bool& STR::IMPLICIT::Generic::IsPredictorState() const
{
  return ispredictor_state_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::Generic::GetDefaultStepLength() const
{
  const Teuchos::ParameterList& p_nox = TimInt().GetDataSDyn().GetNoxParams();
  const std::string& nln_solver = p_nox.get<std::string>("Nonlinear Solver");
  // The pseudo transient implementation holds also a line search object!
  if (nln_solver=="Line Search Based" or
      nln_solver=="Pseudo Transient")
  {
    const Teuchos::ParameterList& p_ls = p_nox.sublist("Line Search");
    const std::string& method = p_ls.get<std::string>("Method");
    const Teuchos::ParameterList& p_method = p_ls.sublist(method);
    if (p_method.isParameter("Default Step"))
      return p_method.get<double>("Default Step");
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
void NOX::NLN::GROUP::PrePostOp::IMPLICIT::Generic::runPreComputeX(
    const NOX::NLN::Group& input_grp,
    const Epetra_Vector& dir,
    const double& step,
    const NOX::NLN::Group& curr_grp)
{
  // set the evaluation parameters
  const Epetra_Vector& xold =
      dynamic_cast<const NOX::Epetra::Vector&>(input_grp.getX()).getEpetraVector();
  Epetra_Vector& dir_mutable = const_cast<Epetra_Vector&>( dir );

  const bool isdefaultstep = (step == default_step_);
  impl_.ModelEval().RunPreComputeX( xold, dir_mutable, step, curr_grp,
      isdefaultstep );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::GROUP::PrePostOp::IMPLICIT::Generic::runPostComputeX(
    const NOX::NLN::Group& input_grp,
    const Epetra_Vector& dir,
    const double& step,
    const NOX::NLN::Group& curr_grp)
{
  // set the evaluation parameters
  const Epetra_Vector& xold =
      dynamic_cast<const NOX::Epetra::Vector&>(input_grp.getX()).getEpetraVector();
  const Epetra_Vector& xnew =
      dynamic_cast<const NOX::Epetra::Vector&>(curr_grp.getX()).getEpetraVector();

  bool isdefaultstep = (step == default_step_);
  impl_.ModelEval().RecoverState( xold, dir, step, xnew, isdefaultstep );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PrePostOp::IMPLICIT::Generic::runPostIterate(
    const NOX::Solver::Generic& solver )
{
  // try to cast the given solver object
  const NOX::NLN::Solver::LineSearchBased* ls_solver =
      dynamic_cast<const NOX::NLN::Solver::LineSearchBased*>( &solver );

  double step = 0.0;
  bool isdefaultstep = false;

  if ( not ls_solver )
  {
    step = default_step_;
    isdefaultstep = true;
  }
  else
  {
    step = ls_solver->getStepSize();
    isdefaultstep = ( step == default_step_ );
  }

  impl_.ModelEval().RunPostIterate( solver, step, isdefaultstep );
}
