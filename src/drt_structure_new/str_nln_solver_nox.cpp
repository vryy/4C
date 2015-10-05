/*
 * str_nln_solver_nox.cpp
 *
 *  Created on: Sep 15, 2015
 *      Author: hiermeier
 */

#include "str_nln_solver_nox.H"       // class header
#include "str_timint_implicit.H"

#include "../solver_nonlin_nox/nox_nln_problem.H"
#include "../solver_nonlin_nox/nox_nln_problem_factory.H"
#include "../solver_nonlin_nox/nox_nln_constraint_interface_required.H"
#include "../solver_nonlin_nox/nox_nln_solver_factory.H"
#include "../linalg/linalg_solver.H"

#include <NOX_Epetra_LinearSystem.H>
#include <NOX_Abstract_Group.H>
#include <NOX_Solver_Generic.H>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::NLN::SOLVER::Nox::Nox()
    : soltypes_(std::vector<enum NOX::NLN::SolutionType>(0)),
      nlnglobaldata_(Teuchos::null),
      problem_(Teuchos::null),
      linsys_(Teuchos::null),
      group_(Teuchos::null),
      ostatus_(Teuchos::null),
      istatus_(Teuchos::null),
      nlnsolver_(Teuchos::null)
{
  //empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::Nox::Setup()
{
  CheckInit();

  // define and initialize the optimization type
  const NOX::NLN::GlobalData::OptimizationProblemType opttype =
      OptimizationType();

  /* Set NOX::Epetra::Interface::Required
   * This interface is necessary for the evaluation of basic things
   * which are evaluated outside of the non-linear solver, but
   * are always necessary. A simple example is the right-hand-side
   * F. (see computeF)
   */
  const Teuchos::RCP<NOX::Epetra::Interface::Required> ireq =
      ImplicitTimIntPtr();

  /* Set NOX::Epetra::Interface::Jacobian
   * This interface is necessary for the evaluation of the jacobian
   * and everything, which is directly related to the jacobian.
   * This interface is optional. You can think of Finite-Differences
   * as one way to circumvent the usage of a evaluation of the jacobian.
   * Nevertheless, we always set this interface ptr in the structural
   * case.
   */
  const Teuchos::RCP<NOX::Epetra::Interface::Jacobian> ijac =
      ImplicitTimIntPtr();

  // convert the INPAR::STR::ModelType to a NOX::NLN::SolType
  ConvertModelType2SolType();

  // set constraint interfaces
  SetConstraintInterfaces();

  // pre-conditioner interface
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iprec = Teuchos::null;

  // build the global data container for the nox_nln_solver
  nlnglobaldata_ =
      Teuchos::rcp(new NOX::NLN::GlobalData(
          DataGlobalState().GetComm(),
          DataSDyn().GetMutableNoxParams(),
          linsolvers_,
          ireq,
          ijac,
          opttype,
          iconstr_,
          iprec));

  // set flag
  issetup_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const enum NOX::NLN::GlobalData::OptimizationProblemType
    STR::NLN::SOLVER::Nox::OptimizationType() const
{
  CheckInit();

  enum NOX::NLN::GlobalData::OptimizationProblemType opttype =
      NOX::NLN::GlobalData::opt_unconstrained;
  // get a const reference to the model types
  const std::set<enum INPAR::STR::ModelType>& modeltypes =
      DataSDyn().GetModelTypes();
  std::set<enum INPAR::STR::ModelType>::const_iterator mt_iter;

  for (mt_iter=modeltypes.begin();mt_iter!=modeltypes.end();++mt_iter)
  {
    switch (*mt_iter)
    {
      // -----------------------------------
      // Inequality constraint
      // active set decision necessary
      // saddle point structure or condensed
      // -----------------------------------
      case INPAR::STR::model_contact:
        return NOX::NLN::GlobalData::opt_inequality_constrained;
        break;
      // -----------------------------------
      // Equality constraint
      // no active set decision necessary
      // saddle point structure or condensed
      // -----------------------------------
      case INPAR::STR::model_meshtying:
      case INPAR::STR::model_lag_pen_constraint:
      case INPAR::STR::model_windkessel:
        opttype = NOX::NLN::GlobalData::opt_equality_constrained;
        break;
      // -----------------------------------
      // Unconstrained problem
      // pure structural problem
      // no saddle point structure
      // -----------------------------------
      case INPAR::STR::model_structure:
      case INPAR::STR::model_springdashpot:
      default:
        break;
    }
  }

  return opttype;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::Nox::ConvertModelType2SolType()
{
  CheckInit();

  // initialize the vector and/or force the length to zero
  if (soltypes_.size()>0)
  {
    soltypes_.clear();
    linsolvers_.clear();
  }

  // get reference to the structural model-types
  const std::set<enum INPAR::STR::ModelType>& modeltypes =
      DataSDyn().GetModelTypes();
  // pre-set the vector size
  soltypes_.reserve(modeltypes.size());

  // The strings of the different enums have to fit!
  std::set<enum INPAR::STR::ModelType>::const_iterator mt_iter;
  for (mt_iter=modeltypes.begin();mt_iter!=modeltypes.end();++mt_iter)
  {
    const std::string name = INPAR::STR::ModelTypeString(*mt_iter);
    const enum NOX::NLN::SolutionType soltype =
        NOX::NLN::String2SolutionType(name);

    // check if the corresponding enum could be found.
    if (soltype == NOX::NLN::sol_unknown)
      dserror("The corresponding solution-type was not found. "
          "Given string: %s", name.c_str());

    soltypes_.push_back(soltype);
    // copy the linsolver pointers into the new map
    linsolvers_[soltype] = DataSDyn().GetMutableLinSolvers()[*mt_iter];
  }

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::Nox::SetConstraintInterfaces()
{
  CheckInit();

  if (iconstr_.size()>0)
    iconstr_.clear();

  std::vector<enum NOX::NLN::SolutionType>::const_iterator st_iter;
  for (st_iter=soltypes_.begin();st_iter!=soltypes_.end();++st_iter)
  {
    switch (*st_iter)
    {
      case NOX::NLN::sol_contact:
//        iconstr_[NOX::NLN::sol_contact] = ImplicitTimInt().GetContactManager();
//        break;
      case NOX::NLN::sol_windkessel:
//        iconstr_[NOX::NLN::sol_windkessel] = ImplicitTimInt().GetWindkesselManager();
//        break;
      case NOX::NLN::sol_lag_pen_constraint:
//        iconstr_[NOX::NLN::sol_windkessel] = ImplicitTimInt().GetLagPenConstrManager();
//        break;
      case NOX::NLN::sol_springdashpot:
//        iconstr_[NOX::NLN::sol_windkessel] = ImplicitTimInt().GetWindkesselManager();
        dserror("Constraint interfaces are not yet considered!");
        break;
      default:
        break;
    }
  }

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::Nox::Reset()
{
  // safety check
  CheckInitSetup();

  // do a hard reset at the beginning to be on the safe side
  problem_   = Teuchos::null;
  linsys_    = Teuchos::null;
  group_     = Teuchos::null;
  ostatus_   = Teuchos::null;
  istatus_   = Teuchos::null;
  nlnsolver_ = Teuchos::null;

  // -------------------------------------------------------------------------
  // reset the parameter list
  // -------------------------------------------------------------------------
  ResetParams();

  // -------------------------------------------------------------------------
  // Create NOX control class: NoxProblem()
  // -------------------------------------------------------------------------
  // FixMe add the solution vector and the complete jacobian here!
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::null;
  Teuchos::RCP<LINALG::SparseOperator> jac = Teuchos::null;
  problem_ = NOX::NLN::Problem::BuildNoxProblem(rhs,jac,nlnglobaldata_);

  // -------------------------------------------------------------------------
  // Create NOX linear system to provide access to Jacobian etc.
  // -------------------------------------------------------------------------
  linsys_ = problem_->CreateLinearSystem();

  // -------------------------------------------------------------------------
  // Create NOX group
  // -------------------------------------------------------------------------
  /* use NOX::NLN::Group to enable access to time integration
   * use NOX::NLN::Constraint::Group to enable access to the constraint data
   */
  group_ = problem_->CreateNoxGroup(linsys_);

  // -------------------------------------------------------------------------
  // Create NOX status test
  // -------------------------------------------------------------------------
  // get the stopping criteria from the nox parameter list
  problem_->CreateStatusTests(ostatus_,istatus_);

  // -------------------------------------------------------------------------
  // Create NOX non-linear solver
  // -------------------------------------------------------------------------
  nlnsolver_ =
     NOX::NLN::Solver::BuildSolver(group_,ostatus_,istatus_,nlnglobaldata_);

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::Nox::ResetParams()
{
  // safety check
  CheckInitSetup();

  const std::string& method = nlnglobaldata_->GetNlnParameterList().
      sublist("Direction",true).get<std::string>("Method");

  if (method == "Newton")
  {
    // get the linear solver sub-sub-sub-list
    Teuchos::ParameterList& lsparams =
        nlnglobaldata_->GetNlnParameterList().sublist("Direction",true).
        sublist("Newton",true).sublist("Linear Solver",true);

    // get current time step and update the parameter list entry
    lsparams.set<int>("Current Time Step",DataGlobalState().GetStepNp());
  }

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::ConvergenceStatus STR::NLN::SOLVER::Nox::Solve()
{
  CheckInitSetup();

  // solve the non-linear step
  NOX::StatusTest::StatusType finalstatus = nlnsolver_->solve();

  // Check if we do something special if the non-linear solver fails,
  // otherwise an error is thrown.
  if (DataSDyn().GetDivergenceAction() == INPAR::STR::divcont_stop)
    problem_->CheckFinalStatus(finalstatus);

  return ConvertFinalStatus(finalstatus);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::ConvergenceStatus STR::NLN::SOLVER::Nox::ConvertFinalStatus(
    const NOX::StatusTest::StatusType& finalstatus) const
{
  CheckInitSetup();

  INPAR::STR::ConvergenceStatus convstatus =
      INPAR::STR::conv_success;

  switch (finalstatus)
  {
    case NOX::StatusTest::Unevaluated:
      convstatus = INPAR::STR::conv_ele_fail;
      break;
    case NOX::StatusTest::Unconverged:
    case NOX::StatusTest::Failed:
      convstatus = INPAR::STR::conv_nonlin_fail;
      break;
    case NOX::StatusTest::Converged:
      convstatus = INPAR::STR::conv_success;
      break;
    default:
      dserror("Conversion of the NOX::StatusTest::StatusType to "
          "a INPAR::STR::ConvergenceStatus is not possible!");
      break;
  }

  return convstatus;
}
