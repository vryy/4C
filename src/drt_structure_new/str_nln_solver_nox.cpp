/*-----------------------------------------------------------*/
/*!
\file str_nln_solver_nox.cpp

\brief Structural non-linear %NOX::NLN solver.

\maintainer Michael Hiermeier

\date Sep 15, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_nln_solver_nox.H"       // class header
#include "str_timint_noxinterface.H"
#include "str_timint_base.H"
#include "str_utils.H"

#include "../solver_nonlin_nox/nox_nln_problem.H"
#include "../solver_nonlin_nox/nox_nln_constraint_interface_required.H"
#include "../solver_nonlin_nox/nox_nln_solver_factory.H"
#include "../solver_nonlin_nox/nox_nln_globaldata.H"
#include "../solver_nonlin_nox/nox_nln_linearsystem.H"

#include "../linalg/linalg_solver.H"

#include <NOX_Epetra_LinearSystem.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Abstract_Group.H>
#include <NOX_Solver_Generic.H>
#include <NOX_Epetra_Scaling.H>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::NLN::SOLVER::Nox::Nox()
    : nlnglobaldata_(Teuchos::null),
      problem_(Teuchos::null),
      linsys_(Teuchos::null),
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

  /* Set NOX::Epetra::Interface::Required
   * This interface is necessary for the evaluation of basic things
   * which are evaluated outside of the non-linear solver, but
   * are always necessary. A simple example is the right-hand-side
   * F. (see computeF) */
  const Teuchos::RCP<NOX::Epetra::Interface::Required> ireq =
      NoxInterfacePtr();

  /* Set NOX::Epetra::Interface::Jacobian
   * This interface is necessary for the evaluation of the jacobian
   * and everything, which is directly related to the jacobian.
   * This interface is optional. You can think of Finite-Differences
   * as one way to circumvent the evaluation of the jacobian.
   * Nevertheless, we always set this interface ptr in the structural
   * case. */
  const Teuchos::RCP<NOX::Epetra::Interface::Jacobian> ijac =
      NoxInterfacePtr();

  // pre-conditioner interface
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iprec =
      NoxInterfacePtr();

  // vector of currently present solution types
  std::vector<enum NOX::NLN::SolutionType> soltypes(0);
  // map of linear solvers, the key is the solution type
  NOX::NLN::LinearSystem::SolverMap linsolvers;
  /* convert the INPAR::STR::ModelType to a NOX::NLN::SolType
   * and fill the linear solver map. */
  STR::NLN::ConvertModelType2SolType(soltypes,linsolvers,
      DataSDyn().GetModelTypes(),DataSDyn().GetLinSolvers());

  // define and initialize the optimization type
  const NOX::NLN::OptimizationProblemType opttype =
      STR::NLN::OptimizationType(soltypes);

  // map of constraint interfaces, the key is the solution type
  NOX::NLN::CONSTRAINT::ReqInterfaceMap iconstr;
  // set constraint interfaces
  STR::NLN::CreateConstraintInterfaces(iconstr,Integrator(),soltypes);

  // preconditioner map for constraint problems
  NOX::NLN::CONSTRAINT::PrecInterfaceMap iconstr_prec;
  STR::NLN::CreateConstraintPreconditioner(iconstr_prec,
      Integrator(),soltypes);

  // create object to scale linear system
  Teuchos::RCP<NOX::Epetra::Scaling> iscale = Teuchos::null;
  STR::NLN::CreateScaling(iscale, DataSDyn(), DataGlobalState());

  // build the global data container for the nox_nln_solver
  nlnglobaldata_ =
      Teuchos::rcp(new NOX::NLN::GlobalData(
          DataGlobalState().GetComm(),
          DataSDyn().GetMutableNoxParams(),
          linsolvers,
          ireq,
          ijac,
          opttype,
          iconstr,
          iprec,
          iconstr_prec,
          iscale));

  // -------------------------------------------------------------------------
  // Create NOX control class: NoxProblem()
  // -------------------------------------------------------------------------
  Teuchos::RCP<NOX::Epetra::Vector> soln =
      DataGlobalState().CreateGlobalVector();
  Teuchos::RCP<LINALG::SparseOperator>& jac =
      DataGlobalState().CreateJacobian();
  problem_ = Teuchos::rcp(new NOX::NLN::Problem(nlnglobaldata_,soln,jac));

  // -------------------------------------------------------------------------
  // Create NOX linear system to provide access to Jacobian etc.
  // -------------------------------------------------------------------------
  linsys_ = problem_->CreateLinearSystem();

  // -------------------------------------------------------------------------
  // Create NOX group
  // -------------------------------------------------------------------------
  /* use NOX::NLN::Group to enable access to time integration
   * use NOX::NLN::Constraint::Group to enable access to the constraint data*/
  GroupPtr() = problem_->CreateGroup(linsys_);

  // -------------------------------------------------------------------------
  // Create NOX status test
  // -------------------------------------------------------------------------
  // get the stopping criteria from the nox parameter list
  problem_->CreateStatusTests(ostatus_,istatus_);

  // set flag
  issetup_ = true;

  return;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::Nox::Reset()
{
  // safety check
  CheckInitSetup();

  // do a hard reset at the beginning to be on the safe side
  nlnsolver_ = Teuchos::null;

  // -------------------------------------------------------------------------
  // reset the parameter list
  // -------------------------------------------------------------------------
  ResetParams();

  // -------------------------------------------------------------------------
  // Create NOX non-linear solver
  // -------------------------------------------------------------------------
  nlnsolver_ =
     NOX::NLN::Solver::BuildSolver(GroupPtr(),ostatus_,istatus_,nlnglobaldata_);

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

  // copy the solution group into the class variable
  Group() = nlnsolver_->getSolutionGroup();

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

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::NLN::SOLVER::Nox::GetNumNlnIterations() const
{
  if (not nlnsolver_.is_null())
    return nlnsolver_->getNumIterations();
  return 0;
}
