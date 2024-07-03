/*-----------------------------------------------------------*/
/*! \file

\brief Structural non-linear %NOX::NLN solver.


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_nln_solver_nox.hpp"  // class header

#include "4C_linear_solver_method_linalg.hpp"
#include "4C_solver_nonlin_nox_constraint_interface_required.hpp"
#include "4C_solver_nonlin_nox_globaldata.hpp"
#include "4C_solver_nonlin_nox_linearsystem.hpp"
#include "4C_solver_nonlin_nox_problem.hpp"
#include "4C_solver_nonlin_nox_solver_factory.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_timint_noxinterface.hpp"
#include "4C_structure_new_utils.hpp"

#include <NOX_Abstract_Group.H>
#include <NOX_Epetra_LinearSystem.H>
#include <NOX_Epetra_Scaling.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Solver_Generic.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Solid::Nln::SOLVER::Nox::Nox()
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::Nln::SOLVER::Nox::setup()
{
  check_init();

  /* Set ::NOX::Epetra::Interface::Required
   * This interface is necessary for the evaluation of basic things
   * which are evaluated outside of the non-linear solver, but
   * are always necessary. A simple example is the right-hand-side
   * F. (see computeF) */
  const Teuchos::RCP<::NOX::Epetra::Interface::Required> ireq = nox_interface_ptr();

  /* Set ::NOX::Epetra::Interface::Jacobian
   * This interface is necessary for the evaluation of the jacobian
   * and everything, which is directly related to the jacobian.
   * This interface is optional. You can think of Finite-Differences
   * as one way to circumvent the evaluation of the jacobian.
   * Nevertheless, we always set this interface ptr in the structural
   * case. */
  const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian> ijac = nox_interface_ptr();

  // pre-conditioner interface
  Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner> iprec = nox_interface_ptr();

  // vector of currently present solution types
  std::vector<enum NOX::Nln::SolutionType> soltypes(0);
  // map of linear solvers, the key is the solution type
  NOX::Nln::LinearSystem::SolverMap linsolvers;
  /* convert the Inpar::Solid::ModelType to a NOX::Nln::SolType
   * and fill the linear solver map. */
  Solid::Nln::ConvertModelType2SolType(
      soltypes, linsolvers, data_sdyn().get_model_types(), data_sdyn().GetLinSolvers());

  // define and initialize the optimization type
  const NOX::Nln::OptimizationProblemType opttype = Solid::Nln::OptimizationType(soltypes);

  // map of constraint interfaces, the key is the solution type
  NOX::Nln::CONSTRAINT::ReqInterfaceMap iconstr;
  // set constraint interfaces
  Solid::Nln::CreateConstraintInterfaces(iconstr, integrator(), soltypes);

  // preconditioner map for constraint problems
  NOX::Nln::CONSTRAINT::PrecInterfaceMap iconstr_prec;
  Solid::Nln::CreateConstraintPreconditioner(iconstr_prec, integrator(), soltypes);

  // create object to scale linear system
  Teuchos::RCP<::NOX::Epetra::Scaling> iscale = Teuchos::null;
  Solid::Nln::CreateScaling(iscale, data_sdyn(), data_global_state());

  // build the global data container for the nox_nln_solver
  nlnglobaldata_ = Teuchos::rcp(
      new NOX::Nln::GlobalData(data_global_state().get_comm(), data_sdyn().get_nox_params(),
          linsolvers, ireq, ijac, opttype, iconstr, iprec, iconstr_prec, iscale));

  // -------------------------------------------------------------------------
  // Create NOX control class: NoxProblem()
  // -------------------------------------------------------------------------
  Teuchos::RCP<::NOX::Epetra::Vector> soln = data_global_state().create_global_vector();
  Teuchos::RCP<Core::LinAlg::SparseOperator>& jac = data_global_state().create_jacobian();
  problem_ = Teuchos::rcp(new NOX::Nln::Problem(nlnglobaldata_, soln, jac));

  // -------------------------------------------------------------------------
  // Create NOX linear system to provide access to Jacobian etc.
  // -------------------------------------------------------------------------
  linsys_ = problem_->create_linear_system();

  // -------------------------------------------------------------------------
  // Create NOX group
  // -------------------------------------------------------------------------
  /* use NOX::Nln::Group to enable access to time integration
   * use NOX::Nln::Constraint::Group to enable access to the constraint data*/
  group_ptr() = problem_->CreateGroup(linsys_);

  // -------------------------------------------------------------------------
  // Create NOX status test
  // -------------------------------------------------------------------------
  // get the stopping criteria from the nox parameter list
  problem_->CreateStatusTests(ostatus_, istatus_);

  // set flag
  issetup_ = true;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::Nln::SOLVER::Nox::reset()
{
  // safety check
  check_init_setup();

  // do a hard reset at the beginning to be on the safe side
  nlnsolver_ = Teuchos::null;

  // -------------------------------------------------------------------------
  // reset the parameter list
  // -------------------------------------------------------------------------
  reset_params();

  // -------------------------------------------------------------------------
  // Create NOX non-linear solver
  // -------------------------------------------------------------------------
  nlnsolver_ = NOX::Nln::Solver::BuildSolver(group_ptr(), ostatus_, istatus_, nlnglobaldata_);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Solid::Nln::SOLVER::Nox::reset_params()
{
  // safety check
  check_init_setup();

  const Teuchos::ParameterList& pdir =
      nlnglobaldata_->GetNlnParameterList().sublist("Direction", true);
  std::string method = pdir.get<std::string>("Method");

  if (method == "User Defined") method = pdir.get<std::string>("User Defined Method");

  if (method == "Newton" or method == "Modified Newton")
  {
    // get the linear solver sub-sub-sub-list
    Teuchos::ParameterList& lsparams = nlnglobaldata_->GetNlnParameterList()
                                           .sublist("Direction", true)
                                           .sublist("Newton", true)
                                           .sublist("Linear Solver", true);

    // get current time step and update the parameter list entry
    lsparams.set<int>("Current Time Step", data_global_state().get_step_np());
  }
  else if (method == "Single Step")
  {
    // get the linear solver sub-sub-list
    Teuchos::ParameterList& lsparams = nlnglobaldata_->GetNlnParameterList()
                                           .sublist("Single Step Solver", true)
                                           .sublist("Linear Solver", true);

    // get current time step and update the parameter list entry
    lsparams.set<int>("Current Time Step", data_global_state().get_step_np());
  }
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum Inpar::Solid::ConvergenceStatus Solid::Nln::SOLVER::Nox::Solve()
{
  check_init_setup();

  // solve the non-linear step
  ::NOX::StatusTest::StatusType finalstatus = nlnsolver_->solve();

  // Check if we do something special if the non-linear solver fails,
  // otherwise an error is thrown.
  if (data_sdyn().get_divergence_action() == Inpar::Solid::divcont_stop)
    problem_->CheckFinalStatus(finalstatus);

  // copy the solution group into the class variable
  group() = nlnsolver_->getSolutionGroup();

  return convert_final_status(finalstatus);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum Inpar::Solid::ConvergenceStatus Solid::Nln::SOLVER::Nox::convert_final_status(
    const ::NOX::StatusTest::StatusType& finalstatus) const
{
  check_init_setup();

  Inpar::Solid::ConvergenceStatus convstatus = Inpar::Solid::conv_success;

  switch (finalstatus)
  {
    case ::NOX::StatusTest::Unevaluated:
      convstatus = Inpar::Solid::conv_ele_fail;
      break;
    case ::NOX::StatusTest::Unconverged:
    case ::NOX::StatusTest::Failed:
      convstatus = Inpar::Solid::conv_nonlin_fail;
      break;
    case ::NOX::StatusTest::Converged:
      convstatus = Inpar::Solid::conv_success;
      break;
    default:
      FOUR_C_THROW(
          "Conversion of the ::NOX::StatusTest::StatusType to "
          "a Inpar::Solid::ConvergenceStatus is not possible!");
      break;
  }

  return convstatus;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int Solid::Nln::SOLVER::Nox::get_num_nln_iterations() const
{
  if (not nlnsolver_.is_null()) return nlnsolver_->getNumIterations();
  return 0;
}

FOUR_C_NAMESPACE_CLOSE
