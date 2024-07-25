/*---------------------------------------------------------------------*/
/*! \file

\brief Main control routine for reduced dimensional airways network
 (time) integration. Includes routines for pure reduced modeling and for
 reduced_airway-tissue coupling.


\level 3

*/
/*---------------------------------------------------------------------*/

#include "4C_red_airways_dyn_drt.hpp"

#include "4C_adapter_str_redairway.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io_control.hpp"
#include "4C_io_pstream.hpp"
#include "4C_red_airways_implicitintegration.hpp"
#include "4C_red_airways_resulttest.hpp"
#include "4C_red_airways_tissue.hpp"
#include "4C_utils_result_test.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <cstdlib>
#include <ctime>
#include <iostream>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 * Main control routine for reduced dimensional airway network including|
 * various solvers                                                      |
 *                                                                      ||
 *----------------------------------------------------------------------*/
void dyn_red_airways_drt() { dyn_red_airways_drt(false); }

Teuchos::RCP<Airway::RedAirwayImplicitTimeInt> dyn_red_airways_drt(bool CoupledTo3D)
{
  if (Global::Problem::instance()->does_exist_dis("red_airway") == false)
  {
    return Teuchos::null;
  }

  // 1. Access the discretization
  Teuchos::RCP<Core::FE::Discretization> actdis = Teuchos::null;
  actdis = Global::Problem::instance()->get_dis("red_airway");

  // Set degrees of freedom in the discretization
  if (!actdis->filled())
  {
    actdis->fill_complete();
  }

  // If discretization is empty, then return empty time integration
  if (actdis->num_global_elements() < 1)
  {
    return Teuchos::null;
  }

  // 2. Context for output and restart
  Teuchos::RCP<Core::IO::DiscretizationWriter> output = actdis->writer();
  output->write_mesh(0, 0.0);

  // 3. Set pointers and variables for ParameterList rawdyn
  const Teuchos::ParameterList& rawdyn =
      Global::Problem::instance()->reduced_d_airway_dynamic_params();

  // 4. Create a linear solver
  // Get the solver number for the LINEAR_SOLVER
  const int linsolvernumber = rawdyn.get<int>("LINEAR_SOLVER");
  // Check if the present solver has a valid solver number
  if (linsolvernumber == (-1))
  {
    FOUR_C_THROW(
        "no linear solver defined. Please set LINEAR_SOLVER in REDUCED DIMENSIONAL AIRWAYS DYNAMIC "
        "to a valid number!");
  }
  // Create the solver
  std::unique_ptr<Core::LinAlg::Solver> solver = std::make_unique<Core::LinAlg::Solver>(
      Global::Problem::instance()->solver_params(linsolvernumber), actdis->get_comm(),
      Global::Problem::instance()->solver_params_callback(),
      Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::instance()->io_params(), "VERBOSITY"));
  actdis->compute_null_space_if_necessary(solver->params());

  // 5. Set parameters in list required for all schemes
  Teuchos::ParameterList airwaystimeparams;

  // Number of degrees of freedom
  const int ndim = Global::Problem::instance()->n_dim();
  airwaystimeparams.set<int>("number of degrees of freedom", 1 * ndim);

  // Time step size
  airwaystimeparams.set<double>("time step size", rawdyn.get<double>("TIMESTEP"));
  // Maximum number of timesteps
  airwaystimeparams.set<int>("max number timesteps", rawdyn.get<int>("NUMSTEP"));

  // Restart
  airwaystimeparams.set("write restart every", rawdyn.get<int>("RESTARTEVRY"));
  // Solution output
  airwaystimeparams.set("write solution every", rawdyn.get<int>("RESULTSEVRY"));

  // Solver type
  airwaystimeparams.set("solver type", rawdyn.get<std::string>("SOLVERTYPE"));
  // Tolerance
  airwaystimeparams.set("tolerance", rawdyn.get<double>("TOLERANCE"));
  // Maximum number of iterations
  airwaystimeparams.set("maximum iteration steps", rawdyn.get<int>("MAXITERATIONS"));
  // SolveScatra
  if (rawdyn.get<std::string>("SOLVESCATRA") == "yes")
    airwaystimeparams.set("SolveScatra", true);
  else
    airwaystimeparams.set("SolveScatra", false);
  // compute Interdependency
  if (rawdyn.get<std::string>("COMPAWACINTER") == "yes")
    airwaystimeparams.set("CompAwAcInter", true);
  else
    airwaystimeparams.set("CompAwAcInter", false);

  // Adjust acini volume with pre-stress condition
  if (rawdyn.get<std::string>("CALCV0PRESTRESS") == "yes")
  {
    airwaystimeparams.set("CalcV0PreStress", true);
    airwaystimeparams.set("transpulmpress", rawdyn.get<double>("TRANSPULMPRESS"));
  }
  else
    airwaystimeparams.set("CalcV0PreStress", false);



  // 6. Create all vectors and variables associated with the time
  // integration (call the constructor);
  // the only parameter from the list required here is the number of
  // velocity degrees of freedom
  Teuchos::RCP<Airway::RedAirwayImplicitTimeInt> airwayimplicit = Teuchos::rcp(
      new Airway::RedAirwayImplicitTimeInt(actdis, std::move(solver), airwaystimeparams, *output));

  // Initialize state save vectors
  if (CoupledTo3D)
  {
    airwayimplicit->init_save_state();
  }

  // Initial field from restart or calculated by given function
  const int restart = Global::Problem::instance()->restart();
  if (restart && !CoupledTo3D)
  {
    // Read the restart information, set vectors and variables
    airwayimplicit->read_restart(restart);
  }

  if (!CoupledTo3D)
  {
    // Call time-integration scheme for 0D problem
    Teuchos::RCP<Teuchos::ParameterList> param_temp;
    airwayimplicit->integrate();

    // Create resulttest
    Teuchos::RCP<Core::UTILS::ResultTest> resulttest =
        Teuchos::rcp(new Airway::RedAirwayResultTest(*airwayimplicit));

    // Resulttest for 0D problem and testing
    Global::Problem::instance()->add_field_test(resulttest);
    Global::Problem::instance()->test_all(actdis->get_comm());

    return airwayimplicit;
  }
  else
  {
    return airwayimplicit;
  }

}  // end of dyn_red_airways_drt()


/*----------------------------------------------------------------------*
 | dyn routine for redairway_tissue coupling                 roth 10/13 |
 *----------------------------------------------------------------------*/
void redairway_tissue_dyn()
{
  const Teuchos::ParameterList& rawdyn =
      Global::Problem::instance()->red_airway_tissue_dynamic_params();
  Teuchos::RCP<Core::FE::Discretization> actdis = Global::Problem::instance()->get_dis("structure");
  Teuchos::RCP<Airway::RedAirwayTissue> redairway_tissue =
      Teuchos::rcp(new Airway::RedAirwayTissue(actdis->get_comm(), rawdyn));

  // Read the restart information, set vectors and variables
  const int restart = Global::Problem::instance()->restart();
  if (restart)
  {
    redairway_tissue->read_restart(restart);
  }

  // Time integration loop for red_airway-tissue coupling
  redairway_tissue->integrate();

  // Print time monitor
  Teuchos::TimeMonitor::summarize(std::cout, false, true, false);

  // Resulttest for red_airway-tissue coupling
  // create result tests for single fields
  Global::Problem::instance()->add_field_test(
      redairway_tissue->red_airway_field()->create_field_test());
  Global::Problem::instance()->add_field_test(
      redairway_tissue->structure_field()->create_field_test());

  // Do the actual testing
  Global::Problem::instance()->test_all(actdis->get_comm());
}

FOUR_C_NAMESPACE_CLOSE
