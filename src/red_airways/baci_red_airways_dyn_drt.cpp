/*---------------------------------------------------------------------*/
/*! \file

\brief Main control routine for reduced dimensional airways network
 (time) integration. Includes routines for pure reduced modeling and for
 reduced_airway-tissue coupling.


\level 3

*/
/*---------------------------------------------------------------------*/

#include "baci_red_airways_dyn_drt.H"

#include "baci_adapter_str_redairway.H"
#include "baci_global_data.H"
#include "baci_inpar_validparameters.H"
#include "baci_io_control.H"
#include "baci_io_pstream.H"
#include "baci_lib_resulttest.H"
#include "baci_red_airways_implicitintegration.H"
#include "baci_red_airways_resulttest.H"
#include "baci_red_airways_tissue.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <cstdlib>
#include <ctime>
#include <iostream>

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 * Main control routine for reduced dimensional airway network including|
 * various solvers                                                      |
 *                                                                      ||
 *----------------------------------------------------------------------*/
void dyn_red_airways_drt() { dyn_red_airways_drt(false); }

Teuchos::RCP<AIRWAY::RedAirwayImplicitTimeInt> dyn_red_airways_drt(bool CoupledTo3D)
{
  if (GLOBAL::Problem::Instance()->DoesExistDis("red_airway") == false)
  {
    return Teuchos::null;
  }

  // 1. Access the discretization
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = GLOBAL::Problem::Instance()->GetDis("red_airway");

  // Set degrees of freedom in the discretization
  if (!actdis->Filled())
  {
    actdis->FillComplete();
  }

  // If discretization is empty, then return empty time integration
  if (actdis->NumGlobalElements() < 1)
  {
    return Teuchos::null;
  }

  // 2. Context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter> output = actdis->Writer();
  output->WriteMesh(0, 0.0);

  // 3. Set pointers and variables for ParameterList rawdyn
  const Teuchos::ParameterList& rawdyn = GLOBAL::Problem::Instance()->ReducedDAirwayDynamicParams();

  // Print default parameters
  if (actdis->Comm().MyPID() == 0)
  {
    INPUT::PrintDefaultParameters(IO::cout, rawdyn);
  }

  // 4. Create a linear solver
  // Get the solver number for the LINEAR_SOLVER
  const int linsolvernumber = rawdyn.get<int>("LINEAR_SOLVER");
  // Check if the present solver has a valid solver number
  if (linsolvernumber == (-1))
  {
    dserror(
        "no linear solver defined. Please set LINEAR_SOLVER in REDUCED DIMENSIONAL AIRWAYS DYNAMIC "
        "to a valid number!");
  }
  // Create the solver
  std::unique_ptr<CORE::LINALG::Solver> solver = std::make_unique<CORE::LINALG::Solver>(
      GLOBAL::Problem::Instance()->SolverParams(linsolvernumber), actdis->Comm());
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // 5. Set parameters in list required for all schemes
  Teuchos::ParameterList airwaystimeparams;

  // Number of degrees of freedom
  const int ndim = GLOBAL::Problem::Instance()->NDim();
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
  Teuchos::RCP<AIRWAY::RedAirwayImplicitTimeInt> airwayimplicit = Teuchos::rcp(
      new AIRWAY::RedAirwayImplicitTimeInt(actdis, std::move(solver), airwaystimeparams, *output));

  // Initialize state save vectors
  if (CoupledTo3D)
  {
    airwayimplicit->InitSaveState();
  }

  // Initial field from restart or calculated by given function
  const int restart = GLOBAL::Problem::Instance()->Restart();
  if (restart && !CoupledTo3D)
  {
    // Read the restart information, set vectors and variables
    airwayimplicit->ReadRestart(restart);
  }

  if (!CoupledTo3D)
  {
    // Call time-integration scheme for 0D problem
    Teuchos::RCP<Teuchos::ParameterList> param_temp;
    airwayimplicit->Integrate();

    // Create resulttest
    Teuchos::RCP<DRT::ResultTest> resulttest =
        Teuchos::rcp(new AIRWAY::RedAirwayResultTest(*airwayimplicit));

    // Resulttest for 0D problem and testing
    GLOBAL::Problem::Instance()->AddFieldTest(resulttest);
    GLOBAL::Problem::Instance()->TestAll(actdis->Comm());

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
      GLOBAL::Problem::Instance()->RedAirwayTissueDynamicParams();
  Teuchos::RCP<DRT::Discretization> actdis = GLOBAL::Problem::Instance()->GetDis("structure");
  Teuchos::RCP<AIRWAY::RedAirwayTissue> redairway_tissue =
      Teuchos::rcp(new AIRWAY::RedAirwayTissue(actdis->Comm(), rawdyn));

  // Read the restart information, set vectors and variables
  const int restart = GLOBAL::Problem::Instance()->Restart();
  if (restart)
  {
    redairway_tissue->ReadRestart(restart);
  }

  // Time integration loop for red_airway-tissue coupling
  redairway_tissue->Integrate();

  // Print time monitor
  Teuchos::TimeMonitor::summarize(std::cout, false, true, false);

  // Resulttest for red_airway-tissue coupling
  // create result tests for single fields
  GLOBAL::Problem::Instance()->AddFieldTest(redairway_tissue->RedAirwayField()->CreateFieldTest());
  GLOBAL::Problem::Instance()->AddFieldTest(redairway_tissue->StructureField()->CreateFieldTest());

  // Do the actual testing
  GLOBAL::Problem::Instance()->TestAll(actdis->Comm());
}

BACI_NAMESPACE_CLOSE
