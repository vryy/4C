/*!----------------------------------------------------------------------
\file red_airways_dyn_drt.cpp
\brief Main control routine for reduced dimensional airways network
 (time) integration.

\maintainer Lena Yoshihara

\level 3

*----------------------------------------------------------------------*/

#include "red_airway_resulttest.H"
#include "red_airways_dyn_drt.H"
#include "airwayimplicitintegration.H"
#include "redairway_tissue.H"
#include "../drt_adapter/ad_str_redairway.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_inpar/drt_validparameters.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <ctime>
#include <cstdlib>
#include <iostream>


/*----------------------------------------------------------------------*
 * Main control routine for reduced dimensional airway network including|
 * various solvers                                                      |
 *                                                                      ||
 *----------------------------------------------------------------------*/
void dyn_red_airways_drt()
{
  dyn_red_airways_drt(false);
}

Teuchos::RCP<AIRWAY::RedAirwayImplicitTimeInt>  dyn_red_airways_drt(bool CoupledTo3D)
{
  if(DRT::Problem::Instance()->DoesExistDis("red_airway")==false)
  {
    return Teuchos::null;
  }

  // 1. Access the discretization
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis("red_airway");

  //Set degrees of freedom in the discretization
  if (!actdis->Filled())
  {
    actdis->FillComplete();
  }

  //If discretization is empty, then return empty time integration
  if (actdis->NumGlobalElements()<1)
  {
    return Teuchos::null;
  }

  // 2. Context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter>  output = actdis->Writer();
  output->WriteMesh(0,0.0);

  // 3. Set pointers and variables for ParameterList rawdyn
  const Teuchos::ParameterList& rawdyn = DRT::Problem::Instance()->ReducedDAirwayDynamicParams();

  //Print default parameters
  if (actdis->Comm().MyPID()==0)
  {
    DRT::INPUT::PrintDefaultParameters(IO::cout, rawdyn);
  }

  // 4. Create a linear solver
  //Get the solver number for the LINEAR_SOLVER
  const int linsolvernumber = rawdyn.get<int>("LINEAR_SOLVER");
  //Check if the present solver has a valid solver number
  if (linsolvernumber == (-1))
  {
    dserror("no linear solver defined. Please set LINEAR_SOLVER in REDUCED DIMENSIONAL AIRWAYS DYNAMIC to a valid number!");
  }
  //Create the solver
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp( new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                                                   actdis->Comm(),
                                                   DRT::Problem::Instance()->ErrorFile()->Handle()),
                                                   false);
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // 5. Set parameters in list required for all schemes
  Teuchos::ParameterList airwaystimeparams;

  //Number of degrees of freedom
  const int ndim = DRT::Problem::Instance()->NDim();
  airwaystimeparams.set<int>              ("number of degrees of freedom" ,1*ndim);

  //Time step size
  airwaystimeparams.set<double>           ("time step size"           ,rawdyn.get<double>("TIMESTEP"));
  //Maximum number of timesteps
  airwaystimeparams.set<int>              ("max number timesteps"     ,rawdyn.get<int>("NUMSTEP"));

  //Restart
  airwaystimeparams.set                  ("write restart every"       ,rawdyn.get<int>("RESTARTEVRY"));
  //Solution output
  airwaystimeparams.set                  ("write solution every"      ,rawdyn.get<int>("RESULTSEVRY"));

  //Solver type
  airwaystimeparams.set                  ("solver type"             ,rawdyn.get<std::string>("SOLVERTYPE"));
  //Tolerance
  airwaystimeparams.set                  ("tolerance"               ,rawdyn.get<double>("TOLERANCE"));
  //Maximum number of iterations
  airwaystimeparams.set                  ("maximum iteration steps" ,rawdyn.get<int>("MAXITERATIONS"));
  //SolveScatra
  if (rawdyn.get<std::string>("SOLVESCATRA")=="yes")
    airwaystimeparams.set                  ("SolveScatra" ,true);
  else
    airwaystimeparams.set                  ("SolveScatra" ,false);
  // compute Interdependency
  if (rawdyn.get<std::string>("COMPAWACINTER")=="yes")
  airwaystimeparams.set ("CompAwAcInter" ,true);
  else
  airwaystimeparams.set ("CompAwAcInter" ,false);

  //Adjust acini volume with pre-stress condition
  if (rawdyn.get<std::string>("CALCV0PRESTRESS")=="yes")
  {
    airwaystimeparams.set                  ("CalcV0PreStress" ,true);
    airwaystimeparams.set                  ("transpulmpress"          ,rawdyn.get<double>("TRANSPULMPRESS"));
  }
  else
    airwaystimeparams.set                  ("CalcV0PreStress" ,false);



  // 6. Create all vectors and variables associated with the time
  // integration (call the constructor);
  // the only parameter from the list required here is the number of
  // velocity degrees of freedom
  Teuchos::RCP<AIRWAY::RedAirwayImplicitTimeInt> airwayimplicit =
    Teuchos::rcp(new AIRWAY::RedAirwayImplicitTimeInt(actdis,*solver,airwaystimeparams,*output));

  // Initialize state save vectors
  if (CoupledTo3D)
  {
    airwayimplicit->InitSaveState();
  }

  //Initial field from restart or calculated by given function
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart && !CoupledTo3D)
  {
    //Read the restart information, set vectors and variables
    airwayimplicit->ReadRestart(restart);
  }

  //Handle errors
  airwaystimeparams.set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

  if (!CoupledTo3D)
  {
    //Call time-integration scheme for 0D problem
    Teuchos::RCP<Teuchos::ParameterList> param_temp;
    airwayimplicit->Integrate();

    //Create resulttest
    Teuchos::RCP<DRT::ResultTest> resulttest
      = Teuchos::rcp(new AIRWAY::RedAirwayResultTest(*airwayimplicit));

    //Resulttest for 0D problem and testing
    DRT::Problem::Instance()->AddFieldTest(resulttest);
    DRT::Problem::Instance()->TestAll(actdis->Comm());

    return airwayimplicit;
  }
  else
  {
    return airwayimplicit;
  }

} // end of dyn_red_airways_drt()


/*----------------------------------------------------------------------*
 | dyn routine for redairway_tissue coupling                 roth 10/13 |
 *----------------------------------------------------------------------*/
void redairway_tissue_dyn()
{
  const Teuchos::ParameterList& rawdyn   = DRT::Problem::Instance()->RedAirwayTissueDynamicParams();
  Teuchos::RCP<DRT::Discretization> actdis = DRT::Problem::Instance()->GetDis("structure");
  Teuchos::RCP<AIRWAY::RedAirwayTissue> redairway_tissue = Teuchos::rcp(new AIRWAY::RedAirwayTissue(actdis->Comm(),rawdyn));

  //Read the restart information, set vectors and variables
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
  {
    redairway_tissue->ReadRestart(restart);
  }

  //Time integration loop for red_airway-tissue coupling
  redairway_tissue->Integrate();

  //Print time monitor
  Teuchos::TimeMonitor::summarize(std::cout, false, true, false);

  //Resulttest for red_airway-tissue coupling
  //create result tests for single fields
  DRT::Problem::Instance()->AddFieldTest(redairway_tissue->RedAirwayField()->CreateFieldTest());
  DRT::Problem::Instance()->AddFieldTest(redairway_tissue->StructureField()->CreateFieldTest());

  //Do the actual testing
  DRT::Problem::Instance()->TestAll(actdis->Comm());

}
