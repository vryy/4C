/*!----------------------------------------------------------------------
\file red_airways_dyn_drt.cpp
\brief Main control routine for reduced dimesional airways network
 (time) integration. Calls


<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/

#include <ctime>
#include <cstdlib>
#include <iostream>

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "red_airway_resulttest.H"
#include "red_airways_dyn_drt.H"
#include "airwayimplicitintegration.H"
#include "redairway_tissue.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_inpar/drt_validparameters.H"


/*----------------------------------------------------------------------*
 * Main control routine for reduced dimensional airway network including
 * various solvers:
 *
 *        o
 *
 *----------------------------------------------------------------------*/

void dyn_red_airways_drt()
{
  dyn_red_airways_drt(false);
}

Teuchos::RCP<AIRWAY::RedAirwayImplicitTimeInt>  dyn_red_airways_drt(bool CoupledTo3D)
{
#ifdef D_RED_AIRWAYS

#if 1

  if(DRT::Problem::Instance()->DoesExistDis("red_airway")==false)
  {
    return Teuchos::null;
  }
#endif
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->GetDis("red_airway");

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled())
  {
    actdis->FillComplete();
  }

  // -------------------------------------------------------------------
  // If discretization is empty, then return empty time integration
  // -------------------------------------------------------------------
  if (actdis->NumGlobalElements()<1)
  {
    return Teuchos::null;
  }

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RCP<IO::DiscretizationWriter>  output = rcp( new IO::DiscretizationWriter(actdis),false);
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  // const Teuchos::ParameterList& probsize = DRT::Problem::Instance()->ProblemSizeParams();
  //  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& rawdyn   = DRT::Problem::Instance()->ReducedDAirwayDynamicParams();

  if (actdis->Comm().MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, rawdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // get the solver number
  const int linsolvernumber = rawdyn.get<int>("LINEAR_SOLVER");
  // check if the present solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined. Please set LINEAR_SOLVER in REDUCED DIMENSIONAL AIRWAYS DYNAMIC to a valid number!");
  RCP<LINALG::Solver> solver = rcp( new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                                                   actdis->Comm(),
                                                   DRT::Problem::Instance()->ErrorFile()->Handle()),
                                    false);
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  ParameterList airwaystimeparams;

  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  const int ndim = DRT::Problem::Instance()->NDim();
  airwaystimeparams.set<int>              ("number of degrees of freedom" ,1*ndim);

  // -------------------------------------------------- time integration
  // the default time step size
  airwaystimeparams.set<double>           ("time step size"           ,rawdyn.get<double>("TIMESTEP"));
  // maximum number of timesteps
  airwaystimeparams.set<int>              ("max number timesteps"     ,rawdyn.get<int>("NUMSTEP"));

  // ----------------------------------------------- restart and output
  // restart
  airwaystimeparams.set                  ("write restart every"       ,rawdyn.get<int>("RESTARTEVRY"));
  // solution output
  airwaystimeparams.set                  ("write solution every"      ,rawdyn.get<int>("UPRES"));

  // ----------------------------------------------- solver parameters
  // solver type
  airwaystimeparams.set                  ("solver type"             ,rawdyn.get<string>("SOLVERTYPE"));
  // tolerance
  airwaystimeparams.set                  ("tolerance"               ,rawdyn.get<double>("TOLERANCE"));
  // maximum number of iterations
  airwaystimeparams.set                  ("maximum iteration steps" ,rawdyn.get<int>("MAXITERATIONS"));

  //------------------------------------------------------------------
  // create all vectors and variables associated with the time
  // integration (call the constructor);
  // the only parameter from the list required here is the number of
  // velocity degrees of freedom
  //------------------------------------------------------------------


  Teuchos::RCP<AIRWAY::RedAirwayImplicitTimeInt> airwayimplicit
    =
    Teuchos::rcp(new AIRWAY::RedAirwayImplicitTimeInt(actdis,*solver,airwaystimeparams,*output));
  // initial field from restart or calculated by given function

  const int restart = DRT::Problem::Instance()->Restart();
  if (restart && !CoupledTo3D)
  {
    // read the restart information, set vectors and variables
    airwayimplicit->ReadRestart(restart);
  }
  else
  {

  }

  airwaystimeparams.set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

  if (!CoupledTo3D)
  {
    // call time-integration (or stationary) scheme
    RCP<ParameterList> param_temp;
    airwayimplicit->Integrate();

    Teuchos::RCP<DRT::ResultTest> resulttest
      = Teuchos::rcp(new AIRWAY::RedAirwayResultTest(*airwayimplicit));
    DRT::Problem::Instance()->AddFieldTest(resulttest);
    DRT::Problem::Instance()->TestAll(actdis->Comm());

    return airwayimplicit;
    //    return  Teuchos::null;
  }
  else
  {
    return airwayimplicit;
  }

#else
  return Teuchos::null;
#endif
} // end of dyn_red_airways_drt()


void redairway_tissue_dyn()
{
  const Teuchos::ParameterList& rawdyn   = DRT::Problem::Instance()->RedAirwayTissueDynamicParams();
  RefCountPtr<DRT::Discretization> actdis = DRT::Problem::Instance()->GetDis("structure");
  Teuchos::RCP<AIRWAY::RedAirwayTissue> myalgo = Teuchos::rcp(new AIRWAY::RedAirwayTissue(actdis->Comm(),rawdyn));

  myalgo->Integrate();
}

