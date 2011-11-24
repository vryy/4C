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
#ifdef CCADISCRET

#include <ctime>
#include <cstdlib>
#include <iostream>

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "red_airway_resulttest.H"
#include "red_airways_dyn_drt.H"
#include "airwayimplicitintegration.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_inpar/drt_validparameters.H"

/*----------------------------------------------------------------------*
  |                                                       ismail 01/10   |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

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
  // -------------------------------------------------------------------
  // check if descretization exits
  // -------------------------------------------------------------------
  if((DRT::Problem::Instance()->NumFields() >= (unsigned) genprob.numartf) && CoupledTo3D)
  {
    return Teuchos::null;
  }
  if(DRT::Problem::Instance()->NumDis(genprob.numawf)<1)
  {
#if 0
    if (actdis->Comm().MyPID()==0)
    {
      cout<<"+--------------------- WARNING ---------------------+"<<endl;
      cout<<"|                                                   |"<<endl;
      cout<<"| Reduced-dimensional airways is compiled, but no   |"<<endl;
      cout<<"| airways elements are defined!                     |"<<endl;
      cout<<"|                                                   |"<<endl;
      cout<<"+---------------------------------------------------+"<<endl;
    }
#endif
    return Teuchos::null;
  }

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numawf,0);

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
  const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& probsize = DRT::Problem::Instance()->ProblemSizeParams();
  //  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& rawdyn   = DRT::Problem::Instance()->ReducedDAirwayDynamicParams();

  if (actdis->Comm().MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, rawdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RCP<LINALG::Solver> solver = rcp( new LINALG::Solver(DRT::Problem::Instance()->ReducedDAirwaySolverParams(),
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
  airwaystimeparams.set<int>              ("number of degrees of freedom" ,1*probsize.get<int>("DIM"));

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

  if (probtype.get<int>("RESTART") && !CoupledTo3D)
  {
    // read the restart information, set vectors and variables
    airwayimplicit->ReadRestart(probtype.get<int>("RESTART"));
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

#endif // #ifdef CCADISCRET
