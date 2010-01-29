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
 * Main control routine for arterial network including various solvers:
 *
 *        o
 *
 *----------------------------------------------------------------------*/

void dyn_red_airways_drt()
{

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;

  actdis = DRT::Problem::Instance()->Dis(genprob.numawf,0);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

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
  // flag for writing the hemodynamic physiological results
  //arterytimeparams.set ("write stresses"  ,Teuchos::getIntegralValue<int>(ioflags,"HEMO_PHYS_RESULTS"));
  //---------------------- A method to initialize the flow inside the 
  //                       arteries.
  //  int init = Teuchos::getIntegralValue<int> (artdyn,"INITIALFIELD");

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

  if (probtype.get<int>("RESTART"))
  {
    // read the restart information, set vectors and variables
    airwayimplicit->ReadRestart(probtype.get<int>("RESTART"));
  }
  else
  {

  }


  airwaystimeparams.set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

  RCP<ParameterList> param_temp;
  airwayimplicit->Integrate();

} // end of dyn_art_net_drt()

#endif // #ifdef CCADISCRET
