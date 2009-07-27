/*!----------------------------------------------------------------------
\file art_net_dyn_drt.cpp
\brief Main control routine for all arterial network  solvers,

     including instationary solvers based on

     o

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

#include "art_net_dyn_drt.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_inpar/drt_validparameters.H"
#include "artnetexplicitintegration.H"

/*----------------------------------------------------------------------*
  |                                                       ismail 01/09   |
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
void dyn_art_net_drt()
{

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numartf,0);
  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  IO::DiscretizationWriter output(actdis);
  output.WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& probsize = DRT::Problem::Instance()->ProblemSizeParams();
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& artdyn   = DRT::Problem::Instance()->ArterialDynamicParams();

  if (actdis->Comm().MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, artdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  LINALG::Solver solver(DRT::Problem::Instance()->ArteryNetworkSolverParams(),
                        actdis->Comm(),
                        DRT::Problem::Instance()->ErrorFile()->Handle());
  actdis->ComputeNullSpaceIfNecessary(solver.Params());

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  ParameterList arterytimeparams;

  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  arterytimeparams.set<int>              ("number of degrees of freedom" ,2*probsize.get<int>("DIM"));

  // -------------------------------------------------- time integration
  // the default time step size
  arterytimeparams.set<double>           ("time step size"           ,artdyn.get<double>("TIMESTEP"));
  // maximum number of timesteps
  arterytimeparams.set<int>              ("max number timesteps"     ,artdyn.get<int>("NUMSTEP"));

  // ----------------------------------------------- restart and output
  // restart
  arterytimeparams.set                  ("write restart every"       ,artdyn.get<int>("RESTARTEVRY"));
  // solution output
  arterytimeparams.set                  ("write solution every"      ,artdyn.get<int>("UPRES"));
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
 ART::ArtNetExplicitTimeInt artnetexplicit(actdis,
                                         solver,
                                         arterytimeparams,
                                         output);

 // initial field from restart or calculated by given function
 if (probtype.get<int>("RESTART"))
 {
   // read the restart information, set vectors and variables
   artnetexplicit.ReadRestart(probtype.get<int>("RESTART"));
 }
 else
 {
   // artnetexplicit.SetInitialData(init,startfuncno);
 }

 arterytimeparams.set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

 // call time-integration (or stationary) scheme
 artnetexplicit.Integrate();


} // end of dyn_art_net_drt()

#endif // #ifdef CCADISCRET
