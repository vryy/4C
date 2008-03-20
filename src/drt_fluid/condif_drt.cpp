/*!----------------------------------------------------------------------
\file condif_drt.cpp
\brief Main control routine for all (in)stationary convect.-diff. solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme 
       (with potential one-step-theta start algorithm)

     o generalized-alpha time-integration scheme

     and stationary solver.

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
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

#include "condif_drt.H"
#include "condifimplicitintegration.H"
#include "condif_genalpha_integration.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"


/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;



/*----------------------------------------------------------------------*
 * Main control routine for convection-diffusion incl. various solvers:
 *
 *        o instationary one-step-theta
 *        o instationary BDF2
 *        o instationary generalized-alpha
 *        o stationary
 *
 *----------------------------------------------------------------------*/
void dyn_condif_drt()
{

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(0,0);

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
  SOLVAR        *actsolv  = &solv[0];

  const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();

  if (actdis->Comm().MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, fdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  LINALG::Solver solver(solveparams,actdis->Comm(),allfiles.out_err);
  solver.TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  ParameterList condiftimeparams;

  // -----------------------------------------------------velocity field
  condiftimeparams.set<int>              ("condif velocity field"     ,Teuchos::getIntegralValue<int>(fdyn,"CD_VELOCITY"));

  // -------------------------------------------------- time integration
  // the default time step size
  condiftimeparams.set<double>           ("time step size"           ,fdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  condiftimeparams.set<double>           ("total time"               ,fdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  condiftimeparams.set<int>              ("max number timesteps"     ,fdyn.get<int>("NUMSTEP"));

  // ----------------------------------------------- restart and output
  // restart
  condiftimeparams.set                  ("write restart every"       ,fdyn.get<int>("RESTARTEVRY"));
  // solution output
  condiftimeparams.set                  ("write solution every"      ,fdyn.get<int>("UPRES"));

  // ---------------------------------(fine-scale) subgrid diffusivity?
  condiftimeparams.set<string>           ("fs subgrid diffusivity"   ,fdyn.get<string>("FSSUGRVISC"));


  // -------------------------------------------------------------------
  // additional parameters and algorithm call depending on respective
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  FLUID_TIMEINTTYPE iop = Teuchos::getIntegralValue<FLUID_TIMEINTTYPE>(fdyn,"TIMEINTEGR");
  if(iop == timeint_stationary or
     iop == timeint_one_step_theta or
     iop == timeint_bdf2
    )
  {
    // -----------------------------------------------------------------
    // set additional parameters in list for OST/BDF2/stationary scheme
    // -----------------------------------------------------------------
    // type of time-integration (or stationary) scheme
    condiftimeparams.set<FLUID_TIMEINTTYPE>("time int algo",iop);
    // parameter theta for time-integration schemes
    condiftimeparams.set<double>           ("theta"                    ,fdyn.get<double>("THETA"));
    // number of steps for potential start algorithm
    condiftimeparams.set<int>              ("number of start steps"    ,fdyn.get<int>("NUMSTASTEPS"));
    // parameter theta for potential start algorithm
    condiftimeparams.set<double>           ("start theta"              ,fdyn.get<double>("START_THETA"));

    //------------------------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor)
    //------------------------------------------------------------------
    CondifImplicitTimeInt condifimplicit(actdis,
                                         solver,
                                         condiftimeparams,
                                         output);

    // initial field from restart
    if (probtype.get<int>("RESTART"))
    {
      // read the restart information, set vectors and variables
      condifimplicit.ReadRestart(probtype.get<int>("RESTART"));
    }

    // call time-integration (or stationary) scheme
    condifimplicit.Integrate();

  }
  else if (iop == timeint_gen_alpha)
  {
    // -------------------------------------------------------------------
    // set additional parameters in list for generalized-alpha scheme
    // -------------------------------------------------------------------
    // parameter alpha_M for for generalized-alpha scheme
    condiftimeparams.set<double>           ("alpha_M"                  ,fdyn.get<double>("ALPHA_M"));
    // parameter alpha_F for for generalized-alpha scheme
    condiftimeparams.set<double>           ("alpha_F"                  ,fdyn.get<double>("ALPHA_F"));

    //------------------------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor)
    //------------------------------------------------------------------
    CondifGenAlphaIntegration genalphaint(actdis,
                                          solver,
                                          condiftimeparams,
                                          output);


    // initial field from restart
    if (probtype.get<int>("RESTART"))
    {
      // read the restart information, set vectors and variables
      genalphaint.ReadRestart(genprob.restart);
    }

    // call generalized-alpha time-integration scheme
    genalphaint.GenAlphaIntegrateTo(fdyn.get<int>("NUMSTEP"),fdyn.get<double>("MAXTIME"));

  }
  else
  {
    dserror("Unknown solver type for drt_condif");
  }

  //---------- this is the end. Beautiful friend. My only friend, The end.
  // thanks to RefCountPtr<> we do not need to delete anything here!

  return;

} // end of dyn_condif_drt()

#endif  // #ifdef CCADISCRET
