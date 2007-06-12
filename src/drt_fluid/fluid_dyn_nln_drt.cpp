/*!----------------------------------------------------------------------
\file
\brief Control routine for fluid time integration. Includes

     o Singele step one-step-theta time integration

     o Two step BDF2 Gear's methode with one-step-theta start step



<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifdef D_FLUID

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "fluid_dyn_nln_drt.H"
#include "fluidimplicitintegration.H"
#include "fluid_genalpha_integration.H"
#include "../drt_lib/drt_resulttest.H"
#include "fluidresulttest.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

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
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;

/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | INT                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern INT            numcurve;
extern struct _CURVE *curve;




/*----------------------------------------------------------------------*
 * Time integration loop for fluid.
 *
 *        o One-step-theta
 *        o BDF2
 *
 *----------------------------------------------------------------------*/
void dyn_fluid_drt()
{

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numff,0);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();


  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  DiscretizationWriter output(actdis);
  output.WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR        *actsolv  = &solv[0];

  FLUID_DYNAMIC *fdyn     = alldyn[0].fdyn;
  fdyn->step              =   0;
  fdyn->acttime           = 0.0;

  // -------------------------------------------------------------------
  // init all applied time curves
  // -------------------------------------------------------------------
  for (int actcurve=0; actcurve<numcurve; actcurve++)
  {
   /* the last three parameters are obsolete!!! */
   dyn_init_curve(actcurve,fdyn->step,fdyn->dt,fdyn->maxtime);
  }

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  LINALG::Solver solver(solveparams,actdis->Comm(),allfiles.out_err);
  solver.TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  if(fdyn->iop == timeint_stationary
     ||
     fdyn->iop == timeint_one_step_theta
     ||
     fdyn->iop == timeint_bdf2
    )
  {
    // -------------------------------------------------------------------
    // create a fluid nonlinear time integrator
    // -------------------------------------------------------------------
    ParameterList fluidtimeparams;
    FluidImplicitTimeInt::SetDefaults(fluidtimeparams);

    // number of degrees of freedom
    fluidtimeparams.set<int>              ("number of velocity degrees of freedom" ,genprob.ndim);
    // the default time step size
    fluidtimeparams.set<double>           ("time step size"           ,fdyn->dt);
    // max. sim. time
    fluidtimeparams.set<double>           ("total time"               ,fdyn->maxtime);
    // parameter for time-integration
    fluidtimeparams.set<double>           ("theta"                    ,fdyn->theta);
    // which kind of time-integration
    fluidtimeparams.set<FLUID_TIMEINTTYPE>("time int algo"            ,fdyn->iop);
    // bound for the number of timesteps
    fluidtimeparams.set<int>              ("max number timesteps"     ,fdyn->nstep);
    // number of steps with start algorithm
    fluidtimeparams.set<int>              ("number of start steps"    ,fdyn->nums);
    // parameter for start algo
    fluidtimeparams.set<double>           ("start theta"              ,fdyn->thetas);


    // ---------------------------------------------- nonlinear iteration
    // maximum number of nonlinear iteration steps
    fluidtimeparams.set<int>             ("max nonlin iter steps"     ,fdyn->itemax);
    // stop nonlinear iteration when both incr-norms are below this bound
    fluidtimeparams.set<double>          ("tolerance for nonlin iter" ,fdyn->ittol);

    // restart
    fluidtimeparams.set                  ("write restart every"       ,fdyn->uprestart);
    // solution output
    fluidtimeparams.set                  ("write solution every"      ,fdyn->upres);    

    //--------------------------------------------------
    // evaluate error for test flows with analytical solutions
    fluidtimeparams.set                  ("eval err for analyt sol"   ,fdyn->init);


    //--------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor)
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom
    FluidImplicitTimeInt fluidimplicit(actdis,
                                       solver,
                                       fluidtimeparams,
                                       output);

    //--------------------------------------------------
    if (genprob.restart)
    {
      // read the restart information, set vectors and variables
      fluidimplicit.ReadRestart(genprob.restart);
    }
    else
    {
      // set initial field for analytical test problems etc
      if(fdyn->init>0)
      {
        fluidimplicit.SetInitialFlowField(fdyn->init);
      }
    }

    //--------------------------------------------------
    // do the time integration (start algo and standard algo)
    fluidimplicit.Integrate();



    //--------------------------------------------------
    // do the result test
#ifdef RESULTTEST
    DRT::ResultTestManager testmanager(actdis->Comm());
    testmanager.AddFieldTest(rcp(new FluidResultTest(fluidimplicit)));
    testmanager.TestAll();
#endif
  }
  else if (fdyn->iop == timeint_gen_alpha)
  {
    // -------------------------------------------------------------------
    // create a generalised alpha time integrator for fluid problems
    // -------------------------------------------------------------------
    // ------------------ set the parameter list
    ParameterList fluidtimeparams;

    // number of degrees of freedom
    fluidtimeparams.set<int>              ("number of velocity degrees of freedom" ,genprob.ndim);
    // the default time step size
    fluidtimeparams.set<double>           ("time step size"           ,fdyn->dt);
    // max. sim. time
    fluidtimeparams.set<double>           ("total time"               ,fdyn->maxtime);
    // parameters for time-integration
    fluidtimeparams.set<double>           ("alpha_M"                  ,fdyn->alpha_m);
    // parameters for time-integration
    fluidtimeparams.set<double>           ("alpha_F"                  ,fdyn->alpha_f);
    // bound for the number of timesteps
    fluidtimeparams.set<int>              ("max number timesteps"     ,fdyn->nstep);

    // ---------------------------------------------- nonlinear iteration
    // maximum number of nonlinear iteration steps
    fluidtimeparams.set<int>             ("max nonlin iter steps"     ,fdyn->itemax);
    // stop nonlinear iteration when both incr-norms are below this bound
    fluidtimeparams.set<double>          ("tolerance for nonlin iter" ,fdyn->ittol);

    // ----------------------------------------------- restart and output
    fluidtimeparams.set                  ("write restart every"       ,fdyn->uprestart);

    //------------evaluate error for test flows with analytical solutions
    fluidtimeparams.set                  ("eval err for analyt sol"   ,fdyn->init);

    //--------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor)
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom
    FluidGenAlphaIntegration genalphaint(actdis,
                                         solver,
                                         fluidtimeparams,
                                         output);


    //------------- initialise the field from input or restart
    if (genprob.restart)
    {
      // read the restart information, set vectors and variables
      genalphaint.ReadRestart(genprob.restart);
    }
    else
    {
      // set initial field for analytical test problems etc
      if(fdyn->init>0)
      {
        genalphaint.SetInitialFlowField(fdyn->init);
      }
    }

    //------------------------- do timeintegration till maxtime
    genalphaint.GenAlphaIntegrateTo(fdyn->nstep,fdyn->maxtime);

    //--------------------------------------------------
    // do the result test
#ifdef RESULTTEST
#if 0
    DRT::ResultTestManager testmanager(actdis->Comm());
    testmanager.AddFieldTest(rcp(new FluidResultTest(genalphaint)));
    testmanager.TestAll();
#endif
#endif

  }
  else
  {
    dserror("Unknown time type for drt fluid");
  }

  //---------- this is the end. Beautiful friend. My only friend, The end.
  // thanks to RefCountPtr<> we do not need to delete anything here!

  return;

} // end of dyn_fluid_drt()

#endif  // #ifdef D_FLUID
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
