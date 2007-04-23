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
#include "drt_resulttest.H"
#include "fluidresulttest.H"


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

  DSTraceHelper dst("dyn_fluid_drt");


  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  {
    vector<RefCountPtr<DRT::Discretization> >* fool =
              (vector<RefCountPtr<DRT::Discretization> >*)field[0].ccadis;
    actdis = (*fool)[0];
  }
  // set degrees of freedom in the discretization
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

  // the only parameter required here is the number of velocity degrees of freedom
  FluidImplicitTimeInt fluidimplicit(*actdis,
                                      solver,
                                      fluidtimeparams,
                                      output);

  fluidimplicit.Integrate();

#ifdef RESULTTEST
  DRT::ResultTestManager testmanager(actdis->Comm());
  testmanager.AddFieldTest(rcp(new FluidResultTest(fluidimplicit)));
  testmanager.TestAll();
#endif

  //---------- this is the end. Beautiful friend. My only friend, The end.
  // thanks to RefCountPtr<> we do not need to delete anything here!

  return;

} // end of dyn_fluid_drt()

#endif  // #ifdef D_FLUID
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
