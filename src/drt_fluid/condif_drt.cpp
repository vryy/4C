/*!----------------------------------------------------------------------
\file
\brief Control routine for convection-diffusion time integration. Includes

     o Singele step one-step-theta time integration

     o Two step BDF2 Gear's methode with one-step-theta start step



<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
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

#include "condif_drt.H"
#include "condifimplicitintegration.H"
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
void dyn_condif_drt()
{

  DSTraceHelper dst("dyn_condif_drt");


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
  ParameterList condiftimeparams;
  CondifImplicitTimeInt::SetDefaults(condiftimeparams);

  // the default time step size
  condiftimeparams.set<double>           ("time step size"           ,fdyn->dt);
  // max. sim. time
  condiftimeparams.set<double>           ("total time"               ,fdyn->maxtime);
  // parameter for time-integration
  condiftimeparams.set<double>           ("theta"                    ,fdyn->theta);
  // which kind of time-integration
  condiftimeparams.set<FLUID_TIMEINTTYPE>("time int algo"            ,fdyn->iop);
  // bound for the number of timesteps
  condiftimeparams.set<int>              ("max number timesteps"     ,fdyn->nstep);
  // number of steps with start algorithm
  condiftimeparams.set<int>              ("number of start steps"    ,fdyn->nums);
  // parameter for start algo
  condiftimeparams.set<double>           ("start theta"              ,fdyn->thetas);
  // restart
  condiftimeparams.set                  ("write restart every"       ,fdyn->uprestart);

  //--------------------------------------------------
  // velocity field
  condiftimeparams.set<int>              ("condif velocity field"     ,fdyn->cdvel);

  //--------------------------------------------------
  // evaluate error for test flows with analytical solutions
  condiftimeparams.set                  ("eval err for analyt sol"   ,fdyn->init);


  //--------------------------------------------------
  // create all vectors and variables associated with the time
  // integration (call the constructor)
  CondifImplicitTimeInt condifimplicit(actdis,
                                      solver,
                                      condiftimeparams,
                                      output);

  //--------------------------------------------------
  // read the restart information, set vectors and variables
  if (genprob.restart) condifimplicit.ReadRestart(genprob.restart);

  //--------------------------------------------------
  // do the time integration (start algo and standard algo)
  condifimplicit.Integrate();

  //---------- this is the end. Beautiful friend. My only friend, The end.
  // thanks to RefCountPtr<> we do not need to delete anything here!

  return;

} // end of dyn_condif_drt()

#endif  // #ifdef D_FLUID
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
