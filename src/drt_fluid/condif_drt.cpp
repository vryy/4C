/*!----------------------------------------------------------------------
\file condif_drt.cpp
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
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

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
 * Time integration loop for convection-diffusion problems
 *
 *        o One-step-theta
 *        o BDF2
 *        o Generalized-alpha
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

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  LINALG::Solver solver(solveparams,actdis->Comm(),allfiles.out_err);
  solver.TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  FLUID_TIMEINTTYPE iop = Teuchos::getIntegralValue<FLUID_TIMEINTTYPE>(fdyn,"TIMEINTEGR");
  if(iop == timeint_stationary or
     iop == timeint_one_step_theta or
     iop == timeint_bdf2
    )
  {
    // -------------------------------------------------------------------
    // create a convection-diffusion one-step-theta/BDF2 time integrator
    // -------------------------------------------------------------------
    ParameterList condiftimeparams;
    CondifImplicitTimeInt::SetDefaults(condiftimeparams);

    // the default time step size
    condiftimeparams.set<double>           ("time step size"           ,fdyn.get<double>("TIMESTEP"));
    // max. sim. time
    condiftimeparams.set<double>           ("total time"               ,fdyn.get<double>("MAXTIME"));
    // parameter for time-integration
    condiftimeparams.set<double>           ("theta"                    ,fdyn.get<double>("THETA"));
    // which kind of time-integration
    condiftimeparams.set<FLUID_TIMEINTTYPE>("time int algo"            ,iop);
    // bound for the number of timesteps
    condiftimeparams.set<int>              ("max number timesteps"     ,fdyn.get<int>("NUMSTEP"));
    // number of steps with start algorithm
    condiftimeparams.set<int>              ("number of start steps"    ,fdyn.get<int>("NUMSTASTEPS"));
    // parameter for start algo
    condiftimeparams.set<double>           ("start theta"              ,fdyn.get<double>("START_THETA"));
    // restart
    condiftimeparams.set                  ("write restart every"       ,fdyn.get<int>("RESTARTEVRY"));
    // solution output
    condiftimeparams.set                  ("write solution every"      ,fdyn.get<int>("UPRES"));

    //--------------------------------------------------
    // velocity field
    condiftimeparams.set<int>              ("condif velocity field"     ,Teuchos::getIntegralValue<int>(fdyn,"CD_VELOCITY"));
    // discontinuity capturing?
    condiftimeparams.set<int>              ("discontinuity capturing"   ,Teuchos::getIntegralValue<int>(fdyn,"DISC_CAPT"));

    //--------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor)
    CondifImplicitTimeInt condifimplicit(actdis,
                                        solver,
                                        condiftimeparams,
                                        output);

    //--------------------------------------------------
    if (probtype.get<int>("RESTART"))
    {
      // read the restart information, set vectors and variables
      condifimplicit.ReadRestart(probtype.get<int>("RESTART"));
    }

    //--------------------------------------------------
    // do the time integration (start algo and standard algo)
    condifimplicit.Integrate();

  }
  else if (iop == timeint_gen_alpha)
  {

    // -------------------------------------------------------------------
    // create a convection-diffusion generalized-alpha time integrator
    // -------------------------------------------------------------------
    // ------------------ set the parameter list
    ParameterList condiftimeparams;

    // the default time step size
    condiftimeparams.set<double>           ("time step size"           ,fdyn.get<double>("TIMESTEP"));
    // max. sim. time
    condiftimeparams.set<double>           ("total time"               ,fdyn.get<double>("MAXTIME"));
    // parameters for time-integration
    condiftimeparams.set<double>           ("alpha_M"                  ,fdyn.get<double>("ALPHA_M"));
    // parameters for time-integration
    condiftimeparams.set<double>           ("alpha_F"                  ,fdyn.get<double>("ALPHA_F"));
    condiftimeparams.set<int>              ("max number timesteps"     ,fdyn.get<int>("NUMSTEP"));
    // restart
    condiftimeparams.set                  ("write restart every"       ,fdyn.get<int>("RESTARTEVRY"));
    // solution output
    condiftimeparams.set                  ("write solution every"      ,fdyn.get<int>("UPRES"));

    //--------------------------------------------------
    // velocity field
    condiftimeparams.set<int>              ("condif velocity field"     ,Teuchos::getIntegralValue<int>(fdyn,"CD_VELOCITY"));
    // discontinuity capturing?
    condiftimeparams.set<int>              ("discontinuity capturing"   ,Teuchos::getIntegralValue<int>(fdyn,"DISC_CAPT"));

    //--------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor)
    CondifGenAlphaIntegration genalphaint(actdis,
                                          solver,
                                          condiftimeparams,
                                          output);


    //------------- initialize the field from input or restart
    if (probtype.get<int>("RESTART"))
    {
      // read the restart information, set vectors and variables
      genalphaint.ReadRestart(genprob.restart);
    }

    //------------------------- do timeintegration till maxtime
    genalphaint.GenAlphaIntegrateTo(fdyn.get<int>("NUMSTEP"),fdyn.get<double>("MAXTIME"));

  }
  else
  {
    dserror("Unknown time type for drt_condif");
  }

  //---------- this is the end. Beautiful friend. My only friend, The end.
  // thanks to RefCountPtr<> we do not need to delete anything here!

  return;

} // end of dyn_condif_drt()



#endif  // #ifdef CCADISCRET
