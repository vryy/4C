/*!----------------------------------------------------------------------
\file fluid_dyn_nln_drt.cpp
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

#include <ctime>
#include <cstdlib>
#include <iostream>

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "fluid_dyn_nln_drt.H"
#include "fluidimplicitintegration.H"
#include "fluid_genalpha_integration.H"
#include "../drt_lib/drt_resulttest.H"
#include "fluidresulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"

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
  IO::DiscretizationWriter output(actdis);
  output.WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR        *actsolv  = &solv[0];

  const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& probsize = DRT::Problem::Instance()->ProblemSizeParams();
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();

  if (actdis->Comm().MyPID()==0)
    DRT::PrintDefaultParameters(std::cout, fdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RCP<ParameterList> solveparams = rcp(new ParameterList());
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
    // create a fluid nonlinear time integrator
    // -------------------------------------------------------------------
    ParameterList fluidtimeparams;
    FluidImplicitTimeInt::SetDefaults(fluidtimeparams);

    // number of degrees of freedom
    fluidtimeparams.set<int>              ("number of velocity degrees of freedom" ,probsize.get<int>("DIM"));
    // the default time step size
    fluidtimeparams.set<double>           ("time step size"           ,fdyn.get<double>("TIMESTEP"));
    // max. sim. time
    fluidtimeparams.set<double>           ("total time"               ,fdyn.get<double>("MAXTIME"));
    // parameter for time-integration
    fluidtimeparams.set<double>           ("theta"                    ,fdyn.get<double>("THETA"));
    // which kind of time-integration
    fluidtimeparams.set<FLUID_TIMEINTTYPE>("time int algo"            ,iop);
    // bound for the number of timesteps
    fluidtimeparams.set<int>              ("max number timesteps"     ,fdyn.get<int>("NUMSTEP"));
    // number of steps with start algorithm
    fluidtimeparams.set<int>              ("number of start steps"    ,fdyn.get<int>("NUMSTASTEPS"));
    // parameter for start algo
    fluidtimeparams.set<double>           ("start theta"              ,fdyn.get<double>("START_THETA"));


    // ---------------------------------------------- nonlinear iteration
    // set linearisation scheme
    fluidtimeparams.set<bool>("Use reaction terms for linearisation",
                              Teuchos::getIntegralValue<int>(fdyn,"NONLINITER")==2);
    // maximum number of nonlinear iteration steps
    fluidtimeparams.set<int>             ("max nonlin iter steps"     ,fdyn.get<int>("ITEMAX"));
    // stop nonlinear iteration when both incr-norms are below this bound
    fluidtimeparams.set<double>          ("tolerance for nonlin iter" ,fdyn.get<double>("CONVTOL"));

    // ----------------------------------------------- restart and output
    // restart
    fluidtimeparams.set                  ("write restart every"       ,fdyn.get<int>("RESTARTEVRY"));
    // solution output
    fluidtimeparams.set                  ("write solution every"      ,fdyn.get<int>("UPRES"));
    // flag for writing stresses
    fluidtimeparams.set                  ("write stresses"            ,Teuchos::getIntegralValue<int>(ioflags,"FLUID_STRESS"));

    //--------------------------------------------------
    // evaluate error for test flows with analytical solutions
    int init = Teuchos::getIntegralValue<int>(fdyn,"INITIALFIELD");
    fluidtimeparams.set                  ("eval err for analyt sol"   ,init);

    // (fine-scale) subgrid viscosity?
    fluidtimeparams.set<int>              ("fs subgrid viscosity"   ,Teuchos::getIntegralValue<int>(fdyn,"SUBGRIDVISC"));

    // hand down the TURBULENCE MODEL parameters to the fluid algorithm
    {
      fluidtimeparams.sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");

      fluidtimeparams.sublist("TURBULENCE MODEL").set<string>("statistics outfile",allfiles.outputfile_kenner);
    }

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
    if (probtype.get<int>("RESTART"))
    {
      // read the restart information, set vectors and variables
      fluidimplicit.ReadRestart(probtype.get<int>("RESTART"));
    }
    else
    {
      // set initial field for analytical test problems etc
      if(init>0)
      {
        int startfuncno = fdyn.get<int>("STARTFUNCNO");
        if (init!=2 and init!=3)
        {
          startfuncno=-1;
        }
        fluidimplicit.SetInitialFlowField(init,startfuncno);
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
  else if (iop == timeint_gen_alpha)
  {
    // -------------------------------------------------------------------
    // create a generalised alpha time integrator for fluid problems
    // -------------------------------------------------------------------
    // ------------------ set the parameter list
    ParameterList fluidtimeparams;

    // number of degrees of freedom
    fluidtimeparams.set<int>              ("number of velocity degrees of freedom" ,probsize.get<int>("DIM"));
    // the default time step size
    fluidtimeparams.set<double>           ("time step size"           ,fdyn.get<double>("TIMESTEP"));
    // max. sim. time
    fluidtimeparams.set<double>           ("total time"               ,fdyn.get<double>("MAXTIME"));
    // parameters for time-integration
    fluidtimeparams.set<double>           ("alpha_M"                  ,fdyn.get<double>("ALPHA_M"));
    // parameters for time-integration
    fluidtimeparams.set<double>           ("alpha_F"                  ,fdyn.get<double>("ALPHA_F"));
    // bound for the number of timesteps
    fluidtimeparams.set<int>              ("max number timesteps"     ,fdyn.get<int>("NUMSTEP"));

    // ---------------------------------------------- nonlinear iteration
    //     // set linearisation scheme
    fluidtimeparams.set<bool>("Use reaction terms for linearisation",
                              Teuchos::getIntegralValue<int>(fdyn,"NONLINITER")==2);
    // maximum number of nonlinear iteration steps
    fluidtimeparams.set<int>             ("max nonlin iter steps"     ,fdyn.get<int>("ITEMAX"));
    // stop nonlinear iteration when both incr-norms are below this bound
    fluidtimeparams.set<double>          ("tolerance for nonlin iter" ,fdyn.get<double>("CONVTOL"));

    // ----------------------------------------------- restart and output
    fluidtimeparams.set                  ("write restart every"       ,fdyn.get<int>("RESTARTEVRY"));
    // solution output
    fluidtimeparams.set                  ("write solution every"      ,fdyn.get<int>("UPRES"));
    // flag for writing stresses
    fluidtimeparams.set                  ("write stresses"            ,Teuchos::getIntegralValue<int>(ioflags,"FLUID_STRESS"));

    //------------evaluate error for test flows with analytical solutions
    int init = Teuchos::getIntegralValue<int>(fdyn,"INITIALFIELD");
    fluidtimeparams.set                  ("eval err for analyt sol"   ,init);

    // hand down the STABILIZATION parameters to the fluid algorithm
    fluidtimeparams.sublist("STABILIZATION")=fdyn.sublist("STABILIZATION");

    // hand down the TURBULENCE MODEL parameters to the fluid algorithm
    {
      fluidtimeparams.sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");

      fluidtimeparams.sublist("TURBULENCE MODEL").set<string>("statistics outfile",allfiles.outputfile_kenner);
    }

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
    if (probtype.get<int>("RESTART"))
    {
      // read the restart information, set vectors and variables
      genalphaint.ReadRestart(probtype.get<int>("RESTART"));
    }
    else
    {
      // set initial field for analytical test problems etc
      if(init>0)
      {
        int startfuncno = fdyn.get<int>("STARTFUNCNO");
        if (init!=2 and init!=3)
        {
          startfuncno=-1;
        }
        genalphaint.SetInitialFlowField(init,startfuncno);
      }
    }

    //------------------------- do timeintegration till maxtime
    genalphaint.GenAlphaIntegrateTo(fdyn.get<int>("NUMSTEP"),fdyn.get<double>("MAXTIME"));

    //--------------------------------------------------
    // do the result test
#ifdef RESULTTEST
    DRT::ResultTestManager testmanager(actdis->Comm());
    testmanager.AddFieldTest(rcp(new FluidResultTest(genalphaint)));
    testmanager.TestAll();
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

#endif  // #ifdef CCADISCRET
