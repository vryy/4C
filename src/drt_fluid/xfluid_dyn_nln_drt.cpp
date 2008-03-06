/*!----------------------------------------------------------------------
\file xfluid_dyn_nln_drt.cpp
\brief Main control routine for all fluid (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     o generalized-alpha time-integration scheme

     and stationary solver.

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
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

#include "xfluid_dyn_nln_drt.H"
#include "xfluidimplicitintegration.H"
#include "../drt_lib/drt_resulttest.H"
#include "xfluidresulttest.H"
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
 * Main control routine for fluid including various solvers:
 *
 *        o instationary one-step-theta
 *        o instationary BDF2
 *        o instationary generalized-alpha
 *        o stationary
 *
 *----------------------------------------------------------------------*/
void xdyn_fluid_drt()
{
  cout << "Hallo, ich bin ein Fluid-XFEM Problem" << endl;

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> fluiddis = null;
  fluiddis = DRT::Problem::Instance()->Dis(genprob.numff,0);
  RCP<DRT::Discretization> soliddis = null;
  soliddis = DRT::Problem::Instance()->Dis(genprob.numsf,0);

  const int fmyrank = fluiddis->Comm().MyPID();
  cout << "FluidProc: " << fmyrank << endl;
  flush(cout);
  const int smyrank = soliddis->Comm().MyPID();
  cout << "SolidProc: " << smyrank << endl;
  flush(cout);
  //cout << *soliddis;

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!fluiddis->Filled()) fluiddis->FillComplete();
  if (!soliddis->Filled()) soliddis->FillComplete();


  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  IO::DiscretizationWriter fluidoutput(fluiddis);
  fluidoutput.WriteMesh(0,0.0);
  IO::DiscretizationWriter solidoutput(soliddis);
  if (soliddis->NumGlobalElements() > 0)
  {
      solidoutput.WriteMesh(0,0.0);
  }

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR        *fluidsolv  = &solv[genprob.numff];
  //SOLVAR        *solidsolv  = &solv[genprob.numsf];

  const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& probsize = DRT::Problem::Instance()->ProblemSizeParams();
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();

  if (fluiddis->Comm().MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, fdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RCP<ParameterList> solveparams = rcp(new ParameterList());
  LINALG::Solver solver(solveparams,fluiddis->Comm(),allfiles.out_err);
  solver.TranslateSolverParameters(*solveparams,fluidsolv);
  fluiddis->ComputeNullSpaceIfNecessary(*solveparams);

  // -------------------------------------------------------------------
  // create a second solver for SIMPLER preconditioner if chosen from input
  // -------------------------------------------------------------------
  if (getIntegralValue<int>(fdyn,"SIMPLER"))
  {
    ParameterList& p = solveparams->sublist("SIMPLER");
    RCP<ParameterList> params = rcp(&p,false);
    LINALG::Solver s(params,fluiddis->Comm(),allfiles.out_err);
    s.TranslateSolverParameters(*params,&solv[genprob.numfld]);
  }

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  ParameterList fluidtimeparams;

  fluidtimeparams.set<int>("Simple Preconditioner",Teuchos::getIntegralValue<int>(fdyn,"SIMPLER"));

  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  fluidtimeparams.set<int>              ("number of velocity degrees of freedom" ,probsize.get<int>("DIM"));

  // ------------------------------------------------ basic scheme, i.e.
  // --------------------- solving nonlinear or linearised flow equation
  fluidtimeparams.set<int>("type of nonlinear solve" ,
					 Teuchos::getIntegralValue<int>(fdyn,"DYNAMICTYP"));

  // -------------------------------------------------- time integration
  // the default time step size
  fluidtimeparams.set<double>           ("time step size"           ,fdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  fluidtimeparams.set<double>           ("total time"               ,fdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  fluidtimeparams.set<int>              ("max number timesteps"     ,fdyn.get<int>("NUMSTEP"));

  // ---------------------------------------------- nonlinear iteration
  // set linearisation scheme
  fluidtimeparams.set<bool>("Use reaction terms for linearisation",
                           Teuchos::getIntegralValue<int>(fdyn,"NONLINITER")==2);
  // maximum number of nonlinear iteration steps
  fluidtimeparams.set<int>             ("max nonlin iter steps"     ,fdyn.get<int>("ITEMAX"));
  // stop nonlinear iteration when both incr-norms are below this bound
  fluidtimeparams.set<double>          ("tolerance for nonlin iter" ,fdyn.get<double>("CONVTOL"));
  // set convergence check
  fluidtimeparams.set<string>          ("CONVCHECK"  ,fdyn.get<string>("CONVCHECK"));
  // set adaptoive linear solver tolerance
  fluidtimeparams.set<bool>            ("ADAPTCONV",getIntegralValue<int>(fdyn,"ADAPTCONV")==1);
  fluidtimeparams.set<double>          ("ADAPTCONV_BETTER",fdyn.get<double>("ADAPTCONV_BETTER"));

  // ----------------------------------------------- restart and output
  // restart
  fluidtimeparams.set                  ("write restart every"       ,fdyn.get<int>("RESTARTEVRY"));
  // solution output
  fluidtimeparams.set                  ("write solution every"      ,fdyn.get<int>("UPRES"));
  // flag for writing stresses
  fluidtimeparams.set                  ("write stresses"            ,Teuchos::getIntegralValue<int>(ioflags,"FLUID_STRESS"));
  // ---------------------------------------------------- lift and drag
  fluidtimeparams.set<int>("liftdrag",Teuchos::getIntegralValue<int>(fdyn,"LIFTDRAG"));

  // -----------evaluate error for test flows with analytical solutions
  int init = Teuchos::getIntegralValue<int>(fdyn,"INITIALFIELD");
  fluidtimeparams.set                  ("eval err for analyt sol"   ,init);

  // ---------------------------- fine-scale subgrid viscosity approach
  fluidtimeparams.set<string>           ("fs subgrid viscosity"   ,fdyn.get<string>("FSSUGRVISC"));

  // -----------------------sublist containing stabilization parameters
  fluidtimeparams.sublist("STABILIZATION")=fdyn.sublist("STABILIZATION");

  // --------------------------sublist containing turbulence parameters
  {
    fluidtimeparams.sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");

    fluidtimeparams.sublist("TURBULENCE MODEL").set<string>("statistics outfile",allfiles.outputfile_kenner);
  }

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
    fluidtimeparams.set<FLUID_TIMEINTTYPE>("time int algo",iop);
    // parameter theta for time-integration schemes
    fluidtimeparams.set<double>           ("theta"                    ,fdyn.get<double>("THETA"));
    // number of steps for potential start algorithm
    fluidtimeparams.set<int>              ("number of start steps"    ,fdyn.get<int>("NUMSTASTEPS"));
    // parameter theta for potential start algorithm
    fluidtimeparams.set<double>           ("start theta"              ,fdyn.get<double>("START_THETA"));

    //------------------------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor);
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom
    //------------------------------------------------------------------
    XFluidImplicitTimeInt fluidimplicit(
    		fluiddis,
    		soliddis,
    		solver,
    		fluidtimeparams,
    		fluidoutput,
    		solidoutput,
    		false);

    // initial field from restart or calculated by given function
    if (probtype.get<int>("RESTART"))
    {
      // read the restart information, set vectors and variables
      fluidimplicit.ReadRestart(probtype.get<int>("RESTART"));
    }
    else
    {
      // set initial field by given function
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

    fluidtimeparams.set<FILE*>("err file",allfiles.out_err);

    // call time-integration (or stationary) scheme
    fluidimplicit.Integrate();

    // do result test if required
    DRT::ResultTestManager testmanager(fluiddis->Comm());
    testmanager.AddFieldTest(rcp(new XFluidResultTest(fluidimplicit)));
    testmanager.TestAll();

  }
  else
  {
    dserror("Unknown solver type for drt_xfluid");
  }

  //---------- this is the end. Beautiful friend. My only friend, The end.
  // thanks to RefCountPtr<> we do not need to delete anything here!

  return;

} // end of xdyn_fluid_drt()

#endif  // #ifdef CCADISCRET
