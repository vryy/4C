/*!----------------------------------------------------------------------
\file xfluid_dyn_nln_drt.cpp
\brief Control routine for fluid time integration. Includes

     o Singele step one-step-theta time integration

     o Two step BDF2 Gear's methode with one-step-theta start step



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
#include <set>
#include <Teuchos_TimeMonitor.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "xfluid_dyn_nln_drt.H"
#include "xfluidimplicitintegration.H"
#include "../drt_lib/drt_resulttest.H"
#include "xfluidresulttest.H"
#include "../drt_lib/drt_globalproblem.H"

extern struct _FIELD    *field;
extern struct _GENPROB  genprob;
extern struct _FILES  	allfiles;
extern struct _IO_FLAGS ioflags;
extern struct _SOLVAR   *solv;
extern         ALLDYNA  *alldyn;
extern struct _CURVE  	*curve;


/*----------------------------------------------------------------------*
 * Time integration loop for fluid.
 *
 *        o One-step-theta
 *        o BDF2
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
  SOLVAR        *solidsolv  = &solv[genprob.numsf];

  dsassert(genprob.timetyp==time_dynamic, "alldyn not allocated!!!\n For stationary computations, choose TIMETYP Dynamic and switch time integration for fluid to stationary instead");
  FLUID_DYNAMIC *fdyn     = alldyn[genprob.numff].fdyn;
  STRUCT_DYNAMIC *sdyn    = alldyn[genprob.numsf].sdyn;
  
  fdyn->step              =   0;
  fdyn->acttime           = 0.0;

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RCP<ParameterList> solveparams = rcp(new ParameterList());
  LINALG::Solver solver(solveparams,fluiddis->Comm(),allfiles.out_err);
  solver.TranslateSolverParameters(*solveparams,fluidsolv);
  fluiddis->ComputeNullSpaceIfNecessary(*solveparams);

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
    XFluidImplicitTimeInt::SetDefaults(fluidtimeparams);

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
    // set linearisation scheme
    if(fdyn->ite==2)
    {
      fluidtimeparams.set<bool>("Use reaction terms for linearisation",true);
    }
    else
    {
      fluidtimeparams.set<bool>("Use reaction terms for linearisation",false);
    }
    // maximum number of nonlinear iteration steps
    fluidtimeparams.set<int>             ("max nonlin iter steps"     ,fdyn->itemax);
    // stop nonlinear iteration when both incr-norms are below this bound
    fluidtimeparams.set<double>          ("tolerance for nonlin iter" ,fdyn->ittol);

    // ----------------------------------------------- restart and output
    // restart
    fluidtimeparams.set                  ("write restart every"       ,fdyn->uprestart);
    // solution output
    fluidtimeparams.set                  ("write solution every"      ,fdyn->upres);
    // flag for writing stresses
    fluidtimeparams.set                  ("write stresses"            ,ioflags.fluid_stress);

    //--------------------------------------------------
    // evaluate error for test flows with analytical solutions
    fluidtimeparams.set                  ("eval err for analyt sol"   ,fdyn->init);


    //--------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor)
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom
    XFluidImplicitTimeInt fluidimplicit(
    		fluiddis,
    		soliddis,
    		solver,
    		fluidtimeparams,
    		fluidoutput,
    		solidoutput,
    		false);

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
        fluidimplicit.SetInitialFlowField(fdyn->init,fdyn->startfuncno);
      }
    }

    //--------------------------------------------------
    // do the time integration (start algo and standard algo)
    fluidimplicit.Integrate();



    //--------------------------------------------------
    // do the result test
#ifdef RESULTTEST
    DRT::ResultTestManager testmanager(fluiddis->Comm());
    testmanager.AddFieldTest(rcp(new XFluidResultTest(fluidimplicit)));
    testmanager.TestAll();
#endif
  }
  else
  {
    dserror("Unknown time type for drt xfluid");
  }

  return;

} // end of xdyn_fluid_drt()

#endif  // #ifdef CCADISCRET
