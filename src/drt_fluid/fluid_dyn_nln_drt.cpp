/*!----------------------------------------------------------------------
\file fluid_dyn_nln_drt.cpp
\brief Main control routine for all fluid (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     o generalized-alpha time-integration scheme

     and stationary solver.

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
#include "fluid_projectionmethod.H"
#include "../drt_lib/drt_resulttest.H"
#include "fluidresulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_inpar/drt_validparameters.H"

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*
 * Main control routine for fluid including various solvers:
 *
 *        o instationary one-step-theta
 *        o instationary BDF2
 *        o instationary generalized-alpha (two versions)
 *        o stationary
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
	if (!actdis->Filled() || !actdis->HaveDofs()) actdis->FillComplete();

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
	const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();

	if (actdis->Comm().MyPID()==0)
		DRT::INPUT::PrintDefaultParameters(std::cout, fdyn);

	// -------------------------------------------------------------------
	// create a solver
	// -------------------------------------------------------------------
	LINALG::Solver solver(DRT::Problem::Instance()->FluidSolverParams(),
			actdis->Comm(),
			DRT::Problem::Instance()->ErrorFile()->Handle());
	actdis->ComputeNullSpaceIfNecessary(solver.Params());

	// -------------------------------------------------------------------
	// create a second solver for SIMPLER preconditioner if chosen from input
	// -------------------------------------------------------------------
	if (getIntegralValue<int>(fdyn,"SIMPLER"))
	{
		solver.PutSolverParamsToSubParams("SIMPLER",
				DRT::Problem::Instance()->FluidPressureSolverParams());
	}

	// -------------------------------------------------------------------
	// set parameters in list required for all schemes
	// -------------------------------------------------------------------
	ParameterList fluidtimeparams;

	fluidtimeparams.set<int>("Simple Preconditioner",Teuchos::getIntegralValue<int>(fdyn,"SIMPLER"));

	// -------------------------------------- number of degrees of freedom
	// number of degrees of freedom
	fluidtimeparams.set<int> ("number of velocity degrees of freedom" ,probsize.get<int>("DIM"));

	// ---------------------------- low-Mach-number or incompressible flow
	fluidtimeparams.set<string> ("low-Mach-number solver"   ,fdyn.get<string>("LOWMACH"));

	// ------------------------------------------------ basic scheme, i.e.
	// --------------------- solving nonlinear or linearised flow equation
	fluidtimeparams.set<int> ("type of nonlinear solve" ,
			Teuchos::getIntegralValue<int>(fdyn,"DYNAMICTYP"));

	// -------------------------------------------------- time integration
	// the default time step size
	fluidtimeparams.set<double> ("time step size"      ,fdyn.get<double>("TIMESTEP"));
	// maximum simulation time
	fluidtimeparams.set<double> ("total time"          ,fdyn.get<double>("MAXTIME"));
	// maximum number of timesteps
	fluidtimeparams.set<int>    ("max number timesteps",fdyn.get<int>("NUMSTEP"));

	// -------- additional parameters in list for generalized-alpha scheme
#if 1
	// parameter alpha_M
	fluidtimeparams.set<double> ("alpha_M", fdyn.get<double>("ALPHA_M"));
	// parameter alpha_F
	fluidtimeparams.set<double> ("alpha_F", fdyn.get<double>("ALPHA_F"));
#else
	// parameter alpha_M
	fluidtimeparams.set<double> ("alpha_M", 1.-fdyn.get<double>("ALPHA_M"));
	// parameter alpha_F
	fluidtimeparams.set<double> ("alpha_F", 1.-fdyn.get<double>("ALPHA_F"));
#endif
	// parameter gamma
	fluidtimeparams.set<double> ("gamma", fdyn.get<double>("GAMMA"));

	// ------------------------------------------------- type of predictor
	fluidtimeparams.set<string> ("predictor", fdyn.get<string>("PREDICTOR"));

	// ---------------------------------------------- nonlinear iteration
	// set linearisation scheme
	fluidtimeparams.set<string>          ("Linearisation", fdyn.get<string>("NONLINITER"));
	// maximum number of nonlinear iteration steps
	fluidtimeparams.set<int>             ("max nonlin iter steps"  ,fdyn.get<int>("ITEMAX"));
	// stop nonlinear iteration when both incr-norms are below this bound
	fluidtimeparams.set<double>          ("tolerance for nonlin iter"  ,fdyn.get<double>("CONVTOL"));
	// set convergence check
	fluidtimeparams.set<string>          ("CONVCHECK" ,fdyn.get<string>("CONVCHECK"));
	// set adaptive linear solver tolerance
	fluidtimeparams.set<bool>            ("ADAPTCONV",getIntegralValue<int>(fdyn,"ADAPTCONV")==1);
	fluidtimeparams.set<double>          ("ADAPTCONV_BETTER", fdyn.get<double>("ADAPTCONV_BETTER"));

	// ----------------------------------------------- restart and output
	// restart
	fluidtimeparams.set ("write restart every", fdyn.get<int>("RESTARTEVRY"));
	// solution output
	fluidtimeparams.set ("write solution every", fdyn.get<int>("UPRES"));
	// flag for writing stresses
	fluidtimeparams.set ("write stresses"  ,Teuchos::getIntegralValue<int>(ioflags,"FLUID_STRESS"));

	// ---------------------------------------------------- lift and drag
	fluidtimeparams.set<int> ("liftdrag", Teuchos::getIntegralValue<int>(fdyn,"LIFTDRAG"));

	// -----------define which initial field u(0)=u_0 to use (by explicit
	//            function, with random perturbation, etc)
	int init = Teuchos::getIntegralValue<int> (fdyn,"INITIALFIELD");

	// -----------for special initial fields (beltrami), we are able to
	//            evaluate the error compared to an analytical solutions
	fluidtimeparams.set ("eval err for analyt sol"   ,init);

	// ------------------------------------------ form of convective term
	fluidtimeparams.set<string> ("form of convective term", fdyn.get<string>("CONVFORM"));

	// ------------------------------------ potential Neumann inflow terms
	fluidtimeparams.set<string> ("Neumann inflow",fdyn.get<string>("NEUMANNINFLOW"));

	// ---------------------------- fine-scale subgrid viscosity approach
	fluidtimeparams.set<string> ("fs subgrid viscosity", fdyn.get<string>("FSSUGRVISC"));

	// -----------------------sublist containing stabilization parameters
	fluidtimeparams.sublist("STABILIZATION")=fdyn.sublist("STABILIZATION");

	// --------------------------sublist containing turbulence parameters
	{
		fluidtimeparams.sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");

		fluidtimeparams.sublist("TURBULENCE MODEL").set<string>("statistics outfile",DRT::Problem::Instance()->OutputControlFile()->FileName());
	}

	// -------------------------------------------------------------------
	// additional parameters and algorithm call depending on respective
	// time-integration (or stationary) scheme
	// -------------------------------------------------------------------
	FLUID_TIMEINTTYPE iop = Teuchos::getIntegralValue<FLUID_TIMEINTTYPE>(fdyn,"TIMEINTEGR");
	if(iop == timeint_stationary or
			iop == timeint_one_step_theta or
			iop == timeint_bdf2 or
			iop == timeint_afgenalpha
	)
	{
		// -----------------------------------------------------------------
		// set additional parameters in list for
		// one-step-theta/BDF2/af-generalized-alpha/stationary scheme
		// -----------------------------------------------------------------
		// type of time-integration (or stationary) scheme
		fluidtimeparams.set<FLUID_TIMEINTTYPE>("time int algo",iop);
		// parameter theta for one-step-theta time-integration scheme
		fluidtimeparams.set<double> ("theta", fdyn.get<double>("THETA"));
		// -----------------------number of steps for potential start algorithm
		// (currently only implemented for af-generalized-alpha with backward-
		//  Euler-type parameter combination for starting algorithm)
		fluidtimeparams.set<int> ("number of start steps", fdyn.get<int>("NUMSTASTEPS"));
		// parameter theta for potential one-step-theta start algorithm
		// (currently not implemented)
		fluidtimeparams.set<double> ("start theta", fdyn.get<double>("START_THETA"));


		// check what fluidsolver is to be used
		int fluidsolver = Teuchos::getIntegralValue<int>(fdyn,"FLUID_SOLVER");
		switch(fluidsolver)
		{
		case fluid_solver_implicit:
		{
			//------------------------------------------------------------------
			// create all vectors and variables associated with the time
			// integration (call the constructor);
			// the only parameter from the list required here is the number of
			// velocity degrees of freedom
			//------------------------------------------------------------------
			FLD::FluidImplicitTimeInt fluidimplicit(actdis,
					solver,
					fluidtimeparams,
					output);

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

			fluidtimeparams.set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

			// call time-integration (or stationary) scheme
			fluidimplicit.Integrate();

			// do result test if required
			DRT::Problem::Instance()->AddFieldTest(rcp(new FLD::FluidResultTest(fluidimplicit)));
			DRT::Problem::Instance()->TestAll(actdis->Comm());
		}
		break;
		case fluid_solver_pressurecorrection:
		case fluid_solver_pressurecorrection_semiimplicit:
		{
			// pressure correction solver

			// check if implicit or semi-implicit projection solver
			fluidtimeparams.set<bool>("PROJ_IMPLICIT",fluidsolver==fluid_solver_pressurecorrection);

			// -------------------------------------------------------------------
			// create a second solver for Projection Methods if chosen from input
			// -------------------------------------------------------------------
			LINALG::Solver psolver(DRT::Problem::Instance()->FluidSolverParams(),
					actdis->Comm(),
					DRT::Problem::Instance()->ErrorFile()->Handle());
			psolver.PutSolverParamsToSubParams("FLUID PRESSURE SOLVER",
					DRT::Problem::Instance()->FluidPressureSolverParams());

			FLD::FluidProjectionMethod fluidprojectionmethod(actdis,
					solver,
					psolver,
					fluidtimeparams,
					output);
			fluidtimeparams.set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

			// initial field from restart or calculated by given function
			if (probtype.get<int>("RESTART"))
			{
				// read the restart information, set vectors and variables
				fluidprojectionmethod.ReadRestart(probtype.get<int>("RESTART"));
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
					//fluidprojectionmethod.SetInitialFlowField(init,startfuncno);
					dserror("no special initial flow field supported");
				}
			}

			// call time-integration scheme
			fluidprojectionmethod.Integrate();

			// do result test if required
			DRT::Problem::Instance()->AddFieldTest(rcp(new FLD::FluidResultTest(fluidprojectionmethod)));
			DRT::Problem::Instance()->TestAll(actdis->Comm());
		}
		break;
		default:
			dserror("fluid solver not known. choose FLUID_SOLVER=Implicit (standard) or Pressure Correction");
		}
	}
	else if (iop == timeint_gen_alpha)
	{
		// -------------------------------------------------------------------
		// no additional parameters required for generalized-alpha scheme
		// -------------------------------------------------------------------
		// create all vectors and variables associated with the time
		// integration (call the constructor);
		// the only parameter from the list required here is the number of
		// velocity degrees of freedom
		//------------------------------------------------------------------
		FLD::FluidGenAlphaIntegration genalphaint(actdis,
				solver,
				fluidtimeparams,
				output,
				false);


		// initial field from restart or calculated by given function
		if (probtype.get<int>("RESTART"))
		{
			// read the restart information, set vectors and variables
			genalphaint.ReadRestart(probtype.get<int>("RESTART"));
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
				genalphaint.SetInitialFlowField(init,startfuncno);
			}
		}

		// call generalized-alpha time-integration scheme
		genalphaint.GenAlphaTimeloop();

		// do result test if required
		DRT::Problem::Instance()->AddFieldTest(rcp(new FLD::FluidResultTest(genalphaint)));
		DRT::Problem::Instance()->TestAll(actdis->Comm());

	}
	else
	{
		dserror("Unknown solver type for drt_fluid");
	}

	return;

} // end of dyn_fluid_drt()


#endif  // #ifdef CCADISCRET
