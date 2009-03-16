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

#include "xfluid_dyn_nln_drt.H"
#include "xfluidimplicitintegration.H"
#include "xfluidresulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/linalg_utils.H"

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
void xdyn_fluid_drt()
{
  // -------------------------------------------------------------------
  // access the discretizations
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(genprob.numff,0);
  Teuchos::RCP<DRT::Discretization> soliddis = DRT::Problem::Instance()->Dis(genprob.numsf,0);
  
  if (fluiddis->Comm().MyPID() == 0)
  {
    std::cout << "Hallo, ich bin ein Fluid-XFEM Problem" << endl;
  }

//  const int fmyrank = fluiddis->Comm().MyPID();
//  std::cout << "FluidProc: " << fmyrank << endl;
//  const int smyrank = soliddis->Comm().MyPID();
//  std::cout << "SolidProc: " << smyrank << endl;

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!soliddis->HaveDofs()) soliddis->FillComplete();
  if (!fluiddis->HaveDofs()) fluiddis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  IO::DiscretizationWriter solidoutput(soliddis);
  if (soliddis->NumGlobalElements() > 0)
  {
      solidoutput.WriteMesh(0,0.0);
      // solid displacement
      Teuchos::RCP<Epetra_Vector> soliddispnp = LINALG::CreateVector(*soliddis->DofRowMap(),true);

      solidoutput.NewStep    (0,0.0);
      soliddispnp->PutScalar(0.0);
      solidoutput.WriteVector("soliddispnp", soliddispnp);
//      solidoutput_.NewStep    (step_,time_);
//      soliddispnp_->PutScalar(0.0);
//      solidoutput_.WriteVector("soliddispnp", soliddispnp_);
  }
  IO::DiscretizationWriter fluidoutput(fluiddis);
  fluidoutput.WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& probsize = DRT::Problem::Instance()->ProblemSizeParams();
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();

  if (fluiddis->Comm().MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, fdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  LINALG::Solver solver(DRT::Problem::Instance()->FluidSolverParams(),
                        fluiddis->Comm(),
                        DRT::Problem::Instance()->ErrorFile()->Handle());

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
  fluidtimeparams.set<bool>("Use reaction terms for linearisation",
                           Teuchos::getIntegralValue<int>(fdyn,"NONLINITER")==2);
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

  // -----------evaluate error for test flows with analytical solutions
  int init = Teuchos::getIntegralValue<int> (fdyn,"INITIALFIELD");
  fluidtimeparams.set ("eval err for analyt sol"   ,init);

  // ---------------------------- fine-scale subgrid viscosity approach
  fluidtimeparams.set<string> ("fs subgrid viscosity", fdyn.get<string>("FSSUGRVISC"));

  // -----------------------sublist containing stabilization parameters
  fluidtimeparams.sublist("STABILIZATION")=fdyn.sublist("STABILIZATION");

  // --------------------------sublist containing turbulence parameters
  {
    fluidtimeparams.sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");

    fluidtimeparams.sublist("TURBULENCE MODEL").set<string>("statistics outfile",DRT::Problem::Instance()->OutputControlFile()->FileName());
  }

  // ----------------------------------------------- XFEM related stuff
  {
    const Teuchos::ParameterList& xdyn = DRT::Problem::Instance()->XFEMGeneralParams();
    fluidtimeparams.sublist("XFEM").set<bool>("DLM_condensation", getIntegralValue<int>(xdyn,"DLM_CONDENSATION")==1 );
    fluidtimeparams.sublist("XFEM").set<bool>("CONDEST", getIntegralValue<int>(xdyn,"CONDEST")==1 );
    fluidtimeparams.sublist("XFEM").set<bool>("EXP_INTERSECTION", getIntegralValue<int>(xdyn,"EXP_INTERSECTION")==1 );
    fluidtimeparams.sublist("XFEM").set<double>("volumeRatioLimit", xdyn.get<double>("volumeRatioLimit"));
    fluidtimeparams.sublist("XFEM").set<double>("boundaryRatioLimit", xdyn.get<double>("boundaryRatioLimit"));
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

    vector<string> conditions_to_copy;
    conditions_to_copy.push_back("XFEMCoupling");
    Teuchos::RCP<DRT::Discretization> boundarydis = DRT::UTILS::CreateDiscretizationFromCondition(soliddis, "XFEMCoupling", "Boundary", "BELE3", conditions_to_copy);
    if (boundarydis->NumGlobalNodes() == 0 and fluiddis->Comm().MyPID() == 0)
      std::cout << "empty discretization detected. XFEMCoupling condition applied?" << endl;
    
    // create node and element distribution with elements and nodes ghosted on all processors
    const Epetra_Map noderowmap = *boundarydis->NodeRowMap();
    const Epetra_Map elemrowmap = *boundarydis->ElementRowMap();
    
    // put all boundary nodes and elements onto all processors
    const Epetra_Map nodecolmap = *LINALG::AllreduceEMap(noderowmap);
    const Epetra_Map elemcolmap = *LINALG::AllreduceEMap(elemrowmap);
    
    // redistribute nodes and elements to column (ghost) map
    boundarydis->ExportColumnNodes(nodecolmap);
    boundarydis->ExportColumnElements(elemcolmap);

    // Now we are done. :)
    const int err = boundarydis->FillComplete();
    if (err) dserror("FillComplete() returned err=%d",err);


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
    
    fluidtimeparams.set<bool>("do explicit predictor",true);

    //------------------------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor);
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom
    //------------------------------------------------------------------
    FLD::XFluidImplicitTimeInt fluidimplicit(
        fluiddis,
        solver,
        fluidtimeparams,
        fluidoutput);

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
        fluidimplicit.SetInitialFlowField(boundarydis,init,startfuncno);
      }
    }

    fluidtimeparams.set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

    // call time-integration (or stationary) scheme
    fluidimplicit.Integrate(boundarydis);

    // do result test if required
    DRT::ResultTestManager testmanager(fluiddis->Comm());
    testmanager.AddFieldTest(rcp(new FLD::XFluidResultTest(fluidimplicit)));
    testmanager.TestAll();

  }
  else
  {
    dserror("Unknown solver type for drt_xfluid");
  }

  return;

} // end of xdyn_fluid_drt()

#endif  // #ifdef CCADISCRET
