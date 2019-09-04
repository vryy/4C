/*----------------------------------------------------------------------*/
/*! \file

\brief base adapter of adjoint fluid equations for topology optimization

\level 2

\maintainer Martin Kronbichler

 *------------------------------------------------------------------------------------------------*/


#include "ad_opt_fluid_adjoint_base.H"
#include "../drt_lib/drt_periodicbc.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_opti/topopt_fluidAdjointImplTimeIntegration.H"
#include "../drt_opti/topopt_utils.H"
#include "../linalg/linalg_solver.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::TopOptFluidAdjointAlgorithm::TopOptFluidAdjointAlgorithm(
    const Teuchos::ParameterList& prbdyn)
{
  SetupAdjointFluid(prbdyn);
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::TopOptFluidAdjointAlgorithm::~TopOptFluidAdjointAlgorithm() {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::TopOptFluidAdjointAlgorithm::ReadRestart(int step)
{
  adjointTimeInt_->ReadRestart(step);
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::TopOptFluidAdjointAlgorithm::SetupAdjointFluid(const Teuchos::ParameterList& prbdyn)
{
  Teuchos::RCP<Teuchos::Time> t =
      Teuchos::TimeMonitor::getNewTimer("ADAPTER::TopOptFluidAdjointAlgorithm::SetupFluid");
  Teuchos::TimeMonitor monitor(*t);

  DRT::Problem* problem = DRT::Problem::Instance();

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> actdis = problem->GetDis("fluid");


  // -------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  // -------------------------------------------------------------------
  PeriodicBoundaryConditions pbc(actdis);
  pbc.UpdateDofsForPeriodicBoundaryConditions();

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->HaveDofs()) dserror("adjoint field solved after fluid field");

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  std::string filename = TOPOPT::modifyFilename(problem->OutputControlFile()->FileName(),
      "xxx_adjoint_", (bool)DRT::Problem::Instance()->Restart(), true);

  // output control for optimization field
  // equal to output for fluid equations except for the filename
  // and the - not necessary - input file name
  Teuchos::RCP<IO::OutputControl> adjointoutput =
      Teuchos::rcp(new IO::OutputControl(actdis->Comm(), problem->ProblemName(),
          problem->SpatialApproximation(), problem->OutputControlFile()->InputFileName(), filename,
          problem->NDim(), problem->Restart(), problem->OutputControlFile()->FileSteps(),
          DRT::INPUT::IntegralValue<int>(problem->IOParams(), "OUTPUT_BIN")));

  Teuchos::RCP<IO::DiscretizationWriter> output =
      Teuchos::rcp(new IO::DiscretizationWriter(actdis));
  output->SetOutput(adjointoutput);
  if (prbdyn.get<int>("NUM_OUTPUT_STEPS") != 1)  // more than one output step
  {
    output->NewResultFile(1);

    // remove not required control file
    filename = filename + ".control";
    std::remove(filename.c_str());
  }
  output->WriteMesh(0, 0.0);


  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  // const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  // const Teuchos::ParameterList& probsize    = DRT::Problem::Instance()->ProblemSizeParams();
  const Teuchos::ParameterList& fdyn = problem->FluidDynamicParams();
  const Teuchos::ParameterList& adjointfdyn =
      problem->OptimizationControlParams().sublist("TOPOLOGY ADJOINT FLUID");
  const Teuchos::ParameterList& opti = problem->OptimizationControlParams();

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // get the solver number used for linear fluid solver
  const int linsolvernumber = fdyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    dserror(
        "no linear solver defined for fluid problem. Please set LINEAR_SOLVER in FLUID DYNAMIC to "
        "a valid number!");
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(new LINALG::Solver(
      problem->SolverParams(linsolvernumber), actdis->Comm(), problem->ErrorFile()->Handle()));

  actdis->ComputeNullSpaceIfNecessary(solver->Params(), true);

  // -------------------------------------------------------------------
  // create a second solver for SIMPLER preconditioner if chosen from input
  // -------------------------------------------------------------------
  if (DRT::INPUT::IntegralValue<int>(fdyn, "SIMPLER"))
  {
    // add Inverse1 block for velocity dofs
    /*Teuchos::ParameterList& inv1 = solver->Params().sublist("Inverse1");
    inv1 = solver->Params();
    inv1.remove("SIMPLER",false); // not necessary
    inv1.remove("Inverse1",false);

    // get the solver number used for SIMPLER SOLVER
    const int linsolvernumber_simpler = fdyn.get<int>("SIMPLER_SOLVER");
    if (linsolvernumber_simpler == (-1))
      dserror("no SIMPLER_SOLVER number set for fluid problem solved with SIMPLER. Please set
    SIMPLER_SOLVER in FLUID DYNAMIC to a valid number!");
    // add Inverse2 block for pressure dofs
    solver->PutSolverParamsToSubParams("Inverse2",
    DRT::Problem::Instance()->SolverParams(linsolvernumber_simpler));
    // use CheapSIMPLE preconditioner (hardwired, change me for others)
    solver->Params().sublist("CheapSIMPLE Parameters").set("Prec Type","CheapSIMPLE");
    solver->Params().set("FLUID",true);*/
    dserror("Fix handling of SIMPLER with top opt adjoint fluid...");
  }

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> fluidadjointtimeparams =
      Teuchos::rcp(new Teuchos::ParameterList());

  fluidadjointtimeparams->set<int>(
      "Simple Preconditioner", DRT::INPUT::IntegralValue<int>(fdyn, "SIMPLER"));
  fluidadjointtimeparams->set<int>("AMG(BS) Preconditioner",
      DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(
          DRT::Problem::Instance()->SolverParams(linsolvernumber), "AZPREC"));

  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  const int ndim = DRT::Problem::Instance()->NDim();
  fluidadjointtimeparams->set<int>("number of velocity degrees of freedom", ndim);

  // -------------------------------------------------- time integration
  // note: here, the values are taken out of the problem-dependent ParameterList prbdyn
  // (which also can be fluiddyn itself!)

  // the default time step size
  fluidadjointtimeparams->set<double>("time step size", 1.0 * prbdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  fluidadjointtimeparams->set<double>("total time", prbdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  fluidadjointtimeparams->set<int>("max number timesteps", prbdyn.get<int>("NUMSTEP"));

  // ---------------------------------------------- nonlinear iteration
  // type of predictor
  fluidadjointtimeparams->set<std::string>("predictor", fdyn.get<std::string>("PREDICTOR"));
  // set linearisation scheme
  fluidadjointtimeparams->set<int>("Linearisation",
      DRT::INPUT::IntegralValue<INPAR::FLUID::LinearisationAction>(fdyn, "NONLINITER"));
  // maximum number of nonlinear iteration steps
  fluidadjointtimeparams->set<int>("max nonlin iter steps", fdyn.get<int>("ITEMAX"));
  // stop nonlinear iteration when both incr-norms are below this bound
  fluidadjointtimeparams->set<double>("tolerance for nonlin iter", fdyn.get<double>("CONVTOL"));
  // set convergence check
  fluidadjointtimeparams->set<std::string>("CONVCHECK", fdyn.get<std::string>("CONVCHECK"));
  // set adaptive linear solver tolerance

  // ---------------------------------------------- objective variables
  // set if objective contains dissipation
  fluidadjointtimeparams->set<INPAR::TOPOPT::ObjectiveDissipation>(
      "OBJECTIVE_DISSIPATION", DRT::INPUT::IntegralValue<INPAR::TOPOPT::ObjectiveDissipation>(
                                   opti, "OBJECTIVE_DISSIPATION"));
  // set if objective contains pressure drop
  fluidadjointtimeparams->set<bool>("OBJECTIVE_PRESSURE_DROP",
      DRT::INPUT::IntegralValue<int>(opti, "OBJECTIVE_PRESSURE_DROP") == 1);
  // set objective's dissipation factor
  fluidadjointtimeparams->set<double>("DISSIPATION_FAC", opti.get<double>("DISSIPATION_FAC"));
  // set objective's pressure drop factor
  fluidadjointtimeparams->set<double>("PRESSURE_DROP_FAC", opti.get<double>("PRESSURE_DROP_FAC"));

  // ----------------------------------------------- restart and output
  //  // restart
  fluidadjointtimeparams->set<int>("write restart every", prbdyn.get<int>("RESTARTEVRY"));
  // solution output
  fluidadjointtimeparams->set<int>("write solution every", prbdyn.get<int>("RESULTSEVRY"));
  // flag for writing fluid field to gmsh
  if (DRT::INPUT::IntegralValue<bool>(DRT::Problem::Instance()->IOParams(), "OUTPUT_GMSH") == false)
    fluidadjointtimeparams->set<bool>("GMSH_OUTPUT", false);
  else
    fluidadjointtimeparams->set<bool>(
        "GMSH_OUTPUT", DRT::INPUT::IntegralValue<bool>(fdyn, "GMSH_OUTPUT"));

  // ----------- initial field for test cases
  INPAR::TOPOPT::InitialAdjointField initfield =
      DRT::INPUT::IntegralValue<INPAR::TOPOPT::InitialAdjointField>(adjointfdyn, "INITIALFIELD");

  // ------------------------------------ special test case for adjoints
  fluidadjointtimeparams->set<INPAR::TOPOPT::AdjointCase>("special test case",
      DRT::INPUT::IntegralValue<INPAR::TOPOPT::AdjointCase>(adjointfdyn, "TESTCASE"));

  // ------------------------------------ test case for optimization routine
  fluidadjointtimeparams->set<double>(
      "theta_obj", prbdyn.sublist("TOPOLOGY OPTIMIZER").get<double>("THETA"));

  // ------------------------------------ test case for optimization routine
  fluidadjointtimeparams->set<int>(
      "opti testcase", DRT::INPUT::IntegralValue<INPAR::TOPOPT::OptiCase>(
                           prbdyn.sublist("TOPOLOGY OPTIMIZER"), "TESTCASE"));

  // -------------------------- potential nonlinear boundary conditions
  fluidadjointtimeparams->set<std::string>(
      "Nonlinear boundary conditions", fdyn.get<std::string>("NONLINEARBC"));

  // -----------------------sublist containing stabilization parameters
  fluidadjointtimeparams->sublist("RESIDUAL-BASED STABILIZATION") =
      fdyn.sublist("RESIDUAL-BASED STABILIZATION");

  // -------------------------------------------------------------------
  // additional parameters and algorithm call depending on respective
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  INPAR::FLUID::TimeIntegrationScheme timeint =
      DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fdyn, "TIMEINTEGR");

  // -------------------------------------------------------------------
  // additional parameters and algorithm call depending on respective
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  if (timeint == INPAR::FLUID::timeint_stationary or
      timeint == INPAR::FLUID::timeint_one_step_theta)
  {
    // -----------------------------------------------------------------
    // set additional parameters in list for
    // one-step-theta/BDF2/af-generalized-alpha/stationary scheme
    // -----------------------------------------------------------------
    // type of adjoint equations
    fluidadjointtimeparams->set<INPAR::TOPOPT::AdjointType>("adjoint type",
        DRT::INPUT::IntegralValue<INPAR::TOPOPT::AdjointType>(adjointfdyn, "ADJOINT_TYPE"));
    // type of time-integration (or stationary) scheme
    fluidadjointtimeparams->set<int>("time int algo", timeint);
    // parameter theta for time-integration schemes
    fluidadjointtimeparams->set<double>("theta", fdyn.get<double>("THETA"));
    // number of steps for potential start algorithm
    fluidadjointtimeparams->set<int>("number of start steps", fdyn.get<int>("NUMSTASTEPS"));
    // parameter theta for potential start algorithm
    fluidadjointtimeparams->set<double>("start theta", fdyn.get<double>("START_THETA"));

    fluidadjointtimeparams->set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());

    // check parameter
    if (DRT::INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn, "PHYSICAL_TYPE") !=
        INPAR::FLUID::topopt)
      dserror("This physical type is not supported for fluid topology optimization");

    if ((fluidadjointtimeparams->get<INPAR::TOPOPT::ObjectiveDissipation>(
             "OBJECTIVE_DISSIPATION") != INPAR::TOPOPT::obj_diss_no) and
        (fluidadjointtimeparams->get<double>("DISSIPATION_FAC") <= 0.0))
      dserror("Optimizing dissipation with negative dissipation coefficient");

    if ((fluidadjointtimeparams->get<bool>("OBJECTIVE_PRESSURE_DROP") == true) and
        (fluidadjointtimeparams->get<double>("PRESSURE_DROP_FAC") <= 0.0))
      dserror("Optimizing pressure drop with negative coefficient");

    //------------------------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor);
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom
    adjointTimeInt_ = Teuchos::rcp(
        new TOPOPT::ADJOINT::ImplicitTimeInt(actdis, solver, fluidadjointtimeparams, output));
  }
  else
  {
    dserror("not implemented for adjoint field");
  }

  // set initial field by given function
  // we do this here, since we have direct access to all necessary parameters
  int startfuncno = adjointfdyn.get<int>("INITFUNCNO");
  if (initfield != INPAR::TOPOPT::initadjointfield_field_by_function) startfuncno = -1;

  adjointTimeInt_->SetInitialAdjointField(initfield, startfuncno);

  if ((timeint != INPAR::FLUID::timeint_stationary) and
      (fluidadjointtimeparams->get<INPAR::TOPOPT::AdjointType>("adjoint type") !=
          INPAR::TOPOPT::discrete_adjoint))
    adjointTimeInt_->Output();

  return;
}
