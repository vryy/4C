/*----------------------------------------------------------------------------*/
/*! \file

\brief General framework for monolithic fsi solution schemes

\maintainer Matthias Mayr

\level 1
*/
/*----------------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Direction_UserDefinedFactory.H>

#include "fsi_monolithic.H"
#include "fsi_debugwriter.H"
#include "fsi_nox_group.H"
#include "fsi_nox_newton.H"
#include "fsi_statustest.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"

#include "../drt_adapter/ad_ale_fsi.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_fld_fluid_fsi.H"
#include "../drt_adapter/ad_ale.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_str_fsi_timint_adaptive.H"

#include "../drt_constraint/constraint_manager.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_structure/stru_aux.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_ale/ale_utils_mapextractor.H"

#include "fsi_overlapprec.H"
#include "fsi_overlapprec_fsiamg.H"
#include "fsi_overlapprec_amgnxn.H"
#include "fsi_overlapprec_hybrid.H"

/*----------------------------------------------------------------------------*/
/* Note: The order of calling the three BaseAlgorithm-constructors is
 * important here! In here control file entries are written. And these
 * entries define the order in which the filters handle the
 * Discretizations, which in turn defines the dof number ordering of the
 * Discretizations.
 */
/*----------------------------------------------------------------------------*/
FSI::MonolithicBase::MonolithicBase(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : AlgorithmBase(comm, timeparams),
      isadastructure_(false),
      isadafluid_(false),
      isadasolver_(false),
      verbosity_(DRT::INPUT::IntegralValue<INPAR::FSI::Verbosity>(
          DRT::Problem::Instance()->FSIDynamicParams(), "VERBOSITY"))
{
  // access the discretizations
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
  Teuchos::RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->GetDis("fluid");
  Teuchos::RCP<DRT::Discretization> aledis = DRT::Problem::Instance()->GetDis("ale");

  CreateStructureTimeIntegrator(timeparams, structdis);
  CreateFluidAndALETimeIntegrator(timeparams, fluiddis, aledis);

  coupsf_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupsa_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupfa_ = Teuchos::rcp(new ADAPTER::Coupling());
  icoupfa_ = Teuchos::rcp(new ADAPTER::Coupling());
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
FSI::MonolithicBase::~MonolithicBase() {}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MonolithicBase::ReadRestart(int step)
{
  StructureField()->ReadRestart(step);
  FluidField()->ReadRestart(step);
  AleField()->ReadRestart(step);

  SetTimeStep(FluidField()->Time(), FluidField()->Step());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MonolithicBase::CreateStructureTimeIntegrator(
    const Teuchos::ParameterList& timeparams, Teuchos::RCP<DRT::Discretization> structdis)
{
  // delete deprecated time integrator
  structure_ = Teuchos::null;

  // access structural dynamic params
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // ask base algorithm for the structural time integrator
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure =
      Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(
          timeparams, const_cast<Teuchos::ParameterList&>(sdyn), structdis));
  structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::FSIStructureWrapper>(structure->StructureField());
  structure_->Setup();

  if (structure_ == Teuchos::null)
    dserror(
        "Cast from ADAPTER::Structure to ADAPTER::FSIStructureWrapper "
        "failed.");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MonolithicBase::CreateFluidAndALETimeIntegrator(const Teuchos::ParameterList& timeparams,
    Teuchos::RCP<DRT::Discretization> fluiddis, Teuchos::RCP<DRT::Discretization> aledis)
{
  // delete deprecated time integrators
  fluid_ = Teuchos::null;
  ale_ = Teuchos::null;

  // ask base algorithm for the fluid time integrator
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid = Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(
      timeparams, DRT::Problem::Instance()->FluidDynamicParams(), "fluid", true));
  fluid_ = Teuchos::rcp_dynamic_cast<ADAPTER::FluidFSI>(fluid->FluidField());

  if (fluid_ == Teuchos::null) dserror("Cast from ADAPTER::Fluid to ADAPTER::FluidFSI failed");

  // ask base algorithm for the ale time integrator
  Teuchos::RCP<ADAPTER::AleBaseAlgorithm> ale =
      Teuchos::rcp(new ADAPTER::AleBaseAlgorithm(timeparams, aledis));
  ale_ = Teuchos::rcp_dynamic_cast<ADAPTER::AleFsiWrapper>(ale->AleField());

  if (ale_ == Teuchos::null) dserror("Cast from ADAPTER::Ale to ADAPTER::AleFsiWrapper failed");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MonolithicBase::PrepareTimeStep()
{
  IncrementTimeAndStep();
  if (verbosity_ >= INPAR::FSI::verbosity_low) PrintHeader();
  PrepareTimeStepPreconditioner();
  PrepareTimeStepFields();

  // Note: it's important to first prepare the single fields and than the fsi problem
  PrepareTimeStepFSI();

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MonolithicBase::PrepareTimeStepFSI()
{
  ddgpred_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->ExtractInterfaceDispnp()));
  ddgpred_->Update(-1.0, *StructureField()->ExtractInterfaceDispn(), 1.0);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MonolithicBase::PrepareTimeStepFields()
{
  StructureField()->PrepareTimeStep();
  FluidField()->PrepareTimeStep();
  AleField()->PrepareTimeStep();

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MonolithicBase::PrepareOutput() { StructureField()->PrepareOutput(); }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MonolithicBase::Output()
{
  /* Note: The order is important here! In here control file entries are
   * written. And these entries define the order in which the filters handle
   * the Discretizations, which in turn defines the dof number ordering of the
   * Discretizations.
   */
  StructureField()->Output();
  FluidField()->Output();
  AleField()->Output();

  if (StructureField()->GetConstraintManager()->HaveMonitor())
  {
    StructureField()->GetConstraintManager()->ComputeMonitorValues(StructureField()->Dispnp());
    if (Comm().MyPID() == 0) StructureField()->GetConstraintManager()->PrintMonitorValues();
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::StructToAle(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsa_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::AleToStruct(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::StructToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::FluidToStruct(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::AleToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupfa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::FluidToAleInterface(
    Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupfa_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::AleToFluidInterface(
    Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupfa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::StructToAle(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::AleToStruct(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::StructToFluid(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::FluidToStruct(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::AleToFluid(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::FluidToAleInterface(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupfa_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::AleToFluidInterface(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupfa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
FSI::Monolithic::Monolithic(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : MonolithicBase(comm, timeparams),
      firstcall_(true),
      noxiter_(0),
      erroraction_(erroraction_stop),
      log_(Teuchos::null),
      logada_(Teuchos::null)
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

  // enable debugging
  if (DRT::INPUT::IntegralValue<int>(fsidyn, "DEBUGOUTPUT") == 1)
  {
    sdbg_ = Teuchos::rcp(new UTILS::DebugWriter(StructureField()->Discretization()));
    // fdbg_ = Teuchos::rcp(new UTILS::DebugWriter(FluidField()->Discretization()));
  }

  // write iterations-file
  std::string fileiter = DRT::Problem::Instance()->OutputControlFile()->FileName();
  fileiter.append(".iteration");
  log_ = Teuchos::rcp(new std::ofstream(fileiter.c_str()));

  // write energy-file
  if (DRT::INPUT::IntegralValue<int>(fsidyn.sublist("MONOLITHIC SOLVER"), "ENERGYFILE") == 1)
  {
    std::string fileiter2 = DRT::Problem::Instance()->OutputControlFile()->FileName();
    fileiter2.append(".fsienergy");
    logenergy_ = Teuchos::rcp(new std::ofstream(fileiter2.c_str()));
  }

  // "Initialize" interface solution increments due to structural predictor
  ddgpred_ = Teuchos::null;

  //-------------------------------------------------------------------------
  // time step size adaptivity
  //-------------------------------------------------------------------------
  const bool timeadapton =
      DRT::INPUT::IntegralValue<bool>(fsidyn.sublist("TIMEADAPTIVITY"), "TIMEADAPTON");

  if (timeadapton)
  {
    InitTimIntAda(fsidyn);
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::SetupSystem()
{
  // right now we use matching meshes at the interface

  const int ndim = DRT::Problem::Instance()->NDim();

  ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  ADAPTER::Coupling& coupsa = StructureAleCoupling();
  ADAPTER::Coupling& coupfa = FluidAleCoupling();
  ADAPTER::Coupling& icoupfa = InterfaceFluidAleCoupling();

  // structure to fluid

  coupsf.SetupConditionCoupling(*StructureField()->Discretization(),
      StructureField()->Interface()->FSICondMap(), *FluidField()->Discretization(),
      FluidField()->Interface()->FSICondMap(), "FSICoupling", ndim);

  // structure to ale

  coupsa.SetupConditionCoupling(*StructureField()->Discretization(),
      StructureField()->Interface()->FSICondMap(), *AleField()->Discretization(),
      AleField()->Interface()->FSICondMap(), "FSICoupling", ndim);

  // fluid to ale at the interface

  icoupfa.SetupConditionCoupling(*FluidField()->Discretization(),
      FluidField()->Interface()->FSICondMap(), *AleField()->Discretization(),
      AleField()->Interface()->FSICondMap(), "FSICoupling", ndim);

  // In the following we assume that both couplings find the same dof
  // map at the structural side. This enables us to use just one
  // interface dof map for all fields and have just one transfer
  // operator from the interface map to the full field map.
  if (not coupsf.MasterDofMap()->SameAs(*coupsa.MasterDofMap()))
    dserror("structure interface dof maps do not match");

  if (coupsf.MasterDofMap()->NumGlobalElements() == 0)
    dserror("No nodes in matching FSI interface. Empty FSI coupling condition?");

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = FluidField()->Discretization()->NodeRowMap();
  const Epetra_Map* alenodemap = AleField()->Discretization()->NodeRowMap();

  coupfa.SetupCoupling(*FluidField()->Discretization(), *AleField()->Discretization(),
      *fluidnodemap, *alenodemap, ndim);

  FluidField()->SetMeshMap(coupfa.MasterDofMap());

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const bool timeadapton =
      DRT::INPUT::IntegralValue<bool>(fsidyn.sublist("TIMEADAPTIVITY"), "TIMEADAPTON");

  // Run time loop with constant or adaptive time step size (depending on the user's will)
  if (not timeadapton)
  {
    // call time loop with constant time step size
    TimeloopConstDt(interface);
  }
  else
  {
    // call time loop with adaptive time step size
    TimeloopAdaDt(interface);
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::TimeloopConstDt(
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
{
  PrepareTimeloop();

  while (NotFinished())
  {
    PrepareTimeStep();
    TimeStep(interface);
    PrepareOutput();
    Update();
    Output();
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::PrepareTimeloop()
{
  // make sure we didn't destroy the maps before we entered the timeloop
  Extractor().CheckForValidMapExtractor();

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = NOXParameterList();

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", Comm().MyPID());

  printParams.set("Output Information",
      NOX::Utils::Error | NOX::Utils::Warning | NOX::Utils::OuterIteration |
          NOX::Utils::InnerIteration |
          // NOX::Utils::Parameters |
          NOX::Utils::Details | NOX::Utils::OuterIterationStatusTest |
          NOX::Utils::LinearSolverDetails | NOX::Utils::TestDetails | NOX::Utils::StepperIteration |
          NOX::Utils::StepperDetails | NOX::Utils::StepperParameters | NOX::Utils::Debug | 0);

  // Create printing utilities
  utils_ = Teuchos::rcp(new NOX::Utils(printParams));

  // write header of log-file
  if (Comm().MyPID() == 0)
  {
    (*log_) << "# num procs      = " << Comm().NumProc() << "\n"
            << "# Method         = " << nlParams.sublist("Direction").get<std::string>("Method")
            << std::endl
            << std::right << std::setw(9) << "# step" << std::right << std::setw(16) << "time"
            << std::right << std::setw(16) << "time/step" << std::right << std::setw(16)
            << "#nliter" << std::right << std::setw(16) << "res-norm" << std::right << std::setw(16)
            << "#liter" << std::right << std::setw(16) << "dt" << std::endl;

    (*log_) << "#\n\n";
  }

  WriteAdaFileHeader();

  // write header of energy-file
  if (Comm().MyPID() == 0 and (not logenergy_.is_null()))
  {
    (*logenergy_) << "# Artificial interface energy due to temporal discretization\n"
                  << "# num procs      = " << Comm().NumProc() << "\n"
                  << "# Method         = "
                  << nlParams.sublist("Direction").get<std::string>("Method") << std::endl
                  << std::right << std::setw(9) << "# step" << std::right << std::setw(16) << "time"
                  << std::right << std::setw(16) << "energy/step" << std::right << std::setw(16)
                  << "sum_of_energy" << std::endl;

    (*logenergy_) << "#\n\n";
  }


  // check for prestressing,
  // do not allow monolithic in the pre-phase
  // allow monolithic in the post-phase
  {
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    INPAR::STR::PreStress pstype =
        Teuchos::getIntegralValue<INPAR::STR::PreStress>(sdyn, "PRESTRESS");
    if (pstype != INPAR::STR::PreStress::none)
    {
      const double pstime = sdyn.get<double>("PRESTRESSTIME");

      if (Time() + Dt() <= pstime)
        dserror("No monolithic FSI in the pre-phase of prestressing, use Aitken!");
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::TimeStep(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::TimeStep");

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = NOXParameterList();

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  // Teuchos::ParameterList& solverOptions = nlParams.sublist("Solver Options");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  // Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", Comm().MyPID());

  switch (verbosity_)
  {
    case INPAR::FSI::verbosity_full:
    {
      printParams.set("Output Information",
          NOX::Utils::Error | NOX::Utils::Warning | NOX::Utils::OuterIteration |
              NOX::Utils::InnerIteration |
              // NOX::Utils::Parameters |
              NOX::Utils::Details |  // weg damit!
              NOX::Utils::OuterIterationStatusTest |
              NOX::Utils::LinearSolverDetails |  // weg damit!
              NOX::Utils::TestDetails | NOX::Utils::StepperIteration | NOX::Utils::StepperDetails |
              NOX::Utils::StepperParameters | NOX::Utils::Debug | 0);
      break;
    }
    case INPAR::FSI::verbosity_medium:
    {
      printParams.set("Output Information",
          NOX::Utils::Error | NOX::Utils::Warning | NOX::Utils::OuterIteration |
              NOX::Utils::InnerIteration |
              // NOX::Utils::Parameters |
              // NOX::Utils::Details | //weg damit!
              NOX::Utils::OuterIterationStatusTest |
              NOX::Utils::LinearSolverDetails |  // weg damit!
              NOX::Utils::TestDetails | NOX::Utils::StepperIteration | NOX::Utils::StepperDetails |
              NOX::Utils::StepperParameters | NOX::Utils::Debug | 0);
      break;
    }
    case INPAR::FSI::verbosity_low:
    case INPAR::FSI::verbosity_subproblem:
    {
      printParams.set(
          "Output Information", NOX::Utils::Error | NOX::Utils::Warning |
                                    //                    NOX::Utils::OuterIteration |
                                    //                    NOX::Utils::InnerIteration |
                                    //                    //NOX::Utils::Parameters |
                                    //  //                  NOX::Utils::Details |
                                    NOX::Utils::OuterIterationStatusTest |
                                    //  //                  NOX::Utils::LinearSolverDetails |
                                    //                    NOX::Utils::TestDetails |
                                    //                    NOX::Utils::StepperIteration |
                                    //                    NOX::Utils::StepperDetails |
                                    //                    NOX::Utils::StepperParameters |
                                    NOX::Utils::Debug | 0);
      break;
    }
    default:
    {
      dserror("Verbosity level not supported!");
      break;
    }
  }

  Teuchos::Time timer("time step timer");

  if (sdbg_ != Teuchos::null) sdbg_->NewTimeStep(Step(), "struct");
  if (fdbg_ != Teuchos::null) fdbg_->NewTimeStep(Step(), "fluid");

  // Single field predictors have been applied, so store the structural
  // interface displacement increment due to predictor or inhomogeneous
  // Dirichlet boundary conditions

  // start time measurement
  Teuchos::RCP<Teuchos::TimeMonitor> timemonitor =
      Teuchos::rcp(new Teuchos::TimeMonitor(timer, true));

  // calculate initial linear system at current position (no increment)
  // This initializes the field algorithms and creates the first linear
  // systems. And this is the reason we know the initial linear system is
  // there when we create the NOX::FSI::Group.
  Evaluate(Teuchos::null);

  // Get initial guess.
  // The initial system is there, so we can happily extract the
  // initial guess. (The Dirichlet conditions are already build in!)
  Teuchos::RCP<Epetra_Vector> initial_guess = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));
  InitialGuess(initial_guess);

  NOX::Epetra::Vector noxSoln(initial_guess, NOX::Epetra::Vector::CreateView);

  // Create the linear system
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys = CreateLinearSystem(nlParams, noxSoln, utils_);

  // Create the Group
  Teuchos::RCP<NOX::FSI::Group> grp =
      Teuchos::rcp(new NOX::FSI::Group(*this, printParams, interface, noxSoln, linSys));

  // Convergence Tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo = CreateStatusTest(nlParams, grp);

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
      NOX::Solver::buildSolver(grp, combo, Teuchos::RCP<Teuchos::ParameterList>(&nlParams, false));

  // we know we already have the first linear system calculated
  grp->CaptureSystemState();

  // solve the whole thing
  noxstatus_ = solver->solve();
  noxiter_ = solver->getNumIterations();

  // recover Lagrange multiplier \lambda_{\Gamma} at the interface at the end of each time step
  // (i.e. condensed traction/forces onto the structure) needed for rhs in next time step
  RecoverLagrangeMultiplier();

  // compute spurious interface energy increment due to temporal discretization
  CalculateInterfaceEnergyIncrement();

  // stop time measurement
  timemonitor = Teuchos::null;

  if (Comm().MyPID() == 0)
  {
    (*log_) << std::right << std::setw(9) << Step() << std::right << std::setw(16) << Time()
            << std::right << std::setw(16) << timer.totalElapsedTime() << std::right
            << std::setw(16) << nlParams.sublist("Output").get<int>("Nonlinear Iterations")
            << std::right << std::setw(16)
            << nlParams.sublist("Output").get<double>("2-Norm of Residual") << std::right
            << std::setw(16)
            << lsParams.sublist("Output").get<int>("Total Number of Linear Iterations")
            << std::right << std::setw(16) << Dt();

    (*log_) << std::endl;
  }
  lsParams.sublist("Output").set("Total Number of Linear Iterations", 0);

  // perform the error check to determine the error action to be performed
  NonLinErrorCheck();

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::Update()
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  bool timeadapton =
      DRT::INPUT::IntegralValue<int>(fsidyn.sublist("TIMEADAPTIVITY"), "TIMEADAPTON");

  if (not timeadapton)
    StructureField()->Update();  // constant dt
  else
  {
    StructureField()->Update(Time());  // variable/adaptive dt

    if (IsAdaStructure())
      Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)
          ->UpdateStepSize(Dt());
  }

  FluidField()->Update();
  AleField()->Update();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::NonLinErrorCheck()
{
  // assume convergence of nonlinear solver
  erroraction_ = erroraction_none;

  // if everything is fine, then return right now
  if (NoxStatus() == NOX::StatusTest::Converged)
  {
    return;
  }

  // The nonlinear solver did not converge. Thus, we have to take some action
  // that depends on the user's will given in the input file

  // get the FSI parameter list
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

  // get the user's will
  const INPAR::FSI::DivContAct divcontype = DRT::INPUT::IntegralValue<INPAR::FSI::DivContAct>(
      fsidyn.sublist("TIMEADAPTIVITY"), ("DIVERCONT"));

  if (NoxStatus() != NOX::StatusTest::Converged)
  {
    switch (divcontype)
    {
      case INPAR::FSI::divcont_stop:
      {
        // set the corresponding error action
        erroraction_ = erroraction_stop;

        // stop the simulation
        dserror("Nonlinear solver did not converge in %i iterations in time step %i.", noxiter_,
            Step());
        break;
      }
      case INPAR::FSI::divcont_continue:
      {
        // set the corresponding error action
        erroraction_ = erroraction_continue;

        // Notify user about non-converged nonlinear solver, but do not abort the simulation
        if (Comm().MyPID() == 0)
        {
          IO::cout << "\n*** Nonlinear solver did not converge in " << noxiter_
                   << " iterations in time step " << Step() << ". Continue ..." << IO::endl;
        }
        break;
      }
      case INPAR::FSI::divcont_halve_step:
      {
        // set the corresponding error action
        erroraction_ = erroraction_halve_step;

        const bool timeadapton =
            DRT::INPUT::IntegralValue<bool>(fsidyn.sublist("TIMEADAPTIVITY"), "TIMEADAPTON");
        if (not timeadapton)
        {
          dserror(
              "Nonlinear solver wants to halve the time step size. This is "
              "not possible in a time integrator with constant Delta t.");
        }

        // Notify user about non-converged nonlinear solver, but do not abort the simulation
        if (Comm().MyPID() == 0)
        {
          IO::cout << IO::endl
                   << "*** Nonlinear solver did not converge in " << noxiter_
                   << " iterations. Halve the time step size." << IO::endl;
        }
        break;
      }
      case INPAR::FSI::divcont_revert_dt:
      {
        // set the corresponding error action
        erroraction_ = erroraction_revert_dt;

        const bool timeadapton =
            DRT::INPUT::IntegralValue<bool>(fsidyn.sublist("TIMEADAPTIVITY"), "TIMEADAPTON");
        if (not timeadapton)
        {
          dserror(
              "Nonlinear solver wants to revert the time step size. This is "
              "not possible in a time integrator with constant Delta t.");
        }

        // Notify user about non-converged nonlinear solver, but do not abort the simulation
        if (Comm().MyPID() == 0)
        {
          IO::cout << IO::endl
                   << "*** Nonlinear solver did not converge in " << noxiter_
                   << " iterations. Revert the time step size to the previous one." << IO::endl;
        }
        break;
      }
      default:
      {
        dserror("Unknown action to cope with non-converged nonlinear solver.");
        break;
      }
    }
  }
  else
  {
    dserror("Unknown NOX::StatusTest::StatusType.");
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::Evaluate(Teuchos::RCP<const Epetra_Vector> x)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::Evaluate");

#ifdef DEBUG
  // check whether all fields have the same time step size
  CheckIfDtsSame();
#endif

  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;
  Teuchos::RCP<const Epetra_Vector> ax;

  if (x != Teuchos::null)
  {
    ExtractFieldVectors(x, sx, fx, ax);

    if (sdbg_ != Teuchos::null)
    {
      sdbg_->NewIteration();
      sdbg_->WriteVector("x", *StructureField()->Interface()->ExtractFSICondVector(sx));
    }
  }

  // Call all elements and assemble rhs and matrices
  // We only need the rhs here because NOX will ask for the rhs
  // only. But the Jacobian is stored internally and will be returned
  // later on without looking at x again!

  if (verbosity_ >= INPAR::FSI::verbosity_medium) Utils()->out() << "\nEvaluate elements\n";

  {
    Epetra_Time ts(Comm());
    StructureField()->Evaluate(sx);
    if (verbosity_ >= INPAR::FSI::verbosity_medium)
      Utils()->out() << "structure: " << ts.ElapsedTime() << " sec\n";
  }

  {
    Epetra_Time ta(Comm());
    AleField()->Evaluate(ax);
    if (verbosity_ >= INPAR::FSI::verbosity_medium)
      Utils()->out() << "ale      : " << ta.ElapsedTime() << " sec\n";
  }

  // transfer the current ale mesh positions to the fluid field
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluid(AleField()->Dispnp());
  FluidField()->ApplyMeshDisplacement(fluiddisp);

  {
    Epetra_Time tf(Comm());
    FluidField()->Evaluate(fx);
    if (verbosity_ >= INPAR::FSI::verbosity_medium)
      Utils()->out() << "fluid    : " << tf.ElapsedTime() << " sec\n";
  }

  if (verbosity_ >= INPAR::FSI::verbosity_medium) Utils()->out() << "\n";
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map>>& maps)
{
  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);
  blockrowdofmap_.Setup(*fullmap, maps);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::SetDefaultParameters(
    const Teuchos::ParameterList& fsidyn, Teuchos::ParameterList& list)
{
  // monolithic solver settings
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = list;

  nlParams.set<std::string>("Nonlinear Solver", "Line Search Based");
  // nlParams.set("Preconditioner", "None");
  // nlParams.set("Norm abs F", fsimono.get<double>("CONVTOL"));

  nlParams.set("Max Iterations", fsimono.get<int>("ITEMAX"));
  // nlParams.set("Max Iterations", 1);

  // ToDo: Remove the CONVTOL-tolerances and replace them by the single field
  //       tolerances given just below also for lung- and constraint-FSI
  //
  // Currently, we have to keep the parameter CONVTOL for lung and constraint FSI
  nlParams.set("Norm abs pres", fsimono.get<double>("CONVTOL"));
  nlParams.set("Norm abs vel", fsimono.get<double>("CONVTOL"));
  nlParams.set("Norm abs disp", fsimono.get<double>("CONVTOL"));

  // set tolerances for nonlinear solver
  nlParams.set("Tol dis res L2", fsimono.get<double>("TOL_DIS_RES_L2"));
  nlParams.set("Tol dis res Inf", fsimono.get<double>("TOL_DIS_RES_INF"));
  nlParams.set("Tol dis inc L2", fsimono.get<double>("TOL_DIS_INC_L2"));
  nlParams.set("Tol dis inc Inf", fsimono.get<double>("TOL_DIS_INC_INF"));
  nlParams.set("Tol fsi res L2", fsimono.get<double>("TOL_FSI_RES_L2"));
  nlParams.set("Tol fsi res Inf", fsimono.get<double>("TOL_FSI_RES_INF"));
  nlParams.set("Tol fsi inc L2", fsimono.get<double>("TOL_FSI_INC_L2"));
  nlParams.set("Tol fsi inc Inf", fsimono.get<double>("TOL_FSI_INC_INF"));
  nlParams.set("Tol pre res L2", fsimono.get<double>("TOL_PRE_RES_L2"));
  nlParams.set("Tol pre res Inf", fsimono.get<double>("TOL_PRE_RES_INF"));
  nlParams.set("Tol pre inc L2", fsimono.get<double>("TOL_PRE_INC_L2"));
  nlParams.set("Tol pre inc Inf", fsimono.get<double>("TOL_PRE_INC_INF"));
  nlParams.set("Tol vel res L2", fsimono.get<double>("TOL_VEL_RES_L2"));
  nlParams.set("Tol vel res Inf", fsimono.get<double>("TOL_VEL_RES_INF"));
  nlParams.set("Tol vel inc L2", fsimono.get<double>("TOL_VEL_INC_L2"));
  nlParams.set("Tol vel inc Inf", fsimono.get<double>("TOL_VEL_INC_INF"));

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& solverOptions = nlParams.sublist("Solver Options");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  // Teuchos::ParameterList& lineSearchParams = nlParams.sublist("Line Search");

  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  dirParams.set<std::string>("Method", "User Defined");
  Teuchos::RCP<NOX::Direction::UserDefinedFactory> newtonfactory = Teuchos::rcp(this, false);
  dirParams.set("User Defined Direction Factory", newtonfactory);


  // status tests are expensive, but instructive
  // solverOptions.set<std::string>("Status Test Check Type","Minimal");
  solverOptions.set<std::string>("Status Test Check Type", "Complete");

  // be explicit about linear solver parameters
  lsParams.set<std::string>("Aztec Solver", "GMRES");
  // lsParams.set<std::string>("BiCGStab","GMRES");
  lsParams.set<std::string>("Orthogonalization", "Modified");

  // "r0", "rhs", "norm", "no scaling", "sol"
  lsParams.set<std::string>("Convergence Test", "r0");

  lsParams.set<int>("Size of Krylov Subspace", fsimono.get<int>("KRYLOV_SIZE"));
  lsParams.set<int>("Max Iterations", fsimono.get<int>("KRYLOV_ITEMAX"));
  lsParams.set<std::string>("Preconditioner", "User Defined");
  lsParams.set<int>("Output Frequency", 10);
  lsParams.set<bool>("Output Solver Details", true);

  // adaptive tolerance settings for linear solver
  lsParams.set<double>("base tolerance", fsimono.get<double>("BASETOL"));  // relative tolerance
  lsParams.set<double>(
      "adaptive distance", fsimono.get<double>("ADAPTIVEDIST"));  // adaptive distance
  lsParams.set<INPAR::FSI::Verbosity>("verbosity", verbosity_);  // verbosity level of FSI algorithm
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Direction::Generic> FSI::Monolithic::buildDirection(
    const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params) const
{
  Teuchos::RCP<NOX::FSI::Newton> newton = Teuchos::rcp(new NOX::FSI::Newton(gd, params));
  for (unsigned i = 0; i < statustests_.size(); ++i)
  {
    statustests_[i]->SetNewton(newton);
  }
  return newton;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool FSI::Monolithic::computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::computeF");
  Evaluate(Teuchos::rcp(&x, false));
  SetupRHS(F);
  return true;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool FSI::Monolithic::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) { return true; }


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool FSI::Monolithic::computePreconditioner(
    const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams)
{
  return true;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::SetupRHS(Epetra_Vector& f, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::SetupRHS");

  firstcall_ = firstcall;

  // We want to add into a zero vector
  f.PutScalar(0.0);

  // contributions of single field residuals
  SetupRHSResidual(f);

  // contributions of Lagrange multiplier from last time step
  SetupRHSLambda(f);

  // contributions of special "first nonlinear iteration" terms
  if (firstcall) SetupRHSFirstiter(f);

  if (dbcmaps_ !=
      Teuchos::null)  // ToDo: Remove if-"Abfrage" after 'dbcmaps_' has been introduced in lung fsi
  {
    // Finally, we take care of Dirichlet boundary conditions
    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(f));
    Teuchos::RCP<const Epetra_Vector> zeros = Teuchos::rcp(new const Epetra_Vector(f.Map(), true));
    LINALG::ApplyDirichlettoSystem(rhs, zeros, *(dbcmaps_->CondMap()));
    f.Update(1.0, *rhs, 0.0);
  }

  // NOX expects the 'positive' residual. The negative sign for the
  // linearized Newton system J*dx=-r is done internally by NOX.
  // Since we assembled the right hand side, we have to invert the sign here.
  f.Scale(-1.);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::InitialGuess");

  CombineFieldVectors(*ig, StructureField()->InitialGuess(), FluidField()->InitialGuess(),
      AleField()->InitialGuess(), true);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::CombineFieldVectors(Epetra_Vector& v, Teuchos::RCP<const Epetra_Vector> sv,
    Teuchos::RCP<const Epetra_Vector> fv, Teuchos::RCP<const Epetra_Vector> av)
{
  Extractor().AddVector(*sv, 0, v);
  Extractor().AddVector(*fv, 1, v);
  Extractor().AddVector(*av, 2, v);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::WriteInterfaceEnergyFile(const double energystep, const double energysum)
{
  // write to energy file
  if (Comm().MyPID() == 0 and (not logenergy_.is_null()))
  {
    (*logenergy_) << std::right << std::setw(9) << Step() << std::right << std::setw(16) << Time()
                  << std::right << std::setw(16) << energystep << std::right << std::setw(16)
                  << energysum;

    (*logenergy_) << std::endl;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
FSI::BlockMonolithic::BlockMonolithic(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : Monolithic(comm, timeparams),
      precondreusecount_(0),
      timeparams_(timeparams),
      interfaceprocs_(0)
{
  // interfaceprocs_.push_back(-1);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool FSI::BlockMonolithic::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::BlockMonolithic::computeJacobian");
  Evaluate(Teuchos::rcp(&x, false));
  LINALG::BlockSparseMatrixBase& mat = Teuchos::dyn_cast<LINALG::BlockSparseMatrixBase>(Jac);
  SetupSystemMatrix(mat);
  return true;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool FSI::BlockMonolithic::computePreconditioner(
    const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::BlockMonolithic::computePreconditioner");

  if (precondreusecount_ <= 0)
  {
    // Create preconditioner operator. The blocks are already there. This is
    // the perfect place to initialize the block preconditioners.
    SystemMatrix()->SetupPreconditioner();

    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
    precondreusecount_ = fsimono.get<int>("PRECONDREUSE");
  }

  precondreusecount_ -= 1;

  if (pcdbg_ != Teuchos::null)
  {
    pcdbg_->NewLinearSystem();
  }

  return true;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::BlockMonolithic::PrepareTimeStepPreconditioner()
{
  const Teuchos::ParameterList& fsimono =
      DRT::Problem::Instance()->FSIDynamicParams().sublist("MONOLITHIC SOLVER");

  if (DRT::INPUT::IntegralValue<int>(fsimono, "REBUILDPRECEVERYSTEP")) precondreusecount_ = 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::BlockMonolithic::CreateSystemMatrix(
    Teuchos::RCP<FSI::OverlappingBlockMatrix>& mat, bool structuresplit)
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

  std::vector<int> pciter;
  std::vector<double> pcomega;
  std::vector<int> spciter;
  std::vector<double> spcomega;
  std::vector<int> fpciter;
  std::vector<double> fpcomega;
  std::vector<int> apciter;
  std::vector<double> apcomega;
  std::vector<std::string> blocksmoother;
  std::vector<double> schuromega;
  {
    std::string word1;
    std::string word2;
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsimono, "PCITER"));
      std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsimono, "PCOMEGA"));
      while (pciterstream >> word1) pciter.push_back(std::atoi(word1.c_str()));
      while (pcomegastream >> word2) pcomega.push_back(std::atof(word2.c_str()));
    }
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsimono, "STRUCTPCITER"));
      std::istringstream pcomegastream(
          Teuchos::getNumericStringParameter(fsimono, "STRUCTPCOMEGA"));
      while (pciterstream >> word1) spciter.push_back(std::atoi(word1.c_str()));
      while (pcomegastream >> word2) spcomega.push_back(std::atof(word2.c_str()));
    }
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsimono, "FLUIDPCITER"));
      std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsimono, "FLUIDPCOMEGA"));
      while (pciterstream >> word1) fpciter.push_back(std::atoi(word1.c_str()));
      while (pcomegastream >> word2) fpcomega.push_back(std::atof(word2.c_str()));
    }
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsimono, "ALEPCITER"));
      std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsimono, "ALEPCOMEGA"));
      while (pciterstream >> word1) apciter.push_back(std::atoi(word1.c_str()));
      while (pcomegastream >> word2) apcomega.push_back(std::atof(word2.c_str()));
    }
    {
      std::string word;
      std::istringstream blocksmootherstream(
          Teuchos::getNumericStringParameter(fsimono, "BLOCKSMOOTHER"));
      while (blocksmootherstream >> word) blocksmoother.push_back(word);
    }
    {
      std::istringstream blocksmootherstream(
          Teuchos::getNumericStringParameter(fsimono, "SCHUROMEGA"));
      while (blocksmootherstream >> word2) schuromega.push_back(std::atof(word2.c_str()));
    }
  }

  INPAR::FSI::LinearBlockSolver linearsolverstrategy =
      DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

  // create block system matrix
  switch (linearsolverstrategy)
  {
    case INPAR::FSI::PreconditionedKrylov:
    case INPAR::FSI::FSIAMG:
    {
      mat = Teuchos::rcp(new OverlappingBlockMatrixFSIAMG(Extractor(), *StructureField(),
          *FluidField(), *AleField(), structuresplit,
          DRT::INPUT::IntegralValue<int>(fsimono, "SYMMETRICPRECOND"), blocksmoother, schuromega,
          pcomega, pciter, spcomega, spciter, fpcomega, fpciter, apcomega, apciter,
          DRT::INPUT::IntegralValue<int>(fsimono, "FSIAMGANALYZE"), linearsolverstrategy,
          verbosity_, DRT::Problem::Instance()->ErrorFile()->Handle()));
      break;
    }
    case INPAR::FSI::AMGnxn:
#ifdef HAVE_MueLu
    {
      // TODO This is a temporary hack to input the xml file without adding a new input parameter in
      // FSI/DYNAMIC MONOLITHIC SOLVER. We assume that the xml file is given in the first position
      // of the BLOCKSMOOTHER list
      std::string amgnxn_xml = "none";
      if ((int)blocksmoother.size() > 0)
        amgnxn_xml = blocksmoother[0];
      else
        dserror("Not found xml file in the first position of the BLOCKSMOOTHER list");
      mat = Teuchos::rcp(new OverlappingBlockMatrixAMGnxn(Extractor(), *StructureField(),
          *FluidField(), *AleField(), structuresplit, amgnxn_xml,
          DRT::Problem::Instance()->ErrorFile()->Handle()));
    }
#else
      dserror("The AMGnxn preconditioner works only with MueLu activated");
#endif  // HAVE_MueLu
    break;
    case INPAR::FSI::HybridSchwarz:
    {
      mat = Teuchos::rcp(new OverlappingBlockMatrixHybridSchwarz(Extractor(), *StructureField(),
          *FluidField(), *AleField(), structuresplit,
          DRT::INPUT::IntegralValue<int>(fsimono, "SYMMETRICPRECOND"), blocksmoother, schuromega,
          pcomega, pciter, spcomega, spciter, fpcomega, fpciter, apcomega, apciter,
          DRT::INPUT::IntegralValue<int>(fsimono, "FSIAMGANALYZE"), linearsolverstrategy,
          interfaceprocs_, verbosity_, DRT::Problem::Instance()->ErrorFile()->Handle()));
      break;
    }
    default:
    {
      dserror("Unsupported type of monolithic solver/preconditioner");
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem> FSI::BlockMonolithic::CreateLinearSystem(
    Teuchos::ParameterList& nlParams, NOX::Epetra::Vector& noxSoln, Teuchos::RCP<NOX::Utils> utils)
{
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys;

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  NOX::Epetra::Interface::Jacobian* iJac = this;
  NOX::Epetra::Interface::Preconditioner* iPrec = this;
  const Teuchos::RCP<Epetra_Operator> J = SystemMatrix();
  const Teuchos::RCP<Epetra_Operator> M = SystemMatrix();

  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  INPAR::FSI::LinearBlockSolver linearsolverstrategy =
      DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

  switch (linearsolverstrategy)
  {
    case INPAR::FSI::PreconditionedKrylov:
    case INPAR::FSI::FSIAMG:
    case INPAR::FSI::AMGnxn:
    {
      linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
          Teuchos::rcp(iJac, false), J, Teuchos::rcp(iPrec, false), M, noxSoln));
      break;
    }
    default:
    {
      dserror("unsupported linear block solver strategy: %d", linearsolverstrategy);
      break;
    }
  }

  return linSys;
}
