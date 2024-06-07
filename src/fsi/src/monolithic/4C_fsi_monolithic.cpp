/*----------------------------------------------------------------------------*/
/*! \file

\brief General framework for monolithic fsi solution schemes


\level 1
*/
/*----------------------------------------------------------------------------*/

#include "4C_fsi_monolithic.hpp"

#include "4C_adapter_ale.hpp"
#include "4C_adapter_ale_fsi.hpp"
#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_adapter_str_fsi_timint_adaptive.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_ale_utils_mapextractor.hpp"
#include "4C_constraint_manager.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fsi_debugwriter.hpp"
#include "4C_fsi_nox_group.hpp"
#include "4C_fsi_nox_linearsystem.hpp"
#include "4C_fsi_nox_newton.hpp"
#include "4C_fsi_overlapprec_fsiamg.hpp"
#include "4C_fsi_overlapprec_hybrid.hpp"
#include "4C_fsi_statustest.hpp"
#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_linear_solver_method_parameters.hpp"
#include "4C_structure_aux.hpp"

#include <NOX_Direction_UserDefinedFactory.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

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
      verbosity_(Core::UTILS::IntegralValue<Inpar::FSI::Verbosity>(
          Global::Problem::Instance()->FSIDynamicParams(), "VERBOSITY"))
{
  // access the discretizations
  Teuchos::RCP<Discret::Discretization> structdis =
      Global::Problem::Instance()->GetDis("structure");
  Teuchos::RCP<Discret::Discretization> fluiddis = Global::Problem::Instance()->GetDis("fluid");
  Teuchos::RCP<Discret::Discretization> aledis = Global::Problem::Instance()->GetDis("ale");

  create_structure_time_integrator(timeparams, structdis);
  create_fluid_and_ale_time_integrator(timeparams, fluiddis, aledis);

  coupsf_ = Teuchos::rcp(new Core::Adapter::Coupling());
  coupsa_ = Teuchos::rcp(new Core::Adapter::Coupling());
  coupfa_ = Teuchos::rcp(new Core::Adapter::Coupling());
  icoupfa_ = Teuchos::rcp(new Core::Adapter::Coupling());
}



/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MonolithicBase::read_restart(int step)
{
  structure_field()->read_restart(step);
  fluid_field()->read_restart(step);
  ale_field()->read_restart(step);

  SetTimeStep(fluid_field()->Time(), fluid_field()->Step());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MonolithicBase::create_structure_time_integrator(
    const Teuchos::ParameterList& timeparams, Teuchos::RCP<Discret::Discretization> structdis)
{
  // delete deprecated time integrator
  structure_ = Teuchos::null;

  // access structural dynamic params
  const Teuchos::ParameterList& sdyn = Global::Problem::Instance()->structural_dynamic_params();

  // ask base algorithm for the structural time integrator
  Teuchos::RCP<Adapter::StructureBaseAlgorithm> structure =
      Teuchos::rcp(new Adapter::StructureBaseAlgorithm(
          timeparams, const_cast<Teuchos::ParameterList&>(sdyn), structdis));
  structure_ =
      Teuchos::rcp_dynamic_cast<Adapter::FSIStructureWrapper>(structure->structure_field());
  structure_->Setup();

  if (structure_ == Teuchos::null)
    FOUR_C_THROW(
        "Cast from Adapter::Structure to Adapter::FSIStructureWrapper "
        "failed.");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MonolithicBase::create_fluid_and_ale_time_integrator(
    const Teuchos::ParameterList& timeparams, Teuchos::RCP<Discret::Discretization> fluiddis,
    Teuchos::RCP<Discret::Discretization> aledis)
{
  // delete deprecated time integrators
  fluid_ = Teuchos::null;
  ale_ = Teuchos::null;

  // ask base algorithm for the fluid time integrator
  Teuchos::RCP<Adapter::FluidBaseAlgorithm> fluid = Teuchos::rcp(new Adapter::FluidBaseAlgorithm(
      timeparams, Global::Problem::Instance()->FluidDynamicParams(), "fluid", true));
  fluid_ = Teuchos::rcp_dynamic_cast<Adapter::FluidFSI>(fluid->fluid_field());

  if (fluid_ == Teuchos::null) FOUR_C_THROW("Cast from Adapter::Fluid to Adapter::FluidFSI failed");

  // ask base algorithm for the ale time integrator
  Teuchos::RCP<Adapter::AleBaseAlgorithm> ale =
      Teuchos::rcp(new Adapter::AleBaseAlgorithm(timeparams, aledis));
  ale_ = Teuchos::rcp_dynamic_cast<Adapter::AleFsiWrapper>(ale->ale_field());

  if (ale_ == Teuchos::null)
    FOUR_C_THROW("Cast from Adapter::Ale to Adapter::AleFsiWrapper failed");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MonolithicBase::prepare_time_step()
{
  increment_time_and_step();
  if (verbosity_ >= Inpar::FSI::verbosity_low) print_header();
  prepare_time_step_preconditioner();
  prepare_time_step_fields();

  // Note: it's important to first prepare the single fields and than the fsi problem
  prepare_time_step_fsi();

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MonolithicBase::prepare_time_step_fsi()
{
  ddgpred_ = Teuchos::rcp(new Epetra_Vector(*structure_field()->extract_interface_dispnp()));
  ddgpred_->Update(-1.0, *structure_field()->extract_interface_dispn(), 1.0);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MonolithicBase::prepare_time_step_fields()
{
  structure_field()->prepare_time_step();
  fluid_field()->prepare_time_step();
  ale_field()->prepare_time_step();

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MonolithicBase::prepare_output(bool force_prepare)
{
  structure_field()->prepare_output(force_prepare);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MonolithicBase::output()
{
  /* Note: The order is important here! In here control file entries are
   * written. And these entries define the order in which the filters handle
   * the Discretizations, which in turn defines the dof number ordering of the
   * Discretizations.
   */
  structure_field()->Output();
  fluid_field()->Output();
  ale_field()->Output();

  if (structure_field()->get_constraint_manager()->HaveMonitor())
  {
    structure_field()->get_constraint_manager()->compute_monitor_values(
        structure_field()->Dispnp());
    if (Comm().MyPID() == 0) structure_field()->get_constraint_manager()->PrintMonitorValues();
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::struct_to_ale(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsa_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::ale_to_struct(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::struct_to_fluid(
    Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::fluid_to_struct(
    Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::ale_to_fluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupfa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::fluid_to_ale_interface(
    Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupfa_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::ale_to_fluid_interface(
    Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupfa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::struct_to_ale(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::ale_to_struct(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::struct_to_fluid(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::fluid_to_struct(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::ale_to_fluid(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::fluid_to_ale_interface(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupfa_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::ale_to_fluid_interface(
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
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();

  // enable debugging
  if (Core::UTILS::IntegralValue<int>(fsidyn, "DEBUGOUTPUT") == 1)
  {
    sdbg_ = Teuchos::rcp(new UTILS::DebugWriter(structure_field()->discretization()));
    // fdbg_ = Teuchos::rcp(new UTILS::DebugWriter(fluid_field()->discretization()));
  }

  // write iterations-file
  std::string fileiter = Global::Problem::Instance()->OutputControlFile()->FileName();
  fileiter.append(".iteration");
  log_ = Teuchos::rcp(new std::ofstream(fileiter.c_str()));

  // write energy-file
  if (Core::UTILS::IntegralValue<int>(fsidyn.sublist("MONOLITHIC SOLVER"), "ENERGYFILE") == 1)
  {
    std::string fileiter2 = Global::Problem::Instance()->OutputControlFile()->FileName();
    fileiter2.append(".fsienergy");
    logenergy_ = Teuchos::rcp(new std::ofstream(fileiter2.c_str()));
  }

  // "Initialize" interface solution increments due to structural predictor
  ddgpred_ = Teuchos::null;

  //-------------------------------------------------------------------------
  // time step size adaptivity
  //-------------------------------------------------------------------------
  const bool timeadapton =
      Core::UTILS::IntegralValue<bool>(fsidyn.sublist("TIMEADAPTIVITY"), "TIMEADAPTON");

  if (timeadapton)
  {
    init_tim_int_ada(fsidyn);
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::SetupSystem()
{
  // right now we use matching meshes at the interface

  const int ndim = Global::Problem::Instance()->NDim();

  Core::Adapter::Coupling& coupsf = structure_fluid_coupling();
  Core::Adapter::Coupling& coupsa = structure_ale_coupling();
  Core::Adapter::Coupling& coupfa = fluid_ale_coupling();
  Core::Adapter::Coupling& icoupfa = interface_fluid_ale_coupling();

  // structure to fluid

  coupsf.setup_condition_coupling(*structure_field()->discretization(),
      structure_field()->Interface()->FSICondMap(), *fluid_field()->discretization(),
      fluid_field()->Interface()->FSICondMap(), "FSICoupling", ndim);

  // structure to ale

  coupsa.setup_condition_coupling(*structure_field()->discretization(),
      structure_field()->Interface()->FSICondMap(), *ale_field()->discretization(),
      ale_field()->Interface()->FSICondMap(), "FSICoupling", ndim);

  // fluid to ale at the interface

  icoupfa.setup_condition_coupling(*fluid_field()->discretization(),
      fluid_field()->Interface()->FSICondMap(), *ale_field()->discretization(),
      ale_field()->Interface()->FSICondMap(), "FSICoupling", ndim);

  // In the following we assume that both couplings find the same dof
  // map at the structural side. This enables us to use just one
  // interface dof map for all fields and have just one transfer
  // operator from the interface map to the full field map.
  if (not coupsf.MasterDofMap()->SameAs(*coupsa.MasterDofMap()))
    FOUR_C_THROW("structure interface dof maps do not match");

  if (coupsf.MasterDofMap()->NumGlobalElements() == 0)
    FOUR_C_THROW("No nodes in matching FSI interface. Empty FSI coupling condition?");

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = fluid_field()->discretization()->NodeRowMap();
  const Epetra_Map* alenodemap = ale_field()->discretization()->NodeRowMap();

  coupfa.setup_coupling(*fluid_field()->discretization(), *ale_field()->discretization(),
      *fluidnodemap, *alenodemap, ndim);

  fluid_field()->SetMeshMap(coupfa.MasterDofMap());

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::Timeloop(const Teuchos::RCP<::NOX::Epetra::Interface::Required>& interface)
{
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  const bool timeadapton =
      Core::UTILS::IntegralValue<bool>(fsidyn.sublist("TIMEADAPTIVITY"), "TIMEADAPTON");

  // Run time loop with constant or adaptive time step size (depending on the user's will)
  if (not timeadapton)
  {
    // call time loop with constant time step size
    timeloop_const_dt(interface);
  }
  else
  {
    // call time loop with adaptive time step size
    timeloop_ada_dt(interface);
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::timeloop_const_dt(
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& interface)
{
  PrepareTimeloop();

  while (NotFinished())
  {
    prepare_time_step();
    TimeStep(interface);
    constexpr bool force_prepare = false;
    prepare_output(force_prepare);
    update();
    output();
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::PrepareTimeloop()
{
  // make sure we didn't destroy the maps before we entered the timeloop
  extractor().check_for_valid_map_extractor();

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = nox_parameter_list();

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", Comm().MyPID());

  printParams.set(
      "Output Information", ::NOX::Utils::Error | ::NOX::Utils::Warning |
                                ::NOX::Utils::OuterIteration | ::NOX::Utils::InnerIteration |
                                // ::NOX::Utils::Parameters |
                                ::NOX::Utils::Details | ::NOX::Utils::OuterIterationStatusTest |
                                ::NOX::Utils::LinearSolverDetails | ::NOX::Utils::TestDetails |
                                ::NOX::Utils::StepperIteration | ::NOX::Utils::StepperDetails |
                                ::NOX::Utils::StepperParameters | ::NOX::Utils::Debug | 0);

  // Create printing utilities
  utils_ = Teuchos::rcp(new ::NOX::Utils(printParams));

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

  write_ada_file_header();

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
  const Inpar::STR::PreStress pstype = Teuchos::getIntegralValue<Inpar::STR::PreStress>(
      Global::Problem::Instance()->structural_dynamic_params(), "PRESTRESS");
  const double pstime =
      Global::Problem::Instance()->structural_dynamic_params().get<double>("PRESTRESSTIME");
  if (pstype != Inpar::STR::PreStress::none && Time() + Dt() <= pstime + 1.0e-15)
    FOUR_C_THROW("No monolithic FSI in the pre-phase of prestressing, use Aitken!");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::TimeStep(const Teuchos::RCP<::NOX::Epetra::Interface::Required>& interface)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::TimeStep");

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = nox_parameter_list();

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
    case Inpar::FSI::verbosity_full:
    {
      printParams.set(
          "Output Information", ::NOX::Utils::Error | ::NOX::Utils::Warning |
                                    ::NOX::Utils::OuterIteration | ::NOX::Utils::InnerIteration |
                                    // ::NOX::Utils::Parameters |
                                    ::NOX::Utils::Details |  // weg damit!
                                    ::NOX::Utils::OuterIterationStatusTest |
                                    ::NOX::Utils::LinearSolverDetails |  // weg damit!
                                    ::NOX::Utils::TestDetails | ::NOX::Utils::StepperIteration |
                                    ::NOX::Utils::StepperDetails | ::NOX::Utils::StepperParameters |
                                    ::NOX::Utils::Debug | 0);
      break;
    }
    case Inpar::FSI::verbosity_medium:
    {
      printParams.set(
          "Output Information", ::NOX::Utils::Error | ::NOX::Utils::Warning |
                                    ::NOX::Utils::OuterIteration | ::NOX::Utils::InnerIteration |
                                    // ::NOX::Utils::Parameters |
                                    // ::NOX::Utils::Details | //weg damit!
                                    ::NOX::Utils::OuterIterationStatusTest |
                                    ::NOX::Utils::LinearSolverDetails |  // weg damit!
                                    ::NOX::Utils::TestDetails | ::NOX::Utils::StepperIteration |
                                    ::NOX::Utils::StepperDetails | ::NOX::Utils::StepperParameters |
                                    ::NOX::Utils::Debug | 0);
      break;
    }
    case Inpar::FSI::verbosity_low:
    case Inpar::FSI::verbosity_subproblem:
    {
      printParams.set(
          "Output Information", ::NOX::Utils::Error | ::NOX::Utils::Warning |
                                    //                    ::NOX::Utils::OuterIteration |
                                    //                    ::NOX::Utils::InnerIteration |
                                    //                    //::NOX::Utils::Parameters |
                                    //  //                  ::NOX::Utils::Details |
                                    ::NOX::Utils::OuterIterationStatusTest |
                                    //  //                  ::NOX::Utils::LinearSolverDetails |
                                    //                    ::NOX::Utils::TestDetails |
                                    //                    ::NOX::Utils::StepperIteration |
                                    //                    ::NOX::Utils::StepperDetails |
                                    //                    ::NOX::Utils::StepperParameters |
                                    ::NOX::Utils::Debug | 0);
      break;
    }
    default:
    {
      FOUR_C_THROW("Verbosity level not supported!");
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
  Teuchos::RCP<Epetra_Vector> initial_guess_v =
      Teuchos::rcp(new Epetra_Vector(*dof_row_map(), true));
  initial_guess(initial_guess_v);

  ::NOX::Epetra::Vector noxSoln(initial_guess_v, ::NOX::Epetra::Vector::CreateView);

  // Create the linear system
  Teuchos::RCP<::NOX::Epetra::LinearSystem> linSys =
      create_linear_system(nlParams, noxSoln, utils_);

  // Create the Group
  Teuchos::RCP<NOX::FSI::Group> grp =
      Teuchos::rcp(new NOX::FSI::Group(*this, printParams, interface, noxSoln, linSys));

  // Convergence Tests
  Teuchos::RCP<::NOX::StatusTest::Combo> combo = create_status_test(nlParams, grp);

  // Create the solver
  Teuchos::RCP<::NOX::Solver::Generic> solver = ::NOX::Solver::buildSolver(
      grp, combo, Teuchos::RCP<Teuchos::ParameterList>(&nlParams, false));

  // we know we already have the first linear system calculated
  grp->CaptureSystemState();

  // solve the whole thing
  noxstatus_ = solver->solve();
  noxiter_ = solver->getNumIterations();

  // recover Lagrange multiplier \lambda_{\Gamma} at the interface at the end of each time step
  // (i.e. condensed traction/forces onto the structure) needed for rhs in next time step
  recover_lagrange_multiplier();

  // compute spurious interface energy increment due to temporal discretization
  calculate_interface_energy_increment();

  // stop time measurement
  timemonitor = Teuchos::null;

  if (Comm().MyPID() == 0)
  {
    (*log_) << std::right << std::setw(9) << Step() << std::right << std::setw(16) << Time()
            << std::right << std::setw(16) << timer.totalElapsedTime(true) << std::right
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
void FSI::Monolithic::update()
{
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  bool timeadapton =
      Core::UTILS::IntegralValue<int>(fsidyn.sublist("TIMEADAPTIVITY"), "TIMEADAPTON");

  if (not timeadapton)
    structure_field()->Update();  // constant dt
  else
  {
    structure_field()->Update(Time());  // variable/adaptive dt

    if (is_ada_structure())
      Teuchos::rcp_dynamic_cast<Adapter::StructureFSITimIntAda>(structure_field(), true)
          ->UpdateStepSize(Dt());
  }

  fluid_field()->Update();
  ale_field()->Update();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::NonLinErrorCheck()
{
  // assume convergence of nonlinear solver
  erroraction_ = erroraction_none;

  // if everything is fine, then return right now
  if (NoxStatus() == ::NOX::StatusTest::Converged)
  {
    return;
  }

  // The nonlinear solver did not converge. Thus, we have to take some action
  // that depends on the user's will given in the input file

  // get the FSI parameter list
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();

  // get the user's will
  const Inpar::FSI::DivContAct divcontype = Core::UTILS::IntegralValue<Inpar::FSI::DivContAct>(
      fsidyn.sublist("TIMEADAPTIVITY"), ("DIVERCONT"));

  if (NoxStatus() != ::NOX::StatusTest::Converged)
  {
    switch (divcontype)
    {
      case Inpar::FSI::divcont_stop:
      {
        // set the corresponding error action
        erroraction_ = erroraction_stop;

        // stop the simulation
        FOUR_C_THROW("Nonlinear solver did not converge in %i iterations in time step %i.",
            noxiter_, Step());
        break;
      }
      case Inpar::FSI::divcont_continue:
      {
        // set the corresponding error action
        erroraction_ = erroraction_continue;

        // Notify user about non-converged nonlinear solver, but do not abort the simulation
        if (Comm().MyPID() == 0)
        {
          Core::IO::cout << "\n*** Nonlinear solver did not converge in " << noxiter_
                         << " iterations in time step " << Step() << ". Continue ..."
                         << Core::IO::endl;
        }
        break;
      }
      case Inpar::FSI::divcont_halve_step:
      {
        // set the corresponding error action
        erroraction_ = erroraction_halve_step;

        const bool timeadapton =
            Core::UTILS::IntegralValue<bool>(fsidyn.sublist("TIMEADAPTIVITY"), "TIMEADAPTON");
        if (not timeadapton)
        {
          FOUR_C_THROW(
              "Nonlinear solver wants to halve the time step size. This is "
              "not possible in a time integrator with constant Delta t.");
        }

        // Notify user about non-converged nonlinear solver, but do not abort the simulation
        if (Comm().MyPID() == 0)
        {
          Core::IO::cout << Core::IO::endl
                         << "*** Nonlinear solver did not converge in " << noxiter_
                         << " iterations. Halve the time step size." << Core::IO::endl;
        }
        break;
      }
      case Inpar::FSI::divcont_revert_dt:
      {
        // set the corresponding error action
        erroraction_ = erroraction_revert_dt;

        const bool timeadapton =
            Core::UTILS::IntegralValue<bool>(fsidyn.sublist("TIMEADAPTIVITY"), "TIMEADAPTON");
        if (not timeadapton)
        {
          FOUR_C_THROW(
              "Nonlinear solver wants to revert the time step size. This is "
              "not possible in a time integrator with constant Delta t.");
        }

        // Notify user about non-converged nonlinear solver, but do not abort the simulation
        if (Comm().MyPID() == 0)
        {
          Core::IO::cout << Core::IO::endl
                         << "*** Nonlinear solver did not converge in " << noxiter_
                         << " iterations. Revert the time step size to the previous one."
                         << Core::IO::endl;
        }
        break;
      }
      default:
      {
        FOUR_C_THROW("Unknown action to cope with non-converged nonlinear solver.");
        break;
      }
    }
  }
  else
  {
    FOUR_C_THROW("Unknown ::NOX::StatusTest::StatusType.");
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::Evaluate(Teuchos::RCP<const Epetra_Vector> step_increment)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::Evaluate");

#ifdef FOUR_C_ENABLE_ASSERTIONS
  // check whether all fields have the same time step size
  check_if_dts_same();
#endif

  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;
  Teuchos::RCP<const Epetra_Vector> ax;

  if (step_increment != Teuchos::null)
  {
    extract_field_vectors(step_increment, sx, fx, ax);

    if (sdbg_ != Teuchos::null)
    {
      sdbg_->NewIteration();
      sdbg_->WriteVector("x", *structure_field()->Interface()->ExtractFSICondVector(sx));
    }
  }

  // Call all elements and assemble rhs and matrices
  // We only need the rhs here because NOX will ask for the rhs
  // only. But the Jacobian is stored internally and will be returned
  // later on without looking at x again!

  if (verbosity_ >= Inpar::FSI::verbosity_medium) utils()->out() << "\nEvaluate elements\n";

  {
    Teuchos::Time ts("structure", true);
    structure_field()->Evaluate(sx);
    if (verbosity_ >= Inpar::FSI::verbosity_medium)
      utils()->out() << "structure: " << ts.totalElapsedTime(true) << " sec\n";
  }

  {
    Teuchos::Time ta("ale", true);
    ale_field()->Evaluate(ax);
    if (verbosity_ >= Inpar::FSI::verbosity_medium)
      utils()->out() << "ale      : " << ta.totalElapsedTime(true) << " sec\n";
  }

  // transfer the current ale mesh positions to the fluid field
  Teuchos::RCP<Epetra_Vector> fluiddisp = ale_to_fluid(ale_field()->Dispnp());
  fluid_field()->apply_mesh_displacement(fluiddisp);

  {
    Teuchos::Time tf("fluid", true);
    fluid_field()->Evaluate(fx);
    if (verbosity_ >= Inpar::FSI::verbosity_medium)
      utils()->out() << "fluid    : " << tf.totalElapsedTime(true) << " sec\n";
  }

  if (verbosity_ >= Inpar::FSI::verbosity_medium) utils()->out() << "\n";
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::set_dof_row_maps(const std::vector<Teuchos::RCP<const Epetra_Map>>& maps)
{
  Teuchos::RCP<Epetra_Map> fullmap = Core::LinAlg::MultiMapExtractor::MergeMaps(maps);
  blockrowdofmap_.Setup(*fullmap, maps);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::set_default_parameters(
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
  Teuchos::RCP<::NOX::Direction::UserDefinedFactory> newtonfactory = Teuchos::rcp(this, false);
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
  lsParams.set<Inpar::FSI::Verbosity>("verbosity", verbosity_);  // verbosity level of FSI algorithm
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Direction::Generic> FSI::Monolithic::buildDirection(
    const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params) const
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
  setup_rhs(F);
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
void FSI::Monolithic::setup_rhs(Epetra_Vector& f, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::setup_rhs");

  firstcall_ = firstcall;

  // We want to add into a zero vector
  f.PutScalar(0.0);

  // contributions of single field residuals
  setup_rhs_residual(f);

  // contributions of Lagrange multiplier from last time step
  setup_rhs_lambda(f);

  // contributions of special "first nonlinear iteration" terms
  if (firstcall) setup_rhs_firstiter(f);

  if (dbcmaps_ !=
      Teuchos::null)  // ToDo: Remove if-"Abfrage" after 'dbcmaps_' has been introduced in lung fsi
  {
    // Finally, we take care of Dirichlet boundary conditions
    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(f));
    Teuchos::RCP<const Epetra_Vector> zeros = Teuchos::rcp(new const Epetra_Vector(f.Map(), true));
    Core::LinAlg::apply_dirichlet_to_system(*rhs, *zeros, *(dbcmaps_->CondMap()));
    f.Update(1.0, *rhs, 0.0);
  }

  // NOX expects the 'positive' residual. The negative sign for the
  // linearized Newton system J*dx=-r is done internally by NOX.
  // Since we assembled the right hand side, we have to invert the sign here.
  f.Scale(-1.);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::initial_guess(Teuchos::RCP<Epetra_Vector> initial_guess)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::initial_guess");

  combine_field_vectors(*initial_guess, structure_field()->initial_guess(),
      fluid_field()->initial_guess(), ale_field()->initial_guess(), true);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::combine_field_vectors(Epetra_Vector& v, Teuchos::RCP<const Epetra_Vector> sv,
    Teuchos::RCP<const Epetra_Vector> fv, Teuchos::RCP<const Epetra_Vector> av)
{
  extractor().AddVector(*sv, 0, v);
  extractor().AddVector(*fv, 1, v);
  extractor().AddVector(*av, 2, v);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::write_interface_energy_file(const double energystep, const double energysum)
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
  Core::LinAlg::BlockSparseMatrixBase& mat =
      Teuchos::dyn_cast<Core::LinAlg::BlockSparseMatrixBase>(Jac);
  setup_system_matrix(mat);
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

    const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
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
void FSI::BlockMonolithic::prepare_time_step_preconditioner()
{
  const Teuchos::ParameterList& fsimono =
      Global::Problem::Instance()->FSIDynamicParams().sublist("MONOLITHIC SOLVER");

  if (Core::UTILS::IntegralValue<int>(fsimono, "REBUILDPRECEVERYSTEP")) precondreusecount_ = 0;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::BlockMonolithic::create_system_matrix(
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>& mat, bool structuresplit)
{
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
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

  Inpar::FSI::LinearBlockSolver linearsolverstrategy =
      Core::UTILS::IntegralValue<Inpar::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

  // create block system matrix
  switch (linearsolverstrategy)
  {
    case Inpar::FSI::PreconditionedKrylov:
    {
      mat = Teuchos::rcp(new OverlappingBlockMatrixFSIAMG(extractor(), *structure_field(),
          *fluid_field(), *ale_field(), structuresplit,
          Core::UTILS::IntegralValue<int>(fsimono, "SYMMETRICPRECOND"), blocksmoother, schuromega,
          pcomega, pciter, spcomega, spciter, fpcomega, fpciter, apcomega, apciter,
          Core::UTILS::IntegralValue<int>(fsimono, "FSIAMGANALYZE"), linearsolverstrategy,
          verbosity_));
      break;
    }
    case Inpar::FSI::HybridSchwarz:
    {
      mat = Teuchos::rcp(new OverlappingBlockMatrixHybridSchwarz(extractor(), *structure_field(),
          *fluid_field(), *ale_field(), structuresplit,
          Core::UTILS::IntegralValue<int>(fsimono, "SYMMETRICPRECOND"), blocksmoother, schuromega,
          pcomega, pciter, spcomega, spciter, fpcomega, fpciter, apcomega, apciter,
          Core::UTILS::IntegralValue<int>(fsimono, "FSIAMGANALYZE"), linearsolverstrategy,
          interfaceprocs_, verbosity_));
      break;
    }
    case Inpar::FSI::LinalgSolver:
    {
      mat = Teuchos::rcp(
          new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
              extractor(), extractor(), 81, false, true));
      break;
    }
    default:
    {
      FOUR_C_THROW("Unsupported type of monolithic solver/preconditioner");
      break;
    }
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Epetra::LinearSystem> FSI::BlockMonolithic::create_linear_system(
    Teuchos::ParameterList& nlParams, ::NOX::Epetra::Vector& noxSoln,
    Teuchos::RCP<::NOX::Utils> utils)
{
  Teuchos::RCP<::NOX::Epetra::LinearSystem> linSys;

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  ::NOX::Epetra::Interface::Jacobian* iJac = this;
  ::NOX::Epetra::Interface::Preconditioner* iPrec = this;
  const Teuchos::RCP<Epetra_Operator> J = SystemMatrix();
  const Teuchos::RCP<Epetra_Operator> M = SystemMatrix();

  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  Inpar::FSI::LinearBlockSolver linearsolverstrategy =
      Core::UTILS::IntegralValue<Inpar::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

  switch (linearsolverstrategy)
  {
    case Inpar::FSI::PreconditionedKrylov:
    case Inpar::FSI::HybridSchwarz:
    {
      linSys = Teuchos::rcp(new ::NOX::Epetra::LinearSystemAztecOO(printParams, lsParams,
          Teuchos::rcp(iJac, false), J, Teuchos::rcp(iPrec, false), M, noxSoln));
      break;
    }
    case Inpar::FSI::LinalgSolver:
    {
      const int linsolvernumber = fsimono.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == -1)
        FOUR_C_THROW(
            "no linear solver defined for monolithic FSI. Please set LINEAR_SOLVER in FSI "
            "DYNAMIC/MONOLITHIC SOLVER to a valid number!");

      const Teuchos::ParameterList& fsisolverparams =
          Global::Problem::Instance()->SolverParams(linsolvernumber);

      auto solver = Teuchos::rcp(new Core::LinAlg::Solver(fsisolverparams, Comm()));

      const auto azprectype = Teuchos::getIntegralValue<Core::LinearSolver::PreconditionerType>(
          fsisolverparams, "AZPREC");

      switch (azprectype)
      {
        case Core::LinearSolver::PreconditionerType::multigrid_muelu_fsi:
        {
          solver->put_solver_params_to_sub_params("Inverse1", fsisolverparams);
          // This might be an alternative to "Core::LinAlg::FixNullspace()", directly calculate
          // nullspace on correct map
          // solver->Params().sublist("Inverse1").set<Teuchos::RCP<Epetra_Map>>("null space: map",
          // Teuchos::rcp(new Epetra_Map(system_matrix()->Matrix(0,0).RowMap())));
          structure_field()->discretization()->compute_null_space_if_necessary(
              solver->Params().sublist("Inverse1"));
          Core::LinearSolver::Parameters::FixNullSpace("Structure",
              *structure_field()->discretization()->dof_row_map(),
              SystemMatrix()->Matrix(0, 0).EpetraMatrix()->RowMap(),
              solver->Params().sublist("Inverse1"));

          solver->put_solver_params_to_sub_params("Inverse2", fsisolverparams);
          // This might be an alternative to "Core::LinAlg::FixNullspace()", directly calculate
          // nullspace on correct map
          // solver->Params().sublist("Inverse2").set<Teuchos::RCP<Epetra_Map>>("null space: map",
          // Teuchos::rcp(new Epetra_Map(system_matrix()->Matrix(1,1).RowMap())));
          fluid_field()->discretization()->compute_null_space_if_necessary(
              solver->Params().sublist("Inverse2"));
          Core::LinearSolver::Parameters::FixNullSpace("Fluid",
              *fluid_field()->discretization()->dof_row_map(),
              SystemMatrix()->Matrix(1, 1).EpetraMatrix()->RowMap(),
              solver->Params().sublist("Inverse2"));

          solver->put_solver_params_to_sub_params("Inverse3", fsisolverparams);
          // This might be an alternative to "Core::LinAlg::FixNullspace()", directly calculate
          // nullspace on correct map
          // solver->Params().sublist("Inverse3").set<Teuchos::RCP<Epetra_Map>>("null space: map",
          // Teuchos::rcp(new Epetra_Map(system_matrix()->Matrix(2,2).RowMap()))); we have to cast
          // the const on the ale discretization away!
          const_cast<Discret::Discretization&>(*(ale_field()->discretization()))
              .compute_null_space_if_necessary(solver->Params().sublist("Inverse3"));
          Core::LinearSolver::Parameters::FixNullSpace("Ale",
              *ale_field()->discretization()->dof_row_map(),
              SystemMatrix()->Matrix(2, 2).EpetraMatrix()->RowMap(),
              solver->Params().sublist("Inverse3"));

          break;
        }
        default:
        {
          std::cout << "\nWARNING:   MueLu FSI block preconditioner expected!" << std::endl;
          break;
        }
      }

      linSys = Teuchos::rcp(new NOX::FSI::LinearSystem(
          printParams, lsParams, Teuchos::rcp(iJac, false), J, noxSoln, solver));

      break;
    }
    default:
    {
      FOUR_C_THROW("unsupported linear block solver strategy: %d", linearsolverstrategy);
      break;
    }
  }

  return linSys;
}

FOUR_C_NAMESPACE_CLOSE
