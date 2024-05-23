/*-----------------------------------------------------------*/
/*! \file

\brief Input/output data container for the structural (time) integration


\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_structure_new_timint_basedataio.hpp"

#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
#include "4C_io_every_iteration_writer.hpp"
#include "4C_solver_nonlin_nox_aux.hpp"
#include "4C_solver_nonlin_nox_linesearch_generic.hpp"
#include "4C_solver_nonlin_nox_linesearch_prepostoperator.hpp"
#include "4C_structure_new_timint_basedataio_monitor_dbc.hpp"
#include "4C_structure_new_timint_basedataio_runtime_vtk_output.hpp"
#include "4C_structure_new_timint_basedataio_runtime_vtp_output.hpp"

#include <NOX_Solver_Generic.H>

FOUR_C_NAMESPACE_OPEN

namespace
{
  inline bool DetermineWriteOutput(int step, int offset, int write_every)
  {
    return (step + offset) % write_every == 0;
  }
}  // namespace

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::BaseDataIO::BaseDataIO()
    : isinit_(false),
      issetup_(false),
      output_(Teuchos::null),
      writer_every_iter_(Teuchos::null),
      params_runtime_vtk_output_(Teuchos::null),
      params_runtime_vtp_output_(Teuchos::null),
      params_monitor_dbc_(Teuchos::null),
      energyfile_(Teuchos::null),
      gmsh_out_(false),
      printlogo_(false),
      printiter_(false),
      outputeveryiter_(false),
      writesurfactant_(false),
      writestate_(false),
      writevelacc_(false),
      writejac2matlab_(false),
      firstoutputofrun_(false),
      printscreen_(-1),
      outputcounter_(-1),
      writerestartevery_(-1),
      writeresultsevery_(-1),
      writeenergyevery_(-1),
      writestress_(INPAR::STR::stress_none),
      writecouplstress_(INPAR::STR::stress_none),
      writestrain_(INPAR::STR::strain_none),
      writeplstrain_(INPAR::STR::strain_none),
      writeoptquantity_(INPAR::STR::optquantity_none),
      conditionnumbertype_(INPAR::STR::ConditionNumber::none)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataIO::Init(const Teuchos::ParameterList& ioparams,
    const Teuchos::ParameterList& sdynparams, const Teuchos::ParameterList& xparams,
    Teuchos::RCP<IO::DiscretizationWriter> output)
{
  // We have to call Setup() after Init()
  issetup_ = false;

  // ---------------------------------------------------------------------------
  // initialize the printing and output parameters
  // ---------------------------------------------------------------------------
  {
    output_ = output;
    printscreen_ = ioparams.get<int>("STDOUTEVRY");
    printlogo_ = printscreen_ > 0;
    gmsh_out_ = (bool)CORE::UTILS::IntegralValue<int>(ioparams, "OUTPUT_GMSH");
    printiter_ = true;
    p_io_every_iteration_ =
        Teuchos::rcp(new Teuchos::ParameterList(ioparams.sublist("EVERY ITERATION")));
    outputeveryiter_ =
        CORE::UTILS::IntegralValue<bool>(*p_io_every_iteration_, "OUTPUT_EVERY_ITER");
    writerestartevery_ = sdynparams.get<int>("RESTARTEVRY");
    writetimestepoffset_ = sdynparams.get<int>("OUTPUT_STEP_OFFSET");
    writestate_ = (bool)CORE::UTILS::IntegralValue<int>(ioparams, "STRUCT_DISP");
    writevelacc_ = (bool)CORE::UTILS::IntegralValue<int>(ioparams, "STRUCT_VEL_ACC");
    writejac2matlab_ = (bool)CORE::UTILS::IntegralValue<int>(ioparams, "STRUCT_JACOBIAN_MATLAB");
    conditionnumbertype_ =
        Teuchos::getIntegralValue<INPAR::STR::ConditionNumber>(ioparams, "STRUCT_CONDITION_NUMBER");
    firstoutputofrun_ = true;
    writeresultsevery_ = sdynparams.get<int>("RESULTSEVRY");
    writecurrentelevolume_ =
        (bool)CORE::UTILS::IntegralValue<int>(ioparams, "STRUCT_CURRENT_VOLUME");
    writestress_ = CORE::UTILS::IntegralValue<INPAR::STR::StressType>(ioparams, "STRUCT_STRESS");
    writecouplstress_ =
        CORE::UTILS::IntegralValue<INPAR::STR::StressType>(ioparams, "STRUCT_COUPLING_STRESS");
    writestrain_ = CORE::UTILS::IntegralValue<INPAR::STR::StrainType>(ioparams, "STRUCT_STRAIN");
    writeplstrain_ =
        CORE::UTILS::IntegralValue<INPAR::STR::StrainType>(ioparams, "STRUCT_PLASTIC_STRAIN");
    writeenergyevery_ = sdynparams.get<int>("RESEVRYERGY");
    writesurfactant_ = (bool)CORE::UTILS::IntegralValue<int>(ioparams, "STRUCT_SURFACTANT");
    writeoptquantity_ = CORE::UTILS::IntegralValue<INPAR::STR::OptQuantityType>(
        ioparams, "STRUCT_OPTIONAL_QUANTITY");

    // build params container for monitoring reaction forces
    params_monitor_dbc_ = Teuchos::rcp(new ParamsMonitorDBC());
    params_monitor_dbc_->Init(ioparams.sublist("MONITOR STRUCTURE DBC"));
    params_monitor_dbc_->Setup();

    // check whether VTK output at runtime is desired
    if (ioparams.sublist("RUNTIME VTK OUTPUT").get<int>("INTERVAL_STEPS") != -1)
    {
      params_runtime_vtk_output_ = Teuchos::rcp(new ParamsRuntimeOutput());

      params_runtime_vtk_output_->Init(ioparams.sublist("RUNTIME VTK OUTPUT"));
      params_runtime_vtk_output_->Setup();
    }

    // check whether VTP output at runtime is desired
    if (ioparams.sublist("RUNTIME VTP OUTPUT STRUCTURE").get<int>("INTERVAL_STEPS") != -1)
    {
      params_runtime_vtp_output_ = Teuchos::rcp(new ParamsRuntimeVtpOutput());

      params_runtime_vtp_output_->Init(ioparams.sublist("RUNTIME VTP OUTPUT STRUCTURE"));
      params_runtime_vtp_output_->Setup();
    }
  }

  isinit_ = true;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataIO::Setup()
{
  // safety check
  FOUR_C_ASSERT(IsInit(), "Init() has not been called, yet!");

  if (outputeveryiter_) writer_every_iter_ = Teuchos::rcp(new IO::EveryIterationWriter());

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataIO::CheckInitSetup() const
{
  FOUR_C_ASSERT(IsInit() and IsSetup(), "Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataIO::init_setup_every_iteration_writer(
    IO::EveryIterationWriterInterface* interface, Teuchos::ParameterList& p_nox)
{
  if (not outputeveryiter_) return;

  writer_every_iter_->Init(output_.get(), interface, *p_io_every_iteration_);
  writer_every_iter_->Setup();

  // insert the every_iter output writer as ppo for the solver object
  Teuchos::ParameterList& p_sol_opt = p_nox.sublist("Solver Options");

  Teuchos::RCP<::NOX::Observer> prepost_solver_ptr = Teuchos::rcp(
      new NOX::NLN::Solver::PrePostOp::TIMINT::WriteOutputEveryIteration(*writer_every_iter_));

  NOX::NLN::AUX::AddToPrePostOpVector(p_sol_opt, prepost_solver_ptr);

  // insert the every_iter output writer as ppo for the linesearch object
  Teuchos::ParameterList& p_linesearch = p_nox.sublist("Line Search");

  // Get the current map. If there is no map, return a new empty one. (reference)
  NOX::NLN::LineSearch::PrePostOperator::map& prepostls_map =
      NOX::NLN::LineSearch::PrePostOperator::GetMap(p_linesearch);

  // insert/replace the old pointer in the map
  prepostls_map[NOX::NLN::LineSearch::prepost_output_every_iter] =
      Teuchos::rcp_dynamic_cast<NOX::NLN::Abstract::PrePostOperator>(prepost_solver_ptr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::BaseDataIO::setup_energy_output_file()
{
  if (energyfile_.is_null())
  {
    std::string energy_file_name =
        GLOBAL::Problem::Instance()->OutputControlFile()->FileName() + "_energy.csv";

    energyfile_ = Teuchos::rcp(new std::ofstream(energy_file_name.c_str()));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::BaseDataIO::write_results_for_this_step(const int step) const
{
  if (step < 0) FOUR_C_THROW("The variable step is not allowed to be negative.");
  return is_write_results_enabled() and
         DetermineWriteOutput(step, get_write_timestep_offset(), get_write_results_every_n_step());
}

bool STR::TIMINT::BaseDataIO::is_write_results_enabled() const
{
  return get_write_results_every_n_step() > 0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::BaseDataIO::write_runtime_vtk_results_for_this_step(const int step) const
{
  if (step < 0) FOUR_C_THROW("The variable step is not allowed to be negative.");
  return (is_runtime_output_enabled() &&
          DetermineWriteOutput(step, get_runtime_output_params()->OutputStepOffset(),
              get_runtime_output_params()->output_interval_in_steps()));
}

bool STR::TIMINT::BaseDataIO::is_runtime_output_enabled() const
{
  return get_runtime_output_params() != Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::BaseDataIO::write_runtime_vtp_results_for_this_step(const int step) const
{
  if (step < 0) FOUR_C_THROW("The variable step is not allowed to be negative.");
  return (get_runtime_vtp_output_params() != Teuchos::null &&
          DetermineWriteOutput(step, get_runtime_output_params()->OutputStepOffset(),
              get_runtime_output_params()->output_interval_in_steps()));
}


bool STR::TIMINT::BaseDataIO::should_write_restart_for_step(const int step) const
{
  return get_write_restart_every_n_step() &&
         DetermineWriteOutput(
             step, get_write_timestep_offset(), get_write_restart_every_n_step()) &&
         step != 0;
}


bool STR::TIMINT::BaseDataIO::should_write_reaction_forces_for_this_step(const int step) const
{
  return GetMonitorDBCParams()->output_interval_in_steps() > 0 &&
         DetermineWriteOutput(
             step, get_write_timestep_offset(), GetMonitorDBCParams()->output_interval_in_steps());
}


bool STR::TIMINT::BaseDataIO::should_write_stress_strain_for_this_step(const int step) const
{
  return write_results_for_this_step(step) &&
         ((GetStressOutputType() != INPAR::STR::stress_none) ||
             (get_coupling_stress_output_type() != INPAR::STR::stress_none) ||
             (GetStrainOutputType() != INPAR::STR::strain_none) ||
             (get_plastic_strain_output_type() != INPAR::STR::strain_none));
}

bool STR::TIMINT::BaseDataIO::should_write_energy_for_this_step(const int step) const
{
  return get_write_energy_every_n_step() > 0 &&
         DetermineWriteOutput(step, get_write_timestep_offset(), get_write_energy_every_n_step());
}

int STR::TIMINT::BaseDataIO::get_last_written_results() const { return lastwrittenresultsstep_; }

void STR::TIMINT::BaseDataIO::set_last_written_results(const int step)
{
  lastwrittenresultsstep_ = step;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Solver::PrePostOp::TIMINT::WriteOutputEveryIteration::WriteOutputEveryIteration(
    IO::EveryIterationWriter& every_iter_writer)
    : every_iter_writer_(every_iter_writer)
{
  /* empty */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PrePostOp::TIMINT::WriteOutputEveryIteration::runPreSolve(
    const ::NOX::Solver::Generic& solver)
{
  every_iter_writer_.InitNewtonIteration();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PrePostOp::TIMINT::WriteOutputEveryIteration::runPostIterate(
    const ::NOX::Solver::Generic& solver)
{
  const int newton_iteration = solver.getNumIterations();
  every_iter_writer_.AddNewtonIteration(newton_iteration);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PrePostOp::TIMINT::WriteOutputEveryIteration::run_pre_modify_step_length(
    const ::NOX::Solver::Generic& solver, const ::NOX::LineSearch::Generic& linesearch)
{
  const int newton_iteration = solver.getNumIterations();
  const int ls_iteration =
      dynamic_cast<const NOX::NLN::LineSearch::Generic&>(linesearch).GetNumIterations();
  every_iter_writer_.add_line_search_iteration(newton_iteration, ls_iteration);
}

FOUR_C_NAMESPACE_CLOSE
