/*-----------------------------------------------------------*/
/*! \file

\brief Base class for all structural time integration strategies.


\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_structure_new_timint_base.hpp"

#include "4C_beaminteraction_str_model_evaluator.hpp"
#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_gmsh.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_structure_new_dbc.hpp"
#include "4C_structure_new_enum_lists.hpp"
#include "4C_structure_new_factory.hpp"
#include "4C_structure_new_integrator.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_model_evaluator_factory.hpp"
#include "4C_structure_new_resulttest.hpp"
#include "4C_structure_new_timint_basedataio_monitor_dbc.hpp"
#include "4C_structure_new_timint_basedataio_runtime_vtk_output.hpp"
#include "4C_structure_new_timint_basedataio_runtime_vtp_output.hpp"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TimeInt::Base::Base()
    : StructureNew(),
      isinit_(false),
      issetup_(false),
      isrestarting_(false),
      state_is_insync_with_noxgroup_(true),
      dataio_(Teuchos::null),
      datasdyn_(Teuchos::null),
      dataglobalstate_(Teuchos::null),
      int_ptr_(Teuchos::null),
      dbc_ptr_(Teuchos::null)
{
  Epetra_Object::SetTracebackMode(1);
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::Init(const Teuchos::RCP<STR::TimeInt::BaseDataIO> dataio,
    const Teuchos::RCP<STR::TimeInt::BaseDataSDyn> datasdyn,
    const Teuchos::RCP<STR::TimeInt::BaseDataGlobalState> dataglobalstate)
{
  // ---------------------------------------------------------------------------
  // We need to call Setup() after Init()
  // ---------------------------------------------------------------------------
  issetup_ = false;

  // ---------------------------------------------------------------------------
  // initialize the data container ptrs
  // ---------------------------------------------------------------------------
  dataio_ = dataio;
  datasdyn_ = datasdyn;
  dataglobalstate_ = dataglobalstate;

  // ---------------------------------------------------------------------------
  // set isInit flag
  // ---------------------------------------------------------------------------
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::Setup()
{
  check_init();

  // ---------------------------------------------------------------------------
  // Create the Dirichlet Boundary Condition handler
  // ---------------------------------------------------------------------------
  dbc_ptr_ = STR::build_dbc(data_sdyn());
  /* FixMe It would be sufficient to use a constant discretization,
   * unfortunately this wasn't considered during the implementation of the
   * discretization routines. Therefore many methods need a slight modification
   * (most times adding a "const" should fix the problem).          hiermeier */
  Teuchos::RCP<Core::FE::Discretization> discret_ptr = data_global_state().get_discret();
  dbc_ptr_->Init(discret_ptr, data_global_state().get_freact_np(), Teuchos::rcp(this, false));
  dbc_ptr_->Setup();

  // ---------------------------------------------------------------------------
  // Create the explicit/implicit integrator
  // ---------------------------------------------------------------------------
  int_ptr_ = STR::build_integrator(data_sdyn());
  int_ptr_->Init(data_s_dyn_ptr(), data_global_state_ptr(), data_io_ptr(), dbc_ptr_,
      Teuchos::rcp(this, false));
  int_ptr_->Setup();
  int_ptr_->post_setup();
  // Initialize and Setup the input/output writer for every Newton iteration
  dataio_->init_setup_every_iteration_writer(this, data_sdyn().get_nox_params());

  // Initialize the output of system energy
  if (dataio_->get_write_energy_every_n_step())
  {
    select_energy_types_to_be_written();

    if (dataglobalstate_->get_my_rank() == 0) initialize_energy_file_stream_and_write_headers();
  }

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::Reset()
{
  FOUR_C_THROW(
      "Reset of all class variables is not yet implemented for "
      "the modelevaluator!");
  // ModelEvaluator().Reset();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::reset_step()
{
  check_init_setup();

  int_ptr_->reset_step_state();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TimeInt::Base::not_finished() const
{
  check_init_setup();

  // check for early stopping
  const bool early_stop = int_ptr_->early_stopping();
  if (early_stop)
  {
    // Simulation is finished regardless of the simulation time or the time step.
    return false;
  }

  // check the current time
  const double& timenp = dataglobalstate_->get_time_np();
  const double& timemax = datasdyn_->GetTimeMax();
  const double& dt = (*dataglobalstate_->get_delta_time())[0];
  // check the step counter
  const int& stepnp = dataglobalstate_->get_step_np();
  const int& stepmax = datasdyn_->GetStepMax();

  return (timenp <= timemax + 1.0e-8 * dt and stepnp <= stepmax);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::set_restart(int stepn, double timen, Teuchos::RCP<Epetra_Vector> disn,
    Teuchos::RCP<Epetra_Vector> veln, Teuchos::RCP<Epetra_Vector> accn,
    Teuchos::RCP<std::vector<char>> elementdata, Teuchos::RCP<std::vector<char>> nodedata)
{
  check_init_setup();

  FOUR_C_THROW("SetRestartState() is deprecated, use the read_restart() routine instead!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map& STR::TimeInt::Base::GetMassDomainMap() const
{
  check_init_setup();
  return dataglobalstate_->get_mass_matrix()->DomainMap();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Core::LinAlg::MapExtractor> STR::TimeInt::Base::GetDBCMapExtractor()
{
  check_init_setup();
  return dbc_ptr_->GetDBCMapExtractor();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const Core::LinAlg::MapExtractor> STR::TimeInt::Base::GetDBCMapExtractor() const
{
  check_init_setup();
  return dbc_ptr_->GetDBCMapExtractor();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::Conditions::LocsysManager> STR::TimeInt::Base::LocsysManager()
{
  check_init_setup();
  return dbc_ptr_->LocSysManagerPtr();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::MODELEVALUATOR::Generic& STR::TimeInt::Base::ModelEvaluator(
    Inpar::STR::ModelType mtype) const
{
  return integrator().model_eval().evaluator(mtype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Generic& STR::TimeInt::Base::ModelEvaluator(Inpar::STR::ModelType mtype)
{
  return integrator().model_eval().evaluator(mtype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TimeInt::Base::TimIntParam() const
{
  check_init_setup();
  return int_ptr_->get_int_param();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::resize_m_step_tim_ada()
{
  check_init_setup();
  // resize time and stepsize fields
  const double& timen = dataglobalstate_->get_time_n();
  dataglobalstate_->get_multi_time()->Resize(-1, 0, timen);
  const double& dtn = (*dataglobalstate_->get_delta_time())[0];
  dataglobalstate_->get_delta_time()->Resize(-1, 0, dtn);

  // resize state vectors, AB2 is a 2-step method, thus we need two
  // past steps at t_{n} and t_{n-1}
  const Epetra_Map* dofrowmap_ptr = dataglobalstate_->dof_row_map_view();
  dataglobalstate_->get_multi_dis()->Resize(-1, 0, dofrowmap_ptr, true);
  dataglobalstate_->get_multi_vel()->Resize(-1, 0, dofrowmap_ptr, true);
  dataglobalstate_->get_multi_acc()->Resize(-1, 0, dofrowmap_ptr, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::Update()
{
  check_init_setup();
  int_ptr_->pre_update();
  int_ptr_->update_structural_energy();
  int_ptr_->update_step_state();
  update_step_time();
  set_number_of_nonlinear_iterations();
  int_ptr_->update_step_element();
  int_ptr_->post_update();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::update_step_time()
{
  check_init_setup();
  double& timenp = dataglobalstate_->get_time_np();
  int& stepnp = dataglobalstate_->get_step_np();
  int& stepn = dataglobalstate_->get_step_n();

  // --------------------------------------------------------------------------
  // update old time and step variables
  // --------------------------------------------------------------------------
  dataglobalstate_->get_multi_time()->UpdateSteps(timenp);
  stepn = stepnp;

  // --------------------------------------------------------------------------
  // update the new time and step variables
  // --------------------------------------------------------------------------
  // get current time step size
  const double& dtn = (*dataglobalstate_->get_delta_time())[0];
  dataglobalstate_->get_delta_time()->UpdateSteps(dtn);
  timenp += dtn;
  stepnp += 1;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::set_number_of_nonlinear_iterations()
{
  int nlniter = 0;

  if (data_sdyn().get_nox_params().isSublist("Output"))
  {
    const Teuchos::ParameterList& nox_output = data_sdyn().get_nox_params().sublist("Output");
    if (nox_output.isParameter("Nonlinear Iterations"))
      nlniter = nox_output.get<int>("Nonlinear Iterations");
  }

  dataglobalstate_->set_nln_iteration_number(nlniter);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::select_energy_types_to_be_written()
{
  STR::MODELEVALUATOR::Data& evaldata = int_ptr_->eval_data();

  // decide which types of energy contributions shall be written separately
  const std::set<enum Inpar::STR::ModelType>& mtypes = datasdyn_->get_model_types();

  std::set<enum Inpar::STR::ModelType>::const_iterator model_iter;
  for (model_iter = mtypes.begin(); model_iter != mtypes.end(); ++model_iter)
  {
    switch (*model_iter)
    {
      case Inpar::STR::model_structure:
      {
        evaldata.insert_energy_type_to_be_considered(STR::internal_energy);
        evaldata.insert_energy_type_to_be_considered(STR::kinetic_energy);
        break;
      }
      case Inpar::STR::model_beaminteraction:
      {
        STR::MODELEVALUATOR::BeamInteraction const beaminteraction_evaluator =
            dynamic_cast<STR::MODELEVALUATOR::BeamInteraction const&>(
                int_ptr_->model_eval_ptr()->evaluator(Inpar::STR::model_beaminteraction));

        if (beaminteraction_evaluator.HaveSubModelType(
                Inpar::BEAMINTERACTION::submodel_beamcontact))
        {
          evaldata.insert_energy_type_to_be_considered(STR::beam_contact_penalty_potential);
        }
        if (beaminteraction_evaluator.HaveSubModelType(Inpar::BEAMINTERACTION::submodel_potential))
        {
          evaldata.insert_energy_type_to_be_considered(STR::beam_interaction_potential);
        }
        if (beaminteraction_evaluator.HaveSubModelType(
                Inpar::BEAMINTERACTION::submodel_crosslinking))
        {
          evaldata.insert_energy_type_to_be_considered(STR::beam_to_beam_link_internal_energy);
          evaldata.insert_energy_type_to_be_considered(STR::beam_to_beam_link_kinetic_energy);
        }
        if (beaminteraction_evaluator.HaveSubModelType(
                Inpar::BEAMINTERACTION::submodel_spherebeamlink))
        {
          evaldata.insert_energy_type_to_be_considered(STR::beam_to_sphere_link_internal_energy);
          evaldata.insert_energy_type_to_be_considered(STR::beam_to_sphere_link_kinetic_energy);
        }
        break;
      }
      default:
      {
        break;
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::initialize_energy_file_stream_and_write_headers()
{
  auto& evaldata = int_ptr_->eval_data();

  dataio_->setup_energy_output_file();

  // write column headers to file
  dataio_->get_energy_output_stream() << std::setw(12) << "#timestep," << std::setw(24) << "time,";

  for (const auto& energy_data : evaldata.get_energy_data())
  {
    dataio_->get_energy_output_stream()
        << std::setw(36) << STR::EnergyType2String(energy_data.first) + ",";
  }

  dataio_->get_energy_output_stream() << std::setw(24) << "total_energy" << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::UTILS::ResultTest> STR::TimeInt::Base::CreateFieldTest()
{
  check_init_setup();
  Teuchos::RCP<STR::ResultTest> resulttest = Teuchos::rcp(new STR::ResultTest());
  resulttest->Init(get_data_global_state(), integrator().eval_data());
  resulttest->Setup();

  return resulttest;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::get_restart_data(Teuchos::RCP<int> step, Teuchos::RCP<double> time,
    Teuchos::RCP<Epetra_Vector> disnp, Teuchos::RCP<Epetra_Vector> velnp,
    Teuchos::RCP<Epetra_Vector> accnp, Teuchos::RCP<std::vector<char>> elementdata,
    Teuchos::RCP<std::vector<char>> nodedata)
{
  check_init_setup();
  // at some point we have to create a copy
  *step = dataglobalstate_->get_step_n();
  *time = dataglobalstate_->get_time_n();
  Teuchos::RCP<const Core::FE::Discretization> discret_ptr =
      Teuchos::rcp_dynamic_cast<const Core::FE::Discretization>(dataglobalstate_->get_discret());
  *elementdata = *(discret_ptr->PackMyElements());
  *nodedata = *(discret_ptr->PackMyNodes());

  // get restart data is only for simple structure problems
  // hence if the model set is larger than one, we throw an error
  if (datasdyn_->get_model_types().size() > 1)
    FOUR_C_THROW("The get_restart_data routine supports the structural model case ONLY!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::prepare_output(bool force_prepare_timestep)
{
  check_init_setup();
  // --- stress, strain and optional quantity calculation ---------------------
  if ((dataio_->is_write_results_enabled() && force_prepare_timestep) ||
      dataio_->write_results_for_this_step(dataglobalstate_->get_step_np()))
  {
    int_ptr_->determine_stress_strain();
    int_ptr_->determine_optional_quantity();

    if (dataio_->is_write_current_ele_volume())
    {
      Teuchos::RCP<Epetra_Vector> elevolumes = Teuchos::null;
      Teuchos::RCP<const Epetra_Vector> disnp = dataglobalstate_->get_dis_np();

      int_ptr_->determine_element_volumes(*disnp, elevolumes);
      int_ptr_->eval_data().set_element_volume_data(elevolumes);
    }
  }
  if ((dataio_->is_runtime_output_enabled() && force_prepare_timestep) ||
      dataio_->write_runtime_vtk_results_for_this_step(dataglobalstate_->get_step_np()) ||
      dataio_->write_runtime_vtp_results_for_this_step(dataglobalstate_->get_step_np()))
  {
    int_ptr_->runtime_pre_output_step_state();
  }
  // --- energy calculation ---------------------------------------------------
  if ((dataio_->get_write_energy_every_n_step() and
          (force_prepare_timestep ||
              dataglobalstate_->get_step_np() % dataio_->get_write_energy_every_n_step() == 0)))
  {
    STR::MODELEVALUATOR::Data& evaldata = int_ptr_->eval_data();
    evaldata.clear_values_for_all_energy_types();

    int_ptr_->determine_energy();

    // sum processor-local values of all separate contributions into global value
    double energy_local = 0.0;
    double energy_global = 0.0;

    for (const auto& energy_data : evaldata.get_energy_data())
    {
      energy_local = energy_data.second;

      dataglobalstate_->get_comm().SumAll(&energy_local, &energy_global, 1);

      evaldata.set_value_for_energy_type(energy_global, energy_data.first);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::Output(bool forced_writerestart)
{
  check_init_setup();
  output_step(forced_writerestart);
  // write Gmsh output
  write_gmsh_struc_output_step();
  int_ptr_->post_output();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::output_step(bool forced_writerestart)
{
  check_init_setup();
  // special treatment is necessary when restart is forced
  if (forced_writerestart)
  {
    // reset possible history data on element level
    reset_step();
    // restart has already been written or simulation has just started
    if (dataio_->should_write_restart_for_step(dataglobalstate_->get_step_n()) or
        dataglobalstate_->get_step_n() == Global::Problem::Instance()->Restart())
      return;
    // if state already exists, add restart information
    if (dataio_->write_results_for_this_step(dataglobalstate_->get_step_n()))
    {
      add_restart_to_output_state();
      return;
    }
  }

  /* This flag indicates whether some form of output has already been written in the current time
   * step. It is passed along subroutines and prevents repeated initialization of output writer,
   * printing of state vectors, or similar.
   */
  bool datawritten = false;

  // output restart (try this first)
  // write restart step
  if (forced_writerestart || dataio_->should_write_restart_for_step(dataglobalstate_->get_step_n()))
  {
    output_restart(datawritten);
    dataio_->set_last_written_results(dataglobalstate_->get_step_n());
  }

  // output results (not necessary if restart in same step)
  if (dataio_->is_write_state() and
      dataio_->write_results_for_this_step(dataglobalstate_->get_step_n()) and (not datawritten))
  {
    new_io_step(datawritten);
    output_state();
    dataio_->set_last_written_results(dataglobalstate_->get_step_n());
  }

  // output results during runtime ( not used for restart so far )
  if (dataio_->write_runtime_vtk_results_for_this_step(dataglobalstate_->get_step_n()) or
      dataio_->write_runtime_vtp_results_for_this_step(dataglobalstate_->get_step_n()))
  {
    runtime_output_state();
  }

  // write reaction forces
  if (dataio_->should_write_reaction_forces_for_this_step(dataglobalstate_->get_step_n()))
  {
    output_reaction_forces();
  }

  // output stress, strain and optional quantity
  if (dataio_->should_write_stress_strain_for_this_step(dataglobalstate_->get_step_n()))
  {
    new_io_step(datawritten);
    output_stress_strain();
    output_optional_quantity();
  }

  if (dataio_->write_results_for_this_step(dataglobalstate_->get_step_n()) and
      dataio_->is_write_current_ele_volume())
  {
    new_io_step(datawritten);
    Core::IO::DiscretizationWriter& iowriter = *(dataio_->get_output_ptr());
    output_element_volume(iowriter);
  }

  // output energy
  if (dataio_->should_write_energy_for_this_step(dataglobalstate_->get_step_n()))
  {
    output_energy();
  }

  //  OutputVolumeMass();

  // ToDo output of nodal positions in current configuration
  //  output_nodal_positions();

  // ToDo write output on micro-scale (multi-scale analysis)
  //  if (HaveMicroMat())
  //    FOUR_C_THROW("OutputMicro() is not yet implemented!"); // OutputMicro();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::new_io_step(bool& datawritten)
{
  if (not datawritten)
  {
    // Make new step
    dataio_->get_output_ptr()->NewStep(
        dataglobalstate_->get_step_n(), dataglobalstate_->get_time_n());

    datawritten = true;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::output_state()
{
  check_init_setup();
  Core::IO::DiscretizationWriter& iowriter = *(dataio_->get_output_ptr());

  output_state(iowriter, dataio_->is_first_output_of_run());

  dataio_->set_first_output_of_run(false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::OutputDebugState(
    Core::IO::DiscretizationWriter& iowriter, bool write_owner) const
{
  output_state(iowriter, write_owner);

  // write element volumes as additional debugging information, if activated
  if (dataio_->is_write_current_ele_volume()) output_element_volume(iowriter);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::output_state(
    Core::IO::DiscretizationWriter& iowriter, bool write_owner) const
{
  // owner of elements is just written once because it does not change during
  // simulation (so far)
  iowriter.WriteElementData(write_owner);
  iowriter.WriteNodeData(write_owner);

  int_ptr_->output_step_state(iowriter);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::runtime_output_state()
{
  check_init_setup();
  int_ptr_->runtime_output_step_state();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::output_reaction_forces()
{
  check_init_setup();
  Core::IO::DiscretizationWriter& iowriter = *(dataio_->get_output_ptr());
  int_ptr_->monitor_dbc(iowriter);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::output_element_volume(Core::IO::DiscretizationWriter& iowriter) const
{
  check_init_setup();

  STR::MODELEVALUATOR::Data& evaldata = int_ptr_->eval_data();

  iowriter.WriteVector("current_ele_volumes",
      Teuchos::rcpFromRef(evaldata.current_element_volume_data()), Core::IO::elementvector);

  evaldata.set_element_volume_data(Teuchos::null);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::output_stress_strain()
{
  check_init_setup();

  STR::MODELEVALUATOR::Data& evaldata = int_ptr_->eval_data();
  Teuchos::RCP<Core::IO::DiscretizationWriter> output_ptr = dataio_->get_output_ptr();

  // ---------------------------------------------------------------------------
  // write stress output
  // ---------------------------------------------------------------------------
  std::string text = "";
  if (dataio_->get_stress_output_type() != Inpar::STR::stress_none)
  {
    switch (dataio_->get_stress_output_type())
    {
      case Inpar::STR::stress_cauchy:
        text = "gauss_cauchy_stresses_xyz";
        break;
      case Inpar::STR::stress_2pk:
        text = "gauss_2PK_stresses_xyz";
        break;
      default:
        FOUR_C_THROW("Requested stress type is not supported!");
        break;
    }
    output_ptr->WriteVector(text, evaldata.stress_data(), *(discretization()->ElementRowMap()));
  }
  // we don't need this anymore
  evaldata.stress_data_ptr() = Teuchos::null;

  // ---------------------------------------------------------------------------
  // write coupling stress output
  // ---------------------------------------------------------------------------
  text.clear();
  if (dataio_->get_coupling_stress_output_type() != Inpar::STR::stress_none)
  {
    switch (dataio_->get_coupling_stress_output_type())
    {
      case Inpar::STR::stress_cauchy:
        text = "gauss_cauchy_coupling_stresses_xyz";
        break;
      case Inpar::STR::stress_2pk:
        text = "gauss_2PK_coupling_stresses_xyz";
        break;
      default:
        FOUR_C_THROW("Requested coupling stress type is not supported!");
        break;
    }
    output_ptr->WriteVector(
        text, evaldata.coupling_stress_data(), *(discretization()->ElementRowMap()));
  }

  // ---------------------------------------------------------------------------
  // write strain output
  // ---------------------------------------------------------------------------
  text.clear();
  if (dataio_->get_strain_output_type() != Inpar::STR::strain_none)
  {
    switch (dataio_->get_strain_output_type())
    {
      case Inpar::STR::strain_none:
        break;
      case Inpar::STR::strain_ea:
        text = "gauss_EA_strains_xyz";
        break;
      case Inpar::STR::strain_gl:
        text = "gauss_GL_strains_xyz";
        break;
      case Inpar::STR::strain_log:
        text = "gauss_LOG_strains_xyz";
        break;
      default:
        FOUR_C_THROW("Requested strain type is not supported!");
        break;
    }
    output_ptr->WriteVector(text, evaldata.strain_data(), *(discretization()->ElementRowMap()));
  }
  // we don't need this anymore
  evaldata.strain_data_ptr() = Teuchos::null;

  // ---------------------------------------------------------------------------
  // write plastic strain output
  // ---------------------------------------------------------------------------
  text.clear();
  if (dataio_->get_plastic_strain_output_type() != Inpar::STR::strain_none)
  {
    switch (dataio_->get_plastic_strain_output_type())
    {
      case Inpar::STR::strain_ea:
        text = "gauss_pl_EA_strains_xyz";
        break;
      case Inpar::STR::strain_gl:
        text = "gauss_pl_GL_strains_xyz";
        break;
      default:
        FOUR_C_THROW("Requested plastic strain type is not supported!");
        break;
    }
    output_ptr->WriteVector(
        text, evaldata.plastic_strain_data(), *(discretization()->ElementRowMap()));
  }
  // we don't need this anymore
  evaldata.plastic_strain_data_ptr() = Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::output_energy() const
{
  check_init_setup();

  if (dataglobalstate_->get_my_rank() == 0)
  {
    std::ostream& energy_output_stream = dataio_->get_energy_output_stream();

    energy_output_stream << std::setw(11) << dataglobalstate_->get_step_n() << std::setw(1) << ","
                         << std::scientific << std::setprecision(14) << std::setw(23)
                         << dataglobalstate_->get_time_n() << std::setw(1) << ",";

    STR::MODELEVALUATOR::Data& evaldata = int_ptr_->eval_data();

    double total_energy = 0.0;

    for (const auto& energy_data : evaldata.get_energy_data())
    {
      energy_output_stream << std::setw(35) << energy_data.second << std::setw(1) << ",";
      total_energy += energy_data.second;
    }

    energy_output_stream << std::setw(24) << total_energy << std::endl;

    Core::IO::cout(Core::IO::verbose) << "\n\nOutput for energy written to file!" << Core::IO::endl;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::output_optional_quantity()
{
  check_init_setup();

  STR::MODELEVALUATOR::Data& evaldata = int_ptr_->eval_data();
  Teuchos::RCP<Core::IO::DiscretizationWriter> output_ptr = dataio_->get_output_ptr();

  // ---------------------------------------------------------------------------
  // write optional quantity output
  // ---------------------------------------------------------------------------
  std::string text = "";
  if (dataio_->get_opt_quantity_output_type() != Inpar::STR::optquantity_none)
  {
    switch (dataio_->get_opt_quantity_output_type())
    {
      case Inpar::STR::optquantity_membranethickness:
        text = "gauss_membrane_thickness";
        break;
      default:
        FOUR_C_THROW("Requested optional quantity type is not supported!");
        break;
    }
    output_ptr->WriteVector(
        text, evaldata.opt_quantity_data(), *(discretization()->ElementRowMap()));
  }
  // we don't need this anymore
  evaldata.opt_quantity_data_ptr() = Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::output_restart(bool& datawritten)
{
  check_init_setup();

  Teuchos::RCP<Core::IO::DiscretizationWriter> output_ptr = dataio_->get_output_ptr();
  // write restart output, please
  if (dataglobalstate_->get_step_n() != 0)
    output_ptr->WriteMesh(dataglobalstate_->get_step_n(), dataglobalstate_->get_time_n());
  new_io_step(datawritten);

  output_ptr->WriteElementData(dataio_->is_first_output_of_run());
  output_ptr->WriteNodeData(dataio_->is_first_output_of_run());
  dataio_->set_first_output_of_run(false);

  // add velocity and acceleration if necessary
  output_ptr->WriteVector("velocity", dataglobalstate_->get_vel_n());
  output_ptr->WriteVector("acceleration", dataglobalstate_->get_acc_n());

  /* Add the restart information of the different time integrators and model
   * evaluators. */
  int_ptr_->write_restart(*output_ptr);

  // info dedicated to user's eyes staring at standard out
  if ((dataglobalstate_->get_my_rank() == 0) and (dataio_->get_print2_screen_every_n_step() > 0) and
      (StepOld() % dataio_->get_print2_screen_every_n_step() == 0))
  {
    Core::IO::cout << "====== Restart for field 'Structure' written in step "
                   << dataglobalstate_->get_step_n() << Core::IO::endl;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::add_restart_to_output_state()
{
  Teuchos::RCP<Core::IO::DiscretizationWriter> output_ptr = dataio_->get_output_ptr();

  // force output of velocity and acceleration in case it is not written previously by the model
  // evaluators
  if (!dataio_->is_write_vel_acc())
  {
    output_ptr->WriteVector("velocity", dataglobalstate_->get_vel_n());
    output_ptr->WriteVector("acceleration", dataglobalstate_->get_acc_n());
  }

  /* Add the restart information of the different time integrators and model
   * evaluators. */
  int_ptr_->write_restart(*output_ptr, true);

  // finally add the missing mesh information, order is important here
  output_ptr->WriteMesh(data_global_state().get_step_n(), data_global_state().get_time_n());

  // info dedicated to user's eyes staring at standard out
  if ((dataglobalstate_->get_my_rank() == 0) and (dataio_->get_print2_screen_every_n_step() > 0) and
      (StepOld() % dataio_->get_print2_screen_every_n_step() == 0))
  {
    Core::IO::cout << "====== Restart for field 'Structure' written in step "
                   << dataglobalstate_->get_step_n() << Core::IO::endl;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::write_gmsh_struc_output_step()
{
  check_init_setup();
  if (!dataio_->is_gmsh()) return;

  const std::string filename =
      Core::IO::Gmsh::GetFileName("struct", disc_writer()->Output()->FileName(),
          dataglobalstate_->get_step_np(), false, dataglobalstate_->get_my_rank());
  std::ofstream gmshfilecontent(filename.c_str());

  // add 'View' to Gmsh postprocessing file
  gmshfilecontent << "View \" "
                  << "struct displacement \" {" << std::endl;
  // draw vector field 'struct displacement' for every element
  Core::IO::Gmsh::VectorFieldDofBasedToGmsh(discretization(), Dispn(), gmshfilecontent, 0, true);
  gmshfilecontent << "};" << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::read_restart(const int stepn)
{
  check_init();
  // that restarting flag
  isrestarting_ = true;

  // create an input/output reader
  Core::IO::DiscretizationReader ioreader(
      discretization(), Global::Problem::Instance()->InputControlFile(), stepn);
  dataglobalstate_->get_step_n() = stepn;
  dataglobalstate_->get_step_np() = stepn + 1;
  dataglobalstate_->get_multi_time() =
      Teuchos::rcp(new TimeStepping::TimIntMStep<double>(0, 0, ioreader.ReadDouble("time")));
  const double& timen = dataglobalstate_->get_time_n();
  const double& dt = (*dataglobalstate_->get_delta_time())[0];
  dataglobalstate_->get_time_np() = timen + dt;

  // ---------------------------------------------------------------------------
  // The order is important at this point!
  // (0) read element and node data --> new discretization
  // (1) Setup() the model evaluator and time integrator
  // (2) read and possibly overwrite the general dynamic state
  // (3) read specific time integrator and model evaluator data
  // ---------------------------------------------------------------------------
  // (0) read element and node data
  ioreader.ReadHistoryData(stepn);

  // (1) Setup() the model evaluator and time integrator
  /* Since we call a redistribution on the structural discretization, we have to
   * setup the structural time integration strategy at this point and not as
   * usually during the adapter call.                         hiermeier 05/16 */
  Setup();

  // (2) read (or overwrite) the general dynamic state
  Teuchos::RCP<Epetra_Vector>& velnp = dataglobalstate_->get_vel_np();
  ioreader.ReadVector(velnp, "velocity");
  dataglobalstate_->get_multi_vel()->UpdateSteps(*velnp);
  Teuchos::RCP<Epetra_Vector>& accnp = dataglobalstate_->get_acc_np();
  ioreader.ReadVector(accnp, "acceleration");
  dataglobalstate_->get_multi_acc()->UpdateSteps(*accnp);

  // (3) read specific time integrator (forces, etc.) and model evaluator data
  int_ptr_->read_restart(ioreader);
  int_ptr_->post_setup();  // compute here the equilibrium system to account for initial
                           // displacement/velocity.

  // short screen output
  if (dataglobalstate_->get_my_rank() == 0)
    Core::IO::cout << "====== Restart of the structural simulation from step " << stepn
                   << Core::IO::endl;

  // end of restarting
  isrestarting_ = false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::FE::Discretization> STR::TimeInt::Base::discretization()
{
  check_init();
  return dataglobalstate_->get_discret();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::SetActionType(const Core::Elements::ActionType& action)
{
  check_init_setup();
  int_ptr_->eval_data().set_action_type(action);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TimeInt::Base::GroupId() const
{
  Teuchos::RCP<Core::Communication::Communicators> group =
      Global::Problem::Instance()->GetCommunicators();
  return group->GroupId();
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TimeInt::Base::post_update() { int_ptr_->post_update(); }

void STR::TimeInt::Base::PostTimeLoop() { int_ptr_->post_time_loop(); }

bool STR::TimeInt::Base::has_final_state_been_written() const
{
  return dataio_->get_last_written_results() == dataglobalstate_->get_step_n();
}

FOUR_C_NAMESPACE_CLOSE
