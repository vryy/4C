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
STR::TIMINT::Base::Base()
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
void STR::TIMINT::Base::Init(const Teuchos::RCP<STR::TIMINT::BaseDataIO> dataio,
    const Teuchos::RCP<STR::TIMINT::BaseDataSDyn> datasdyn,
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> dataglobalstate)
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
void STR::TIMINT::Base::Setup()
{
  CheckInit();

  // ---------------------------------------------------------------------------
  // Create the Dirichlet Boundary Condition handler
  // ---------------------------------------------------------------------------
  dbc_ptr_ = STR::BuildDbc(DataSDyn());
  /* FixMe It would be sufficient to use a constant discretization,
   * unfortunately this wasn't considered during the implementation of the
   * discretization routines. Therefore many methods need a slight modification
   * (most times adding a "const" should fix the problem).          hiermeier */
  Teuchos::RCP<DRT::Discretization> discret_ptr = DataGlobalState().GetDiscret();
  dbc_ptr_->Init(discret_ptr, DataGlobalState().GetFreactNp(), Teuchos::rcp(this, false));
  dbc_ptr_->Setup();

  // ---------------------------------------------------------------------------
  // Create the explicit/implicit integrator
  // ---------------------------------------------------------------------------
  int_ptr_ = STR::BuildIntegrator(DataSDyn());
  int_ptr_->Init(
      DataSDynPtr(), DataGlobalStatePtr(), DataIOPtr(), dbc_ptr_, Teuchos::rcp(this, false));
  int_ptr_->Setup();
  int_ptr_->PostSetup();
  // Initialize and Setup the input/output writer for every Newton iteration
  dataio_->init_setup_every_iteration_writer(this, DataSDyn().GetNoxParams());

  // Initialize the output of system energy
  if (dataio_->get_write_energy_every_n_step())
  {
    select_energy_types_to_be_written();

    if (dataglobalstate_->GetMyRank() == 0) initialize_energy_file_stream_and_write_headers();
  }

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::Reset()
{
  FOUR_C_THROW(
      "Reset of all class variables is not yet implemented for "
      "the modelevaluator!");
  // ModelEvaluator().Reset();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::ResetStep()
{
  CheckInitSetup();

  int_ptr_->ResetStepState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::Base::NotFinished() const
{
  CheckInitSetup();

  // check for early stopping
  const bool early_stop = int_ptr_->EarlyStopping();
  if (early_stop)
  {
    // Simulation is finished regardless of the simulation time or the time step.
    return false;
  }

  // check the current time
  const double& timenp = dataglobalstate_->GetTimeNp();
  const double& timemax = datasdyn_->GetTimeMax();
  const double& dt = (*dataglobalstate_->GetDeltaTime())[0];
  // check the step counter
  const int& stepnp = dataglobalstate_->GetStepNp();
  const int& stepmax = datasdyn_->GetStepMax();

  return (timenp <= timemax + 1.0e-8 * dt and stepnp <= stepmax);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::SetRestart(int stepn, double timen, Teuchos::RCP<Epetra_Vector> disn,
    Teuchos::RCP<Epetra_Vector> veln, Teuchos::RCP<Epetra_Vector> accn,
    Teuchos::RCP<std::vector<char>> elementdata, Teuchos::RCP<std::vector<char>> nodedata)
{
  CheckInitSetup();

  FOUR_C_THROW("SetRestartState() is deprecated, use the ReadRestart() routine instead!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Map& STR::TIMINT::Base::GetMassDomainMap() const
{
  CheckInitSetup();
  return dataglobalstate_->GetMassMatrix()->DomainMap();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::MapExtractor> STR::TIMINT::Base::GetDBCMapExtractor()
{
  CheckInitSetup();
  return dbc_ptr_->GetDBCMapExtractor();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::MapExtractor> STR::TIMINT::Base::GetDBCMapExtractor() const
{
  CheckInitSetup();
  return dbc_ptr_->GetDBCMapExtractor();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::UTILS::LocsysManager> STR::TIMINT::Base::LocsysManager()
{
  CheckInitSetup();
  return dbc_ptr_->LocSysManagerPtr();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::MODELEVALUATOR::Generic& STR::TIMINT::Base::ModelEvaluator(
    INPAR::STR::ModelType mtype) const
{
  return Integrator().ModelEval().Evaluator(mtype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Generic& STR::TIMINT::Base::ModelEvaluator(INPAR::STR::ModelType mtype)
{
  return Integrator().ModelEval().Evaluator(mtype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::TIMINT::Base::TimIntParam() const
{
  CheckInitSetup();
  return int_ptr_->GetIntParam();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::ResizeMStepTimAda()
{
  CheckInitSetup();
  // resize time and stepsize fields
  const double& timen = dataglobalstate_->GetTimeN();
  dataglobalstate_->GetMultiTime()->Resize(-1, 0, timen);
  const double& dtn = (*dataglobalstate_->GetDeltaTime())[0];
  dataglobalstate_->GetDeltaTime()->Resize(-1, 0, dtn);

  // resize state vectors, AB2 is a 2-step method, thus we need two
  // past steps at t_{n} and t_{n-1}
  const Epetra_Map* dofrowmap_ptr = dataglobalstate_->DofRowMapView();
  dataglobalstate_->GetMultiDis()->Resize(-1, 0, dofrowmap_ptr, true);
  dataglobalstate_->GetMultiVel()->Resize(-1, 0, dofrowmap_ptr, true);
  dataglobalstate_->GetMultiAcc()->Resize(-1, 0, dofrowmap_ptr, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::Update()
{
  CheckInitSetup();
  int_ptr_->PreUpdate();
  int_ptr_->update_structural_energy();
  int_ptr_->UpdateStepState();
  UpdateStepTime();
  set_number_of_nonlinear_iterations();
  int_ptr_->UpdateStepElement();
  int_ptr_->PostUpdate();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::UpdateStepTime()
{
  CheckInitSetup();
  double& timenp = dataglobalstate_->GetTimeNp();
  int& stepnp = dataglobalstate_->GetStepNp();
  int& stepn = dataglobalstate_->GetStepN();

  // --------------------------------------------------------------------------
  // update old time and step variables
  // --------------------------------------------------------------------------
  dataglobalstate_->GetMultiTime()->UpdateSteps(timenp);
  stepn = stepnp;

  // --------------------------------------------------------------------------
  // update the new time and step variables
  // --------------------------------------------------------------------------
  // get current time step size
  const double& dtn = (*dataglobalstate_->GetDeltaTime())[0];
  dataglobalstate_->GetDeltaTime()->UpdateSteps(dtn);
  timenp += dtn;
  stepnp += 1;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::set_number_of_nonlinear_iterations()
{
  int nlniter = 0;

  if (DataSDyn().GetNoxParams().isSublist("Output"))
  {
    const Teuchos::ParameterList& nox_output = DataSDyn().GetNoxParams().sublist("Output");
    if (nox_output.isParameter("Nonlinear Iterations"))
      nlniter = nox_output.get<int>("Nonlinear Iterations");
  }

  dataglobalstate_->set_nln_iteration_number(nlniter);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::select_energy_types_to_be_written()
{
  STR::MODELEVALUATOR::Data& evaldata = int_ptr_->EvalData();

  // decide which types of energy contributions shall be written separately
  const std::set<enum INPAR::STR::ModelType>& mtypes = datasdyn_->GetModelTypes();

  std::set<enum INPAR::STR::ModelType>::const_iterator model_iter;
  for (model_iter = mtypes.begin(); model_iter != mtypes.end(); ++model_iter)
  {
    switch (*model_iter)
    {
      case INPAR::STR::model_structure:
      {
        evaldata.insert_energy_type_to_be_considered(STR::internal_energy);
        evaldata.insert_energy_type_to_be_considered(STR::kinetic_energy);
        break;
      }
      case INPAR::STR::model_beaminteraction:
      {
        STR::MODELEVALUATOR::BeamInteraction const beaminteraction_evaluator =
            dynamic_cast<STR::MODELEVALUATOR::BeamInteraction const&>(
                int_ptr_->ModelEvalPtr()->Evaluator(INPAR::STR::model_beaminteraction));

        if (beaminteraction_evaluator.HaveSubModelType(
                INPAR::BEAMINTERACTION::submodel_beamcontact))
        {
          evaldata.insert_energy_type_to_be_considered(STR::beam_contact_penalty_potential);
        }
        if (beaminteraction_evaluator.HaveSubModelType(INPAR::BEAMINTERACTION::submodel_potential))
        {
          evaldata.insert_energy_type_to_be_considered(STR::beam_interaction_potential);
        }
        if (beaminteraction_evaluator.HaveSubModelType(
                INPAR::BEAMINTERACTION::submodel_crosslinking))
        {
          evaldata.insert_energy_type_to_be_considered(STR::beam_to_beam_link_internal_energy);
          evaldata.insert_energy_type_to_be_considered(STR::beam_to_beam_link_kinetic_energy);
        }
        if (beaminteraction_evaluator.HaveSubModelType(
                INPAR::BEAMINTERACTION::submodel_spherebeamlink))
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
void STR::TIMINT::Base::initialize_energy_file_stream_and_write_headers()
{
  auto& evaldata = int_ptr_->EvalData();

  dataio_->setup_energy_output_file();

  // write column headers to file
  dataio_->get_energy_output_stream() << std::setw(12) << "#timestep," << std::setw(24) << "time,";

  for (const auto& energy_data : evaldata.GetEnergyData())
  {
    dataio_->get_energy_output_stream()
        << std::setw(36) << STR::EnergyType2String(energy_data.first) + ",";
  }

  dataio_->get_energy_output_stream() << std::setw(24) << "total_energy" << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::UTILS::ResultTest> STR::TIMINT::Base::CreateFieldTest()
{
  CheckInitSetup();
  Teuchos::RCP<STR::ResultTest> resulttest = Teuchos::rcp(new STR::ResultTest());
  resulttest->Init(GetDataGlobalState(), Integrator().EvalData());
  resulttest->Setup();

  return resulttest;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::GetRestartData(Teuchos::RCP<int> step, Teuchos::RCP<double> time,
    Teuchos::RCP<Epetra_Vector> disnp, Teuchos::RCP<Epetra_Vector> velnp,
    Teuchos::RCP<Epetra_Vector> accnp, Teuchos::RCP<std::vector<char>> elementdata,
    Teuchos::RCP<std::vector<char>> nodedata)
{
  CheckInitSetup();
  // at some point we have to create a copy
  *step = dataglobalstate_->GetStepN();
  *time = dataglobalstate_->GetTimeN();
  Teuchos::RCP<const DRT::Discretization> discret_ptr =
      Teuchos::rcp_dynamic_cast<const DRT::Discretization>(dataglobalstate_->GetDiscret());
  *elementdata = *(discret_ptr->PackMyElements());
  *nodedata = *(discret_ptr->PackMyNodes());

  // get restart data is only for simple structure problems
  // hence if the model set is larger than one, we throw an error
  if (datasdyn_->GetModelTypes().size() > 1)
    FOUR_C_THROW("The GetRestartData routine supports the structural model case ONLY!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::PrepareOutput(bool force_prepare_timestep)
{
  CheckInitSetup();
  // --- stress, strain and optional quantity calculation ---------------------
  if ((dataio_->is_write_results_enabled() && force_prepare_timestep) ||
      dataio_->write_results_for_this_step(dataglobalstate_->GetStepNp()))
  {
    int_ptr_->determine_stress_strain();
    int_ptr_->determine_optional_quantity();

    if (dataio_->is_write_current_ele_volume())
    {
      Teuchos::RCP<Epetra_Vector> elevolumes = Teuchos::null;
      Teuchos::RCP<const Epetra_Vector> disnp = dataglobalstate_->GetDisNp();

      int_ptr_->determine_element_volumes(*disnp, elevolumes);
      int_ptr_->EvalData().set_element_volume_data(elevolumes);
    }
  }
  if ((dataio_->is_runtime_output_enabled() && force_prepare_timestep) ||
      dataio_->write_runtime_vtk_results_for_this_step(dataglobalstate_->GetStepNp()) ||
      dataio_->write_runtime_vtp_results_for_this_step(dataglobalstate_->GetStepNp()))
  {
    int_ptr_->runtime_pre_output_step_state();
  }
  // --- energy calculation ---------------------------------------------------
  if ((dataio_->get_write_energy_every_n_step() and
          (force_prepare_timestep ||
              dataglobalstate_->GetStepNp() % dataio_->get_write_energy_every_n_step() == 0)))
  {
    STR::MODELEVALUATOR::Data& evaldata = int_ptr_->EvalData();
    evaldata.clear_values_for_all_energy_types();

    int_ptr_->DetermineEnergy();

    // sum processor-local values of all separate contributions into global value
    double energy_local = 0.0;
    double energy_global = 0.0;

    for (const auto& energy_data : evaldata.GetEnergyData())
    {
      energy_local = energy_data.second;

      dataglobalstate_->GetComm().SumAll(&energy_local, &energy_global, 1);

      evaldata.set_value_for_energy_type(energy_global, energy_data.first);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::Output(bool forced_writerestart)
{
  CheckInitSetup();
  OutputStep(forced_writerestart);
  // write Gmsh output
  write_gmsh_struc_output_step();
  int_ptr_->PostOutput();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputStep(bool forced_writerestart)
{
  CheckInitSetup();
  // special treatment is necessary when restart is forced
  if (forced_writerestart)
  {
    // reset possible history data on element level
    ResetStep();
    // restart has already been written or simulation has just started
    if (dataio_->should_write_restart_for_step(dataglobalstate_->GetStepN()) or
        dataglobalstate_->GetStepN() == GLOBAL::Problem::Instance()->Restart())
      return;
    // if state already exists, add restart information
    if (dataio_->write_results_for_this_step(dataglobalstate_->GetStepN()))
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
  if (forced_writerestart || dataio_->should_write_restart_for_step(dataglobalstate_->GetStepN()))
  {
    OutputRestart(datawritten);
    dataio_->set_last_written_results(dataglobalstate_->GetStepN());
  }

  // output results (not necessary if restart in same step)
  if (dataio_->IsWriteState() and
      dataio_->write_results_for_this_step(dataglobalstate_->GetStepN()) and (not datawritten))
  {
    NewIOStep(datawritten);
    OutputState();
    dataio_->set_last_written_results(dataglobalstate_->GetStepN());
  }

  // output results during runtime ( not used for restart so far )
  if (dataio_->write_runtime_vtk_results_for_this_step(dataglobalstate_->GetStepN()) or
      dataio_->write_runtime_vtp_results_for_this_step(dataglobalstate_->GetStepN()))
  {
    RuntimeOutputState();
  }

  // write reaction forces
  if (dataio_->should_write_reaction_forces_for_this_step(dataglobalstate_->GetStepN()))
  {
    output_reaction_forces();
  }

  // output stress, strain and optional quantity
  if (dataio_->should_write_stress_strain_for_this_step(dataglobalstate_->GetStepN()))
  {
    NewIOStep(datawritten);
    OutputStressStrain();
    output_optional_quantity();
  }

  if (dataio_->write_results_for_this_step(dataglobalstate_->GetStepN()) and
      dataio_->is_write_current_ele_volume())
  {
    NewIOStep(datawritten);
    IO::DiscretizationWriter& iowriter = *(dataio_->GetOutputPtr());
    OutputElementVolume(iowriter);
  }

  // output energy
  if (dataio_->should_write_energy_for_this_step(dataglobalstate_->GetStepN()))
  {
    OutputEnergy();
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
void STR::TIMINT::Base::NewIOStep(bool& datawritten)
{
  if (not datawritten)
  {
    // Make new step
    dataio_->GetOutputPtr()->NewStep(dataglobalstate_->GetStepN(), dataglobalstate_->GetTimeN());

    datawritten = true;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputState()
{
  CheckInitSetup();
  IO::DiscretizationWriter& iowriter = *(dataio_->GetOutputPtr());

  OutputState(iowriter, dataio_->IsFirstOutputOfRun());

  dataio_->SetFirstOutputOfRun(false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputDebugState(IO::DiscretizationWriter& iowriter, bool write_owner) const
{
  OutputState(iowriter, write_owner);

  // write element volumes as additional debugging information, if activated
  if (dataio_->is_write_current_ele_volume()) OutputElementVolume(iowriter);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputState(IO::DiscretizationWriter& iowriter, bool write_owner) const
{
  // owner of elements is just written once because it does not change during
  // simulation (so far)
  iowriter.WriteElementData(write_owner);
  iowriter.WriteNodeData(write_owner);

  int_ptr_->OutputStepState(iowriter);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::RuntimeOutputState()
{
  CheckInitSetup();
  int_ptr_->runtime_output_step_state();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::output_reaction_forces()
{
  CheckInitSetup();
  IO::DiscretizationWriter& iowriter = *(dataio_->GetOutputPtr());
  int_ptr_->MonitorDbc(iowriter);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputElementVolume(IO::DiscretizationWriter& iowriter) const
{
  CheckInitSetup();

  STR::MODELEVALUATOR::Data& evaldata = int_ptr_->EvalData();

  iowriter.WriteVector("current_ele_volumes",
      Teuchos::rcpFromRef(evaldata.current_element_volume_data()), IO::elementvector);

  evaldata.set_element_volume_data(Teuchos::null);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputStressStrain()
{
  CheckInitSetup();

  STR::MODELEVALUATOR::Data& evaldata = int_ptr_->EvalData();
  Teuchos::RCP<IO::DiscretizationWriter> output_ptr = dataio_->GetOutputPtr();

  // ---------------------------------------------------------------------------
  // write stress output
  // ---------------------------------------------------------------------------
  std::string text = "";
  if (dataio_->GetStressOutputType() != INPAR::STR::stress_none)
  {
    switch (dataio_->GetStressOutputType())
    {
      case INPAR::STR::stress_cauchy:
        text = "gauss_cauchy_stresses_xyz";
        break;
      case INPAR::STR::stress_2pk:
        text = "gauss_2PK_stresses_xyz";
        break;
      default:
        FOUR_C_THROW("Requested stress type is not supported!");
        break;
    }
    output_ptr->WriteVector(text, evaldata.StressData(), *(Discretization()->ElementRowMap()));
  }
  // we don't need this anymore
  evaldata.StressDataPtr() = Teuchos::null;

  // ---------------------------------------------------------------------------
  // write coupling stress output
  // ---------------------------------------------------------------------------
  text.clear();
  if (dataio_->get_coupling_stress_output_type() != INPAR::STR::stress_none)
  {
    switch (dataio_->get_coupling_stress_output_type())
    {
      case INPAR::STR::stress_cauchy:
        text = "gauss_cauchy_coupling_stresses_xyz";
        break;
      case INPAR::STR::stress_2pk:
        text = "gauss_2PK_coupling_stresses_xyz";
        break;
      default:
        FOUR_C_THROW("Requested coupling stress type is not supported!");
        break;
    }
    output_ptr->WriteVector(
        text, evaldata.CouplingStressData(), *(Discretization()->ElementRowMap()));
  }

  // ---------------------------------------------------------------------------
  // write strain output
  // ---------------------------------------------------------------------------
  text.clear();
  if (dataio_->GetStrainOutputType() != INPAR::STR::strain_none)
  {
    switch (dataio_->GetStrainOutputType())
    {
      case INPAR::STR::strain_none:
        break;
      case INPAR::STR::strain_ea:
        text = "gauss_EA_strains_xyz";
        break;
      case INPAR::STR::strain_gl:
        text = "gauss_GL_strains_xyz";
        break;
      case INPAR::STR::strain_log:
        text = "gauss_LOG_strains_xyz";
        break;
      default:
        FOUR_C_THROW("Requested strain type is not supported!");
        break;
    }
    output_ptr->WriteVector(text, evaldata.StrainData(), *(Discretization()->ElementRowMap()));
  }
  // we don't need this anymore
  evaldata.StrainDataPtr() = Teuchos::null;

  // ---------------------------------------------------------------------------
  // write plastic strain output
  // ---------------------------------------------------------------------------
  text.clear();
  if (dataio_->get_plastic_strain_output_type() != INPAR::STR::strain_none)
  {
    switch (dataio_->get_plastic_strain_output_type())
    {
      case INPAR::STR::strain_ea:
        text = "gauss_pl_EA_strains_xyz";
        break;
      case INPAR::STR::strain_gl:
        text = "gauss_pl_GL_strains_xyz";
        break;
      default:
        FOUR_C_THROW("Requested plastic strain type is not supported!");
        break;
    }
    output_ptr->WriteVector(
        text, evaldata.PlasticStrainData(), *(Discretization()->ElementRowMap()));
  }
  // we don't need this anymore
  evaldata.plastic_strain_data_ptr() = Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputEnergy() const
{
  CheckInitSetup();

  if (dataglobalstate_->GetMyRank() == 0)
  {
    std::ostream& energy_output_stream = dataio_->get_energy_output_stream();

    energy_output_stream << std::setw(11) << dataglobalstate_->GetStepN() << std::setw(1) << ","
                         << std::scientific << std::setprecision(14) << std::setw(23)
                         << dataglobalstate_->GetTimeN() << std::setw(1) << ",";

    STR::MODELEVALUATOR::Data& evaldata = int_ptr_->EvalData();

    double total_energy = 0.0;

    for (const auto& energy_data : evaldata.GetEnergyData())
    {
      energy_output_stream << std::setw(35) << energy_data.second << std::setw(1) << ",";
      total_energy += energy_data.second;
    }

    energy_output_stream << std::setw(24) << total_energy << std::endl;

    IO::cout(IO::verbose) << "\n\nOutput for energy written to file!" << IO::endl;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::output_optional_quantity()
{
  CheckInitSetup();

  STR::MODELEVALUATOR::Data& evaldata = int_ptr_->EvalData();
  Teuchos::RCP<IO::DiscretizationWriter> output_ptr = dataio_->GetOutputPtr();

  // ---------------------------------------------------------------------------
  // write optional quantity output
  // ---------------------------------------------------------------------------
  std::string text = "";
  if (dataio_->get_opt_quantity_output_type() != INPAR::STR::optquantity_none)
  {
    switch (dataio_->get_opt_quantity_output_type())
    {
      case INPAR::STR::optquantity_membranethickness:
        text = "gauss_membrane_thickness";
        break;
      default:
        FOUR_C_THROW("Requested optional quantity type is not supported!");
        break;
    }
    output_ptr->WriteVector(text, evaldata.OptQuantityData(), *(Discretization()->ElementRowMap()));
  }
  // we don't need this anymore
  evaldata.OptQuantityDataPtr() = Teuchos::null;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::OutputRestart(bool& datawritten)
{
  CheckInitSetup();

  Teuchos::RCP<IO::DiscretizationWriter> output_ptr = dataio_->GetOutputPtr();
  // write restart output, please
  if (dataglobalstate_->GetStepN() != 0)
    output_ptr->WriteMesh(dataglobalstate_->GetStepN(), dataglobalstate_->GetTimeN());
  NewIOStep(datawritten);

  output_ptr->WriteElementData(dataio_->IsFirstOutputOfRun());
  output_ptr->WriteNodeData(dataio_->IsFirstOutputOfRun());
  dataio_->SetFirstOutputOfRun(false);

  // add velocity and acceleration if necessary
  output_ptr->WriteVector("velocity", dataglobalstate_->GetVelN());
  output_ptr->WriteVector("acceleration", dataglobalstate_->GetAccN());

  /* Add the restart information of the different time integrators and model
   * evaluators. */
  int_ptr_->WriteRestart(*output_ptr);

  // info dedicated to user's eyes staring at standard out
  if ((dataglobalstate_->GetMyRank() == 0) and (dataio_->get_print2_screen_every_n_step() > 0) and
      (StepOld() % dataio_->get_print2_screen_every_n_step() == 0))
  {
    IO::cout << "====== Restart for field 'Structure' written in step "
             << dataglobalstate_->GetStepN() << IO::endl;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::add_restart_to_output_state()
{
  Teuchos::RCP<IO::DiscretizationWriter> output_ptr = dataio_->GetOutputPtr();

  // force output of velocity and acceleration in case it is not written previously by the model
  // evaluators
  if (!dataio_->IsWriteVelAcc())
  {
    output_ptr->WriteVector("velocity", dataglobalstate_->GetVelN());
    output_ptr->WriteVector("acceleration", dataglobalstate_->GetAccN());
  }

  /* Add the restart information of the different time integrators and model
   * evaluators. */
  int_ptr_->WriteRestart(*output_ptr, true);

  // finally add the missing mesh information, order is important here
  output_ptr->WriteMesh(DataGlobalState().GetStepN(), DataGlobalState().GetTimeN());

  // info dedicated to user's eyes staring at standard out
  if ((dataglobalstate_->GetMyRank() == 0) and (dataio_->get_print2_screen_every_n_step() > 0) and
      (StepOld() % dataio_->get_print2_screen_every_n_step() == 0))
  {
    IO::cout << "====== Restart for field 'Structure' written in step "
             << dataglobalstate_->GetStepN() << IO::endl;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::write_gmsh_struc_output_step()
{
  CheckInitSetup();
  if (!dataio_->IsGmsh()) return;

  const std::string filename = IO::GMSH::GetFileName("struct", DiscWriter()->Output()->FileName(),
      dataglobalstate_->GetStepNp(), false, dataglobalstate_->GetMyRank());
  std::ofstream gmshfilecontent(filename.c_str());

  // add 'View' to Gmsh postprocessing file
  gmshfilecontent << "View \" "
                  << "struct displacement \" {" << std::endl;
  // draw vector field 'struct displacement' for every element
  IO::GMSH::VectorFieldDofBasedToGmsh(Discretization(), Dispn(), gmshfilecontent, 0, true);
  gmshfilecontent << "};" << std::endl;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::ReadRestart(const int stepn)
{
  CheckInit();
  // that restarting flag
  isrestarting_ = true;

  // create an input/output reader
  IO::DiscretizationReader ioreader(
      Discretization(), GLOBAL::Problem::Instance()->InputControlFile(), stepn);
  dataglobalstate_->GetStepN() = stepn;
  dataglobalstate_->GetStepNp() = stepn + 1;
  dataglobalstate_->GetMultiTime() =
      Teuchos::rcp(new TIMESTEPPING::TimIntMStep<double>(0, 0, ioreader.ReadDouble("time")));
  const double& timen = dataglobalstate_->GetTimeN();
  const double& dt = (*dataglobalstate_->GetDeltaTime())[0];
  dataglobalstate_->GetTimeNp() = timen + dt;

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
  Teuchos::RCP<Epetra_Vector>& velnp = dataglobalstate_->GetVelNp();
  ioreader.ReadVector(velnp, "velocity");
  dataglobalstate_->GetMultiVel()->UpdateSteps(*velnp);
  Teuchos::RCP<Epetra_Vector>& accnp = dataglobalstate_->GetAccNp();
  ioreader.ReadVector(accnp, "acceleration");
  dataglobalstate_->GetMultiAcc()->UpdateSteps(*accnp);

  // (3) read specific time integrator (forces, etc.) and model evaluator data
  int_ptr_->ReadRestart(ioreader);
  int_ptr_->PostSetup();  // compute here the equilibrium system to account for initial
                          // displacement/velocity.

  // short screen output
  if (dataglobalstate_->GetMyRank() == 0)
    IO::cout << "====== Restart of the structural simulation from step " << stepn << IO::endl;

  // end of restarting
  isrestarting_ = false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> STR::TIMINT::Base::Discretization()
{
  CheckInit();
  return dataglobalstate_->GetDiscret();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::SetActionType(const DRT::ELEMENTS::ActionType& action)
{
  CheckInitSetup();
  int_ptr_->EvalData().SetActionType(action);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::TIMINT::Base::GroupId() const
{
  Teuchos::RCP<CORE::COMM::Communicators> group = GLOBAL::Problem::Instance()->GetCommunicators();
  return group->GroupId();
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::Base::PostUpdate() { int_ptr_->PostUpdate(); }

void STR::TIMINT::Base::PostTimeLoop() { int_ptr_->PostTimeLoop(); }

bool STR::TIMINT::Base::has_final_state_been_written() const
{
  return dataio_->get_last_written_results() == dataglobalstate_->GetStepN();
}

FOUR_C_NAMESPACE_CLOSE
