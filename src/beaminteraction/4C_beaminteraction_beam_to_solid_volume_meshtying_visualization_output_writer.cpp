/*----------------------------------------------------------------------*/
/*! \file

\brief Object to handle beam to solid volume meshtying output creation.

\level 3

*/


#include "4C_beaminteraction_beam_to_solid_volume_meshtying_visualization_output_writer.hpp"

#include "4C_beaminteraction_beam_to_solid_mortar_manager.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_base.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_utils.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"
#include "4C_beaminteraction_beam_to_solid_volume_meshtying_visualization_output_params.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_beaminteraction_str_model_evaluator_datastate.hpp"
#include "4C_beaminteraction_submodel_evaluator_beamcontact.hpp"
#include "4C_beaminteraction_submodel_evaluator_beamcontact_assembly_manager_direct.hpp"
#include "4C_beaminteraction_submodel_evaluator_beamcontact_assembly_manager_indirect.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"

#include <unordered_set>
#include <utility>

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BEAMINTERACTION::BeamToSolidVolumeMeshtyingVisualizationOutputWriter::
    BeamToSolidVolumeMeshtyingVisualizationOutputWriter(
        Core::IO::VisualizationParameters visualization_params)
    : isinit_(false),
      issetup_(false),
      output_params_ptr_(Teuchos::null),
      output_writer_base_ptr_(Teuchos::null),
      visualization_params_(std::move(visualization_params))
{
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingVisualizationOutputWriter::init()
{
  issetup_ = false;
  isinit_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingVisualizationOutputWriter::setup(
    Teuchos::RCP<const Solid::TimeInt::ParamsRuntimeOutput> visualization_output_params,
    Teuchos::RCP<const BEAMINTERACTION::BeamToSolidVolumeMeshtyingVisualizationOutputParams>
        output_params_ptr)
{
  check_init();

  // Set beam to solid volume mesh tying output parameters.
  output_params_ptr_ = output_params_ptr;

  // Initialize the writer base object and add the desired visualizations.
  output_writer_base_ptr_ = Teuchos::rcp<BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase>(
      new BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase(
          "beam-to-solid-volume", visualization_output_params, visualization_params_));

  // Whether or not to write unique cell and node IDs.
  const bool write_unique_ids = output_params_ptr_->get_write_unique_i_ds_flag();

  // Depending on the selected input parameters, create the needed writers. All node / cell data
  // fields that should be output eventually have to be defined here. This helps to prevent issues
  // with ranks that do not contribute to a certain writer.
  {
    if (output_params_ptr_->get_nodal_force_output_flag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->add_visualization_writer("nodal-forces", "btsvc-nodal-forces");
      auto& visualization_data = visualization_writer->get_visualization_data();
      visualization_data.register_point_data<double>("displacement", 3);
      visualization_data.register_point_data<double>("force_beam", 3);
      visualization_data.register_point_data<double>("force_solid", 3);
      if (write_unique_ids) visualization_data.register_point_data<double>("uid_0_node_id", 1);
    }

    if (output_params_ptr_->get_mortar_lambda_discret_output_flag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->add_visualization_writer("mortar", "btsvc-mortar");
      auto& visualization_data = visualization_writer->get_visualization_data();
      visualization_data.register_point_data<double>("displacement", 3);
      visualization_data.register_point_data<double>("lambda", 3);
      if (write_unique_ids)
      {
        visualization_data.register_point_data<double>("uid_0_pair_beam_id", 1);
        visualization_data.register_point_data<double>("uid_1_pair_solid_id", 1);
      }
    }

    if (output_params_ptr_->get_mortar_lambda_continuous_output_flag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->add_visualization_writer(
              "mortar-continuous", "btsvc-mortar-continuous");
      auto& visualization_data = visualization_writer->get_visualization_data();
      visualization_data.register_point_data<double>("displacement", 3);
      visualization_data.register_point_data<double>("lambda", 3);
      if (write_unique_ids)
      {
        visualization_data.register_point_data<double>("uid_0_pair_beam_id", 1);
        visualization_data.register_point_data<double>("uid_1_pair_solid_id", 1);
        visualization_data.register_cell_data<double>("uid_0_pair_beam_id", 1);
        visualization_data.register_cell_data<double>("uid_1_pair_solid_id", 1);
      }
    }

    if (output_params_ptr_->get_integration_points_output_flag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->add_visualization_writer(
              "integration-points", "btsvc-integration-points");
      auto& visualization_data = visualization_writer->get_visualization_data();
      visualization_data.register_point_data<double>("displacement", 3);
      visualization_data.register_point_data<double>("force", 3);
      if (write_unique_ids)
      {
        visualization_data.register_point_data<double>("uid_0_pair_beam_id", 1);
        visualization_data.register_point_data<double>("uid_1_pair_solid_id", 1);
      }
    }

    if (output_params_ptr_->get_segmentation_output_flag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->add_visualization_writer("segmentation", "btsvc-segmentation");
      auto& visualization_data = visualization_writer->get_visualization_data();
      visualization_data.register_point_data<double>("displacement", 3);
      if (write_unique_ids)
      {
        visualization_data.register_point_data<double>("uid_0_pair_beam_id", 1);
        visualization_data.register_point_data<double>("uid_1_pair_solid_id", 1);
      }
    }
  }

  issetup_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingVisualizationOutputWriter::write_output_runtime(
    const BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact* beam_contact) const
{
  check_init_setup();

  // Get the time step and time for the output file. If output is desired at every iteration, the
  // values are padded. The runtime output is written when the time step is already set to the
  // next step.
  auto [output_time, output_step] =
      Core::IO::GetTimeAndTimeStepIndexForOutput(visualization_params_,
          beam_contact->g_state().get_time_n(), beam_contact->g_state().get_step_n());
  write_output_beam_to_solid_volume_mesh_tying(beam_contact, output_step, output_time);
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingVisualizationOutputWriter::
    write_output_runtime_iteration(
        const BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact* beam_contact, int i_iteration) const
{
  check_init_setup();

  if (output_params_ptr_->get_output_every_iteration())
  {
    auto [output_time, output_step] = Core::IO::GetTimeAndTimeStepIndexForOutput(
        visualization_params_, beam_contact->g_state().get_time_n(),
        beam_contact->g_state().get_step_n(), i_iteration);
    write_output_beam_to_solid_volume_mesh_tying(beam_contact, output_step, output_time);
  }
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingVisualizationOutputWriter::
    write_output_beam_to_solid_volume_mesh_tying(
        const BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact* beam_contact, int i_step,
        double time) const
{
  // Parameter list that will be passed to all contact pairs when they create their visualization.
  Teuchos::ParameterList visualization_params;
  visualization_params.set<Teuchos::RCP<const BeamToSolidVolumeMeshtyingVisualizationOutputParams>>(
      "btsvc-output_params_ptr", output_params_ptr_);


  // Add the nodal forces resulting from beam contact. The forces are split up into beam and solid
  // nodes.
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization =
      output_writer_base_ptr_->get_visualization_writer("btsvc-nodal-forces");
  if (visualization != Teuchos::null)
    AddBeamInteractionNodalForces(visualization, beam_contact->discret_ptr(),
        beam_contact->beam_interaction_data_state().get_dis_np(),
        beam_contact->beam_interaction_data_state().get_force_np(),
        output_params_ptr_->get_write_unique_i_ds_flag());


  // Loop over the assembly managers and add the visualization for the pairs contained in the
  // assembly managers.
  for (auto& assembly_manager : beam_contact->get_assembly_managers())
  {
    // Add pair specific output for direct assembly managers.
    auto direct_assembly_manager = Teuchos::rcp_dynamic_cast<
        BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerDirect>(assembly_manager);
    if (not(direct_assembly_manager == Teuchos::null))
    {
      for (const auto& pair : direct_assembly_manager->get_contact_pairs())
        pair->get_pair_visualization(output_writer_base_ptr_, visualization_params);
    }

    // Add pair specific output for indirect assembly managers.
    auto indirect_assembly_manager = Teuchos::rcp_dynamic_cast<
        BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerInDirect>(assembly_manager);
    if (not(indirect_assembly_manager == Teuchos::null))
    {
      // Get the global vector with the Lagrange Multiplier values and add it to the parameter
      // list that will be passed to the pairs.
      Teuchos::RCP<Epetra_Vector> lambda =
          indirect_assembly_manager->get_mortar_manager()->get_global_lambda_col();
      visualization_params.set<Teuchos::RCP<Epetra_Vector>>("lambda", lambda);

      // The pairs will need the mortar manager to extract their Lambda DOFs.
      visualization_params.set<Teuchos::RCP<const BEAMINTERACTION::BeamToSolidMortarManager>>(
          "mortar_manager", indirect_assembly_manager->get_mortar_manager());

      // This map is used to ensure, that each discrete Lagrange multiplier is only written once
      // per beam element.
      Teuchos::RCP<std::unordered_set<int>> beam_tracker =
          Teuchos::rcp(new std::unordered_set<int>());
      visualization_params.set<Teuchos::RCP<std::unordered_set<int>>>("beam_tracker", beam_tracker);

      // Add the pair specific output.
      for (const auto& pair : indirect_assembly_manager->get_mortar_manager()->get_contact_pairs())
        pair->get_pair_visualization(output_writer_base_ptr_, visualization_params);

      // Reset assembly manager specific values in the parameter list passed to the individual
      // pairs.
      visualization_params.remove("lambda");
      visualization_params.remove("mortar_manager");
      visualization_params.remove("beam_tracker");
    }
  }

  // Write the data to disc. The data will be cleared in this method.
  output_writer_base_ptr_->write(i_step, time);
}


/**
 * \brief Checks the init and setup status.
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingVisualizationOutputWriter::check_init_setup() const
{
  if (!isinit_ or !issetup_) FOUR_C_THROW("Call init() and setup() first!");
}

/**
 * \brief Checks the init status.
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingVisualizationOutputWriter::check_init() const
{
  if (!isinit_) FOUR_C_THROW("init() has not been called, yet!");
}

FOUR_C_NAMESPACE_CLOSE
