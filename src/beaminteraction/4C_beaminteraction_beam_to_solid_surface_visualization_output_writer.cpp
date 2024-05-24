/*----------------------------------------------------------------------*/
/*! \file

\brief Object to handle beam to solid surface output creation.

\level 3

*/


#include "4C_beaminteraction_beam_to_solid_surface_visualization_output_writer.hpp"

#include "4C_beaminteraction_beam_to_solid_conditions.hpp"
#include "4C_beaminteraction_beam_to_solid_mortar_manager.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_visualization_output_params.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_base.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_utils.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_beaminteraction_str_model_evaluator_datastate.hpp"
#include "4C_beaminteraction_submodel_evaluator_beamcontact.hpp"
#include "4C_beaminteraction_submodel_evaluator_beamcontact_assembly_manager_direct.hpp"
#include "4C_beaminteraction_submodel_evaluator_beamcontact_assembly_manager_indirect.hpp"
#include "4C_geometry_pair_line_to_surface_evaluation_data.hpp"
#include "4C_io_visualization_parameters.hpp"
#include "4C_lib_discret.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"

#include <Epetra_FEVector.h>
#include <Epetra_MpiComm.h>

#include <unordered_set>
#include <utility>

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputWriter::
    BeamToSolidSurfaceVisualizationOutputWriter(IO::VisualizationParameters visualization_params)
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
void BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputWriter::Init()
{
  issetup_ = false;
  isinit_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputWriter::Setup(
    Teuchos::RCP<const STR::TIMINT::ParamsRuntimeOutput> visualization_output_params,
    Teuchos::RCP<const BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputParams>
        output_params_ptr)
{
  check_init();

  // Set beam to solid surface interactions output parameters.
  output_params_ptr_ = output_params_ptr;

  // Initialize the writer base object and add the desired visualizations.
  output_writer_base_ptr_ = Teuchos::rcp<BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase>(
      new BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase(
          "beam-to-solid-surface", visualization_output_params, visualization_params_));

  // Whether or not to write unique cell and node IDs.
  const bool write_unique_ids = output_params_ptr_->get_write_unique_i_ds_flag();

  // Depending on the selected input parameters, create the needed writers. All node / cell data
  // fields that should be output eventually have to be defined here. This helps to prevent issues
  // with ranks that do not contribute to a certain writer.
  {
    if (output_params_ptr_->get_nodal_force_output_flag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->add_visualization_writer("nodal-forces", "btssc-nodal-forces");
      auto& visualization_data = visualization_writer->get_visualization_data();
      visualization_data.RegisterPointData<double>("displacement", 3);
      visualization_data.RegisterPointData<double>("force_beam", 3);
      visualization_data.RegisterPointData<double>("force_solid", 3);
      visualization_data.RegisterFieldData<double>("sum_coupling_force_beam");
      visualization_data.RegisterFieldData<double>("sum_coupling_moment_beam");
      visualization_data.RegisterFieldData<double>("sum_coupling_force_solid");
      visualization_data.RegisterFieldData<double>("sum_coupling_moment_solid");
      if (write_unique_ids) visualization_data.RegisterPointData<double>("uid_0_node_id", 1);
    }

    if (output_params_ptr_->get_averaged_normals_output_flag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->add_visualization_writer(
              "averaged-normals", "btssc-averaged-normals");
      auto& visualization_data = visualization_writer->get_visualization_data();
      visualization_data.RegisterPointData<double>("displacement", 3);
      visualization_data.RegisterPointData<double>("normal_averaged", 3);
      visualization_data.RegisterPointData<double>("normal_element", 3);
      visualization_data.RegisterPointData<double>("coupling_id", 1);
      if (write_unique_ids) visualization_data.RegisterPointData<double>("uid_0_face_id", 1);
    }

    if (output_params_ptr_->get_mortar_lambda_discret_output_flag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->add_visualization_writer("mortar", "btssc-mortar");
      auto& visualization_data = visualization_writer->get_visualization_data();
      visualization_data.RegisterPointData<double>("displacement", 3);
      visualization_data.RegisterPointData<double>("lambda", 3);
      if (write_unique_ids)
      {
        visualization_data.RegisterPointData<double>("uid_0_pair_beam_id", 1);
        visualization_data.RegisterPointData<double>("uid_1_pair_solid_id", 1);
      }
    }

    if (output_params_ptr_->get_mortar_lambda_continuous_output_flag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->add_visualization_writer(
              "mortar-continuous", "btssc-mortar-continuous");
      auto& visualization_data = visualization_writer->get_visualization_data();
      visualization_data.RegisterPointData<double>("displacement", 3);
      visualization_data.RegisterPointData<double>("lambda", 3);
      if (write_unique_ids)
      {
        visualization_data.RegisterPointData<double>("uid_0_pair_beam_id", 1);
        visualization_data.RegisterPointData<double>("uid_1_pair_solid_id", 1);
        visualization_data.RegisterCellData<double>("uid_0_pair_beam_id", 1);
        visualization_data.RegisterCellData<double>("uid_1_pair_solid_id", 1);
      }
    }

    if (output_params_ptr_->get_integration_points_output_flag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->add_visualization_writer(
              "integration-points", "btssc-integration-points");
      auto& visualization_data = visualization_writer->get_visualization_data();
      visualization_data.RegisterPointData<double>("displacement", 3);
      visualization_data.RegisterPointData<double>("projection_direction", 3);
      if (write_unique_ids)
      {
        visualization_data.RegisterPointData<double>("uid_0_pair_beam_id", 1);
        visualization_data.RegisterPointData<double>("uid_1_pair_solid_id", 1);
      }
    }

    if (output_params_ptr_->get_segmentation_output_flag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->add_visualization_writer("segmentation", "btssc-segmentation");
      auto& visualization_data = visualization_writer->get_visualization_data();
      visualization_data.RegisterPointData<double>("displacement", 3);
      visualization_data.RegisterPointData<double>("projection_direction", 3);
      if (write_unique_ids)
      {
        visualization_data.RegisterPointData<double>("uid_0_pair_beam_id", 1);
        visualization_data.RegisterPointData<double>("uid_1_pair_solid_id", 1);
      }
    }
  }

  issetup_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputWriter::write_output_runtime(
    const BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact* beam_contact) const
{
  check_init_setup();

  // Get the time step and time for the output file. If output is desired at every iteration, the
  // values are padded. The runtime output is written when the time step is already set to the next
  // step.
  auto [output_time, output_step] = IO::GetTimeAndTimeStepIndexForOutput(
      visualization_params_, beam_contact->GState().GetTimeN(), beam_contact->GState().GetStepN());
  write_output_beam_to_solid_surface(beam_contact, output_step, output_time);
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputWriter::write_output_runtime_iteration(
    const BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact* beam_contact, int i_iteration) const
{
  check_init_setup();

  if (output_params_ptr_->get_output_every_iteration())
  {
    auto [output_time, output_step] = IO::GetTimeAndTimeStepIndexForOutput(visualization_params_,
        beam_contact->GState().GetTimeN(), beam_contact->GState().GetStepN(), i_iteration);
    write_output_beam_to_solid_surface(beam_contact, output_step, output_time);
  }
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputWriter::
    write_output_beam_to_solid_surface(
        const BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact* beam_contact, int i_step,
        double time) const
{
  // Parameter list that will be passed to all contact pairs when they create their visualization.
  Teuchos::ParameterList visualization_params;
  visualization_params.set<Teuchos::RCP<const BeamToSolidSurfaceVisualizationOutputParams>>(
      "btssc-output_params_ptr", output_params_ptr_);


  // Add the averaged nodal normal output.
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization>
      visualization_averaged_normals =
          output_writer_base_ptr_->get_visualization_writer("btssc-averaged-normals");
  if (visualization_averaged_normals != Teuchos::null)
  {
    const std::vector<Teuchos::RCP<BeamInteractionConditionBase>>& surface_condition_vector =
        beam_contact->GetConditions()->GetConditionMap().at(
            INPAR::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_surface_meshtying);
    for (const auto& condition : surface_condition_vector)
    {
      // Get the line-to-surface evaluation data for the current condition.
      auto beam_to_surface_condition =
          Teuchos::rcp_dynamic_cast<const BeamToSolidConditionSurface>(condition, true);
      Teuchos::RCP<const GEOMETRYPAIR::LineToSurfaceEvaluationData> surface_evaluation_data =
          Teuchos::rcp_dynamic_cast<const GEOMETRYPAIR::LineToSurfaceEvaluationData>(
              beam_to_surface_condition->get_geometry_evaluation_data(), true);

      // Get the coupling ID for the current condition.
      const int coupling_id =
          beam_to_surface_condition->GetOtherCondition()->parameters().Get<int>("COUPLING_ID");

      // Create the output for the averaged normal field.
      AddAveragedNodalNormals(visualization_averaged_normals,
          surface_evaluation_data->GetFaceElements(), coupling_id,
          output_params_ptr_->get_write_unique_i_ds_flag());
    }
  }


  // Add the nodal forces resulting from beam contact. The forces are split up into beam and
  // solid nodes.
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> nodal_force_visualization =
      output_writer_base_ptr_->get_visualization_writer("btssc-nodal-forces");
  if (nodal_force_visualization != Teuchos::null)
    AddBeamInteractionNodalForces(nodal_force_visualization, beam_contact->DiscretPtr(),
        beam_contact->beam_interaction_data_state().GetDisNp(),
        beam_contact->beam_interaction_data_state().GetForceNp(),
        output_params_ptr_->get_write_unique_i_ds_flag());


  // Loop over the assembly managers and add the visualization for the pairs contained in the
  // assembly managers.
  for (auto& assembly_manager : beam_contact->GetAssemblyManagers())
  {
    // Add pair specific output for direct assembly managers.
    auto direct_assembly_manager = Teuchos::rcp_dynamic_cast<
        BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerDirect>(assembly_manager);
    if (not(direct_assembly_manager == Teuchos::null))
    {
      for (const auto& pair : direct_assembly_manager->GetContactPairs())
        pair->get_pair_visualization(output_writer_base_ptr_, visualization_params);
    }

    // Add pair specific output for indirect assembly managers.
    auto indirect_assembly_manager = Teuchos::rcp_dynamic_cast<
        BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerInDirect>(assembly_manager);
    if (not(indirect_assembly_manager == Teuchos::null))
    {
      // If needed, setup the vector with the global moments around the origin.
      if (nodal_force_visualization != Teuchos::null)
      {
        // This array will hold the global coupling moment around the origin.
        auto global_coupling_moment_origin =
            Teuchos::rcp(new CORE::LINALG::Matrix<3, 1, double>(true));
        visualization_params.set("global_coupling_moment_origin", global_coupling_moment_origin);
      }

      // Get the global vector with the Lagrange Multiplier values and add it to the parameter list
      // that will be passed to the pairs.
      Teuchos::RCP<Epetra_Vector> lambda =
          indirect_assembly_manager->GetMortarManager()->GetGlobalLambdaCol();
      visualization_params.set<Teuchos::RCP<Epetra_Vector>>("lambda", lambda);

      // The pairs will need the mortar manager to extract their Lambda DOFs.
      visualization_params.set<Teuchos::RCP<const BEAMINTERACTION::BeamToSolidMortarManager>>(
          "mortar_manager", indirect_assembly_manager->GetMortarManager());

      // This map is used to ensure, that each discrete Lagrange multiplier is only written once per
      // beam element.
      Teuchos::RCP<std::unordered_set<int>> beam_tracker =
          Teuchos::rcp(new std::unordered_set<int>());
      visualization_params.set<Teuchos::RCP<std::unordered_set<int>>>("beam_tracker", beam_tracker);

      // Add the pair specific output.
      for (const auto& pair : indirect_assembly_manager->GetMortarManager()->GetContactPairs())
        pair->get_pair_visualization(output_writer_base_ptr_, visualization_params);

      if (nodal_force_visualization != Teuchos::null)
      {
        // Get the global force and moment resultants from the nodal forces. The first column
        // represents the forces, the second one the moments.
        CORE::LINALG::Matrix<3, 2, double> beam_resultant(true);
        CORE::LINALG::Matrix<3, 2, double> solid_resultant(true);
        GetGlobalCouplingForceResultants(beam_contact->Discret(),
            *(beam_contact->beam_interaction_data_state().GetForceNp()),
            *(beam_contact->beam_interaction_data_state().GetDisNp()), beam_resultant,
            solid_resultant);

        // The beam coupling moments are calculated in each pair and replace the ones evaluated in
        // the previous function.
        auto global_coupling_moment_origin =
            visualization_params.get<Teuchos::RCP<CORE::LINALG::Matrix<3, 1, double>>>(
                "global_coupling_moment_origin");
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          beam_resultant(i_dim, 1) = (*global_coupling_moment_origin)(i_dim);

        // Sum the values over all ranks.
        CORE::LINALG::Matrix<3, 2, double> beam_resultant_global(true);
        CORE::LINALG::Matrix<3, 2, double> solid_resultant_global(true);
        MPI_Allreduce(beam_resultant.A(), beam_resultant_global.A(),
            beam_resultant.numRows() * beam_resultant.numCols(), MPI_DOUBLE, MPI_SUM,
            dynamic_cast<const Epetra_MpiComm*>(&(beam_contact->Discret().Comm()))->Comm());
        MPI_Allreduce(solid_resultant.A(), solid_resultant_global.A(),
            solid_resultant.numRows() * solid_resultant.numCols(), MPI_DOUBLE, MPI_SUM,
            dynamic_cast<const Epetra_MpiComm*>(&(beam_contact->Discret().Comm()))->Comm());

        // Add to the visualization output writer.
        auto& visualization_data = nodal_force_visualization->get_visualization_data();
        std::vector<double>& field_data_beam_force =
            visualization_data.GetFieldData<double>("sum_coupling_force_beam");
        std::vector<double>& field_data_beam_moment =
            visualization_data.GetFieldData<double>("sum_coupling_moment_beam");
        field_data_beam_force.clear();
        field_data_beam_force.resize(3);
        field_data_beam_moment.clear();
        field_data_beam_moment.resize(3);
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        {
          field_data_beam_force[i_dim] = beam_resultant_global(i_dim, 0);
          field_data_beam_moment[i_dim] = beam_resultant_global(i_dim, 1);
        }
        std::vector<double>& field_data_solid_force =
            visualization_data.GetFieldData<double>("sum_coupling_force_solid");
        std::vector<double>& field_data_solid_moment =
            visualization_data.GetFieldData<double>("sum_coupling_moment_solid");
        field_data_solid_force.clear();
        field_data_solid_force.resize(3);
        field_data_solid_moment.clear();
        field_data_solid_moment.resize(3);
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        {
          field_data_solid_force[i_dim] = solid_resultant_global(i_dim, 0);
          field_data_solid_moment[i_dim] = solid_resultant_global(i_dim, 1);
        }

        visualization_params.remove("global_coupling_moment_origin", false);
      }

      // Reset assembly manager specific values in the parameter list passed to the individual
      // pairs.
      visualization_params.remove("lambda");
      visualization_params.remove("mortar_manager");
      visualization_params.remove("beam_tracker");
    }
  }


  // Write the data to disc. The data will be cleared in this method.
  output_writer_base_ptr_->Write(i_step, time);
}


/**
 * \brief Checks the init and setup status.
 */
void BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputWriter::check_init_setup() const
{
  if (!isinit_ or !issetup_) FOUR_C_THROW("Call Init() and Setup() first!");
}

/**
 * \brief Checks the init status.
 */
void BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputWriter::check_init() const
{
  if (!isinit_) FOUR_C_THROW("Init() has not been called, yet!");
}

FOUR_C_NAMESPACE_CLOSE
