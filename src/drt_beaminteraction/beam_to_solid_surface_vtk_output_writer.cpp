/*----------------------------------------------------------------------*/
/*! \file

\brief Object to handle beam to solid surface output creation.

\level 3

\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_surface_vtk_output_writer.H"

#include "beam_contact_pair.H"
#include "beam_to_solid_surface_vtk_output_params.H"
#include "beam_to_solid_vtu_output_writer_base.H"
#include "beam_to_solid_vtu_output_writer_visualization.H"
#include "beam_to_solid_vtu_output_writer_utils.H"
#include "beam_to_solid_mortar_manager.H"
#include "beaminteraction_submodel_evaluator_beamcontact.H"
#include "beaminteraction_submodel_evaluator_beamcontact_assembly_manager_direct.H"
#include "beaminteraction_submodel_evaluator_beamcontact_assembly_manager_indirect.H"
#include "beam_to_solid_conditions.H"
#include "../drt_geometry_pair/geometry_pair_line_to_surface_evaluation_data.H"
#include "str_model_evaluator_beaminteraction_datastate.H"
#include "../drt_structure_new/str_timint_basedataglobalstate.H"

#include <Epetra_FEVector.h>
#include <unordered_set>


/**
 *
 */
BEAMINTERACTION::BeamToSolidSurfaceVtkOutputWriter::BeamToSolidSurfaceVtkOutputWriter()
    : isinit_(false),
      issetup_(false),
      output_params_ptr_(Teuchos::null),
      output_writer_base_ptr_(Teuchos::null)
{
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceVtkOutputWriter::Init()
{
  issetup_ = false;
  isinit_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceVtkOutputWriter::Setup(
    Teuchos::RCP<const STR::TIMINT::ParamsRuntimeVtkOutput> vtk_params,
    Teuchos::RCP<const BEAMINTERACTION::BeamToSolidSurfaceVtkOutputParams> output_params_ptr,
    double restart_time)
{
  CheckInit();

  // Set beam to solid surface interactions output parameters.
  output_params_ptr_ = output_params_ptr;

  // Initialize the writer base object and add the desired visualizations.
  output_writer_base_ptr_ = Teuchos::rcp<BEAMINTERACTION::BeamToSolidVtuOutputWriterBase>(
      new BEAMINTERACTION::BeamToSolidVtuOutputWriterBase(
          "beam-to-solid-surface", vtk_params, restart_time));

  // Depending on the selected input parameters, create the needed writers. All node / cell data
  // fields that should be output eventually have to be defined here. This helps to prevent issues
  // with ranks that do not contribute to a certain writer.
  {
    if (output_params_ptr_->GetNodalForceOutputFlag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->AddVisualizationWriter("nodal-forces", "btssc-nodal-forces");
      visualization_writer->AddPointDataVector("displacement", 3);
      visualization_writer->AddPointDataVector("force_beam", 3);
      visualization_writer->AddPointDataVector("force_solid", 3);
      visualization_writer->AddFieldDataVector("sum_coupling_force_beam");
      visualization_writer->AddFieldDataVector("sum_coupling_moment_beam");
      visualization_writer->AddFieldDataVector("sum_coupling_force_solid");
      visualization_writer->AddFieldDataVector("sum_coupling_moment_solid");
    }

    if (output_params_ptr_->GetAveragedNormalsOutputFlag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->AddVisualizationWriter(
              "averaged-normals", "btssc-averaged-normals");
      visualization_writer->AddPointDataVector("displacement", 3);
      visualization_writer->AddPointDataVector("normal_averaged", 3);
      visualization_writer->AddPointDataVector("normal_element", 3);
      visualization_writer->AddPointDataVector("coupling_id", 1);
    }

    if (output_params_ptr_->GetMortarLambdaDiscretOutputFlag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->AddVisualizationWriter("mortar", "btssc-mortar");
      visualization_writer->AddPointDataVector("displacement", 3);
      visualization_writer->AddPointDataVector("lambda", 3);
    }

    if (output_params_ptr_->GetMortarLambdaContinuousOutputFlag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->AddVisualizationWriter(
              "mortar-continuous", "btssc-mortar-continuous");
      visualization_writer->AddPointDataVector("displacement", 3);
      visualization_writer->AddPointDataVector("lambda", 3);
    }

    if (output_params_ptr_->GetIntegrationPointsOutputFlag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->AddVisualizationWriter(
              "integration-points", "btssc-integration-points");
      visualization_writer->AddPointDataVector("displacement", 3);
      visualization_writer->AddPointDataVector("projection_direction", 3);
    }

    if (output_params_ptr_->GetSegmentationOutputFlag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->AddVisualizationWriter("segmentation", "btssc-segmentation");
      visualization_writer->AddPointDataVector("displacement", 3);
      visualization_writer->AddPointDataVector("projection_direction", 3);
    }
  }

  issetup_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceVtkOutputWriter::WriteOutputRuntime(
    const BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact* beam_contact) const
{
  CheckInitSetup();

  // Get the time step and time for the output file. If output is desired at every iteration, the
  // values are padded. The runtime output is written when the time step is already set to the next
  // step.
  int i_step = beam_contact->GState().GetStepN();
  double time = beam_contact->GState().GetTimeN();
  if (output_params_ptr_->GetOutputEveryIteration()) i_step *= 10000;

  WriteOutputBeamToSolidSurface(beam_contact, i_step, time);
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceVtkOutputWriter::WriteOutputRuntimeIteration(
    const BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact* beam_contact, int i_iteration) const
{
  CheckInitSetup();

  if (output_params_ptr_->GetOutputEveryIteration())
  {
    // Get the time step and time for the output file. If output is desired at every iteration, the
    // values are padded.
    int i_step = 10000 * beam_contact->GState().GetStepN() + i_iteration;
    double time = beam_contact->GState().GetTimeN() + 1e-8 * i_iteration;

    WriteOutputBeamToSolidSurface(beam_contact, i_step, time);
  }
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceVtkOutputWriter::WriteOutputBeamToSolidSurface(
    const BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact* beam_contact, int i_step,
    double time) const
{
  // Parameter list that will be passed to all contact pairs when they create their visualization.
  Teuchos::ParameterList visualization_params;
  visualization_params.set<Teuchos::RCP<const BeamToSolidSurfaceVtkOutputParams>>(
      "btssc-output_params_ptr", output_params_ptr_);


  // Add the averaged nodal normal output.
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization>
      visualization_averaged_normals =
          output_writer_base_ptr_->GetVisualizationWriter("btssc-averaged-normals");
  if (visualization_averaged_normals != Teuchos::null)
  {
    const std::vector<Teuchos::RCP<BeamInteractionConditionBase>>& surface_condition_vector =
        beam_contact->GetConditions()->GetConditionMap().at(
            INPAR::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_surface_meshtying);
    for (const auto& condition : surface_condition_vector)
    {
      // Get the line-to-surface evaluation data for the current condition.
      Teuchos::RCP<const GEOMETRYPAIR::LineToSurfaceEvaluationData> surface_evaluation_data =
          Teuchos::rcp_dynamic_cast<const GEOMETRYPAIR::LineToSurfaceEvaluationData>(
              condition->GetGeometryEvaluationData(), true);

      // Get the coupling ID for the current condition.
      auto beam_to_surface_condition =
          Teuchos::rcp_dynamic_cast<const BeamToSolidConditionSurfaceMeshtying>(condition, true);
      const int coupling_id = beam_to_surface_condition->GetOtherCondition()->GetInt("COUPLING_ID");

      // Create the output for the averaged normal field.
      AddAveragedNodalNormals(
          visualization_averaged_normals, surface_evaluation_data->GetFaceElements(), coupling_id);
    }
  }


  // Add the nodal forces resulting from beam contact. The forces are split up into beam and
  // solid nodes.
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> nodal_force_visualization =
      output_writer_base_ptr_->GetVisualizationWriter("btssc-nodal-forces");
  if (nodal_force_visualization != Teuchos::null)
    AddBeamInteractionNodalForces(nodal_force_visualization, beam_contact->DiscretPtr(),
        beam_contact->BeamInteractionDataState().GetDisNp(),
        beam_contact->BeamInteractionDataState().GetForceNp());


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
        pair->GetPairVisualization(output_writer_base_ptr_, visualization_params);
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
        auto global_coupling_moment_origin = Teuchos::rcp(new LINALG::Matrix<3, 1, double>(true));
        visualization_params.set("global_coupling_moment_origin", global_coupling_moment_origin);
      }

      // Get the global vector with the Lagrange Multiplier values and add it to the parameter list
      // that will be passed to the pairs.
      Teuchos::RCP<Epetra_Vector> lambda =
          indirect_assembly_manager->GetMortarManager()->GetGlobalLambdaCol(
              beam_contact->GState().GetDisNp());
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
        pair->GetPairVisualization(output_writer_base_ptr_, visualization_params);

      if (nodal_force_visualization != Teuchos::null)
      {
        // Get the global force and moment resultants from the nodal forces. The first column
        // represents the forces, the second one the moments.
        LINALG::Matrix<3, 2, double> beam_resultant(true);
        LINALG::Matrix<3, 2, double> solid_resultant(true);
        GetGlobalCouplingForceResultants(beam_contact->Discret(),
            *(beam_contact->BeamInteractionDataState().GetForceNp()),
            *(beam_contact->BeamInteractionDataState().GetDisNp()), beam_resultant,
            solid_resultant);

        // The beam coupling moments are calculated in each pair and replace the ones evaluated in
        // the previous function.
        auto global_coupling_moment_origin =
            visualization_params.get<Teuchos::RCP<LINALG::Matrix<3, 1, double>>>(
                "global_coupling_moment_origin");
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          beam_resultant(i_dim, 1) = (*global_coupling_moment_origin)(i_dim);

        // Sum the values over all ranks.
        LINALG::Matrix<3, 2, double> beam_resultant_global(true);
        LINALG::Matrix<3, 2, double> solid_resultant_global(true);
        MPI_Allreduce(beam_resultant.A(), beam_resultant_global.A(),
            beam_resultant.M() * beam_resultant.N(), MPI_DOUBLE, MPI_SUM,
            dynamic_cast<const Epetra_MpiComm*>(&(beam_contact->Discret().Comm()))->Comm());
        MPI_Allreduce(solid_resultant.A(), solid_resultant_global.A(),
            solid_resultant.M() * solid_resultant.N(), MPI_DOUBLE, MPI_SUM,
            dynamic_cast<const Epetra_MpiComm*>(&(beam_contact->Discret().Comm()))->Comm());

        // Add to the vtk output writer.
        std::vector<double>& field_data_beam_force =
            nodal_force_visualization->GetMutableFieldDataVector("sum_coupling_force_beam");
        std::vector<double>& field_data_beam_moment =
            nodal_force_visualization->GetMutableFieldDataVector("sum_coupling_moment_beam");
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
            nodal_force_visualization->GetMutableFieldDataVector("sum_coupling_force_solid");
        std::vector<double>& field_data_solid_moment =
            nodal_force_visualization->GetMutableFieldDataVector("sum_coupling_moment_solid");
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
void BEAMINTERACTION::BeamToSolidSurfaceVtkOutputWriter::CheckInitSetup() const
{
  if (!isinit_ or !issetup_) dserror("Call Init() and Setup() first!");
}

/**
 * \brief Checks the init status.
 */
void BEAMINTERACTION::BeamToSolidSurfaceVtkOutputWriter::CheckInit() const
{
  if (!isinit_) dserror("Init() has not been called, yet!");
}
