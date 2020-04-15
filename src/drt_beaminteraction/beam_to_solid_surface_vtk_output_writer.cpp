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
#include "beaminteraction_submodel_evaluator_beamcontact.H"
#include "beam_to_solid_conditions.H"
#include "../drt_geometry_pair/geometry_pair_line_to_surface_evaluation_data.H"
#include "str_model_evaluator_beaminteraction_datastate.H"
#include "../drt_structure_new/str_timint_basedataglobalstate.H"

#include <Epetra_FEVector.h>


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
          output_writer_base_ptr_->AddVisualizationWriter("nodal-forces");
      visualization_writer->AddPointDataVector("displacement", 3);
      visualization_writer->AddPointDataVector("force_beam", 3);
      visualization_writer->AddPointDataVector("force_solid", 3);
    }

    if (output_params_ptr_->GetAveragedNormalsOutputFlag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->AddVisualizationWriter("averaged-normals");
      visualization_writer->AddPointDataVector("displacement", 3);
      visualization_writer->AddPointDataVector("normal_averaged", 3);
      visualization_writer->AddPointDataVector("normal_element", 3);
      visualization_writer->AddPointDataVector("coupling_id", 1);
    }

    if (output_params_ptr_->GetMortarLambdaDiscretOutputFlag())
    {
      dserror("Mortar output not yet implemented.");
    }

    if (output_params_ptr_->GetMortarLambdaContinuousOutputFlag())
    {
      dserror("Mortar output not yet implemented.");
    }

    if (output_params_ptr_->GetIntegrationPointsOutputFlag())
    {
      dserror("Integration point output not yet implemented.");
    }

    if (output_params_ptr_->GetSegmentationOutputFlag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->AddVisualizationWriter("segmentation");
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
      "output_params_ptr", output_params_ptr_);


  // Add the averaged nodal normal output.
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization>
      visualization_averaged_normals =
          output_writer_base_ptr_->GetVisualizationWriter("averaged-normals");
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
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization =
      output_writer_base_ptr_->GetVisualizationWriter("nodal-forces");
  if (visualization != Teuchos::null)
    AddBeamInteractionNodalForces(visualization, beam_contact->DiscretPtr(),
        beam_contact->BeamInteractionDataState().GetDisNp(),
        beam_contact->BeamInteractionDataState().GetForceNp());


  // Add the pair specific visualization by looping over the individual contact pairs.
  for (const auto& pair : beam_contact->GetContactPairs())
    pair->GetPairVisualization(output_writer_base_ptr_, visualization_params);


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
