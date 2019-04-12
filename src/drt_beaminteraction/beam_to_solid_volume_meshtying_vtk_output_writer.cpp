/*!
\file beam_to_solid_volume_meshtying_vtk_output_writer.cpp

\brief Object to handle beam to solid volume meshtying output creation.

\level 3

\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_volume_meshtying_vtk_output_writer.H"

#include "beam_to_solid_volume_meshtying_vtk_output_params.H"
#include "beam_to_solid_vtu_output_writer_base.H"
#include "beam_to_solid_vtu_output_writer_visualization.H"
#include "beaminteraction_submodel_evaluator_beamcontact.H"

#include "../drt_structure_new/str_timint_basedataglobalstate.H"


/**
 *
 */
BEAMINTERACTION::BeamToSolidVolumeMeshtyingVtkOutputWriter::
    BeamToSolidVolumeMeshtyingVtkOutputWriter()
    : isinit_(false),
      issetup_(false),
      output_params_ptr_(Teuchos::null),
      output_writer_base_ptr_(Teuchos::null)
{
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingVtkOutputWriter::Init()
{
  issetup_ = false;
  isinit_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingVtkOutputWriter::Setup(
    Teuchos::RCP<const STR::TIMINT::ParamsRuntimeVtkOutput> vtk_params,
    Teuchos::RCP<const BEAMINTERACTION::BeamToSolidVolumeMeshtyingVtkOutputParams>
        output_params_ptr)
{
  CheckInit();

  // Set beam to solid volume mesh tying output parameters.
  output_params_ptr_ = output_params_ptr;

  // Initialize the writer base object and add the desired visualizations.
  output_writer_base_ptr_ = Teuchos::rcp<BEAMINTERACTION::BeamToSolidVtuOutputWriterBase>(
      new BEAMINTERACTION::BeamToSolidVtuOutputWriterBase("beam-to-solid-volume", vtk_params));

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

    if (output_params_ptr_->GetMortarLambdaDiscretOutputFlag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->AddVisualizationWriter("mortar");
      visualization_writer->AddPointDataVector("displacement", 3);
      visualization_writer->AddPointDataVector("lambda", 3);
    }

    if (output_params_ptr_->GetIntegrationPointsOutputFlag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->AddVisualizationWriter("integration-points");
      visualization_writer->AddPointDataVector("displacement", 3);
      visualization_writer->AddPointDataVector("force", 3);
    }

    if (output_params_ptr_->GetSegmentationOutputFlag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->AddVisualizationWriter("segmentation");
      visualization_writer->AddPointDataVector("displacement", 3);
    }
  }

  issetup_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingVtkOutputWriter::WriteOutputRuntime(
    const BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact* beam_contact) const
{
  CheckInitSetup();

  // Get the time step and time for the output file. If output is desired at every iteration, the
  // values are padded.
  int i_step = beam_contact->GState().GetStepNp();
  double time = beam_contact->GState().GetTimeNp();
  if (output_params_ptr_->GetOutputEveryIteration()) i_step *= 10000;

  WriteOutputBeamToSolidVolumeMeshTying(beam_contact, i_step, time);
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingVtkOutputWriter::WriteOutputRuntimeIteration(
    const BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact* beam_contact, int i_iteration) const
{
  CheckInitSetup();

  if (output_params_ptr_->GetOutputEveryIteration())
  {
    // Get the time step and time for the output file. If output is desired at every iteration, the
    // values are padded.
    int i_step = 10000 * beam_contact->GState().GetStepN() + i_iteration;
    double time = beam_contact->GState().GetTimeN() + 1e-8 * i_iteration;

    WriteOutputBeamToSolidVolumeMeshTying(beam_contact, i_step, time);
  }
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingVtkOutputWriter::
    WriteOutputBeamToSolidVolumeMeshTying(
        const BEAMINTERACTION::SUBMODELEVALUATOR::BeamContact* beam_contact, int i_step,
        double time) const
{
}


/**
 * \brief Checks the init and setup status.
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingVtkOutputWriter::CheckInitSetup() const
{
  if (!isinit_ or !issetup_) dserror("Call Init() and Setup() first!");
}

/**
 * \brief Checks the init status.
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingVtkOutputWriter::CheckInit() const
{
  if (!isinit_) dserror("Init() has not been called, yet!");
}
