/*----------------------------------------------------------------------*/
/*! \file

\brief Object to handle beam to solid volume meshtying output creation.

\level 3

\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_volume_meshtying_vtk_output_writer.H"

#include "beam_contact_pair.H"
#include "beam_to_solid_volume_meshtying_vtk_output_params.H"
#include "beam_to_solid_vtu_output_writer_base.H"
#include "beam_to_solid_vtu_output_writer_visualization.H"
#include "beaminteraction_submodel_evaluator_beamcontact.H"
#include "beam_to_solid_mortar_manager.H"
#include "beaminteraction_calc_utils.H"
#include "str_model_evaluator_beaminteraction_datastate.H"
#include "beaminteraction_submodel_evaluator_beamcontact_assembly_manager_indirect.H"
#include "../drt_structure_new/str_timint_basedataglobalstate.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"

#include <Epetra_FEVector.h>


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
        output_params_ptr,
    double restart_time)
{
  CheckInit();

  // Set beam to solid volume mesh tying output parameters.
  output_params_ptr_ = output_params_ptr;

  // Initialize the writer base object and add the desired visualizations.
  output_writer_base_ptr_ = Teuchos::rcp<BEAMINTERACTION::BeamToSolidVtuOutputWriterBase>(
      new BEAMINTERACTION::BeamToSolidVtuOutputWriterBase(
          "beam-to-solid-volume", vtk_params, restart_time));

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

    if (output_params_ptr_->GetMortarLambdaContinuousOutputFlag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->AddVisualizationWriter("mortar-continuous");
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
  // values are padded. The runtime output is written when the time step is already set to the next
  // step.
  int i_step = beam_contact->GState().GetStepN();
  double time = beam_contact->GState().GetTimeN();
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
  // Parameter list that will be passed to all contact pairs when they create their visualization.
  Teuchos::ParameterList visualization_params;
  visualization_params.set<Teuchos::RCP<const BeamToSolidVolumeMeshtyingVtkOutputParams>>(
      "output_params_ptr", output_params_ptr_);


  // Add the nodal forces resulting from beam contact. The forces are split up into beam and solid
  // nodes.
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization =
      output_writer_base_ptr_->GetVisualizationWriter("nodal-forces");
  if (visualization != Teuchos::null)
  {
    // Add the reference geometry and displacement to the visualization.
    visualization->AddDiscretizationNodalReferencePosition(beam_contact->DiscretPtr());
    visualization->AddDiscretizationNodalData(
        "displacement", beam_contact->BeamInteractionDataState().GetDisNp());

    // Create maps with the GIDs of beam and solid nodes.
    std::vector<int> gid_beam_dof;
    std::vector<int> gid_solid_dof;
    std::vector<int> gid_node;
    for (int i_lid = 0; i_lid < beam_contact->DiscretPtr()->NumMyRowNodes(); i_lid++)
    {
      gid_node.clear();
      DRT::Node* current_node = beam_contact->DiscretPtr()->lRowNode(i_lid);
      beam_contact->DiscretPtr()->Dof(current_node, gid_node);
      if (BEAMINTERACTION::UTILS::IsBeamNode(*current_node))
        for (unsigned int dim = 0; dim < 3; ++dim) gid_beam_dof.push_back(gid_node[dim]);
      else
        for (unsigned int dim = 0; dim < 3; ++dim) gid_solid_dof.push_back(gid_node[dim]);
    }
    Epetra_Map beam_dof_map(
        -1, gid_beam_dof.size(), &gid_beam_dof[0], 0, beam_contact->DiscretPtr()->Comm());
    Epetra_Map solid_dof_map(
        -1, gid_solid_dof.size(), &gid_solid_dof[0], 0, beam_contact->DiscretPtr()->Comm());

    // Extract the forces and add them to the discretization.
    Teuchos::RCP<Epetra_Vector> force_beam =
        Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(beam_dof_map, true));
    LINALG::Export(*beam_contact->BeamInteractionDataState().GetForceNp(), *force_beam);
    Teuchos::RCP<Epetra_Vector> force_solid =
        Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(solid_dof_map, true));
    LINALG::Export(*beam_contact->BeamInteractionDataState().GetForceNp(), *force_solid);
    visualization->AddDiscretizationNodalData("force_beam", force_beam);
    visualization->AddDiscretizationNodalData("force_solid", force_solid);
  }


  // Add the discrete Lagrange multiplicator values at the nodes of the Lagrange multiplicator
  // shape function. To do this we need to calculate the global lambda vector. It will be added to
  // the parameter list and each pair can get the values it needs and generate the visualization.
  visualization = output_writer_base_ptr_->GetVisualizationWriter("mortar");
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_continuous =
      output_writer_base_ptr_->GetVisualizationWriter("mortar-continuous");
  if (visualization != Teuchos::null || visualization_continuous != Teuchos::null)
  {
    // This output only works if there is an indirect assembly manager in the beam contact submodel
    // evaluator.
    Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerInDirect>
        indirect_assembly_manager = Teuchos::null;
    for (auto& assembly_manager : beam_contact->GetAssemblyManagers())
    {
      indirect_assembly_manager = Teuchos::rcp_dynamic_cast<
          BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManagerInDirect>(assembly_manager);
      if (indirect_assembly_manager != Teuchos::null) break;
    }

    if (indirect_assembly_manager != Teuchos::null)
    {
      // Get the global vector with the Lagrange Multiplier values and add it to the parameter list
      // that will be passed to the pairs.
      Teuchos::RCP<Epetra_Vector> lambda =
          indirect_assembly_manager->GetMortarManager()->GetGlobalLambdaCol(
              beam_contact->GState().GetDisNp());
      visualization_params.set<Teuchos::RCP<Epetra_Vector>>("lambda", lambda);

      // The pairs will need the mortar manager to extract their Lambda DOFs.
      visualization_params.set<Teuchos::RCP<const BEAMINTERACTION::BeamToSolidMortarManager>>(
          "mortar_manager", indirect_assembly_manager->GetMortarManager());
    }
  }


  // Add the pair specific visualization by looping over the individual contact pairs.
  for (const auto& pair : beam_contact->GetContactPairs())
    pair->GetPairVisualization(output_writer_base_ptr_, visualization_params);


  // Write the data to disc. The data will be cleared in this method.
  output_writer_base_ptr_->Write(i_step, time);
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
