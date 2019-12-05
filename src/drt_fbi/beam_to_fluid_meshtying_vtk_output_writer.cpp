/*----------------------------------------------------------------------*/
/*! \file

\brief Object to handle beam to solid volume meshtying output creation.

\level 2

\maintainer Nora Hagmeyer
*/


#include "beam_to_fluid_meshtying_vtk_output_writer.H"

#include "constraintenforcer_fbi.H"
#include "ad_fbi_constraintbridge_penalty.H"
#include "beam_to_fluid_meshtying_params.H"

#include "../drt_adapter/ad_str_fsiwrapper.H"

#include "../drt_beaminteraction/beam_contact_pair.H"
#include "../drt_fbi/beam_to_fluid_meshtying_vtk_output_params.H"
#include "../drt_beaminteraction/beam_to_solid_vtu_output_writer_base.H"
#include "../drt_beaminteraction/beam_to_solid_vtu_output_writer_visualization.H"
#include "../drt_beaminteraction/beaminteraction_calc_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

#include <Epetra_FEVector.h>


/**
 *
 */
BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter::BeamToFluidMeshtyingVtkOutputWriter()
    : isinit_(false),
      issetup_(false),
      output_params_ptr_(Teuchos::null),
      output_writer_base_ptr_(Teuchos::null)
{
}

/**
 *
 */
void BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter::Init()
{
  issetup_ = false;
  isinit_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter::Setup(
    Teuchos::RCP<const STR::TIMINT::ParamsRuntimeVtkOutput> vtk_params,
    Teuchos::RCP<const FBI::BeamToFluidMeshtyingVtkOutputParams> output_params_ptr,
    double restart_time)
{
  CheckInit();

  // Set beam to solid volume mesh tying output parameters.
  output_params_ptr_ = output_params_ptr;

  // Initialize the writer base object and add the desired visualizations.
  output_writer_base_ptr_ = Teuchos::rcp<BEAMINTERACTION::BeamToSolidVtuOutputWriterBase>(
      new BEAMINTERACTION::BeamToSolidVtuOutputWriterBase(
          "beam-to-fluid", vtk_params, restart_time));

  // Depending on the selected input parameters, create the needed writers. All node / cell data
  // fields that should be output eventually have to be defined here. This helps to prevent issues
  // with ranks that do not contribute to a certain writer.
  {
    if (output_params_ptr_->GetNodalForceOutputFlag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->AddVisualizationWriter("nodal-forces");
      visualization_writer->AddPointDataVector("velocity", 3);
      visualization_writer->AddPointDataVector("displacement", 3);
      visualization_writer->AddPointDataVector("force", 3);
    }

    if (output_params_ptr_->GetMortarLambdaDiscretOutputFlag())
    {
      dserror("Mortar discretization not implemented for beam-to-fluid meshtying!\n");
    }

    if (output_params_ptr_->GetMortarLambdaContinuousOutputFlag())
    {
      dserror("Mortar discretization not implemented for beam-to-fluid meshtying!\n");
    }

    if (output_params_ptr_->GetIntegrationPointsOutputFlag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->AddVisualizationWriter("integration-points");
      visualization_writer->AddPointDataVector("velocity", 3);
      visualization_writer->AddPointDataVector("displacement", 3);
      visualization_writer->AddPointDataVector("force", 3);
    }
  }

  issetup_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter::WriteOutputRuntime(
    const Teuchos::RCP<ADAPTER::FBIConstraintenforcer>& couplingenforcer, int i_step,
    double time) const
{
  CheckInitSetup();

  // If output is desired at every iteration, the
  // values are padded. The runtime output is written when the time step is already set to the next
  // step.
  printf("We are in the output writer\n");
  if (output_params_ptr_->GetOutputEveryIteration()) i_step *= 10000;

  WriteOutputBeamToFluidMeshTying(couplingenforcer, i_step, time);
}

/**
 *
 */
void BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter::WriteOutputBeamToFluidMeshTying(
    const Teuchos::RCP<ADAPTER::FBIConstraintenforcer>& couplingenforcer, int i_step,
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
    visualization->AddDiscretizationNodalReferencePosition(
        couplingenforcer->GetStructure()->GetDiscretization());
    visualization->AddDiscretizationNodalData(
        "velocity", couplingenforcer->GetStructure()->Velnp());
    visualization->AddDiscretizationNodalData(
        "displacement", couplingenforcer->GetStructure()->Dispnp());

    // Create maps with the GIDs of beam nodes
    std::vector<int> gid_beam_dof;
    std::vector<int> gid_node;
    for (int i_lid = 0;
         i_lid < couplingenforcer->GetStructure()->GetDiscretization()->NumMyRowNodes(); i_lid++)
    {
      gid_node.clear();
      DRT::Node* current_node =
          couplingenforcer->GetStructure()->GetDiscretization()->lRowNode(i_lid);
      couplingenforcer->GetStructure()->GetDiscretization()->Dof(current_node, gid_node);
      if (BEAMINTERACTION::UTILS::IsBeamNode(*current_node))
        for (unsigned int dim = 0; dim < 3; ++dim) gid_beam_dof.push_back(gid_node[dim]);
    }
    Epetra_Map beam_dof_map(-1, gid_beam_dof.size(), &gid_beam_dof[0], 0,
        couplingenforcer->GetStructure()->GetDiscretization()->Comm());

    // Extract the forces and add them to the discretization.
    Teuchos::RCP<Epetra_Vector> force_beam =
        Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(beam_dof_map, true));
    LINALG::Export(*Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(
                       couplingenforcer->GetBridge(), true)
                        ->GetFs(),
        *force_beam);

    if (force_beam->Scale((-1) * couplingenforcer->GetBridge()->GetParams()->GetPenaltyParameter()))
      dserror("Scaling of the penalty force was unsuccessful!\n");

    visualization->AddDiscretizationNodalData("force", force_beam);
  }

  // Add the pair specific visualization by looping over the individual contact pairs.
  for (const auto& pair : *couplingenforcer->GetBridge()->GetPairs())
    pair->GetPairVisualization(output_writer_base_ptr_, visualization_params);


  // Write the data to disc. The data will be cleared in this method.
  output_writer_base_ptr_->Write(i_step, time);
}


/**
 * \brief Checks the init and setup status.
 */
void BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter::CheckInitSetup() const
{
  if (!isinit_ or !issetup_) dserror("Call Init() and Setup() first!");
}

/**
 * \brief Checks the init status.
 */
void BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter::CheckInit() const
{
  if (!isinit_) dserror("Init() has not been called, yet!");
}
