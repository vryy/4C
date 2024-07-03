/*----------------------------------------------------------------------*/
/*! \file

\brief Object to handle beam to solid volume meshtying output creation.

\level 2

*/


#include "4C_fbi_beam_to_fluid_meshtying_output_writer.hpp"

#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_base.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_fbi_adapter_constraintbridge_penalty.hpp"
#include "4C_fbi_beam_to_fluid_meshtying_output_params.hpp"
#include "4C_fbi_beam_to_fluid_meshtying_params.hpp"
#include "4C_fbi_constraintenforcer.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN


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
void BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter::init()
{
  issetup_ = false;
  isinit_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter::setup(
    const Core::IO::VisualizationParameters& visualization_params,
    Teuchos::RCP<const Solid::TimeInt::ParamsRuntimeOutput> visualization_output_params,
    Teuchos::RCP<const FBI::BeamToFluidMeshtyingVtkOutputParams> output_params_ptr)
{
  check_init();

  // Set beam to solid volume mesh tying output parameters.
  output_params_ptr_ = output_params_ptr;

  visualization_params_ = visualization_params;

  // Initialize the writer base object and add the desired visualizations.
  output_writer_base_ptr_ = Teuchos::rcp<BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase>(
      new BEAMINTERACTION::BeamToSolidVisualizationOutputWriterBase(
          "beam-to-fluid", visualization_output_params, visualization_params));

  // Depending on the selected input parameters, create the needed writers. All node / cell data
  // fields that should be output eventually have to be defined here. This helps to prevent issues
  // with ranks that do not contribute to a certain writer.
  {
    if (output_params_ptr_->get_nodal_force_output_flag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->add_visualization_writer("nodal-forces");
      auto& visualization_data = visualization_writer->get_visualization_data();
      visualization_data.RegisterPointData<double>("velocity", 3);
      visualization_data.RegisterPointData<double>("displacement", 3);
      visualization_data.RegisterPointData<double>("force", 3);
    }

    if (output_params_ptr_->get_mortar_lambda_discret_output_flag())
    {
      FOUR_C_THROW("Mortar discretization not implemented for beam-to-fluid meshtying!\n");
    }

    if (output_params_ptr_->get_mortar_lambda_continuous_output_flag())
    {
      FOUR_C_THROW("Mortar discretization not implemented for beam-to-fluid meshtying!\n");
    }

    if (output_params_ptr_->get_integration_points_output_flag())
    {
      Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_writer =
          output_writer_base_ptr_->add_visualization_writer("integration-points");
      auto& visualization_data = visualization_writer->get_visualization_data();
      visualization_data.RegisterPointData<double>("displacement", 3);
    }
  }

  issetup_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter::write_output_runtime(
    const Teuchos::RCP<Adapter::FBIConstraintenforcer>& couplingenforcer, int i_step,
    double time) const
{
  check_init_setup();

  auto [output_time, output_step] =
      Core::IO::GetTimeAndTimeStepIndexForOutput(visualization_params_, time, i_step);
  write_output_beam_to_fluid_mesh_tying(couplingenforcer, output_step, output_time);
}

/**
 *
 */
void BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter::write_output_beam_to_fluid_mesh_tying(
    const Teuchos::RCP<Adapter::FBIConstraintenforcer>& couplingenforcer, int i_step,
    double time) const
{
  // Parameter list that will be passed to all contact pairs when they create their visualization.
  Teuchos::ParameterList visualization_params;
  visualization_params.set<Teuchos::RCP<const BeamToSolidVolumeMeshtyingVisualizationOutputParams>>(
      "output_params_ptr", output_params_ptr_);


  // Add the nodal forces resulting from beam contact. The forces are split up into beam and solid
  // nodes.
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization =
      output_writer_base_ptr_->get_visualization_writer("nodal-forces");
  if (visualization != Teuchos::null)
  {
    // Add the reference geometry and displacement to the visualization.
    visualization->add_discretization_nodal_reference_position(
        couplingenforcer->GetStructure()->get_discretization());
    visualization->add_discretization_nodal_data(
        "velocity", couplingenforcer->GetStructure()->Velnp());
    visualization->add_discretization_nodal_data(
        "displacement", couplingenforcer->GetStructure()->Dispnp());

    // Create maps with the GIDs of beam nodes
    std::vector<int> gid_beam_dof;
    std::vector<int> gid_node;
    for (int i_lid = 0;
         i_lid < couplingenforcer->GetStructure()->get_discretization()->NumMyRowNodes(); i_lid++)
    {
      gid_node.clear();
      Core::Nodes::Node* current_node =
          couplingenforcer->GetStructure()->get_discretization()->lRowNode(i_lid);
      couplingenforcer->GetStructure()->get_discretization()->Dof(current_node, gid_node);
      if (BEAMINTERACTION::UTILS::IsBeamNode(*current_node))
        for (unsigned int dim = 0; dim < 3; ++dim) gid_beam_dof.push_back(gid_node[dim]);
    }
    Epetra_Map beam_dof_map(-1, gid_beam_dof.size(), gid_beam_dof.data(), 0,
        couplingenforcer->GetStructure()->get_discretization()->Comm());

    // Extract the forces and add them to the discretization.
    Teuchos::RCP<Epetra_Vector> force_beam =
        Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(beam_dof_map, true));
    Core::LinAlg::Export(*couplingenforcer->assemble_structure_coupling_residual(), *force_beam);


    visualization->add_discretization_nodal_data("force", force_beam);
  }

  // Add the pair specific visualization by looping over the individual contact pairs.
  for (const auto& pair : *couplingenforcer->GetBridge()->GetPairs())
    pair->get_pair_visualization(output_writer_base_ptr_, visualization_params);


  // Write the data to disc. The data will be cleared in this method.
  output_writer_base_ptr_->Write(i_step, time);
}


/**
 * \brief Checks the init and setup status.
 */
void BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter::check_init_setup() const
{
  if (!isinit_ or !issetup_) FOUR_C_THROW("Call init() and setup() first!");
}

/**
 * \brief Checks the init status.
 */
void BEAMINTERACTION::BeamToFluidMeshtyingVtkOutputWriter::check_init() const
{
  if (!isinit_) FOUR_C_THROW("init() has not been called, yet!");
}

FOUR_C_NAMESPACE_CLOSE
