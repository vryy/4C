/*----------------------------------------------------------------------*/
/*! \file

\brief Object that stores the relevant data for a single output file.

\level 3

*/


#include "baci_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"

#include "baci_comm_utils.hpp"
#include "baci_global_data.hpp"
#include "baci_io_visualization_manager.hpp"
#include "baci_lib_discret.hpp"
#include "baci_linalg_utils_sparse_algebra_manipulation.hpp"
#include "baci_structure_new_timint_basedataio_runtime_vtk_output.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
BEAMINTERACTION::BeamToSolidOutputWriterVisualization::BeamToSolidOutputWriterVisualization(
    const std::string& writer_full_name, IO::VisualizationParameters visualization_params,
    Teuchos::RCP<const STR::TIMINT::ParamsRuntimeOutput> visualization_output_params)
    : IO::VisualizationManager(std::move(visualization_params),
          *(GLOBAL::Problem::Instance()->GetCommunicators()->GlobalComm()), writer_full_name),
      visualization_output_params_(visualization_output_params),
      writer_full_name_(writer_full_name),
      discret_(Teuchos::null),
      node_gid_map_(Teuchos::null)
{
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidOutputWriterVisualization::AddDiscretizationNodalReferencePosition(
    const Teuchos::RCP<const DRT::Discretization>& discret)
{
  auto& visualization_data = GetVisualizationData();

  // Check that the discretization is not already set, and that all data in the writer is empty.
  if (discret_ != Teuchos::null)
    FOUR_C_THROW(
        "When calling AddDiscretizationNodalReferencePosition, the discretization can not be "
        "already set. Did you forget to reset the writer?");
  if (visualization_data.GetPointCoordinates().size() != 0)
    FOUR_C_THROW("Point coordinate vector is not empty!");
  for (const auto& point_data_name : visualization_data.GetPointDataNames())
    if (visualization_data.GetPointDataSize(point_data_name) != 0)
      FOUR_C_THROW("Point data for '%s' is not empty!", point_data_name.c_str());
  if (visualization_data.GetCellTypes().size() != 0)
    FOUR_C_THROW("Cell types vector is not empty!");
  if (visualization_data.GetCellOffsets().size() != 0)
    FOUR_C_THROW("Cell offsets vector is not empty!");
  for (const auto& cell_data_name : visualization_data.GetCellDataNames())
    if (visualization_data.GetCellDataSize(cell_data_name) != 0)
      FOUR_C_THROW("Cell data for '%s' is not empty!", cell_data_name.c_str());

  // Set the discretization for this writer.
  discret_ = discret;

  // Setup variables for the position and map.
  unsigned int num_my_nodes = discret_->NumMyRowNodes();
  std::vector<int> my_global_dof_ids;
  std::vector<int> node_global_dof_ids;
  std::vector<double>& point_coordinates = visualization_data.GetPointCoordinates(3 * num_my_nodes);

  // Check that the position vector is empty.
  if (point_coordinates.size() != 0)
    FOUR_C_THROW("The position vector has to be empty when adding nodal reference data!");

  // Loop over the nodes on this rank.
  for (unsigned int i_node = 0; i_node < num_my_nodes; i_node++)
  {
    const DRT::Node* current_node = discret_->lRowNode(i_node);
    node_global_dof_ids.clear();
    discret_->Dof(current_node, node_global_dof_ids);
    for (unsigned int dim = 0; dim < 3; ++dim)
    {
      my_global_dof_ids.push_back(node_global_dof_ids[dim]);
      point_coordinates.push_back(current_node->X()[dim]);
    }
  }
  node_gid_map_ = Teuchos::rcp<Epetra_Map>(
      new Epetra_Map(-1, my_global_dof_ids.size(), my_global_dof_ids.data(), 0, discret_->Comm()));
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidOutputWriterVisualization::AddDiscretizationNodalData(
    const std::string& data_name, const Teuchos::RCP<const Epetra_MultiVector>& vector)
{
  if (discret_ == Teuchos::null || node_gid_map_ == Teuchos::null)
    FOUR_C_THROW("Discretization or node GID map is not set!");

  // Extract the vector according to the GIDs needed on this rank.
  Teuchos::RCP<Epetra_Vector> vector_extract =
      Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(*node_gid_map_, true));
  CORE::LINALG::Export(*vector, *vector_extract);

  // Add the values form the vector to the writer data.
  const int num_my_gid = node_gid_map_->NumMyElements();
  std::vector<double>& data_vector =
      GetVisualizationData().GetPointData<double>(data_name, 3 * num_my_gid);
  data_vector.reserve(3 * num_my_gid);
  for (int i_lid = 0; i_lid < num_my_gid; i_lid++) data_vector.push_back((*vector_extract)[i_lid]);
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidOutputWriterVisualization::Write(
    const unsigned int timestep_number, const double time)
{
  // Finalize everything and write all required vtk files to filesystem.
  WriteToDisk(time, timestep_number);

  // Reset the data.
  discret_ = Teuchos::null;
  node_gid_map_ = Teuchos::null;
  ClearData();
}

FOUR_C_NAMESPACE_CLOSE
