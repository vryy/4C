/*----------------------------------------------------------------------*/
/*! \file

\brief Object that stores the relevant data for a single output file.

\level 3

\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_vtu_output_writer_visualization.H"

#include "../drt_structure_new/str_timint_basedataio_runtime_vtk_output.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

#include <mpi.h>


/**
 *
 */
BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization::BeamToSolidVtuOutputWriterVisualization(
    const std::string& writer_full_name,
    Teuchos::RCP<const STR::TIMINT::ParamsRuntimeVtkOutput> vtk_params, double restart_time)
    : RuntimeVtuWriter(),
      vtk_params_(vtk_params),
      writer_full_name_(writer_full_name),
      discret_(Teuchos::null),
      node_gid_map_(Teuchos::null)
{
  int my_rank;
  int total_ranks;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &total_ranks);

  // Todo: we need a better upper bound for total number of time steps here
  // however, this 'only' affects the number of leading zeros in the vtk file names
  unsigned int num_timesteps_in_simulation_upper_bound = 1000000;

  if (vtk_params_->OutputEveryIteration()) num_timesteps_in_simulation_upper_bound *= 10000;

  // Determine path of output directory.
  const std::string outputfilename(DRT::Problem::Instance()->OutputControlFile()->FileName());
  size_t pos = outputfilename.find_last_of("/");
  if (pos == outputfilename.npos)
    pos = 0ul;
  else
    pos++;
  const std::string output_directory_path(outputfilename.substr(0ul, pos));

  Initialize(my_rank, total_ranks, num_timesteps_in_simulation_upper_bound, output_directory_path,
      DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix(), writer_full_name_,
      DRT::Problem::Instance()->OutputControlFile()->RestartName(), restart_time,
      vtk_params_->WriteBinaryOutput());
}

/**
 *
 */
std::vector<double>&
BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization::GetMutablePointCoordinateVector(
    int additional_reserve)
{
  std::vector<double>& return_value = RuntimeVtuWriter::GetMutablePointCoordinateVector();
  if (additional_reserve > 0) return_value.reserve(return_value.size() + additional_reserve);
  return return_value;
}

/**
 *
 */
std::vector<double>&
BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization::GetMutablePointDataVector(
    std::string data_name, int additional_reserve)
{
  // Check if the point data vector already exists.
  std::map<std::string, std::pair<std::vector<double>, unsigned int>>& point_data_map =
      GetMutablePointDataMap();
  const auto& it = point_data_map.find(data_name);
  if (it != point_data_map.end())
  {
    if (additional_reserve > 0)
      it->second.first.reserve(it->second.first.size() + additional_reserve);
    return it->second.first;
  }
  else
  {
    dserror("The data vector with the name '%s' does not exist.", data_name.c_str());
    return it->second.first;
  }
}

/**
 *
 */
std::vector<double>&
BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization::GetMutableCellDataVector(
    std::string data_name, int additional_reserve)
{
  // Check if the cell data vector already exists.
  std::map<std::string, std::pair<std::vector<double>, unsigned int>>& cell_data_map =
      GetMutableCellDataMap();
  const auto& it = cell_data_map.find(data_name);
  if (it != cell_data_map.end())
  {
    if (additional_reserve > 0)
      it->second.first.reserve(it->second.first.size() + additional_reserve);
    return it->second.first;
  }
  else
  {
    dserror("The cell data vector with the name '%s' does not exist.", data_name.c_str());
    return it->second.first;
  }
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization::AddPointDataVector(
    std::string data_name, unsigned int n_dim, int reserve)
{
  // Check if the point data vector already exists.
  std::map<std::string, std::pair<std::vector<double>, unsigned int>>& point_data_map =
      GetMutablePointDataMap();
  const auto& it = point_data_map.find(data_name);
  if (it != point_data_map.end())
    dserror("The point data vector with the name '%s' you want to add already exists.",
        data_name.c_str());
  else
  {
    std::pair<std::vector<double>, unsigned int> new_pair;
    point_data_map[data_name] = new_pair;
    std::pair<std::vector<double>, unsigned int>& added_pair = point_data_map[data_name];
    added_pair.second = n_dim;
    if (reserve > 0) added_pair.first.reserve(reserve);
  }
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization::AddCellDataVector(
    std::string data_name, unsigned int n_dim, int reserve)
{
  // Check if the cell data vector already exists.
  std::map<std::string, std::pair<std::vector<double>, unsigned int>>& cell_data_map =
      GetMutableCellDataMap();
  const auto& it = cell_data_map.find(data_name);
  if (it != cell_data_map.end())
    dserror("The cell data vector with the name '%s' you want to add already exists.",
        data_name.c_str());
  else
  {
    std::pair<std::vector<double>, unsigned int> new_pair;
    cell_data_map[data_name] = new_pair;
    std::pair<std::vector<double>, unsigned int>& added_pair = cell_data_map[data_name];
    added_pair.second = n_dim;
    if (reserve > 0) added_pair.first.reserve(reserve);
  }
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization::
    AddDiscretizationNodalReferencePosition(const Teuchos::RCP<const DRT::Discretization>& discret)
{
  // Check that the discretization is not already set, and that all data in the writer is empty.
  if (discret_ != Teuchos::null)
    dserror(
        "When calling AddDiscretizationNodalReferencePosition, the discretization can not be "
        "already set. Did you forget to reset the writer?");
  if (GetMutablePointCoordinateVector().size() != 0)
    dserror("Point coordinate vector is not empty!");
  for (const auto& point_data : GetMutablePointDataMap())
    if (point_data.second.first.size() != 0)
      dserror("Point data for '%s' is not empty!", point_data.first.c_str());
  if (GetMutableCellTypeVector().size() != 0) dserror("Cell type vector is not empty!");
  if (GetMutableCellOffsetVector().size() != 0) dserror("Cell offset vector is not empty!");
  for (const auto& cell_data : GetMutableCellDataMap())
    if (cell_data.second.first.size() != 0)
      dserror("Cell data for '%s' is not empty!", cell_data.first.c_str());

  // Set the discretization for this writer.
  discret_ = discret;

  // Setup variables for the position and map.
  unsigned int num_my_nodes = discret_->NumMyRowNodes();
  std::vector<int> my_global_dof_ids;
  std::vector<int> node_global_dof_ids;
  std::vector<double>& point_coordinates = GetMutablePointCoordinateVector(3 * num_my_nodes);

  // Check that the position vector is empty.
  if (point_coordinates.size() != 0)
    dserror("The position vector has to be empty when adding nodal reference data!");

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
      new Epetra_Map(-1, my_global_dof_ids.size(), &my_global_dof_ids[0], 0, discret_->Comm()));
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization::AddDiscretizationNodalData(
    const std::string& data_name, const Teuchos::RCP<const Epetra_MultiVector>& vector)
{
  if (discret_ == Teuchos::null || node_gid_map_ == Teuchos::null)
    dserror("Discretization or node GID map is not set!");

  // Extract the vector according to the GIDs needed on this rank.
  Teuchos::RCP<Epetra_Vector> vector_extract =
      Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(*node_gid_map_, true));
  LINALG::Export(*vector, *vector_extract);

  // Add the values form the vector to the writer data.
  const int num_my_gid = node_gid_map_->NumMyElements();
  std::vector<double>& data_vector = GetMutablePointDataVector(data_name, 3 * num_my_gid);
  data_vector.reserve(3 * num_my_gid);
  for (int i_lid = 0; i_lid < num_my_gid; i_lid++) data_vector.push_back((*vector_extract)[i_lid]);
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization::Write(
    const unsigned int timestep_number, const double time)
{
  // Do some sanity checks on the data that will be written.
  const unsigned int n_points = GetMutablePointCoordinateVector().size() / 3;
  if (GetMutablePointCoordinateVector().size() % 3 != 0)
    dserror("The size of the coordinate vector (%d) is not a multiple of 3!",
        GetMutablePointCoordinateVector().size());
  for (const auto& point_data : GetMutablePointDataMap())
    if (point_data.second.first.size() / point_data.second.second != n_points ||
        point_data.second.first.size() % point_data.second.second != 0)
      dserror(
          "The size of the point data vector (%d) '%s' does not match the number of points (%d).",
          point_data.second.first.size(), point_data.first.c_str(), n_points);

  // Reset time and time step and geometry name in the writer object.
  SetupForNewTimeStepAndGeometry(time, timestep_number, writer_full_name_);
  // Finalize everything and write all required vtk files to filesystem.
  WriteFiles();
  // Write a collection file summarizing all previously written files.
  WriteCollectionFileOfAllWrittenFiles(
      DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix() + "-" +
      writer_full_name_);

  // Reset the data.
  discret_ = Teuchos::null;
  node_gid_map_ = Teuchos::null;
  GetMutablePointCoordinateVector().clear();
  GetMutableCellTypeVector().clear();
  GetMutableCellOffsetVector().clear();
  for (auto& point_data : GetMutablePointDataMap()) point_data.second.first.clear();
  for (auto& cell_data : GetMutableCellDataMap()) cell_data.second.first.clear();
}
