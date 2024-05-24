/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Define a data container that stores visualization data in the VTU data format

\level 0

*/
/*-----------------------------------------------------------------------------------------------*/


#include "4C_io_visualization_data.hpp"

#include <algorithm>

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
size_t IO::VisualizationData::get_point_coordinates_number_of_points() const
{
  // Check if number of entries is consistent
  if (point_coordinates_.size() % n_dim_ != 0)
  {
    FOUR_C_THROW(
        "The size of the point coordinate vector (%d) is not a multiple of the spatial "
        "dimension (%d)!",
        point_coordinates_.size(), n_dim_);
  }
  return point_coordinates_.size() / n_dim_;
}

/**
 *
 */
void IO::VisualizationData::clear_data()
{
  point_coordinates_.clear();
  cell_types_.clear();
  cell_connectivity_.clear();
  cell_offsets_.clear();
  face_connectivity_.clear();
  face_offsets_.clear();
  for (const auto& data_name : GetPointDataNames()) clear_data(point_data_, data_name);
  for (const auto& data_name : GetCellDataNames()) clear_data(cell_data_, data_name);
  for (const auto& data_name : GetFieldDataNames()) clear_data(field_data_, data_name);
}

/**
 *
 */
void IO::VisualizationData::ResetContainer()
{
  clear_data();
  point_data_.clear();
  cell_data_.clear();
  field_data_.clear();
}

/**
 *
 */
void IO::VisualizationData::consistency_check_and_complete_data()
{
  complete_cell_connectivity();
  complete_face_connectivity();
  ConsistencyCheck();
}

/**
 *
 */
void IO::VisualizationData::ConsistencyCheck() const
{
  // Check basic point coordinates consistency (is done with get_point_coordinates_number_of_points)
  const auto n_points = get_point_coordinates_number_of_points();

  // Check basic cell consistency
  const auto n_cells = cell_types_.size();
  if (cell_offsets_.size() != n_cells)
    FOUR_C_THROW("The cell offset data array length (%d) does not match the number of cells (%d)",
        cell_offsets_.size(), n_cells);
  if (cell_connectivity_.size() < n_cells)
    FOUR_C_THROW("The number of cells is %d, however, there are fewer connectivity entries (%d)",
        n_cells, cell_connectivity_.size());
  if (n_cells == 0 && cell_connectivity_.size() > 0)
    FOUR_C_THROW(
        "There are no cells, but there are %d connectivity entries", cell_connectivity_.size());

  // Check basic face consistency. Only check if the number of entries makes sense, a more
  // sophisticated check (but also more expensive) is performed in CompleteFaceData
  if (face_offsets_.size() != 0 && face_offsets_.size() != cell_offsets_.size())
  {
    FOUR_C_THROW(
        "Error in the consistency check for the face connectivity and offsets. Did you call "
        "complete_face_connectivity?");
  }

  // Check consistency of point and cell data
  const auto check_data_dimensions =
      [&](const auto& data, const std::string& data_name, const auto& n_items)
  {
    const auto map_item = get_data_map_item(data, data_name);
    const auto size = get_data_vector_size(get_data_vector_from_map_item(map_item));
    const auto n_dim = get_data_dimension(map_item);
    if (size % n_dim != 0)
      FOUR_C_THROW("The size of the data vector %s (%d) is not a multiple of the dimension (%d)!",
          data_name.c_str(), size, n_dim);
    if (size / n_dim != n_items)
    {
      FOUR_C_THROW(
          "The size of the data vector %s (%d) does not match the total number of point/cells "
          "(%d with dimension = %d)!",
          data_name.c_str(), size, n_items, n_dim);
    }
  };

  for (const auto& data_name : GetPointDataNames())
    check_data_dimensions(point_data_, data_name, n_points);
  for (const auto& data_name : GetCellDataNames())
    check_data_dimensions(cell_data_, data_name, n_cells);
}

/**
 *
 */
void IO::VisualizationData::complete_cell_connectivity()
{
  if (cell_connectivity_.size() == 0 && cell_offsets_.size() > 0)
  {
    const auto n_points = get_point_coordinates_number_of_points();
    cell_connectivity_.reserve(n_points);
    for (size_t i_point = 0; i_point < n_points; ++i_point) cell_connectivity_.push_back(i_point);
  }
}

/**
 *
 */
void IO::VisualizationData::complete_face_connectivity()
{
  const int vtk_polyhedron_cell_id = polyhedron_cell_type;
  const size_t n_polyhedron =
      std::count(cell_types_.begin(), cell_types_.end(), vtk_polyhedron_cell_id);

  if (face_offsets_.size() == cell_offsets_.size())
  {
    // Basic face offsets consistency is fulfilled, since we have the same size in the cell and
    // face offsets
  }
  else if (n_polyhedron == 0 && face_offsets_.size() == 0 && face_connectivity_.size() == 0)
  {
    // There are no polyhedrons and the face data is empty, we can continue
  }
  else if (n_polyhedron == face_offsets_.size() && face_connectivity_.size() != 0)
  {
    // We have offsets for each polyhedron, now we need to fill the entries for all non
    // polyhedron cells
    size_t face_offset_count = 0;
    const auto face_offsets_old = face_offsets_;
    face_offsets_.clear();
    face_offsets_.reserve(cell_types_.size());
    for (const auto& cell_type : cell_types_)
    {
      if (cell_type == vtk_polyhedron_cell_id)
      {
        face_offsets_.push_back(face_offsets_old[face_offset_count]);
        face_offset_count++;
      }
      else
      {
        face_offsets_.push_back(-1);
      }
    }
  }
  else
  {
    FOUR_C_THROW(
        "Error in the consistency of the face connectivity and offsets. Number of polyhedrons "
        "(%d), cell size (%d), face "
        "connectivity size (%d), face offsets size (%d)",
        n_polyhedron, cell_offsets_.size(), face_connectivity_.size(), face_offsets_.size());
  }
}

/**
 *
 */
std::string IO::VisualizationData::get_data_type(
    const std::map<std::string, std::pair<visualization_vector_type_variant, uint8_t>>& data) const
{
  if (&point_data_ == &data) return "point data";
  if (&cell_data_ == &data) return "cell data";
  FOUR_C_THROW("Could not determine the data type");
}

/**
 *
 */
std::string IO::VisualizationData::get_data_type(
    const std::map<std::string, visualization_vector_type_variant>& data) const
{
  if (&field_data_ == &data) return "field data";
  FOUR_C_THROW("Could not determine the data type");
}

FOUR_C_NAMESPACE_CLOSE
