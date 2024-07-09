/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Utility functions for the visualization output

\level 1

*/
/*-----------------------------------------------------------------------------------------------*/


#include "4C_io_visualization_utils.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Core::IO
{

  /**
   *
   */
  template <unsigned int n_dim>
  void AppendPolyhedronToVisualizationData(VisualizationData& visualization_data,
      const std::vector<Core::LinAlg::Matrix<n_dim, 1>>& point_coordinates,
      const std::vector<std::vector<int>>& face_connectivity)
  {
    const std::size_t n_points_data_old =
        visualization_data.get_point_coordinates_number_of_points();

    // Add the points, connectivity ids and offsets
    for (const auto& point : point_coordinates)
    {
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      {
        visualization_data.get_point_coordinates().push_back(point(i_dim));
      }
    }
    for (unsigned int i_point = 0; i_point < point_coordinates.size(); i_point++)
    {
      visualization_data.get_cell_connectivity().push_back(i_point + n_points_data_old);
    }
    visualization_data.get_cell_offsets().push_back(
        visualization_data.get_cell_connectivity().size());

    // Add the number of polygons for this cell
    visualization_data.get_face_connectivity().push_back(face_connectivity.size());

    // Add the polygon faces
    for (const auto& polygon : face_connectivity)
    {
      // Set the number of points for this face
      visualization_data.get_face_connectivity().push_back(polygon.size());

      // Set the connectivity
      for (const auto& id : polygon)
        visualization_data.get_face_connectivity().push_back(id + n_points_data_old);
    }

    // Set the face offset
    visualization_data.get_face_offsets().push_back(
        visualization_data.get_face_connectivity().size());

    // Set the cell id
    visualization_data.get_cell_types().push_back(polyhedron_cell_type);
  }

  //! Explicit template definition
  template void AppendPolyhedronToVisualizationData(VisualizationData&,
      const std::vector<Core::LinAlg::Matrix<3, 1>>&, const std::vector<std::vector<int>>&);

}  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE
