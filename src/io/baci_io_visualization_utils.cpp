/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Utility functions for the visualization output

\level 1

*/
/*-----------------------------------------------------------------------------------------------*/


#include "baci_io_visualization_utils.H"


namespace IO
{

  /**
   *
   */
  template <unsigned int n_dim>
  void AppendPolyhedronToVisualizationData(VisualizationData& visualization_data,
      const std::vector<CORE::LINALG::Matrix<n_dim, 1>>& point_coordinates,
      const std::vector<std::vector<int>>& face_connectivity)
  {
    const std::size_t n_points_data_old = visualization_data.GetPointCoordinatesNumberOfPoints();

    // Add the points, connectivity ids and offsets
    for (const auto& point : point_coordinates)
    {
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      {
        visualization_data.GetPointCoordinates().push_back(point(i_dim));
      }
    }
    for (unsigned int i_point = 0; i_point < point_coordinates.size(); i_point++)
    {
      visualization_data.GetCellConnectivity().push_back(i_point + n_points_data_old);
    }
    visualization_data.GetCellOffsets().push_back(visualization_data.GetCellConnectivity().size());

    // Add the number of polygons for this cell
    visualization_data.GetFaceConnectivity().push_back(face_connectivity.size());

    // Add the polygon faces
    for (const auto& polygon : face_connectivity)
    {
      // Set the number of points for this face
      visualization_data.GetFaceConnectivity().push_back(polygon.size());

      // Set the connectivity
      for (const auto& id : polygon)
        visualization_data.GetFaceConnectivity().push_back(id + n_points_data_old);
    }

    // Set the face offset
    visualization_data.GetFaceOffsets().push_back(visualization_data.GetFaceConnectivity().size());

    // Set the cell id
    visualization_data.GetCellTypes().push_back(42);
  }

  //! Explicit template definition
  template void AppendPolyhedronToVisualizationData(VisualizationData&,
      const std::vector<CORE::LINALG::Matrix<3, 1>>&, const std::vector<std::vector<int>>&);

}  // namespace IO
