/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Utility functions for the visualization output

\level 1

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BACI_IO_VISUALIZATION_UTILS_H
#define BACI_IO_VISUALIZATION_UTILS_H

#include "baci_config.H"

#include "baci_io_visualization_data.H"
#include "baci_linalg_fixedsizematrix.H"

#include <vector>

BACI_NAMESPACE_OPEN

namespace IO
{
  /**
   * @brief Add a polyhedron cell to a visualization data object
   *
   * @tparam n_dim Spatial dimension of the problem
   * @param visualization_data (in/out) Visualization data that the polyhedron is added to
   * @param point_coordinates (in) Points of the polyhedron
   * @param face_connectivity (in) Connectivity for each face of the polyhedron
   */
  template <unsigned int n_dim>
  void AppendPolyhedronToVisualizationData(VisualizationData& visualization_data,
      const std::vector<CORE::LINALG::Matrix<n_dim, 1>>& point_coordinates,
      const std::vector<std::vector<int>>& face_connectivity);
}  // namespace IO


BACI_NAMESPACE_CLOSE

#endif
