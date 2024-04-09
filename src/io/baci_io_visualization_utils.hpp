/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Utility functions for the visualization output

\level 1

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_IO_VISUALIZATION_UTILS_HPP
#define FOUR_C_IO_VISUALIZATION_UTILS_HPP

#include "baci_config.hpp"

#include "baci_io_visualization_data.hpp"
#include "baci_linalg_fixedsizematrix.hpp"

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
