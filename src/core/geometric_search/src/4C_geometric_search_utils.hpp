// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_GEOMETRIC_SEARCH_UTILS_HPP
#define FOUR_C_GEOMETRIC_SEARCH_UTILS_HPP

#include "4C_config.hpp"

#include "4C_geometric_search_bounding_volume.hpp"
#include "4C_io_pstream.hpp"

#include <mpi.h>

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::GeometricSearch
{
  /*! \brief Storing information on the geometric search
   */
  struct GeometricSearchInfo
  {
    size_t primitive_size;
    size_t predicate_size;
    size_t coupling_pair_size;
  };

  /*! \brief Options to steer the output of the geometric search
   */
  struct GeometricSearchOutput
  {
    Core::IO::Verbositylevel verbosity_Level;
    std::string search_type_label;
  };

  /*! \brief Prints details on the geometric search algorithm
   */
  void print_geometric_search_details(
      const MPI_Comm comm, const GeometricSearchInfo info, const GeometricSearchOutput options);

  /*! \brief Get the polyhedron representation of a k-DOP
   *
   * @param boundingVolume Bounding volume enclosing the respective element (as k-DOP)
   * @return Points of the polyhedron and connecting polygons
   */
  std::pair<std::vector<LinAlg::Matrix<3, 1>>, std::vector<std::vector<int>>>
  get_k_dop_polyhedron_representation(const BoundingVolume boundingVolume);

}  // namespace Core::GeometricSearch

FOUR_C_NAMESPACE_CLOSE

#endif
