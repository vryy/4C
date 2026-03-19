// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_geometric_search_bvh.hpp"

#include "4C_geometric_search_access_traits.hpp"
#include "4C_geometric_search_bounding_volume.hpp"
#include "4C_geometric_search_utils.hpp"
#include "4C_io_pstream.hpp"

#include <Teuchos_TimeMonitor.hpp>

#ifdef FOUR_C_WITH_ARBORX
#include <ArborX.hpp>
#endif

FOUR_C_NAMESPACE_OPEN

namespace Core::GeometricSearch
{
  std::vector<CollisionSearchResult> collision_search(
      const std::vector<std::pair<int, BoundingVolume>>& primitives,
      const std::vector<std::pair<int, BoundingVolume>>& predicates)
  {
#ifndef FOUR_C_WITH_ARBORX
    FOUR_C_THROW(
        "Core::GeometricSearch::CollisionSearch can only be used with ArborX."
        "To use it, enable ArborX during the configure process.");
    return {};
#else

    TEUCHOS_FUNC_TIME_MONITOR("Core::GeometricSearch::CollisionSearch");

    std::vector<CollisionSearchResult> pairs;

    if (primitives.size() == 0 or predicates.size() == 0)
    {
      // This is the special case, where we can a priori say that there are no collisions. ArborX
      // produces a floating point exception if primitives.size() == 0, therefore, we take care of
      // this special case here.
      pairs.resize(0);
    }
    else
    {
      // Build tree structure containing all primitives.
      BoundingVolumeHierarchy bounding_volume_hierarchy(
          BoundingVolumeVectorPlaceholder<PrimitivesTag>{primitives});

      // Search for collisions between predicates and primitives.
      const auto [indices_full, offset_full] = bounding_volume_hierarchy.query(predicates);

      // Create the vector with the pairs.
      pairs.reserve(indices_full.size());
      for (size_t i_offset = 0; i_offset < offset_full.size() - 1; i_offset++)
      {
        const int gid_predicate = predicates[i_offset].first;
        for (int j = offset_full[i_offset]; j < offset_full[i_offset + 1]; j++)
        {
          pairs.emplace_back(CollisionSearchResult{
              .gid_predicate = gid_predicate,
              .gid_primitive = indices_full[j],
          });
        }
      }
    }

    return pairs;

#endif
  }

  std::vector<CollisionSearchResult> collision_search_print_results(
      const std::vector<std::pair<int, BoundingVolume>>& primitives,
      const std::vector<std::pair<int, BoundingVolume>>& predicates, const MPI_Comm comm,
      const Core::IO::Verbositylevel verbosity)
  {
    auto pairs = collision_search(primitives, predicates);
    print_geometric_search_details(comm,
        {.primitive_size = primitives.size(),
            .predicate_size = predicates.size(),
            .coupling_pair_size = pairs.size()},
        {.verbosity_Level = verbosity, .search_type_label = "local"});
    return pairs;
  }
}  // namespace Core::GeometricSearch

FOUR_C_NAMESPACE_CLOSE
