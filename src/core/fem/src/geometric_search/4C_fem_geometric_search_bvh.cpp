/*-----------------------------------------------------------*/
/*! \file

\brief Contains a local bounding volume hierarchy search implementation
       based on bounding volumes.

\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_fem_geometric_search_bvh.hpp"

#include "4C_fem_geometric_search_access_traits.hpp"
#include "4C_fem_geometric_search_bounding_volume.hpp"
#include "4C_fem_geometric_search_utils.hpp"
#include "4C_io_pstream.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_TimeMonitor.hpp>

#ifdef FOUR_C_WITH_ARBORX
#include <ArborX.hpp>
#endif

FOUR_C_NAMESPACE_OPEN

namespace Core::GeometricSearch
{
  std::pair<std::vector<int>, std::vector<int>> CollisionSearch(
      const std::vector<std::pair<int, BoundingVolume>>& primitives,
      const std::vector<std::pair<int, BoundingVolume>>& predicates, const Epetra_Comm& comm,
      const Core::IO::Verbositylevel verbosity)
  {
#ifndef FOUR_C_WITH_ARBORX
    FOUR_C_THROW(
        "Core::GeometricSearch::CollisionSearch can only be used with ArborX."
        "To use it, enable ArborX during the configure process.");
    return {};
#else

    TEUCHOS_FUNC_TIME_MONITOR("Core::GeometricSearch::CollisionSearch");

    std::vector<int> indices_final;
    std::vector<int> offsets_final;

    if (primitives.size() == 0 or predicates.size() == 0)
    {
      // This is the special case, where we can a priori say that there are no collisions. ArborX
      // produces a floating point exception if primitves.size() == 0, therefore, we take care of
      // this special case here.

      // Indices stays empty, as there are no collisions.
      indices_final.resize(0);
      // Offsets have to be filled so put all indices to 0, to align with the expected ArborX
      // output.
      offsets_final.resize(predicates.size() + 1, 0);
    }
    else
    {
      using memory_space = Kokkos::HostSpace;

      // Build tree structure containting all primitives.
      ArborX::BoundingVolumeHierarchy<memory_space> bounding_volume_hierarchy(
          Kokkos::DefaultExecutionSpace{}, primitives);

      Kokkos::View<int*, Kokkos::HostSpace> indices_full("indices_full", 0);
      Kokkos::View<int*, Kokkos::HostSpace> offset_full("offset_full", 0);

      // The currently used ArborX version only supports tree structures for points and
      // axis-aligned-boundary-boxes (AABB). We convert the k-DOPs to AABB in the creation of the
      // tree. Therefore, the standard query would only give the intersections of the AABB of the
      // primitives with the predicates. With this callback we perform the query in the sense, that
      // we check each k-DOP of the primitives with each k-DOP of the predicates. This can lead to a
      // drastic reduction of the found pairs.
      auto IntersectActualVolumeType =
          KOKKOS_LAMBDA(const auto predicate, const int primitive_index, const auto& out)->void
      {
        const auto& primitive_geometry = primitives[primitive_index].second.bounding_volume_;

        if (predicate(primitive_geometry)) out(primitive_index);
      };

      // Perform the collision check.
      bounding_volume_hierarchy.query(Kokkos::DefaultExecutionSpace{}, predicates,
          IntersectActualVolumeType, indices_full, offset_full);

      // Copy kokkos view to std::vector
      indices_final.insert(
          indices_final.begin(), indices_full.data(), indices_full.data() + indices_full.extent(0));
      offsets_final.insert(
          offsets_final.begin(), offset_full.data(), offset_full.data() + offset_full.extent(0));
    }

    if (verbosity == Core::IO::verbose)
    {
      Core::GeometricSearch::GeometricSearchInfo info = {static_cast<int>(primitives.size()),
          static_cast<int>(predicates.size()), static_cast<int>(indices_final.size())};
      Core::GeometricSearch::PrintGeometricSearchDetails(comm, info);
    }

    return {indices_final, offsets_final};
#endif
  }
}  // namespace Core::GeometricSearch

FOUR_C_NAMESPACE_CLOSE
