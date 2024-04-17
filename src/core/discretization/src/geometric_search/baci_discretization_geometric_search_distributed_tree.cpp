/*-----------------------------------------------------------*/
/*! \file

\brief Contains a global distributed tree search implementation
       based on bounding volumes.

\level 3

*/
/*-----------------------------------------------------------*/

#include "baci_discretization_geometric_search_distributed_tree.hpp"

#include "baci_discretization_geometric_search_access_traits.hpp"
#include "baci_discretization_geometric_search_bounding_volume.hpp"
#include "baci_discretization_geometric_search_utils.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_TimeMonitor.hpp>

#ifdef BACI_WITH_ARBORX
#include <ArborX.hpp>
#endif

FOUR_C_NAMESPACE_OPEN

namespace CORE::GEOMETRICSEARCH
{
  std::vector<std::tuple<int, int, int, int, int>> GlobalCollisionSearch(
      const std::vector<std::pair<int, BoundingVolume>>& primitives,
      const std::vector<std::pair<int, BoundingVolume>>& predicates, const Epetra_Comm& comm,
      const IO::verbositylevel verbosity)
  {
#ifndef BACI_WITH_ARBORX
    dserror(
        "GEOMETRICSEARCH::GlobalCollisionSearch can only be used with ArborX."
        "To use it, enable ArborX during the configure process.");
    return {};
#else

    TEUCHOS_FUNC_TIME_MONITOR("CORE::GEOMETRICSEARCH::GlobalCollisionSearch");

    int myrank = comm.MyPID();

    using memory_space = Kokkos::HostSpace;

    // Build tree structure containting all primitives.
    ArborX::DistributedTree<memory_space> distributed_tree(
        dynamic_cast<const Epetra_MpiComm*>(&comm)->Comm(), Kokkos::DefaultExecutionSpace{},
        primitives);

    // TODO: Check for better data structure in Kokkos (something like a tuple)
    Kokkos::View<Kokkos::pair<int, Kokkos::pair<int, int>>*, Kokkos::HostSpace> indices_ranks_full(
        "indices_ranks_full", 0);
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
      const auto& primitive_vector_entry = primitives[primitive_index];
      const auto& primitive_geometry = primitive_vector_entry.second.bounding_volume_;
      if (predicate(primitive_geometry))
        out({primitive_index, {primitive_vector_entry.first, myrank}});
    };

    // Perform the collision check.
    distributed_tree.query(Kokkos::DefaultExecutionSpace{}, predicates, IntersectActualVolumeType,
        indices_ranks_full, offset_full);

    // Create the vector with the pairs.
    // {lid_predicate, gid_predicate, lid_primitive, gid_primitive, pid_primitive}
    std::vector<std::tuple<int, int, int, int, int>> pairs;
    for (size_t i_offset = 0; i_offset < offset_full.size() - 1; i_offset++)
    {
      const int gid_predicate = predicates[i_offset].first;
      for (int j = offset_full[i_offset]; j < offset_full[i_offset + 1]; j++)
      {
        pairs.emplace_back(std::tuple{i_offset, gid_predicate, indices_ranks_full[j].first,
            indices_ranks_full[j].second.first, indices_ranks_full[j].second.second});
      }
    }

    if (verbosity == IO::verbose)
    {
      CORE::GEOMETRICSEARCH::GeometricSearchInfo info = {static_cast<int>(primitives.size()),
          static_cast<int>(predicates.size()), static_cast<int>(indices_ranks_full.size())};
      CORE::GEOMETRICSEARCH::PrintGeometricSearchDetails(comm, info);
    }

    return pairs;
#endif
  }
}  // namespace CORE::GEOMETRICSEARCH

FOUR_C_NAMESPACE_CLOSE
