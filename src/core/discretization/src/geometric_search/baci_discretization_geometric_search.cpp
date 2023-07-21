/*-----------------------------------------------------------*/
/*! \file

\brief Contains a brute force search implementation based on
       bounding boxes.

\level 3

*/
/*-----------------------------------------------------------*/

#include <Epetra_MpiComm.h>

#include "baci_io_pstream.H"

#include <Teuchos_TimeMonitor.hpp>

#ifdef BACI_WITH_ARBORX
#include <ArborX.hpp>
#endif

#include "baci_discretization_geometric_search.H"
#include "baci_discretization_geometric_search_utils.H"

#ifdef BACI_WITH_ARBORX
namespace ArborX
{
  template <>
  struct AccessTraits<std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>>,
      PrimitivesTag>
  {
    using memory_space = Kokkos::HostSpace;

    static std::size_t size(
        const std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>>& vector)
    {
      return vector.size();
    }

    static auto get(
        const std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>>& vector,
        std::size_t i)
    {
      return ArborX::Box{vector[i].second.bounding_volume_};
    }
  };

  template <>
  struct AccessTraits<std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>>,
      PredicatesTag>
  {
    using memory_space = Kokkos::HostSpace;

    static std::size_t size(
        const std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>>& vector)
    {
      return vector.size();
    }

    static auto get(
        const std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>>& vector,
        std::size_t i)
    {
      return intersects(vector[i].second.bounding_volume_);
    }
  };
}  // namespace ArborX

#endif

namespace CORE::GEOMETRICSEARCH
{
  std::pair<std::vector<int>, std::vector<int>> CollisionSearch(
      const std::vector<std::pair<int, BoundingVolume>>& primitives,
      const std::vector<std::pair<int, BoundingVolume>>& predicates, const Epetra_Comm& comm,
      const IO::verbositylevel verbosity)
  {
#ifndef BACI_WITH_ARBORX
    dserror(
        "CORE::GEOMETRICSEARCH::CollisionSearch can only be used with ArborX."
        "To use it, enable ArborX during the configure process.");
    return {};
#else

    TEUCHOS_FUNC_TIME_MONITOR("CORE::GEOMETRICSEARCH::CollisionSearch");

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

    if (verbosity == IO::verbose)
    {
      CORE::GEOMETRICSEARCH::GeometricSearchInfo info = {static_cast<int>(primitives.size()),
          static_cast<int>(predicates.size()), static_cast<int>(indices_final.size())};
      CORE::GEOMETRICSEARCH::PrintGeometricSearchDetails(comm, info);
    }

    return {indices_final, offsets_final};
#endif
  }
}  // namespace CORE::GEOMETRICSEARCH
