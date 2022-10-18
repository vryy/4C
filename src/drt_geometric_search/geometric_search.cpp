/*-----------------------------------------------------------*/
/*! \file

\brief Contains a brute force search implementation based on
       bounding boxes.

\level 3

*/
/*-----------------------------------------------------------*/

#include <Epetra_MpiComm.h>

#include "../drt_io/io_pstream.H"

#include <Teuchos_TimeMonitor.hpp>

#ifdef HAVE_ARBORX
#include <ArborX.hpp>
#endif


#include "geometric_search.H"
#include "geometric_search_utils.H"

#ifdef HAVE_ARBORX
namespace ArborX
{
  template <>
  struct AccessTraits<std::vector<std::pair<int, GEOMETRICSEARCH::BoundingVolume>>, PrimitivesTag>
  {
    using memory_space = Kokkos::HostSpace;

    static std::size_t size(
        const std::vector<std::pair<int, GEOMETRICSEARCH::BoundingVolume>>& vector)
    {
      return vector.size();
    }

    static auto get(
        const std::vector<std::pair<int, GEOMETRICSEARCH::BoundingVolume>>& vector, std::size_t i)
    {
      return vector[i].second.bounding_volume_;
    }
  };

  template <>
  struct AccessTraits<std::vector<std::pair<int, GEOMETRICSEARCH::BoundingVolume>>, PredicatesTag>
  {
    using memory_space = Kokkos::HostSpace;

    static std::size_t size(
        const std::vector<std::pair<int, GEOMETRICSEARCH::BoundingVolume>>& vector)
    {
      return vector.size();
    }

    static auto get(
        const std::vector<std::pair<int, GEOMETRICSEARCH::BoundingVolume>>& vector, std::size_t i)
    {
      return intersects(vector[i].second.bounding_volume_);
    }
  };
}  // namespace ArborX
#endif

namespace GEOMETRICSEARCH
{
  std::pair<std::vector<int>, std::vector<int>> CollisionSearch(
      const std::vector<std::pair<int, BoundingVolume>>& primitives,
      const std::vector<std::pair<int, BoundingVolume>>& predicates, const Epetra_Comm& comm,
      const IO::verbositylevel verbosity)
  {
#ifndef HAVE_ARBORX
    dserror(
        "GEOMETRICSEARCH::CollisionSearch can only be used with ArborX."
        "To use it, enable ArborX during the configure process.");
    return {};
#else

    TEUCHOS_FUNC_TIME_MONITOR("GEOMETRICSEARCH::CollisionSearch");

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


      Kokkos::View<int*, memory_space> indices("indices", 0);
      Kokkos::View<int*, memory_space> offset("offset", 0);

      // Perform the collision check.
      bounding_volume_hierarchy.query(Kokkos::DefaultExecutionSpace{}, predicates, indices, offset);

      // Convert the Kokkos data arrays to std vectors.
      indices_final.insert(
          indices_final.begin(), indices.data(), indices.data() + indices.extent(0));
      offsets_final.insert(offsets_final.begin(), offset.data(), offset.data() + offset.extent(0));
    }

    if (verbosity == IO::verbose)
    {
      UTILS::GeometricSearchInfo info = {static_cast<int>(primitives.size()),
          static_cast<int>(predicates.size()), static_cast<int>(indices_final.size())};
      UTILS::PrintGeometricSearchDetails(comm, info);
    }

    return {indices_final, offsets_final};
#endif
  }
}  // namespace GEOMETRICSEARCH
