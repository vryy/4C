/*-----------------------------------------------------------*/
/*! \file

\brief Contains a brute force search implementation based on
       bounding boxes.

\level 3

*/
/*-----------------------------------------------------------*/

#include <Epetra_MpiComm.h>

#include "io_pstream.H"

#include <Teuchos_TimeMonitor.hpp>

#ifdef HAVE_ARBORX
#include <ArborX.hpp>
#endif


#include "discretization_geometric_search.H"
#include "discretization_geometric_search_utils.H"

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
      return ArborX::Box{vector[i].second.bounding_volume_};
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

      Kokkos::View<int*, Kokkos::HostSpace> indices_full("indices_full", 0);
      Kokkos::View<int*, Kokkos::HostSpace> offset_full("offset_full", 0);

      // Perform the collision check.
      bounding_volume_hierarchy.query(
          Kokkos::DefaultExecutionSpace{}, predicates, indices_full, offset_full);

      // ArborX only supports tree structures for points and axis-aligned-boundary-boxes (AABB). We
      // convert the k-DOPs to AABB in the creation of the tree. Therefore, the query only gives the
      // intersections of the AABB of the primitives with the predicates. In this loop we perform a
      // post-processing of the query results in the sense, that we check each k-DOP of the
      // primitives with each k-DOP of the predicates. This only marginally effects the performance
      // of this function, but can lead to a drastic reduction of the found pairs.
      offsets_final.push_back(0);
      for (size_t i_predicate = 0; i_predicate < predicates.size(); i_predicate++)
      {
        int collisions_for_predicate_i = 0;
        const auto& predicate = predicates[i_predicate].second.bounding_volume_;
        for (int j = offset_full[i_predicate]; j < offset_full[i_predicate + 1]; j++)
        {
          // Check for actual intersection.
          int i_primitive = indices_full[j];
          if (ArborX::Experimental::intersects(
                  predicate, primitives[i_primitive].second.bounding_volume_))
          {
            indices_final.push_back(i_primitive);
            collisions_for_predicate_i += 1;
          }
        }

        offsets_final.push_back(offsets_final.back() + collisions_for_predicate_i);
      }
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
