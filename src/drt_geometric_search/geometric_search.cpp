/*-----------------------------------------------------------*/
/*! \file

\brief Contains a brute force search implementation based on
       bounding boxes.

\level 3

*/
/*-----------------------------------------------------------*/

#include <Epetra_MpiComm.h>

#include "geometric_search.H"
#include "geometric_search_utils.H"

namespace GEOMETRICSEARCH
{
  std::pair<std::vector<int>, std::vector<int>> CollisionSearch(
      const std::vector<std::pair<int, BoundingVolume>>& primitives,
      const std::vector<std::pair<int, BoundingVolume>>& predicates, const Epetra_Comm& comm,
      const IO::verbositylevel verbosity)
  {
    std::vector<int> indices, offsets;
    offsets.push_back(0);

    for (const auto& predicate : predicates)
    {
      int predicate_collisions = 0;
      for (size_t i_primitive = 0; i_primitive < primitives.size(); i_primitive++)
      {
        const auto& primitive = primitives[i_primitive];
        if (BoxesIntersect(primitive.second, predicate.second))
        {
          indices.push_back(i_primitive);
          predicate_collisions += 1;
        }
      }
      offsets.push_back(offsets.back() + predicate_collisions);
    }

    if (verbosity == IO::verbose)
    {
      UTILS::GeometricSearchInfo info = {static_cast<int>(primitives.size()),
          static_cast<int>(predicates.size()), static_cast<int>(indices.size())};
      UTILS::PrintGeometricSearchDetails(comm, info);
    }

    return {indices, offsets};
  }
}  // namespace GEOMETRICSEARCH