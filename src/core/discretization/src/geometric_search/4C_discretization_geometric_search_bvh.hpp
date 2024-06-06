/*-----------------------------------------------------------*/
/*! \file

\brief Contains a local bounding volume hierarchy search implementation
       based on bounding volumes.

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_GEOMETRIC_SEARCH_BVH_HPP
#define FOUR_C_DISCRETIZATION_GEOMETRIC_SEARCH_BVH_HPP

#include "4C_config.hpp"

#include "4C_io_pstream.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::GeometricSearch
{
  struct BoundingVolume;

  /*! \brief Finds all primitives meeting the predicates and record results in {indices, offsets}.
   *
   * (From ArborX documentation)
   * Finds all primitives meeting the predicates and record results in {indices, offsets}. indices
   * stores the indices of the objects that satisfy the predicates. offsets stores the locations in
   * indices that start a predicate, that is, predicates(i) is satisfied by primitives(indices(j))
   * for offsets(i) <= j < offsets(i+1). Following the usual convention, offsets(n) ==
   * indices.size(), where n is the number of queries that were performed and indices.size() is the
   * total number of collisions. Be aware, that this function will return all in primitives vs. all
   * in predicates!
   *
   * @param primitives Bounding volumes to search for intersections
   * @param predicates Bounding volumes to intersect with
   * @param comm Communicator object of the discretization
   * @param verbosity Enabeling printout of the geometric search information
   * @return {indices, offsets} stores indices of the objects that satisfy the predicates and
   *          offsets stores the locations in the indices view that start a predicate
   */
  std::pair<std::vector<int>, std::vector<int>> CollisionSearch(
      const std::vector<std::pair<int, BoundingVolume>>& primitives,
      const std::vector<std::pair<int, BoundingVolume>>& predicates, const Epetra_Comm& comm,
      const Core::IO::Verbositylevel verbosity);

}  // namespace Core::GeometricSearch

FOUR_C_NAMESPACE_CLOSE

#endif
