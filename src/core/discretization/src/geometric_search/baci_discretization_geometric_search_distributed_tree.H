/*-----------------------------------------------------------*/
/*! \file

\brief Contains a global distributed tree search implementation
       based on bounding volumes.

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_DISCRETIZATION_GEOMETRIC_SEARCH_DISTRIBUTED_TREE_H
#define BACI_DISCRETIZATION_GEOMETRIC_SEARCH_DISTRIBUTED_TREE_H

#include "baci_config.H"

#include "baci_discretization_geometric_search_bounding_volume.H"
#include "baci_io_pstream.H"

#include <vector>

class Epetra_Comm;

BACI_NAMESPACE_OPEN

namespace CORE::GEOMETRICSEARCH
{
  struct BoundingVolume;

  /*! \brief Finds all primitives on different ranks meeting the locally owned predicates and
   * records results in {lid_predicate, gid_predicate, lid_primitive, gid_primitive, pid_primitive}.
   *
   * Hereby the local/global index of the given predicate and the local/global index as well as
   * the MPI processor rank of the found primitive is returned.
   *
   * @param primitives Bounding volumes to search for intersections
   * @param predicates Bounding volumes to intersect with
   * @param comm Communicator object of the discretization
   * @param verbosity Enabeling printout of the geometric search information
   * @return {lid_predicate, gid_predicate, lid_primitive, gid_primitive, pid_primitive}
   */
  std::vector<std::tuple<int, int, int, int, int>> GlobalCollisionSearch(
      const std::vector<std::pair<int, BoundingVolume>>& primitives,
      const std::vector<std::pair<int, BoundingVolume>>& predicates, const Epetra_Comm& comm,
      const IO::verbositylevel verbosity);

}  // namespace CORE::GEOMETRICSEARCH

BACI_NAMESPACE_CLOSE

#endif
