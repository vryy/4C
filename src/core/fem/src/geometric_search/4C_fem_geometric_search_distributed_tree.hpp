/*-----------------------------------------------------------*/
/*! \file

\brief Contains a global distributed tree search implementation
       based on bounding volumes.

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FEM_GEOMETRIC_SEARCH_DISTRIBUTED_TREE_HPP
#define FOUR_C_FEM_GEOMETRIC_SEARCH_DISTRIBUTED_TREE_HPP

#include "4C_config.hpp"

#include "4C_fem_geometric_search_bounding_volume.hpp"
#include "4C_io_pstream.hpp"

#include <vector>

class Epetra_Comm;

FOUR_C_NAMESPACE_OPEN

namespace Core::GeometricSearch
{
  struct BoundingVolume;

  /*! \brief Structure to hold the resulting global/local ID pairings that are found during a global
   * collision search
   */
  struct GlobalCollisionSearchResult
  {
    //! Local ID of the predicate (on this rank)
    int lid_predicate;
    //! Global ID of the predicate
    int gid_predicate;
    //! Local ID of the primitives (on the primitives rank)
    int lid_primitive;
    //! Global ID of the primitives
    int gid_primitive;
    //! Processor ID owning the primitive
    int pid_primitive;
  };

  /*! \brief Finds all primitives on different ranks meeting the locally owned predicates and
   * return results.
   *
   * Hereby the local/global index of the given predicate and the local/global index as well as
   * the MPI processor rank of the found primitive is returned.
   *
   * @param primitives Bounding volumes to search for intersections
   * @param predicates Bounding volumes to intersect with
   * @param comm Communicator object of the discretization
   * @param verbosity Enabeling printout of the geometric search information
   * @return Collision pairs found with their global and local IDs
   */
  std::vector<GlobalCollisionSearchResult> GlobalCollisionSearch(
      const std::vector<std::pair<int, BoundingVolume>>& primitives,
      const std::vector<std::pair<int, BoundingVolume>>& predicates, const Epetra_Comm& comm,
      const Core::IO::Verbositylevel verbosity);

}  // namespace Core::GeometricSearch

FOUR_C_NAMESPACE_CLOSE

#endif
