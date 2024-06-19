/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for std vector with nodal GIDs

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_COMM_UTILS_GID_VECTOR_HPP
#define FOUR_C_COMM_UTILS_GID_VECTOR_HPP

#include "4C_config.hpp"

#include <set>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE
namespace Core::Communication
{
  /*!
   * \brief Add nodal GID on this processor to existing list of GIDs
   *
   * @param[in] dis                   discretization, that holds nodes with GIDs
   * @param[in] nodegid               nodal GID
   * @param[out] my_gid_vec           vector/set with my node GIDs
   */
  void AddOwnedNodeGID(
      const Core::FE::Discretization& dis, int nodegid, std::vector<int>& my_gid_vec);
  void AddOwnedNodeGID(const Core::FE::Discretization& dis, int nodegid, std::set<int>& my_gid_set);

  /*!
   * \brief Add nodal GIDs on this processor to existing list from list with global GIDs
   *
   * @param[in] dis                   discretization, that holds nodes with GIDs
   * @param[in] global_node_gid_vec   vector/set with all node GIDs
   * @param[out] my_gid_vec           vector/set with my node GIDs
   */
  template <typename T, typename U>
  void AddOwnedNodeGIDFromList(
      const Core::FE::Discretization& dis, const T& global_node_gid_vec, U& my_gid_list)
  {
    for (const int nodegid : global_node_gid_vec) AddOwnedNodeGID(dis, nodegid, my_gid_list);
  }

  /*!
   * \brief check, whether node with GID is owned by this processor
   *
   * @param[in] dis                   discretization, that holds nodes with GIDs
   * @param[in] node_gid              GID of node to be checked
   * @return                          indicates, whether node is owned by this processor
   */
  bool IsNodeGIDOnThisProc(const Core::FE::Discretization& dis, int node_gid);
}  // namespace Core::Communication

FOUR_C_NAMESPACE_CLOSE

#endif
