// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

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
  void add_owned_node_gid(
      const Core::FE::Discretization& dis, int nodegid, std::vector<int>& my_gid_vec);
  void add_owned_node_gid(
      const Core::FE::Discretization& dis, int nodegid, std::set<int>& my_gid_set);

  /*!
   * \brief Add nodal GIDs on this processor to existing list from list with global GIDs
   *
   * @param[in] dis                   discretization, that holds nodes with GIDs
   * @param[in] global_node_gid_vec   vector/set with all node GIDs
   * @param[out] my_gid_vec           vector/set with my node GIDs
   */
  template <typename T, typename U>
  void add_owned_node_gid_from_list(
      const Core::FE::Discretization& dis, const T& global_node_gid_vec, U& my_gid_list)
  {
    for (const int nodegid : global_node_gid_vec) add_owned_node_gid(dis, nodegid, my_gid_list);
  }

  /*!
   * \brief check, whether node with GID is owned by this processor
   *
   * @param[in] dis                   discretization, that holds nodes with GIDs
   * @param[in] node_gid              GID of node to be checked
   * @return                          indicates, whether node is owned by this processor
   */
  bool is_node_gid_on_this_proc(const Core::FE::Discretization& dis, int node_gid);
}  // namespace Core::Communication

FOUR_C_NAMESPACE_CLOSE

#endif
