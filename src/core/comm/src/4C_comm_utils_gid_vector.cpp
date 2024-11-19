// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_comm_utils_gid_vector.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_node.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::Communication::add_owned_node_gid(
    const Core::FE::Discretization& dis, const int nodegid, std::vector<int>& my_gid_vec)
{
  if (is_node_gid_on_this_proc(dis, nodegid)) my_gid_vec.push_back(nodegid);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::Communication::add_owned_node_gid(
    const Core::FE::Discretization& dis, const int nodegid, std::set<int>& my_gid_set)
{
  if (is_node_gid_on_this_proc(dis, nodegid)) my_gid_set.emplace(nodegid);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Core::Communication::is_node_gid_on_this_proc(
    const Core::FE::Discretization& dis, const int node_gid)
{
  return (dis.have_global_node(node_gid) and
          dis.g_node(node_gid)->owner() == Core::Communication::my_mpi_rank(dis.get_comm()));
}

FOUR_C_NAMESPACE_CLOSE
