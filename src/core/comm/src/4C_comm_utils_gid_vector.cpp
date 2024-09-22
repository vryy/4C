/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for std vector with nodal GIDs

\level 0


*/
/*---------------------------------------------------------------------*/

#include "4C_comm_utils_gid_vector.hpp"

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
  return (
      dis.have_global_node(node_gid) and dis.g_node(node_gid)->owner() == dis.get_comm().MyPID());
}

FOUR_C_NAMESPACE_CLOSE
