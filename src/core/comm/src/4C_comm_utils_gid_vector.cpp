/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for std vector with nodal GIDs

\level 0


*/
/*---------------------------------------------------------------------*/

#include "4C_comm_utils_gid_vector.hpp"

#include "4C_fem_discretization.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::Communication::AddOwnedNodeGID(
    const Core::FE::Discretization& dis, const int nodegid, std::vector<int>& my_gid_vec)
{
  if (IsNodeGIDOnThisProc(dis, nodegid)) my_gid_vec.push_back(nodegid);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Core::Communication::AddOwnedNodeGID(
    const Core::FE::Discretization& dis, const int nodegid, std::set<int>& my_gid_set)
{
  if (IsNodeGIDOnThisProc(dis, nodegid)) my_gid_set.emplace(nodegid);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Core::Communication::IsNodeGIDOnThisProc(
    const Core::FE::Discretization& dis, const int node_gid)
{
  return (dis.HaveGlobalNode(node_gid) and dis.gNode(node_gid)->Owner() == dis.Comm().MyPID());
}

FOUR_C_NAMESPACE_CLOSE
