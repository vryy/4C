/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for std vector with nodal GIDs

\level 0


*/
/*---------------------------------------------------------------------*/

#include "baci_lib_utils_gid_vector.hpp"

#include "baci_lib_discret.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::AddOwnedNodeGID(
    const Discretization& dis, const int nodegid, std::vector<int>& my_gid_vec)
{
  if (IsNodeGIDOnThisProc(dis, nodegid)) my_gid_vec.push_back(nodegid);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::AddOwnedNodeGID(
    const Discretization& dis, const int nodegid, std::set<int>& my_gid_set)
{
  if (IsNodeGIDOnThisProc(dis, nodegid)) my_gid_set.emplace(nodegid);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::UTILS::IsNodeGIDOnThisProc(const DRT::Discretization& dis, const int node_gid)
{
  return (dis.HaveGlobalNode(node_gid) and dis.gNode(node_gid)->Owner() == dis.Comm().MyPID());
}

BACI_NAMESPACE_CLOSE
