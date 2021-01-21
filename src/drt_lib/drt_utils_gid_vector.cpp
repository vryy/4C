/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for std vector with nodal GIDs

\level 0


*/
/*---------------------------------------------------------------------*/

#include "drt_discret.H"
#include "drt_utils_gid_vector.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::AddOwnedNodeGIDVector(Teuchos::RCP<DRT::Discretization> dis,
    const std::vector<int>& global_node_gid_vec, std::vector<int>& my_gid_vec)
{
  for (const int nodegid : global_node_gid_vec)
    if (IsNodeGIDOnThisProc(dis, nodegid)) my_gid_vec.push_back(nodegid);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::UTILS::RemoveNodeGIDsFromVector(Teuchos::RCP<DRT::Discretization> dis,
    const std::vector<int>& node_gids_to_remove, std::vector<int>& node_gid_vec)
{
  for (const int node_gid_to_remove : node_gids_to_remove)
  {
    if (IsNodeGIDOnThisProc(dis, node_gid_to_remove))
    {
      node_gid_vec.erase(std::remove(node_gid_vec.begin(), node_gid_vec.end(), node_gid_to_remove),
          node_gid_vec.end());
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::UTILS::IsNodeGIDOnThisProc(Teuchos::RCP<DRT::Discretization> dis, int node_gid)
{
  return DRT::UTILS::IsNodeGIDOnThisProc(*dis, node_gid);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::UTILS::IsNodeGIDOnThisProc(DRT::DiscretizationInterface& dis, int node_gid)
{
  return (dis.HaveGlobalNode(node_gid) and dis.gNode(node_gid)->Owner() == dis.Comm().MyPID());
}
