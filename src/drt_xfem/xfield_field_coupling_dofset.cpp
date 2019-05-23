/*----------------------------------------------------------------------------*/
/**

\brief DoF set for coupling a xfield and a field discretization at a common
       interface

\maintainer Matthias Mayr

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "xfield_field_coupling_dofset.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XFEM::XFieldField::CouplingDofSet::CouplingDofSet(const int& my_num_reserve_dof_per_node,
    const int& g_node_index_range, const int& g_num_std_dof_per_node,
    const std::map<int, int>& my_num_dofs_per_node)
    : DRT::FixedSizeDofSet(my_num_reserve_dof_per_node, g_node_index_range),
      my_num_dof_per_node_(my_num_dofs_per_node),
      g_num_std_dof_per_node_(g_num_std_dof_per_node)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::XFieldField::CouplingDofSet::Dof(
    std::vector<int>& dofs, const DRT::Node* node, unsigned nodal_dofset_id) const
{
  const int lid = node->LID();
  if (lid == -1) return;
  const int num_dof = NumStandardDofPerNode();
  const int idx = (*idxcolnodes_)[lid] + nodal_dofset_id * num_dof;
  dofs.resize(num_dof, 0);
  for (int i = 0; i < num_dof; ++i) dofs[i] = idx + i;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XFEM::XFieldField::CouplingDofSet::NumDofPerNode(const DRT::Node& node) const
{
  return MyNumDofPerNode(node.Id());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XFEM::XFieldField::CouplingDofSet::MyNumDofPerNode(const int& node_gid) const
{
  std::map<int, int>::const_iterator pos = my_num_dof_per_node_.find(node_gid);
  if (pos == my_num_dof_per_node_.end())
    dserror("The given node GID %d is no coupling interface node!", node_gid);

  return pos->second;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XFEM::XFieldField::CouplingDofSet::NumStandardDofPerNode() const
{
  return g_num_std_dof_per_node_;
}
