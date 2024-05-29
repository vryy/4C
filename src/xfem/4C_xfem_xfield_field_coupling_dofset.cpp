/*----------------------------------------------------------------------------*/
/** \file

\brief DoF set for coupling a xfield and a field discretization at a common
       interface


\level 3

*/
/*----------------------------------------------------------------------------*/

#include "4C_xfem_xfield_field_coupling_dofset.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XFEM::XFieldField::CouplingDofSet::CouplingDofSet(const int& my_num_reserve_dof_per_node,
    const int& g_node_index_range, const int& g_num_std_dof_per_node,
    const std::map<int, int>& my_num_dofs_per_node)
    : CORE::Dofsets::FixedSizeDofSet(my_num_reserve_dof_per_node, g_node_index_range),
      my_num_dof_per_node_(my_num_dofs_per_node),
      g_num_std_dof_per_node_(g_num_std_dof_per_node)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XFEM::XFieldField::CouplingDofSet::Dof(
    std::vector<int>& dofs, const CORE::Nodes::Node* node, unsigned nodal_dofset_id) const
{
  const int lid = node->LID();
  if (lid == -1) return;
  const int num_dof = num_standard_dof_per_node();
  const int idx = (*idxcolnodes_)[lid] + nodal_dofset_id * num_dof;
  dofs.resize(num_dof, 0);
  for (int i = 0; i < num_dof; ++i) dofs[i] = idx + i;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XFEM::XFieldField::CouplingDofSet::NumDofPerNode(const CORE::Nodes::Node& node) const
{
  return my_num_dof_per_node(node.Id());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XFEM::XFieldField::CouplingDofSet::my_num_dof_per_node(const int& node_gid) const
{
  std::map<int, int>::const_iterator pos = my_num_dof_per_node_.find(node_gid);
  if (pos == my_num_dof_per_node_.end())
    FOUR_C_THROW("The given node GID %d is no coupling interface node!", node_gid);

  return pos->second;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XFEM::XFieldField::CouplingDofSet::num_standard_dof_per_node() const
{
  return g_num_std_dof_per_node_;
}

FOUR_C_NAMESPACE_CLOSE
