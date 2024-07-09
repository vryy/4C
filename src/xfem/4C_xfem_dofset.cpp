/*----------------------------------------------------------------------*/
/*! \file

\brief provides a general XFEM dofset which uses the information from the cut-library to determine
the number of dofs per node when multiple sets of degrees of freedom per node have to be used


\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_xfem_dofset.hpp"

#include "4C_cut_cutwizard.hpp"
#include "4C_cut_node.hpp"

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*
 |  Get the gid of all dofs of a node                      schott 12/14 |
 *----------------------------------------------------------------------*/
void XFEM::XFEMDofSet::dof(
    std::vector<int>& dofs, const Core::Nodes::Node* node, unsigned nodal_dofset_id) const
{
  const int lid = node->lid();
  if (lid == -1) return;
  int numdf = Core::DOFSets::DofSet::num_dof_per_node(*node);
  const int idx = (*idxcolnodes_)[lid] + nodal_dofset_id * numdf;
  dofs.reserve(numdf);
  for (int i = 0; i < numdf; ++i)
  {
    dofs.push_back(idx + i);
  }
}

/*----------------------------------------------------------------------*
 |  Get the gid of all dofs of a node                      schott 12/14 |
 *----------------------------------------------------------------------*/
int XFEM::XFEMDofSet::num_dof_per_node(const Core::Nodes::Node& node) const
{
  Core::Geo::Cut::Node* n = wizard_.get_node(node.id());
  if (n != nullptr)
  {
    int numdofpernode = Core::DOFSets::DofSet::num_dof_per_node(node);
    return numdofpernode * n->num_dof_sets();
  }
  return Core::DOFSets::DofSet::num_dof_per_node(node);
}

FOUR_C_NAMESPACE_CLOSE
