/*----------------------------------------------------------------------*/
/*! \file

\brief transparent independent dofset

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_xfem_dofset_transparent_independent.hpp"

#include "4C_cut_cutwizard.hpp"
#include "4C_cut_node.hpp"

FOUR_C_NAMESPACE_OPEN


XFEM::XFEMTransparentIndependentDofSet::XFEMTransparentIndependentDofSet(
    Teuchos::RCP<Core::FE::Discretization> sourcedis, bool parallel,
    Teuchos::RCP<Core::Geo::CutWizard> wizard)
    : Core::DOFSets::TransparentIndependentDofSet(sourcedis, parallel), wizard_(wizard)
{
  return;
}

int XFEM::XFEMTransparentIndependentDofSet::num_dof_per_node(const Core::Nodes::Node &node) const
{
  if (wizard_ != Teuchos::null)
  {
    Core::Geo::Cut::Node *n = wizard_->get_node(node.id());
    if (n != nullptr)
    {
      int numdofpernode = Core::DOFSets::DofSet::num_dof_per_node(node);
      return numdofpernode * n->num_dof_sets();
    }
  }
  return Core::DOFSets::DofSet::num_dof_per_node(node);
}

FOUR_C_NAMESPACE_CLOSE
