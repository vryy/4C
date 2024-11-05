// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_xfem_dofset_transparent_independent.hpp"

#include "4C_cut_cutwizard.hpp"
#include "4C_cut_node.hpp"

FOUR_C_NAMESPACE_OPEN


XFEM::XFEMTransparentIndependentDofSet::XFEMTransparentIndependentDofSet(
    std::shared_ptr<Core::FE::Discretization> sourcedis, bool parallel,
    std::shared_ptr<Cut::CutWizard> wizard)
    : Core::DOFSets::TransparentIndependentDofSet(sourcedis, parallel), wizard_(wizard)
{
  return;
}

int XFEM::XFEMTransparentIndependentDofSet::num_dof_per_node(const Core::Nodes::Node &node) const
{
  if (wizard_ != nullptr)
  {
    Cut::Node *n = wizard_->get_node(node.id());
    if (n != nullptr)
    {
      int numdofpernode = Core::DOFSets::DofSet::num_dof_per_node(node);
      return numdofpernode * n->num_dof_sets();
    }
  }
  return Core::DOFSets::DofSet::num_dof_per_node(node);
}

FOUR_C_NAMESPACE_CLOSE
