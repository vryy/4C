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
    Teuchos::RCP<DRT::Discretization> sourcedis, bool parallel,
    Teuchos::RCP<CORE::GEO::CutWizard> wizard)
    : CORE::Dofsets::TransparentIndependentDofSet(sourcedis, parallel), wizard_(wizard)
{
  return;
}

int XFEM::XFEMTransparentIndependentDofSet::NumDofPerNode(const DRT::Node &node) const
{
  if (wizard_ != Teuchos::null)
  {
    CORE::GEO::CUT::Node *n = wizard_->GetNode(node.Id());
    if (n != nullptr)
    {
      int numdofpernode = CORE::Dofsets::DofSet::NumDofPerNode(node);
      return numdofpernode * n->NumDofSets();
    }
  }
  return CORE::Dofsets::DofSet::NumDofPerNode(node);
}

FOUR_C_NAMESPACE_CLOSE
