/*----------------------------------------------------------------------*/
/*! \file

\brief transparent independent dofset

\level 1

*/
/*----------------------------------------------------------------------*/

#include "baci_xfem_dofset_transparent_independent.hpp"

#include "baci_cut_cutwizard.hpp"
#include "baci_cut_node.hpp"

BACI_NAMESPACE_OPEN


XFEM::XFEMTransparentIndependentDofSet::XFEMTransparentIndependentDofSet(
    Teuchos::RCP<DRT::Discretization> sourcedis, bool parallel,
    Teuchos::RCP<CORE::GEO::CutWizard> wizard)
    : DRT::TransparentIndependentDofSet(sourcedis, parallel), wizard_(wizard)
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
      int numdofpernode = DRT::DofSet::NumDofPerNode(node);
      return numdofpernode * n->NumDofSets();
    }
  }
  return DRT::DofSet::NumDofPerNode(node);
}

BACI_NAMESPACE_CLOSE
