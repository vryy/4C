/*----------------------------------------------------------------------*/
/*!
\file xfem_dofset_transparent_independent.cpp

\brief transparent independent dofset

\level 1

<pre>
\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
*/
/*----------------------------------------------------------------------*/

#include "xfem_dofset_transparent_independent.H"

#include "../drt_cut/cut_node.H"
#include "../drt_cut/cut_cutwizard.H"


XFEM::XFEMTransparentIndependentDofSet::XFEMTransparentIndependentDofSet(
  Teuchos::RCP<DRT::Discretization> sourcedis,
  bool parallel,
  Teuchos::RCP<GEO::CutWizard> wizard)
  : DRT::TransparentIndependentDofSet(sourcedis, parallel),
    wizard_(wizard)
{
  return;
}

int XFEM::XFEMTransparentIndependentDofSet::NumDofPerNode( const DRT::Node & node ) const
{
  if (wizard_ != Teuchos::null)
  {
    GEO::CUT::Node *n = wizard_->GetNode( node.Id() );
    if ( n!=NULL )
    {
      int numdofpernode = DRT::DofSet::NumDofPerNode( node );
      return numdofpernode * n->NumDofSets();
    }
  }
  return DRT::DofSet::NumDofPerNode( node );
}
