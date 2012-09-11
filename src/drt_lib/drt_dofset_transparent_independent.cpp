/*!----------------------------------------------------------------------
\file drt_dofset_transparent_independent.cpp

\brief transparent independent dofset

<pre>
Maintainer: Shadan Shahmiri
            Shahmiri@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/

#include "drt_dofset_transparent_independent.H"
#include "drt_dofset.H"
#include "../linalg/linalg_utils.H"
#include "../drt_cut/cut_node.H"
#include "../drt_xfem/xfem_fluidwizard.H"
#include "../drt_geometry/geo_intersection.H"


DRT::TransparentIndependentDofSet::TransparentIndependentDofSet(
  RCP<DRT::Discretization> sourcedis,
  bool parallel,
  Teuchos::RCP<XFEM::FluidWizard> wizard = Teuchos::null)
  : DRT::TransparentDofSet(sourcedis, parallel),
    wizard_(wizard)
{
  return;
}

int DRT::TransparentIndependentDofSet::AssignDegreesOfFreedom(const DRT::Discretization& dis, const unsigned dspos, const int start)
{

  // first, we call the standard AssignDegreesOfFreedom from the base class
  int count = DRT::IndependentDofSet::AssignDegreesOfFreedom(dis,dspos,start);

  if(!parallel_)
  {
    TransferDegreesOfFreedom(*sourcedis_, dis, start);
  }
  else
  {
    ParallelTransferDegreesOfFreedom(*sourcedis_, dis, start);
  }

  // tell all proxies (again!)
  NotifyAssigned();

  return count;
}

int DRT::TransparentIndependentDofSet::NumDofPerNode( const DRT::Node & node, unsigned dspos ) const
{
  if (wizard_ != Teuchos::null)
  {
    GEO::CUT::Node *n = wizard_->CutWizard().GetNode( node.Id() );
    if ( n!=NULL )
    {
      int numdofpernode = DRT::DofSet::NumDofPerNode( node, dspos );
      return numdofpernode * n->NumDofSets();
    }
  }
  return DRT::DofSet::NumDofPerNode( node, dspos );
}


