
#include "xfem_fluiddofset.H"
#include "xfem_fluidwizard.H"

#include "../drt_geometry/geo_intersection.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_cut/cut_node.H"


void XFEM::FluidDofSet::Dof( const DRT::Node * node, int nodaldofset, std::vector<int> & dofs ) const
{
  const int lid = node->LID();
  if (lid==-1)
    return;
  int numdf = DRT::DofSet::NumDofPerNode( *node, dspos_ );
  const int idx = (*idxcolnodes_)[lid] + nodaldofset*numdf;
  dofs.reserve( numdf );
  for ( int i=0; i<numdf; ++i )
  {
    dofs.push_back( idx+i );
  }
}

/// Get the gid of all dofs of a node
void XFEM::FluidDofSet::Dof(std::vector<int>& dofs, const DRT::Node* node,unsigned nodaldofset) const
{
  Dof(node,nodaldofset,dofs);
}

int XFEM::FluidDofSet::NumDofPerNode( const DRT::Node & node, unsigned dspos ) const
{
  GEO::CUT::Node * n = wizard_->CutWizard().GetNode( node.Id() );
  if ( n!=NULL )
  {
    int numdofpernode = DRT::DofSet::NumDofPerNode( node, dspos );
    return numdofpernode * n->NumDofSets();
  }
  return DRT::DofSet::NumDofPerNode( node, dspos );
}

int XFEM::FluidDofSet::AssignDegreesOfFreedom(const DRT::Discretization& dis, const unsigned dspos, const int start)
{
  int count = DRT::DofSet::AssignDegreesOfFreedom( dis, dspos, start );
  return count;
}

void XFEM::FluidDofSet::MinGID(int mingid)
{
  // set the minimal GID of the fixed-size-dofset
  minGID_ = mingid;
}
