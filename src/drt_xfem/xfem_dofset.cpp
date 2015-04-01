/*!-----------------------------------------------------------------------------------------------*
\file xfem_dofset.cpp

\brief provides a general XFEM dofset which uses the information from the cut-library to determine
the number of dofs per node when multiple sets of degrees of freedom per node have to be used

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*------------------------------------------------------------------------------------------------*/

#include "xfem_dofset.H"

#include "../drt_cut/cut_node.H"
#include "../drt_cut/cut_cutwizard.H"



/*----------------------------------------------------------------------*
 |  Get the gid of all dofs of a node                      schott 12/14 |
 *----------------------------------------------------------------------*/
void XFEM::XFEMDofSet::Dof(std::vector<int>& dofs, const DRT::Node* node,unsigned nds) const
{
  const int lid = node->LID();
  if (lid==-1)
    return;
  int numdf = DRT::DofSet::NumDofPerNode( *node );
  const int idx = (*idxcolnodes_)[lid] + nds*numdf;
  dofs.reserve( numdf );
  for ( int i=0; i<numdf; ++i )
  {
    dofs.push_back( idx+i );
  }
}

/*----------------------------------------------------------------------*
 |  Get the gid of all dofs of a node                      schott 12/14 |
 *----------------------------------------------------------------------*/
int XFEM::XFEMDofSet::NumDofPerNode( const DRT::Node & node ) const
{
  GEO::CUT::Node * n = wizard_.GetNode( node.Id() );
  if ( n!=NULL )
  {
    int numdofpernode = DRT::DofSet::NumDofPerNode( node );
    return numdofpernode * n->NumDofSets();
  }
  return DRT::DofSet::NumDofPerNode( node );
}
