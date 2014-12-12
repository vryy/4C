/*!-----------------------------------------------------------------------------------------------*
\file xfem_dofset.H

\brief provides a general XFEM dofset which uses the information from the cut-library to determine
the number of dofs per node when multiple sets of degrees of freedom per node have to be used

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*------------------------------------------------------------------------------------------------*/


#ifndef XFEM_DOFSET_H
#define XFEM_DOFSET_H


#include "../drt_lib/drt_dofset_fixed_size.H"
#include "../drt_lib/drt_discret.H"

namespace DRT
{
  class Discretization;
}

namespace GEO
{
  class CutWizard;
}

namespace XFEM
{

  class XFEMDofSet : public DRT::FixedSizeDofSet
  {
  public:

    /// constructor
    XFEMDofSet( const Teuchos::RCP<GEO::CutWizard>& wizard , int numMyReservedDofsperNode, const Teuchos::RCP<DRT::Discretization>& dis)
      : FixedSizeDofSet(numMyReservedDofsperNode,dis->NumGlobalNodes()),
        wizard_( wizard ),
        dis_(dis)
    {
    }

    /// equality relational operator for two XFEM::DofSets based on the number of nodal dofsets for all nodes
    bool operator==(XFEMDofSet const& other)
    {
      const int numnode = dis_->NumMyRowNodes();
      for(int lid=0; lid<numnode; lid++)
      {
        int gid = dis_->NodeRowMap()->GID(lid);
        DRT::Node* node = dis_->gNode(gid);
        if(NumDofPerNode( *node ) != other.NumDofPerNode( *node ))
            return false; // dofsets not equal if at least one node has a different number of nodal dofsets
      }

      return true;
    }

    /// equality relational operator for two XFEM::DofSets based on the number of dofsets for all nodes
    bool operator!=(XFEMDofSet const& other)
    {
      return !(*this==other);
    }

    /*!
    \brief Get the gid of all dofs of a node.

    Ask the current DofSet for the gids of the dofs of this node. The
    required vector is created and filled on the fly. So better keep it
    if you need more than one dof gid.
    - HaveDofs()==true prerequisite (produced by call to AssignDegreesOfFreedom()))

    Additional input nodal dof set: If the node contains more than one set of dofs, which can be
    evaluated, the number of the set needs to be given. Currently only the case for XFEM.

    \param dof (in)         : vector of dof gids (to be filled)
    \param nds (in)         : number of dofset
    \param nodaldofset (in) : number of nodal dofset
    \param node (in)        : the node
    */
    virtual void Dof(std::vector<int>& dofs, const DRT::Node* node,unsigned nodaldofset) const;

  protected:

    /// get number of nodal dofs for this element at this node
    virtual int NumDofPerNode( const DRT::Node & node ) const;

  private:

    /// the cut wizard, holds information about the number of XFEM dofsets per node
    Teuchos::RCP<GEO::CutWizard> wizard_;

    /// background discretization
    Teuchos::RCP<DRT::Discretization> dis_;
  };
}

#endif
