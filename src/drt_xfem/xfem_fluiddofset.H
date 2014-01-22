#ifndef XFEM_FLUIDDOFSET_H
#define XFEM_FLUIDDOFSET_H


#include "../drt_lib/drt_dofset_fixed_size.H"


namespace XFEM
{
  class FluidWizard;

  class FluidDofSet : public DRT::FixedSizeDofSet
  {
  public:

    FluidDofSet( FluidWizard * wizard , int numMyReservedDofs , const DRT::Discretization& dis)
      : FixedSizeDofSet(numMyReservedDofs), wizard_( wizard ), dis_ (dis)
    {

    }

    void Dof( const DRT::Node * node, int nodaldofset, std::vector<int> & dofs ) const;

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

    virtual void MinGID(int mingid);

    /// Assign dof numbers using all elements and nodes of the discretization.
    virtual int AssignDegreesOfFreedom(const DRT::Discretization& dis, const unsigned dspos, const int start);

  protected:

    /// get number of nodal dofs for this element at this node
    virtual int NumDofPerNode( const DRT::Node & node, unsigned dspos ) const;

  private:

    FluidWizard * wizard_;

    const DRT::Discretization& dis_;
  };
}

#endif
