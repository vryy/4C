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

    void Dof( DRT::Node & node, int nodaldofset, std::vector<int> & dofs ) const;

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
