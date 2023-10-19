/*----------------------------------------------------------------------*/
/*! \file

\brief A set of degrees of freedom on two discretizations

\level 2


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_LIB_DOFSET_TRANSPARENT_INDEPENDENT_H
#define BACI_LIB_DOFSET_TRANSPARENT_INDEPENDENT_H


#include "baci_lib_dofset_independent.H"
#include "baci_lib_dofset_transparent.H"
#include "baci_lib_exporter.H"

namespace DRT
{
  class Discretization;
}

namespace CORE::GEO
{
  class CutWizard;
}

namespace DRT
{
  /// Alias dofset that shares dof numbers with another dofset
  /*!
  A special set of degrees of freedom, implemented in order to assign the same degrees of freedom to
  nodes belonging to two discretizations. This way two discretizations can assemble into the same
  position of the system matrix. As internal variable it holds a source discretization
  (Constructor). If such a nodeset is assigned to a sub-discretization, its dofs are assigned
  according to the dofs of the source.

  */
  class TransparentIndependentDofSet : public IndependentDofSet, public TransparentDofSet
  {
   public:
    /*!
    \brief Standard Constructor
    */
    explicit TransparentIndependentDofSet(
        Teuchos::RCP<DRT::Discretization> sourcedis, bool parallel);

    /*!
    \brief Destructor
    */
    ~TransparentIndependentDofSet() override { return; }

    /// create a copy of this object
    Teuchos::RCP<DofSet> Clone() override { return Teuchos::rcp(new IndependentDofSet(*this)); }

    /// Assign dof numbers to all elements and nodes of the discretization.
    int AssignDegreesOfFreedom(
        const DRT::Discretization& dis, const unsigned dspos, const int start) override;

   protected:
    int NumDofPerNode(const DRT::Node& node) const override;


  };  // class TransparentDofSet
}  // namespace DRT

#endif  // LIB_DOFSET_TRANSPARENT_INDEPENDENT_H
