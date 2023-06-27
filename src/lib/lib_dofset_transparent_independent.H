/*----------------------------------------------------------------------*/
/*! \file

\brief A set of degrees of freedom on two discretizations

\level 2


*/
/*----------------------------------------------------------------------*/

#ifndef LIB_DOFSET_TRANSPARENT_INDEPENDENT_H
#define LIB_DOFSET_TRANSPARENT_INDEPENDENT_H


#include "lib_exporter.H"
#include "lib_dofset_independent.H"
#include "lib_dofset_transparent.H"

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
    virtual ~TransparentIndependentDofSet() { return; }

    /// create a copy of this object
    virtual Teuchos::RCP<DofSet> Clone() { return Teuchos::rcp(new IndependentDofSet(*this)); }

    /// Assign dof numbers to all elements and nodes of the discretization.
    virtual int AssignDegreesOfFreedom(
        const DRT::Discretization& dis, const unsigned dspos, const int start);

   protected:
    virtual int NumDofPerNode(const DRT::Node& node) const;


  };  // class TransparentDofSet
}  // namespace DRT

#endif  // LIB_DOFSET_TRANSPARENT_INDEPENDENT_H
