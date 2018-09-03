/*----------------------------------------------------------------------*/
/*!
\file drt_dofset_transparent_independent.H

\brief A set of degrees of freedom special for ???

\level 2

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/

#ifndef TRANSPARENTINDEPENDENTDOFSET_H
#define TRANSPARENTINDEPENDENTDOFSET_H


#include "drt_exporter.H"
#include "drt_dofset_independent.H"
#include "drt_dofset_transparent.H"

namespace DRT
{
  class Discretization;
}

namespace GEO
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

#endif  // #ifndef  TRANSPARENTINDEPENDENTDOFSET_H
