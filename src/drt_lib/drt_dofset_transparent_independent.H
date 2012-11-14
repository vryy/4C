/*!----------------------------------------------------------------------
\file drt_dofset_transparent_independent.H

\brief A set of degrees of freedom special for contact

<pre>
Maintainer: Shadan Shahmiri
            shahmiri@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifndef TRANSPARENTINDEPENDENTDOFSET_H
#define TRANSPARENTINDEPENDENTDOFSET_H


#include "drt_exporter.H"
#include "drt_dofset_independent.H"
#include "drt_dofset_transparent.H"
#include "../drt_lib/drt_discret.H"


namespace XFEM
{
  class FluidWizard;
}

namespace DRT
{

  /// Alias dofset that shares dof numbers with another dofset
  /*!
  A special set of degrees of freedom, implemented in order to assign the same degrees of freedom to nodes
  belonging to two discretizations. This way two discretizations can assemble into the same position of the
  system matrix.
  As internal variable it holds a source discretization (Constructor). If such a nodeset is assigned to a
  sub-discretization, its dofs are assigned according to the dofs of the
  source.
  The source discretization can be a xfem discretization. In this case this  should be called with a
  Fluidwizard not equal to zero to determine the  number of xfem dofs.

  */
  class TransparentIndependentDofSet : public virtual IndependentDofSet, public virtual TransparentDofSet
  {
  public:

    /*!
    \brief Standard Constructor

    */

    explicit TransparentIndependentDofSet(RCP<DRT::Discretization> sourcedis,
                                          bool parallel,
                                          Teuchos::RCP<XFEM::FluidWizard> wizard);

    /*!
    \brief Destructor

    */
    virtual ~TransparentIndependentDofSet() { return; }

    /// create a copy of this object
    virtual Teuchos::RCP<DofSet> Clone() { return Teuchos::rcp(new IndependentDofSet(*this)); }

    /// Assign dof numbers to all elements and nodes of the discretization.
    virtual int AssignDegreesOfFreedom(
            const DRT::Discretization& dis,
            const unsigned dspos,
            const int start);

    virtual int NumDofPerNode( const DRT::Node & node, unsigned dspos ) const;


  private:

    Teuchos::RCP<XFEM::FluidWizard> wizard_;

  }; // class TransparentDofSet
} //namespace

#endif  // #ifndef  TRANSPARENTINDEPENDENTDOFSET_H
