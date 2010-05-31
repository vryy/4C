/*!----------------------------------------------------------------------
\file drt_dofset_transparent.H

\brief A set of degrees of freedom special for contact

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef TRANSPARENTDOFSET_H
#define TRANSPARENTDOFSET_H


#include "../drt_lib/drt_dofset.H"
#include "../drt_lib/drt_discret.H"
using namespace std;


namespace DRT
{
  /*!
  A special set of degrees of freedom, implemented in order to assign the same degrees of freedom to nodes 
  belonging to two discretizations. This way two discretizations can assemble into the same position of the 
  system matrix. 
  As internal variable it holds a source discretization (Constructor). If such a nodeset is assigned to a 
  sub-discretization, its dofs are assigned according to the dofs of the source.
  
  */
  class TransparentDofSet : public DRT::DofSet
  {
  public:

    /*!
    \brief Standard Constructor

    */
    explicit TransparentDofSet(RCP<DRT::Discretization> sourcedis);

    /*!
    \brief Destructor

    */
    virtual ~TransparentDofSet() { return; }

    /// create a copy of this object
    virtual Teuchos::RCP<DofSet> Clone() { return Teuchos::rcp(new TransparentDofSet(*this)); }

    /// Assign dof numbers to all elements and nodes of the discretization.
    virtual int AssignDegreesOfFreedom(
            const DRT::Discretization& dis,
            const unsigned dspos,
            const int start);

    /// Assign dof numbers for new discretization using dof numbering from source discretization.
    void TransferDegreesOfFreedom(
            const DRT::Discretization& sourcedis,    ///< source discret
            const DRT::Discretization& newdis,       ///< discretization that gets dof numbering from source discret
            const int                  start         ///< offset for dof numbering (obsolete)
            );

  private:

    RCP<DRT::Discretization> sourcedis_; ///< source discretization


  }; // class TransparentDofSet
} //namespace

#endif  // #ifndef DRT_TRANSPARENDOFSET_H
#endif  // #ifdef CCADISCRET
