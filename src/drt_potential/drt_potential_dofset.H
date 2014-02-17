/*!----------------------------------------------------------------------
\file potential_dofset.H

\brief A set of degrees of freedom for potential discretization

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifndef POTENTIAL_DOFSET_H_
#define POTENTIAL_DOFSET_H_


#include "../drt_lib/drt_dofset.H"
#include "../drt_lib/drt_discret.H"


namespace POTENTIAL
{
  /*!
  \brief A set of degrees of freedom for potential discretization
  */
  class PotentialDofSet : public DRT::DofSet
  {
  public:

    /*!
    \brief Standard Constructor

    */
    explicit PotentialDofSet(Teuchos::RCP<DRT::Discretization> sourcedis);

    /*!
    \brief Destructor

    */
    virtual ~PotentialDofSet() { return; }

    /// create a copy of this object
    virtual Teuchos::RCP<DofSet> Clone() { dserror("not supported"); return Teuchos::null; }

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

    // no copy constructor and assignment operator
    PotentialDofSet operator = (const PotentialDofSet& old);
    PotentialDofSet(const PotentialDofSet& old);

    Teuchos::RCP<DRT::Discretization> sourcedis_; ///< source discretization


  }; // class PotentialDofSet
} //namespace

#endif  /*POTENTIAL_DOFSET_H_*/






