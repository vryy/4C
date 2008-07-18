/*!----------------------------------------------------------------------
\file mpcdofset.H

\brief A set of degrees of freedom special for contact

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef MPCDOFSET_H
#define MPCDOFSET_H


#include "../drt_lib/drt_dofset.H"
#include "../drt_lib/drt_discret.H"
using namespace std;


namespace UTILS
{
  /*!
  \brief A set of degrees of freedom special for contact
  */
  class MPCDofSet : public DRT::DofSet
  {
  public:
  
    /*!
    \brief Standard Constructor
  
    */
    explicit MPCDofSet(RCP<DRT::Discretization> sourcedis);
  
    /*!
    \brief Destructor
  
    */
    virtual ~MPCDofSet() { return; }
  
    /// Assign dof numbers to all elements and nodes of the discretization.
    virtual int AssignDegreesOfFreedom(
            const DRT::Discretization& dis,
            const int start);
  
    /// Assign dof numbers for new discretization using dof numbering from source discretization.
    void TransferDegreesOfFreedom(
            const DRT::Discretization& sourcedis,    ///< source discret
            const DRT::Discretization& newdis,       ///< discretization that gets dof numbering from source discret
            const int                  start         ///< offset for dof numbering (obsolete)
            );
  
  private:
    
    RCP<DRT::Discretization> sourcedis_; ///< source discretization
  
  
  }; // class MPCDofSet
} //namespace

#endif  // #ifndef DRT_CDOFSET_H
#endif  // #ifdef CCADISCRET
