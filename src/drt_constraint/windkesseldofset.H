/*!----------------------------------------------------------------------
\file windkesseldofset.H

\brief A set of degrees of freedom

\level 2

<pre>
\maintainer Marc Hirschvogel
            hirschvogel@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10363
</pre>

*----------------------------------------------------------------------*/
#ifndef WINDKESSELDOFSET_H
#define WINDKESSELDOFSET_H

#include <vector>
#include <list>

#include <Epetra_IntVector.h>
#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>

#include "../drt_lib/drt_dofset_base.H"
#include "../drt_lib/drt_discret.H"

namespace UTILS
{

/*!
\brief A set of degrees of freedom

\note This is an internal class of the Windkessel manager that one
should not need to touch on an ordinary day. It is here to support the
Windkessel manager class. And does all the degree of freedom assignmets
for the Windkessels.

<h3>Purpose</h3>

This class represents one set of degrees of freedom for the
Windkessels in the usual parallel fashion. That is there is a
DofRowMap() and a DofColMap() that return the maps of the global FE
system of equation in row and column setting respectively. These maps
are used by the algorithm's Epetra_Vector classes amoung others.

It is not connected to elements or nodes.
<h3>Invariants</h3>

There are two possible states in this class: Reset and setup. To
change back and forth use AssignDegreesOfFreedom() and Reset().


\author tk
*/
class WindkesselDofSet: public DRT::DofSetBase
{
public:

  /*!
  \brief Standard Constructor

  */
  WindkesselDofSet();

  /*!
  \brief Destructor

  */
  virtual ~WindkesselDofSet();


  //! @name Access methods

  virtual int FirstGID()
  {
    int lmin = dofrowmap_->MinMyGID();
    if(dofrowmap_->NumMyElements() == 0)
      lmin = INT_MAX;
    int gmin = INT_MAX;
    dofrowmap_->Comm().MinAll(&lmin,&gmin,1);
    return gmin;
  };

  //@}

  //! @name Construction

  /// Assign dof numbers using all elements and nodes of the discretization.
  virtual int AssignDegreesOfFreedom(const Teuchos::RCP<DRT::Discretization> dis, const int numdofs, const int start);

  /// reset all internal variables
  virtual void Reset();

  //@}

protected:

}; // class WindkesselDofSet
} // namespace UTILS


// << operator
std::ostream& operator << (std::ostream& os, const UTILS::WindkesselDofSet& dofset);


#endif  // #ifndef DOFSET_H
