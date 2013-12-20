/*!----------------------------------------------------------------------
\file windkesseldofset.H

\brief A set of degrees of freedom

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Prof. Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Marc Hirschvogel
            hirschvogel@lnm.mw.tum.de
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
