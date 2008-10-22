/*!----------------------------------------------------------------------
\file drt_cdofset.H

\brief A set of degrees of freedom special for contact

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich
              
Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is propietary software. It must not be published, distributed, 
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de) 
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de                   

-------------------------------------------------------------------------
<\pre>

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DRT_CDOFSET_H
#define DRT_CDOFSET_H


#include "../drt_lib/drt_dofset.H"
#include "../drt_lib/drt_discret.H"
using namespace std;

/*!
\brief CONTACT: namespace of the baci contact algorithms

*/
namespace CONTACT
{
/*!
\brief A set of degrees of freedom special for contact

*/
class CDofSet : public DRT::DofSet
{
public:

  //! @name Constructors and destructors and related methods
  
  /*!
  \brief Standard Constructor

  */
  explicit CDofSet();

  /*!
  \brief Destructor

  */
  virtual ~CDofSet() { return; }

  //@}

  //! @name Construction

  /// Assign dof numbers to all elements and nodes of the discretization.
  virtual int AssignDegreesOfFreedom(const DRT::Discretization& dis, const int start);

  //@}

protected:


}; // class CDofSet
} // namespace CONTACT




#endif  // #ifndef DRT_CDOFSET_H
#endif  // #ifdef CCADISCRET
