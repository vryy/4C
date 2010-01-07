/*!----------------------------------------------------------------------
\file mortar_dofset.H

\brief A set of degrees of freedom special for mortar coupling

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef MORTAR_DOFSET_H
#define MORTAR_DOFSET_H

#include "../drt_lib/drt_dofset.H"
#include "../drt_lib/drt_discret.H"

using namespace std;

/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
namespace MORTAR
{
/*!
\brief A set of degrees of freedom special for mortar coupling

*/
class MortarDofSet : public DRT::DofSet
{
public:

  //! @name Constructors and destructors and related methods

  /*!
  \brief Standard Constructor

  */
  explicit MortarDofSet();

  /*!
  \brief Destructor

  */
  virtual ~MortarDofSet() { return; }

  /// create a copy of this object
  virtual Teuchos::RCP<DofSet> Clone() { return Teuchos::rcp(new MortarDofSet(*this)); }

  //@}

  //! @name Construction

  /// Assign dof numbers to all elements and nodes of the discretization.
  virtual int AssignDegreesOfFreedom(const DRT::Discretization& dis, const unsigned dspos, const int start);

  //@}

protected:


}; // class MortarDofSet
} // namespace MORTAR

#endif  // #ifndef MORTAR_DOFSET_H
#endif  // #ifdef CCADISCRET
