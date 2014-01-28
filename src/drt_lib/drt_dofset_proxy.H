/*!
\file drt_dofset_proxy.H

\brief Proxy to a set of degrees of freedom

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

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>

*/

#ifndef DRT_DOFSET_PROXY_H
#define DRT_DOFSET_PROXY_H

#include "drt_dofset.H"


namespace DRT
{


/// Proxy to a DofSet that does not own dofs itself
/*!

  As the Discretization handles DofSets a little implicit, a DofSetProxy is
  needed to change the DofSet behaviour. We need a DofSet that

  - returns dof numbers owned by a different DofSet
  - does not reset those dof numbers
  - does not assign dof numbers itself

  The DofSetProxy is meant to be used as secondary DofSet in a Discretization
  if there are two volume coupled Discretizations with matching nodes. Think
  of Structure-Thermo coupling.

  \author u.kue
  \date 12/09
 */
class DofSetProxy : public DofSet
{
public:

  explicit DofSetProxy(DofSet* dofset);

  virtual ~DofSetProxy();

  /// create a copy of this object
  virtual Teuchos::RCP<DofSet> Clone() { return Teuchos::rcp(new DofSetProxy(*this)); }

  //! @name Construction

  /// Add Dof Set to list #static_dofsets_
  virtual void AddDofSettoList();

  /// Assign dof numbers using all elements and nodes of the discretization.
  virtual int AssignDegreesOfFreedom(const Discretization& dis, const unsigned dspos, const int start);

  virtual bool Filled() const;

  //@}

  //! @name Communication
  /// The original DofSet sends notifications if it changes.

  /// original DofSet has new dofs
  void NotifyAssigned();

  /// original DofSet has been reset
  virtual void NotifyReset();

  /// our original DofSet dies
  void Disconnect(DofSet* dofset);

  //@}

private:

  DofSet* dofset_;
};

}


#endif
