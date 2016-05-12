/*!
\file drt_dofset_proxy.H

\brief Proxy to a set of degrees of freedom

<pre>
\level 1

\maintainer Martin Kronbichler
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
