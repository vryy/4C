/*!----------------------------------------------------------------------
\file drt_dofset_base.H

\brief Base class implementing common functionality and dofset registration

<pre>
\level 0

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef DOFSETBASE_H
#define DOFSETBASE_H

#include <list>

#include "drt_dofset_interface.H"

namespace DRT
{

class Discretization;


/*! \brief Base class set of degrees of freedom

  This base class manages the static list, all DofSets are
  written into and the list of DofSets the current DofSet is connected to.

  \author tk    */
class DofSetBase : public DofSetInterface
{
public:

  /// Standard Constructor
  DofSetBase();

  /// Destructor
  virtual ~DofSetBase();


  //! @name Utility Methods

  /// Print this class
  virtual void Print(std::ostream& os) const;

  /// Get Max of all GID assigned in the DofSets in front of current one in the list #static_dofsets_
  virtual int MaxGIDinList(const Epetra_Comm& comm) const;

  //@}

  //! @name Construction

  /// Add Dof Set to list #static_dofsets_
  virtual void AddDofSettoList();

  /// Replace a Dof Set in list #static_dofsets_ with this
  virtual void ReplaceInStaticDofsets(Teuchos::RCP<DofSetInterface> olddofset);

  //@}

  /// Print the dofsets in the static_dofsets_ list
  virtual void PrintAllDofsets(const Epetra_Comm& comm) const;


  //! @name DofSet management
  /// Registered DofSets need to know about changes to the DofSet.

  /// Notify proxies of new dofs
  virtual void NotifyAssigned();

  /// Notify proxies of reset
  virtual void NotifyReset();

  /// Register new dofset to notify
  virtual void Register(DofSetInterface* dofset);

  /// Remove dofset from list
  virtual void Unregister(DofSetInterface* dofset);

  //@}

private:

  /*! \brief list of registered dof sets

    Whenever you request a proxy of any specific DofSet implementation, this proxy
    has to be registered in the list registered_dofsets_ . See also \ref DofSetProxy.
    Also other, special implementations of a DofSet, which are linked to another
    DofSet in some way may register here, see e.g. \ref DofSetDefinedMappingWrapper and
    \ref DofSetGIDBasedWrapper. This way the registered DofSet will be notified
    of any state changes of the original DofSet by calling \ref NotifyAssigned and \ref NotifyReset .*/
  std::list<DofSetInterface*> registered_dofsets_;

  /*! \brief store dofset in static list, if derived class chooses so using AddDofSettoList()

      This is hack to get unique dof numbers on all dof sets. With this in place we
      can combine the maps from any dof sets to form block systems and the like.    */
  static std::list<DofSetInterface*> static_dofsets_;


}; // class DofSetBase
} // namespace DRT

#endif  // #ifndef DofSetBase_H
