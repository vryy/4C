/*---------------------------------------------------------------------*/
/*! \file

\brief Base class implementing common functionality and dofset registration

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_LIB_DOFSET_BASE_HPP
#define FOUR_C_LIB_DOFSET_BASE_HPP

#include "4C_config.hpp"

#include "4C_lib_dofset_interface.hpp"

#include <list>

FOUR_C_NAMESPACE_OPEN

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
    ~DofSetBase() override;


    //! @name Utility Methods

    /// Print this class
    void Print(std::ostream& os) const override;

    /// Get Max of all GID assigned in the DofSets in front of current one in the list
    /// #static_dofsets_
    int MaxGIDinList(const Epetra_Comm& comm) const override;

    //@}

    //! @name Construction

    /// Add Dof Set to list #static_dofsets_
    void AddDofSettoList() override;

    /// Replace a Dof Set in list #static_dofsets_ with this
    void ReplaceInStaticDofsets(Teuchos::RCP<DofSetInterface> olddofset) override;

    //@}

    /// Print the dofsets in the static_dofsets_ list
    void PrintAllDofsets(const Epetra_Comm& comm) const override;


    //! @name DofSet management
    /// Registered DofSets need to know about changes to the DofSet.

    /// Notify proxies of new dofs
    void NotifyAssigned() override;

    /// Notify proxies of reset
    void NotifyReset() override;

    /// Register new dofset to notify
    void Register(DofSetInterface* dofset) override;

    /// Remove dofset from list
    void Unregister(DofSetInterface* dofset) override;

    //@}

   private:
    /*! \brief list of registered dof sets

      Whenever you request a proxy of any specific DofSet implementation, this proxy
      has to be registered in the list registered_dofsets_ . See also \ref DofSetProxy.
      Also other, special implementations of a DofSet, which are linked to another
      DofSet in some way may register here, see e.g. \ref DofSetDefinedMappingWrapper and
      \ref DofSetGIDBasedWrapper. This way the registered DofSet will be notified
      of any state changes of the original DofSet by calling \ref NotifyAssigned and \ref
      NotifyReset .*/
    std::list<DofSetInterface*> registered_dofsets_;

    /*! \brief store dofset in static list, if derived class chooses so using AddDofSettoList()

        This is hack to get unique dof numbers on all dof sets. With this in place we
        can combine the maps from any dof sets to form block systems and the like.    */
    static std::list<DofSetInterface*> static_dofsets_;


  };  // class DofSetBase
}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
