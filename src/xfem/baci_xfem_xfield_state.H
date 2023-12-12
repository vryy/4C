/*---------------------------------------------------------------------------*/
/*! \file

\brief state management for one XFEM discretization

\level 3


*/
/*---------------------------------------------------------------------------*/

#ifndef BACI_XFEM_XFIELD_STATE_H
#define BACI_XFEM_XFIELD_STATE_H

#include "baci_utils_exceptions.H"

#include <Teuchos_RCP.hpp>

class Epetra_Map;

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}  // namespace DRT
namespace CORE::GEO
{
  class CutWizard;
}  // namespace CORE::GEO

namespace XFEM
{
  class ConditionManager;
  class XFEMDofSet;

  class XFieldState
  {
   public:
    /// constructor
    XFieldState();

    /// destructor
    virtual ~XFieldState() = default;

    /** \brief initialize member variables for xfield<-->field couplings
     *
     *  An examples is the XFluidFluid problem. */
    void Init(const Teuchos::RCP<XFEM::ConditionManager>& condition_manager,
        const Teuchos::RCP<CORE::GEO::CutWizard>& wizard,
        const Teuchos::RCP<XFEM::XFEMDofSet>& xdofset,
        const Teuchos::RCP<DRT::Discretization>& xfielddiscret,
        const Teuchos::RCP<DRT::Discretization>& fielddiscret);

    /// setup the stored state objects
    virtual void Setup() = 0;

    /// destroy the stored objects
    virtual bool Destroy() = 0;

    /// transfer the old to a new state object
    virtual void TransferToNewState(
        const DRT::Discretization& new_discret, XFEM::XFieldState& new_xstate) const = 0;

    virtual void ResetNonStandardDofs(const DRT::Discretization& full_discret) = 0;

    virtual void SetNewState(const XFEM::XFieldState& xstate);

    /// @name Accessors
    /// @{
    /// Get cut wizard
    CORE::GEO::CutWizard& CutWizard()
    {
      CheckInit();
      if (wizard_.is_null()) dserror("The CutWizard was not initialized! (Teuchos::null)");
      return *wizard_;
    }

    /// Get condition manager
    XFEM::ConditionManager& ConditionManager()
    {
      CheckInit();
      if (condition_manager_.is_null())
        dserror("The condition_manager was not initialized! (Teuchos::null)");
      return *condition_manager_;
    }

    /// Get dofset of the cut discretization
    XFEM::XFEMDofSet& XDofSet()
    {
      CheckInit();
      if (xdofset_.is_null()) dserror("The xDoF set was not initialized! (Teuchos::null)");
      return *xdofset_;
    }

   protected:
    /// Get cut wizard pointer
    Teuchos::RCP<CORE::GEO::CutWizard>& CutWizardPtr() { return wizard_; }

    /// Get condition manager pointer
    Teuchos::RCP<XFEM::ConditionManager>& ConditionManagerPtr() { return condition_manager_; }

    /// Get pointer to the dofset of the cut discretization
    Teuchos::RCP<XFEM::XFEMDofSet>& XDofSetPtr() { return xdofset_; }

    /// Returns the xFEM field discretizaton
    DRT::Discretization& XFieldDiscret()
    {
      if (xfield_discret_ptr_.is_null()) dserror("xfield_discret_ptr_ is nullptr!");

      return *xfield_discret_ptr_;
    }

    /// Returns a pointer to the xFEM discretization
    Teuchos::RCP<DRT::Discretization>& XFieldDiscretPtr() { return xfield_discret_ptr_; }

    /// Returns the standard field discretizaton
    DRT::Discretization& FieldDiscret()
    {
      if (field_discret_ptr_.is_null()) dserror("field_discret_ptr_ is nullptr!");

      return *field_discret_ptr_;
    }

    /// Returns a pointer to the standard discretization
    Teuchos::RCP<DRT::Discretization>& FieldDiscretPtr() { return field_discret_ptr_; }

    /// @}


   protected:
    //! check the initialization indicator
    inline void CheckInit() const
    {
      if (not isinit_) dserror("Call XFEM::XFieldState::Init() first!");
    }

    //! check the initialization and setup indicators
    inline void CheckInitSetup() const
    {
      if (not issetup_ or not isinit_) dserror("Call Init() and Setup() first!");
    }

   protected:
    /// init indicator
    bool isinit_;

    /// setup indicator
    bool issetup_;

   private:
    /// cut wizard
    Teuchos::RCP<CORE::GEO::CutWizard> wizard_;

    /// condition manager
    Teuchos::RCP<XFEM::ConditionManager> condition_manager_;

    /// XFEM dofset
    Teuchos::RCP<XFEM::XFEMDofSet> xdofset_;

    /// XFEM field discretization pointer
    Teuchos::RCP<DRT::Discretization> xfield_discret_ptr_;

    /// field discretization pointer
    Teuchos::RCP<DRT::Discretization> field_discret_ptr_;
  };  // class XFieldState
}  // namespace XFEM


BACI_NAMESPACE_CLOSE

#endif  // XFEM_XFIELD_STATE_H
