/*---------------------------------------------------------------------------*/
/*! \file

\brief state management for one XFEM discretization

\level 3


*/
/*---------------------------------------------------------------------------*/

#ifndef FOUR_C_XFEM_XFIELD_STATE_HPP
#define FOUR_C_XFEM_XFIELD_STATE_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Geo
{
  class CutWizard;
}  // namespace Core::Geo

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
        const Teuchos::RCP<Core::Geo::CutWizard>& wizard,
        const Teuchos::RCP<XFEM::XFEMDofSet>& xdofset,
        const Teuchos::RCP<Core::FE::Discretization>& xfielddiscret,
        const Teuchos::RCP<Core::FE::Discretization>& fielddiscret);

    /// setup the stored state objects
    virtual void Setup() = 0;

    /// destroy the stored objects
    virtual bool Destroy() = 0;

    /// transfer the old to a new state object
    virtual void TransferToNewState(
        const Core::FE::Discretization& new_discret, XFEM::XFieldState& new_xstate) const = 0;

    virtual void reset_non_standard_dofs(const Core::FE::Discretization& full_discret) = 0;

    virtual void SetNewState(const XFEM::XFieldState& xstate);

    /// @name Accessors
    /// @{
    /// Get cut wizard
    Core::Geo::CutWizard& CutWizard()
    {
      check_init();
      if (wizard_.is_null()) FOUR_C_THROW("The CutWizard was not initialized! (Teuchos::null)");
      return *wizard_;
    }

    /// Get condition manager
    XFEM::ConditionManager& ConditionManager()
    {
      check_init();
      if (condition_manager_.is_null())
        FOUR_C_THROW("The condition_manager was not initialized! (Teuchos::null)");
      return *condition_manager_;
    }

    /// Get dofset of the cut discretization
    XFEM::XFEMDofSet& XDofSet()
    {
      check_init();
      if (xdofset_.is_null()) FOUR_C_THROW("The xDoF set was not initialized! (Teuchos::null)");
      return *xdofset_;
    }

   protected:
    /// Get cut wizard pointer
    Teuchos::RCP<Core::Geo::CutWizard>& cut_wizard_ptr() { return wizard_; }

    /// Get condition manager pointer
    Teuchos::RCP<XFEM::ConditionManager>& condition_manager_ptr() { return condition_manager_; }

    /// Get pointer to the dofset of the cut discretization
    Teuchos::RCP<XFEM::XFEMDofSet>& x_dof_set_ptr() { return xdofset_; }

    /// Returns the xFEM field discretizaton
    Core::FE::Discretization& x_field_discret()
    {
      if (xfield_discret_ptr_.is_null()) FOUR_C_THROW("xfield_discret_ptr_ is nullptr!");

      return *xfield_discret_ptr_;
    }

    /// Returns a pointer to the xFEM discretization
    Teuchos::RCP<Core::FE::Discretization>& x_field_discret_ptr() { return xfield_discret_ptr_; }

    /// Returns the standard field discretizaton
    Core::FE::Discretization& field_discret()
    {
      if (field_discret_ptr_.is_null()) FOUR_C_THROW("field_discret_ptr_ is nullptr!");

      return *field_discret_ptr_;
    }

    /// Returns a pointer to the standard discretization
    Teuchos::RCP<Core::FE::Discretization>& field_discret_ptr() { return field_discret_ptr_; }

    /// @}


   protected:
    //! check the initialization indicator
    inline void check_init() const
    {
      if (not isinit_) FOUR_C_THROW("Call XFEM::XFieldState::Init() first!");
    }

    //! check the initialization and setup indicators
    inline void check_init_setup() const
    {
      if (not issetup_ or not isinit_) FOUR_C_THROW("Call Init() and Setup() first!");
    }

   protected:
    /// init indicator
    bool isinit_;

    /// setup indicator
    bool issetup_;

   private:
    /// cut wizard
    Teuchos::RCP<Core::Geo::CutWizard> wizard_;

    /// condition manager
    Teuchos::RCP<XFEM::ConditionManager> condition_manager_;

    /// XFEM dofset
    Teuchos::RCP<XFEM::XFEMDofSet> xdofset_;

    /// XFEM field discretization pointer
    Teuchos::RCP<Core::FE::Discretization> xfield_discret_ptr_;

    /// field discretization pointer
    Teuchos::RCP<Core::FE::Discretization> field_discret_ptr_;
  };  // class XFieldState
}  // namespace XFEM


FOUR_C_NAMESPACE_CLOSE

#endif
