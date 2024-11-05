// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_XFEM_XFIELD_STATE_HPP
#define FOUR_C_XFEM_XFIELD_STATE_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE


namespace Cut
{
  class CutWizard;
}


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
    void init(const std::shared_ptr<XFEM::ConditionManager>& condition_manager,
        const std::shared_ptr<Cut::CutWizard>& wizard,
        const std::shared_ptr<XFEM::XFEMDofSet>& xdofset,
        const std::shared_ptr<Core::FE::Discretization>& xfielddiscret,
        const std::shared_ptr<Core::FE::Discretization>& fielddiscret);

    /// setup the stored state objects
    virtual void setup() = 0;

    /// destroy the stored objects
    virtual bool destroy() = 0;

    /// transfer the old to a new state object
    virtual void transfer_to_new_state(
        const Core::FE::Discretization& new_discret, XFEM::XFieldState& new_xstate) const = 0;

    virtual void reset_non_standard_dofs(const Core::FE::Discretization& full_discret) = 0;

    virtual void set_new_state(const XFEM::XFieldState& xstate);

    /// @name Accessors
    /// @{
    /// Get cut wizard
    Cut::CutWizard& cut_wizard()
    {
      check_init();
      if (!wizard_) FOUR_C_THROW("The CutWizard was not initialized! (nullptr)");
      return *wizard_;
    }

    /// Get condition manager
    XFEM::ConditionManager& condition_manager()
    {
      check_init();
      if (!condition_manager_) FOUR_C_THROW("The condition_manager was not initialized! (nullptr)");
      return *condition_manager_;
    }

    /// Get dofset of the cut discretization
    XFEM::XFEMDofSet& x_dof_set()
    {
      check_init();
      if (!xdofset_) FOUR_C_THROW("The xDoF set was not initialized! (nullptr)");
      return *xdofset_;
    }

   protected:
    /// Get cut wizard pointer
    std::shared_ptr<Cut::CutWizard>& cut_wizard_ptr() { return wizard_; }

    /// Get condition manager pointer
    std::shared_ptr<XFEM::ConditionManager>& condition_manager_ptr() { return condition_manager_; }

    /// Get pointer to the dofset of the cut discretization
    std::shared_ptr<XFEM::XFEMDofSet>& x_dof_set_ptr() { return xdofset_; }

    /// Returns the xFEM field discretizaton
    Core::FE::Discretization& x_field_discret()
    {
      if (!xfield_discret_ptr_) FOUR_C_THROW("xfield_discret_ptr_ is nullptr!");

      return *xfield_discret_ptr_;
    }

    /// Returns a pointer to the xFEM discretization
    std::shared_ptr<Core::FE::Discretization>& x_field_discret_ptr() { return xfield_discret_ptr_; }

    /// Returns the standard field discretizaton
    Core::FE::Discretization& field_discret()
    {
      if (!field_discret_ptr_) FOUR_C_THROW("field_discret_ptr_ is nullptr!");

      return *field_discret_ptr_;
    }

    /// Returns a pointer to the standard discretization
    std::shared_ptr<Core::FE::Discretization>& field_discret_ptr() { return field_discret_ptr_; }

    /// @}


   protected:
    //! check the initialization indicator
    inline void check_init() const
    {
      if (not isinit_) FOUR_C_THROW("Call XFEM::XFieldState::init() first!");
    }

    //! check the initialization and setup indicators
    inline void check_init_setup() const
    {
      if (not issetup_ or not isinit_) FOUR_C_THROW("Call init() and setup() first!");
    }

   protected:
    /// init indicator
    bool isinit_;

    /// setup indicator
    bool issetup_;

   private:
    /// cut wizard
    std::shared_ptr<Cut::CutWizard> wizard_;

    /// condition manager
    std::shared_ptr<XFEM::ConditionManager> condition_manager_;

    /// XFEM dofset
    std::shared_ptr<XFEM::XFEMDofSet> xdofset_;

    /// XFEM field discretization pointer
    std::shared_ptr<Core::FE::Discretization> xfield_discret_ptr_;

    /// field discretization pointer
    std::shared_ptr<Core::FE::Discretization> field_discret_ptr_;
  };  // class XFieldState
}  // namespace XFEM


FOUR_C_NAMESPACE_CLOSE

#endif
