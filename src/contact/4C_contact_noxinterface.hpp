/*---------------------------------------------------------------------*/
/*! \file
\brief Concrete mplementation of all the %NOX::NLN::CONSTRAINT::Interface::Required
       (pure) virtual routines.

\level 3


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_NOXINTERFACE_HPP
#define FOUR_C_CONTACT_NOXINTERFACE_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_constraint_interface_required.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  class AbstractStrategy;

  class NoxInterface : public NOX::NLN::CONSTRAINT::Interface::Required
  {
   public:
    /// constructor
    NoxInterface();

    /// initialize important member variables
    void Init(const Teuchos::RCP<CONTACT::AbstractStrategy>& strategy_ptr);

    /** \brief Setup important new member variables
     *
     *  Supposed to be overloaded by derived classes. */
    virtual void Setup();

    /// @name Supported basic interface functions
    /// @{
    //! Returns the constraint right-hand-side norms [derived]
    double get_constraint_rhs_norms(const Epetra_Vector& F,
        NOX::NLN::StatusTest::QuantityType checkQuantity, ::NOX::Abstract::Vector::NormType type,
        bool isScaled) const override;

    /// Returns the root mean square (abbr.: RMS) of the Lagrange multiplier updates [derived]
    double get_lagrange_multiplier_update_rms(const Epetra_Vector& xNew, const Epetra_Vector& xOld,
        double aTol, double rTol, NOX::NLN::StatusTest::QuantityType checkQuantity,
        bool disable_implicit_weighting) const override;

    /// Returns the increment norm of the largange multiplier DoFs
    double get_lagrange_multiplier_update_norms(const Epetra_Vector& xNew,
        const Epetra_Vector& xOld, NOX::NLN::StatusTest::QuantityType checkQuantity,
        ::NOX::Abstract::Vector::NormType type, bool isScaled) const override;

    /// Returns the previous solution norm of the largange multiplier DoFs
    double get_previous_lagrange_multiplier_norms(const Epetra_Vector& xOld,
        NOX::NLN::StatusTest::QuantityType checkQuantity, ::NOX::Abstract::Vector::NormType type,
        bool isScaled) const override;

    /// Returns the active set info [derived]
    enum ::NOX::StatusTest::StatusType GetActiveSetInfo(
        enum NOX::NLN::StatusTest::QuantityType checkQuantity, int& activesetsize) const override;

    /// Returns the current active set map
    Teuchos::RCP<const Epetra_Map> get_current_active_set_map(
        enum NOX::NLN::StatusTest::QuantityType checkQuantity) const override;

    /// Returns the old active set map of the previous Newton step
    Teuchos::RCP<const Epetra_Map> GetOldActiveSetMap(
        enum NOX::NLN::StatusTest::QuantityType checkQuantity) const override;
    /// @}

    /// @name Merit function support functions
    /// @{

    double GetModelValue(NOX::NLN::MeritFunction::MeritFctName name) const override;

    double get_linearized_model_terms(const Epetra_Vector& dir,
        const enum NOX::NLN::MeritFunction::MeritFctName name,
        const enum NOX::NLN::MeritFunction::LinOrder linorder,
        const enum NOX::NLN::MeritFunction::LinType lintype) const override;

    /// @}

   protected:
    /// get the init indicator state
    inline const bool& is_init() const { return isinit_; };

    /// get the setup indicator state
    inline const bool& is_setup() const { return issetup_; };

    /// Check if Init() has been called
    inline void check_init() const
    {
      if (not is_init()) FOUR_C_THROW("Call Init() first!");
    };

    /// Check if Init() and Setup() have been called, yet.
    inline void check_init_setup() const
    {
      if (not is_init() or not is_setup()) FOUR_C_THROW("Call Init() and Setup() first!");
    };

    /// Access the underlying strategy
    const CONTACT::AbstractStrategy& Strategy() const
    {
      check_init();
      return *strategy_ptr_;
    };

   protected:
    /// flag indicating if Init() has been called
    bool isinit_;

    /// flag indicating if Setup() has been called
    bool issetup_;

   private:
    Teuchos::RCP<CONTACT::AbstractStrategy> strategy_ptr_;

    std::vector<Teuchos::RCP<Epetra_Map>> cycling_maps_;
  };
}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
