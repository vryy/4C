/*---------------------------------------------------------------------*/
/*! \file
\brief Concrete mplementation of all the %NOX::Nln::CONSTRAINT::Interface::Required
       (pure) virtual routines.

\level 3


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_MESHTYING_NOXINTERFACE_HPP
#define FOUR_C_CONTACT_MESHTYING_NOXINTERFACE_HPP


#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_constraint_interface_preconditioner.hpp"
#include "4C_solver_nonlin_nox_constraint_interface_required.hpp"
#include "4C_structure_new_timint_basedataglobalstate.hpp"

FOUR_C_NAMESPACE_OPEN


namespace CONTACT
{
  class MtNoxInterface : public NOX::Nln::CONSTRAINT::Interface::Required
  {
   public:
    /// constructor
    MtNoxInterface();

    /// initialize important member variables
    void Init(const Teuchos::RCP<STR::TimeInt::BaseDataGlobalState>& gstate_ptr);

    /** \brief Setup important new member variables
     *
     *  Supposed to be overloaded by derived classes. */
    virtual void Setup();

    /// @name Supported basic interface functions
    /// @{
    //! Returns the constraint right-hand-side norms [derived]
    double get_constraint_rhs_norms(const Epetra_Vector& F,
        NOX::Nln::StatusTest::QuantityType checkQuantity, ::NOX::Abstract::Vector::NormType type,
        bool isScaled) const override;

    /// Returns the root mean square (abbr.: RMS) of the Lagrange multiplier updates [derived]
    double get_lagrange_multiplier_update_rms(const Epetra_Vector& xNew, const Epetra_Vector& xOld,
        double aTol, double rTol, NOX::Nln::StatusTest::QuantityType checkQuantity,
        bool disable_implicit_weighting) const override;

    /// Returns the increment norm of the largange multiplier DoFs
    double get_lagrange_multiplier_update_norms(const Epetra_Vector& xNew,
        const Epetra_Vector& xOld, NOX::Nln::StatusTest::QuantityType checkQuantity,
        ::NOX::Abstract::Vector::NormType type, bool isScaled) const override;

    /// Returns the previous solution norm of the largange multiplier DoFs
    double get_previous_lagrange_multiplier_norms(const Epetra_Vector& xOld,
        NOX::Nln::StatusTest::QuantityType checkQuantity, ::NOX::Abstract::Vector::NormType type,
        bool isScaled) const override;
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


   protected:
    /// flag indicating if Init() has been called
    bool isinit_;

    /// flag indicating if Setup() has been called
    bool issetup_;

   private:
    //! global state data container
    Teuchos::RCP<STR::TimeInt::BaseDataGlobalState> gstate_ptr_;
  };

}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
