/*---------------------------------------------------------------------*/
/*! \file
\brief Concrete mplementation of all the %NOX::NLN::CONSTRAINT::Interface::Required
       (pure) virtual routines.

\level 3


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_CONTACT_NOXINTERFACE_HPP
#define BACI_CONTACT_NOXINTERFACE_HPP

#include "baci_config.hpp"

#include "baci_solver_nonlin_nox_constraint_interface_required.hpp"

#include <vector>

BACI_NAMESPACE_OPEN

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
    double GetConstraintRHSNorms(const Epetra_Vector& F,
        NOX::NLN::StatusTest::QuantityType checkQuantity, ::NOX::Abstract::Vector::NormType type,
        bool isScaled) const override;

    /// Returns the root mean square (abbr.: RMS) of the Lagrange multiplier updates [derived]
    double GetLagrangeMultiplierUpdateRMS(const Epetra_Vector& xNew, const Epetra_Vector& xOld,
        double aTol, double rTol, NOX::NLN::StatusTest::QuantityType checkQuantity,
        bool disable_implicit_weighting) const override;

    /// Returns the increment norm of the largange multiplier DoFs
    double GetLagrangeMultiplierUpdateNorms(const Epetra_Vector& xNew, const Epetra_Vector& xOld,
        NOX::NLN::StatusTest::QuantityType checkQuantity, ::NOX::Abstract::Vector::NormType type,
        bool isScaled) const override;

    /// Returns the previous solution norm of the largange multiplier DoFs
    double GetPreviousLagrangeMultiplierNorms(const Epetra_Vector& xOld,
        NOX::NLN::StatusTest::QuantityType checkQuantity, ::NOX::Abstract::Vector::NormType type,
        bool isScaled) const override;

    /// Returns the active set info [derived]
    enum ::NOX::StatusTest::StatusType GetActiveSetInfo(
        enum NOX::NLN::StatusTest::QuantityType checkQuantity, int& activesetsize) const override;

    /// Returns the current active set map
    Teuchos::RCP<const Epetra_Map> GetCurrentActiveSetMap(
        enum NOX::NLN::StatusTest::QuantityType checkQuantity) const override;

    /// Returns the old active set map of the previous Newton step
    Teuchos::RCP<const Epetra_Map> GetOldActiveSetMap(
        enum NOX::NLN::StatusTest::QuantityType checkQuantity) const override;
    /// @}

    /// @name Merit function support functions
    /// @{

    double GetModelValue(NOX::NLN::MeritFunction::MeritFctName name) const override;

    double GetLinearizedModelTerms(const Epetra_Vector& dir,
        const enum NOX::NLN::MeritFunction::MeritFctName name,
        const enum NOX::NLN::MeritFunction::LinOrder linorder,
        const enum NOX::NLN::MeritFunction::LinType lintype) const override;

    /// @}

   protected:
    /// get the init indicator state
    inline const bool& IsInit() const { return isinit_; };

    /// get the setup indicator state
    inline const bool& IsSetup() const { return issetup_; };

    /// Check if Init() has been called
    inline void CheckInit() const
    {
      if (not IsInit()) dserror("Call Init() first!");
    };

    /// Check if Init() and Setup() have been called, yet.
    inline void CheckInitSetup() const
    {
      if (not IsInit() or not IsSetup()) dserror("Call Init() and Setup() first!");
    };

    /// Access the underlying strategy
    const CONTACT::AbstractStrategy& Strategy() const
    {
      CheckInit();
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


BACI_NAMESPACE_CLOSE

#endif  // CONTACT_NOXINTERFACE_H
