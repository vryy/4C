/*---------------------------------------------------------------------*/
/*! \file

\brief Concrete mplementation of all the %NOX::NLN::CONSTRAINT::Interface::Required
       (pure) virtual routines.

\level 3


\date July 29, 2016

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_CONSTRAINT_LAGPENCONSTRAINT_NOXINTERFACE_HPP
#define FOUR_C_CONSTRAINT_LAGPENCONSTRAINT_NOXINTERFACE_HPP


#include "baci_config.hpp"

#include "baci_solver_nonlin_nox_constraint_interface_preconditioner.hpp"
#include "baci_solver_nonlin_nox_constraint_interface_required.hpp"
#include "baci_structure_new_timint_basedataglobalstate.hpp"

BACI_NAMESPACE_OPEN


namespace LAGPENCONSTRAINT
{
  class NoxInterface : public NOX::NLN::CONSTRAINT::Interface::Required
  {
   public:
    /// constructor
    NoxInterface();

    /// initialize important member variables
    void Init(const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr);

    /** \brief Setup important new member variables
     *
     *  Supposed to be overloaded by derived classes. */
    virtual void Setup();

    /// @name Supported basic interface functions
    /// @{
    //! Returns the constraint right-hand-side norms [derived]
    double GetConstraintRHSNorms(const Epetra_Vector& F, NOX::NLN::StatusTest::QuantityType chQ,
        ::NOX::Abstract::Vector::NormType type, bool isScaled) const override;

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


   protected:
    /// flag indicating if Init() has been called
    bool isinit_;

    /// flag indicating if Setup() has been called
    bool issetup_;

   private:
    //! global state data container
    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate_ptr_;
  };


  class NoxInterfacePrec : public NOX::NLN::CONSTRAINT::Interface::Preconditioner
  {
   public:
    /// constructor
    NoxInterfacePrec();

    /// initialize important member variables
    void Init(const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr);

    /** \brief Setup important new member variables
     *
     *  Supposed to be overloaded by derived classes. */
    virtual void Setup();


    bool IsSaddlePointSystem() const override;

    bool IsCondensedSystem() const override;

    void FillMapsForPreconditioner(std::vector<Teuchos::RCP<Epetra_Map>>& maps) const override;

    bool computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M,
        Teuchos::ParameterList* precParams = nullptr) override;
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

   protected:
    /// flag indicating if Init() has been called
    bool isinit_;

    /// flag indicating if Setup() has been called
    bool issetup_;

   private:
    //! global state data container
    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate_ptr_;
  };

}  // namespace LAGPENCONSTRAINT


BACI_NAMESPACE_CLOSE

#endif  // CONSTRAINT_LAGPENCONSTRAINT_NOXINTERFACE_H
