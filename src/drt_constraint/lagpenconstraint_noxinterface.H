/*---------------------------------------------------------------------*/
/*!

\brief Concrete mplementation of all the %NOX::NLN::CONSTRAINT::Interface::Required
       (pure) virtual routines.

\level 3

\maintainer Amadeus Gebauer

\date July 29, 2016

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_LAGPENCONSTRAINT_LAGPENCONSTRAINT_NOXINTERFACE_H_
#define SRC_DRT_LAGPENCONSTRAINT_LAGPENCONSTRAINT_NOXINTERFACE_H_


#include "../solver_nonlin_nox/nox_nln_constraint_interface_required.H"
#include "../solver_nonlin_nox/nox_nln_constraint_interface_preconditioner.H"

#include "../drt_structure_new/str_timint_basedataglobalstate.H"


namespace LAGPENCONSTRAINT
{
  class NoxInterface : public NOX::NLN::CONSTRAINT::Interface::Required
  {
   public:
    /// constructor
    NoxInterface();

    /// destructor
    virtual ~NoxInterface(){};

    /// initialize important member variables
    void Init(const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr);

    /** \brief Setup important new member variables
     *
     *  Supposed to be overloaded by derived classes. */
    virtual void Setup();

    /// @name Supported basic interface functions
    /// @{
    //! Returns the constraint right-hand-side norms [derived]
    virtual double GetConstraintRHSNorms(const Epetra_Vector& F,
        NOX::NLN::StatusTest::QuantityType chQ, NOX::Abstract::Vector::NormType type,
        bool isScaled) const;

    /// Returns the root mean square (abbr.: RMS) of the Lagrange multiplier updates [derived]
    virtual double GetLagrangeMultiplierUpdateRMS(const Epetra_Vector& xNew,
        const Epetra_Vector& xOld, double aTol, double rTol,
        NOX::NLN::StatusTest::QuantityType checkQuantity, bool disable_implicit_weighting) const;

    /// Returns the increment norm of the largange multiplier DoFs
    virtual double GetLagrangeMultiplierUpdateNorms(const Epetra_Vector& xNew,
        const Epetra_Vector& xOld, NOX::NLN::StatusTest::QuantityType checkQuantity,
        NOX::Abstract::Vector::NormType type, bool isScaled) const;

    /// Returns the previous solution norm of the largange multiplier DoFs
    virtual double GetPreviousLagrangeMultiplierNorms(const Epetra_Vector& xOld,
        NOX::NLN::StatusTest::QuantityType checkQuantity, NOX::Abstract::Vector::NormType type,
        bool isScaled) const;
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

    /// destructor
    virtual ~NoxInterfacePrec(){};

    /// initialize important member variables
    void Init(const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr);

    /** \brief Setup important new member variables
     *
     *  Supposed to be overloaded by derived classes. */
    virtual void Setup();


    virtual bool IsSaddlePointSystem() const;

    virtual bool IsCondensedSystem() const;

    virtual void FillMapsForPreconditioner(std::vector<Teuchos::RCP<Epetra_Map>>& maps) const;

    virtual bool computePreconditioner(
        const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams = 0);
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


#endif /* SRC_DRT_LAGPENCONSTRAINT_LAGPENCONSTRAINT_NOXINTERFACE_H_ */
