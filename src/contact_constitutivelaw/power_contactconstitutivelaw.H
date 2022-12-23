/*----------------------------------------------------------------------*/
/*! \file

\brief implements a simple power law as contact constitutive law

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef POWER_COCONSTLAW_H_
#define POWER_COCONSTLAW_H_


#include "contactconstitutivelaw.H"
#include "contactconstitutivelaw_parameter.H"


namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    /*----------------------------------------------------------------------*/
    /** \brief constitutive law parameters for a power contact law \f$ Ax^B \f$ relating the gap to
     * the contact pressure
     *
     */
    class PowerConstitutiveLawParams : public Parameter
    {
     public:
      /** \brief standard constructor
       * \param[in] container containing the law parameter from the input file
       */
      PowerConstitutiveLawParams(
          const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> container);

      /// destructor
      virtual ~PowerConstitutiveLawParams() = default;

      /// create constitutive law instance of matching type with my parameters
      Teuchos::RCP<ConstitutiveLaw> CreateConstitutiveLaw() override;

      /// @name get-functions for the Constitutive Law parameters of a power law function
      //@{
      /// Get the scaling factor
      double GetA() const { return a_; };
      /// Get the power coefficient
      double GetB() const { return b_; };
      //@}

     private:
      /// @name Constitutive Law parameters of a power function
      //@{
      /// scaling factor
      const double a_;
      /// power coefficient
      const double b_;
      //@}
    };  // class

    /*----------------------------------------------------------------------*/
    /** \brief implements a power contact constitutive law \f$ Ax^B \f$ relating the gap to the
     * contact pressure
     *
     */
    class PowerConstitutiveLaw : public ConstitutiveLaw
    {
     public:
      /// construct the constitutive law object given a set of parameters
      explicit PowerConstitutiveLaw(CONTACT::CONSTITUTIVELAW::PowerConstitutiveLawParams* params);

      //! @name Access methods

      /// contact constitutive law type
      INPAR::CONTACT::ConstitutiveLawType GetConstitutiveLawType() const override
      {
        return INPAR::CONTACT::ConstitutiveLawType::colaw_power;
      }

      /// Get scaling factor of power law
      double GetA() { return params_->GetA(); }
      /// Get power coefficient of power law
      double GetB() { return params_->GetB(); }

      /// Return quick accessible contact constitutive law parameter data
      CONTACT::CONSTITUTIVELAW::Parameter* Parameter() const override { return params_; }

      //! @name Evaluation methods
      //@{
      /// evaluate the constitutive law
      double Evaluate(double gap) override;
      /// Evaluate derivative of the constitutive law
      double EvaluateDeriv(double gap) override;
      //@}

     private:
      /// my constitutive law parameters
      CONTACT::CONSTITUTIVELAW::PowerConstitutiveLawParams* params_;
    };
  }  // namespace CONSTITUTIVELAW
}  // namespace CONTACT

#endif /*POWER_COCONSTLAW_H_*/
