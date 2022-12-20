/*----------------------------------------------------------------------*/
/*! \file

\brief implements a Cubic contact constitutive law

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef Cubic_ConstitutiveLaw_H_
#define Cubic_ConstitutiveLaw_H_


#include "contactconstitutivelaw.H"
#include "contactconstitutivelaw_parameter.H"

namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    /*----------------------------------------------------------------------*/
    /**
     * \brief Contact constitutive parameters for a cubic constitutive law relating the gap to the
     * contact pressure
     *
     * This law has the coefficients \f$ Ax^3+Bx^2+Cx+D \f$
     */
    class CubicConstitutiveLawParams : public Parameter
    {
     public:
      /// standard constructor
      CubicConstitutiveLawParams(
          const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> container);

      /// destructor
      virtual ~CubicConstitutiveLawParams() = default;

      /// create ConstitutiveLaw instance of matching type with my parameters
      Teuchos::RCP<ConstitutiveLaw> CreateConstitutiveLaw() override;

      /// @name get-functions for the parameters of a cubic polynomial
      double GetA() { return a_; }
      double GetB() { return b_; }
      double GetC() { return c_; }
      double GetD() { return d_; }
      //@}

     private:
      /// @name Constitutive Law parameters of a cubic polynomial
      //@{
      const double a_;
      const double b_;
      const double c_;
      const double d_;
      //@}
    };  // class

    /*----------------------------------------------------------------------*/
    /** \brief implements a cubic contact constitutive law \f$ Ax^3+Bx^2+Cx+D \f$ relating the gap
     * to the contact pressure
     *
     */
    class CubicConstitutiveLaw : public ConstitutiveLaw
    {
     public:
      /// construct the ConstitutiveLaw object given a set of parameters
      explicit CubicConstitutiveLaw(CONTACT::CONSTITUTIVELAW::CubicConstitutiveLawParams* params);

      //! @name Access methods

      /// contact constitutive law type
      INPAR::CONTACT::ConstitutiveLawType GetConstitutiveLawType() const override
      {
        return INPAR::CONTACT::ConstitutiveLawType::colaw_cubic;
      }

      double GetA() { return params_->GetA(); }
      double GetB() { return params_->GetB(); }
      double GetC() { return params_->GetC(); }
      double GetD() { return params_->GetD(); }

      /// Return quick accessible constitutive law parameter data
      virtual CONTACT::CONSTITUTIVELAW::Parameter* Parameter() const override { return params_; }

      //@}

      //! @name Evaluation methods

      /// evaluate contact constitutive law
      double Evaluate(double gap) override;
      /// evaluate derivative of the contact constitutive law
      double EvaluateDeriv(double gap) override;
      //@}

     private:
      /// my material parameters
      CONTACT::CONSTITUTIVELAW::CubicConstitutiveLawParams* params_;
    };
  }  // namespace CONSTITUTIVELAW
}  // namespace CONTACT

#endif /*CUBIC_ConstitutiveLaw_H_*/
