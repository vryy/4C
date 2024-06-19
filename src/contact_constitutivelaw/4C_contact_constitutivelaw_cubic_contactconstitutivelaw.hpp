/*----------------------------------------------------------------------*/
/*! \file

\brief implements a Cubic contact constitutive law

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_CONTACT_CONSTITUTIVELAW_CUBIC_CONTACTCONSTITUTIVELAW_HPP
#define FOUR_C_CONTACT_CONSTITUTIVELAW_CUBIC_CONTACTCONSTITUTIVELAW_HPP


#include "4C_config.hpp"

#include "4C_contact_constitutivelaw_contactconstitutivelaw.hpp"
#include "4C_contact_constitutivelaw_contactconstitutivelaw_parameter.hpp"

FOUR_C_NAMESPACE_OPEN

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


      /// create ConstitutiveLaw instance of matching type with my parameters
      Teuchos::RCP<ConstitutiveLaw> create_constitutive_law() override;

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
      Inpar::CONTACT::ConstitutiveLawType get_constitutive_law_type() const override
      {
        return Inpar::CONTACT::ConstitutiveLawType::colaw_cubic;
      }

      double GetA() { return params_->GetA(); }
      double GetB() { return params_->GetB(); }
      double GetC() { return params_->GetC(); }
      double GetD() { return params_->GetD(); }

      /// Return quick accessible constitutive law parameter data
      CONTACT::CONSTITUTIVELAW::Parameter* Parameter() const override { return params_; }

      //@}

      //! @name Evaluation methods

      /// evaluate contact constitutive law
      double evaluate(double gap, CONTACT::Node* cnode) override;
      /// evaluate derivative of the contact constitutive law
      double EvaluateDeriv(double gap, CONTACT::Node* cnode) override;
      //@}

     private:
      /// my material parameters
      CONTACT::CONSTITUTIVELAW::CubicConstitutiveLawParams* params_;
    };
  }  // namespace CONSTITUTIVELAW
}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif
