/*----------------------------------------------------------------------*/
/*! \file
\brief Implements a default linear penalty function as contact constitutive law

\level 3

\maintainer Nora Hagmeyer
*/
/*----------------------------------------------------------------------*/
#ifndef LINEAR_COCONSTLAW_H_
#define LINEAR_COCONSTLAW_H_


#include "coconstlaw_parameter.H"
#include "coconstitutivelaw.H"


namespace CONTACT
{
  namespace CONSTITUTIVELAW
  {
    /*----------------------------------------------------------------------*/
    /// constitutive law parameters for a linear contact law Ax+B
    class LinearConstitutiveLawParams : public Parameter
    {
     public:
      /** \brief standard constructor
       * \params[in] container containing the law parameter from the input file
       */
      LinearConstitutiveLawParams(
          const Teuchos::RCP<const CONTACT::CONSTITUTIVELAW::Container> container);

      /// destructor
      virtual ~LinearConstitutiveLawParams() = default;

      /// create constitutive law instance of matching type with my parameters
      Teuchos::RCP<CONTACT::CONSTITUTIVELAW::ConstitutiveLaw> CreateConstitutiveLaw() override;

      /// @name get-functions for the Constitutive Law parameters of a broken rational function
      //@{
      /// Get the slope
      double GetA() const { return a_; };
      /// Get the y intercept
      double GetB() const { return b_; };
      //@}

     private:
      /// @name Constitutive Law parameters of a linear function
      //@{
      /// slope
      double a_;
      /// y intercept
      const double b_;
      //@}
    };  // class

    /// For every DRT::ParObject a corresponding DRT::ParObjectType has to be implemented
    class LinearConstitutiveLawType : public ConstitutiveLawType
    {
     public:
      /// return name of the contact constitutive law
      std::string Name() const override { return "LinearConstitutiveLaw"; }

      static LinearConstitutiveLawType& Instance() { return instance_; };

     private:
      static LinearConstitutiveLawType instance_;
    };

    /*----------------------------------------------------------------------*/
    /// Wrapper for linear contact constitutive law
    class LinearConstitutiveLaw : public ConstitutiveLaw
    {
     public:
      /// construct the constitutive law object given parameters
      explicit LinearConstitutiveLaw(CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams* params);

      /**
       *  \brief Return unique ParObject id
       *
       * every class implementing ParObject needs a unique id defined at the
       *  top of drt_parobject.H (this file) and should return it in this method.
       */
      int UniqueParObjectId() const override
      {
        return LinearConstitutiveLawType::Instance().UniqueParObjectId();
      }

      /**
       *  \brief Pack this class so it can be communicated
       *
       *  Resizes the vector data and stores all information of a class in it.
       *  The first information to be stored in data has to be the
       *  unique parobject id delivered by UniqueParObjectId() which will then
       *  identify the exact class on the receiving processor.
       *
       *  \param data (in/out): char vector to store class information
       */
      void Pack(DRT::PackBuffer& data) const override;

      /**
       *  \brief Unpack data from a char vector into this class
       *
       *  The vector data contains all information to rebuild the
       *  exact copy of an instance of a class on a different processor.
       *  The first entry in data has to be an integer which is the unique
       *  parobject id defined at the top of this file and delivered by
       *  UniqueParObjectId().
       *
       *  \param data (in) : vector storing all data to be unpacked into this
       *  instance.
       */
      void Unpack(const std::vector<char>& data) override;

      //! @name Access methods

      /// contact constitutive law type
      INPAR::CONTACT::ConstitutiveLawType GetConstitutiveLawType() const override
      {
        return INPAR::CONTACT::ConstitutiveLawType::colaw_linear;
      }

      /// return copy of this contact constitutive law object
      Teuchos::RCP<ConstitutiveLaw> Clone() const override
      {
        return Teuchos::rcp(new LinearConstitutiveLaw(*this));
      }

      /// Get slope of linear polynomial
      double GetA() { return params_->GetA(); }
      /// Get y intercept of linear polynomial
      double GetB() { return params_->GetB(); }

      /// Return quick accessible contact constitutive law parameter data
      virtual CONTACT::CONSTITUTIVELAW::Parameter* Parameter() const { return params_; }

      //! @name Evaluation methods
      //@{
      /// evaluate the constitutive law
      double Evaluate(double gap) override;

      /// Evaluate derivative of the constitutive law
      double EvaluateDeriv(double gap) override;
      //@}

     private:
      /// my constitutive law parameters
      CONTACT::CONSTITUTIVELAW::LinearConstitutiveLawParams* params_;
    };
  }  // namespace CONSTITUTIVELAW
}  // namespace CONTACT

#endif /*LINEAR_COCONSTLAW_H_*/
