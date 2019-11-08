/*----------------------------------------------------------------------*/
/*! \file

\brief implements a Cubic contact constitutive law

\level 1

\maintainer Nora Hagmeyer
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
     * \brief contact constitutive parameters for a cubic constitutive law
     *
     * This law has the coefficients Ax^3+Bx^2+Cx+D
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

    /// For every DRT::ParObject a corresponding DRT::ParObjectType has to be implemented
    class CubicConstitutiveLawType : public ConstitutiveLawType
    {
     public:
      std::string Name() const { return "CubicConstitutiveLaw"; }

      static CubicConstitutiveLawType& Instance() { return instance_; };

     private:
      static CubicConstitutiveLawType instance_;
    };

    /*----------------------------------------------------------------------*/
    /// Wrapper for cubic contact constitutive law
    class CubicConstitutiveLaw : public ConstitutiveLaw
    {
     public:
      /// construct the ConstitutiveLaw object given the parameters
      explicit CubicConstitutiveLaw(CONTACT::CONSTITUTIVELAW::CubicConstitutiveLawParams* params);

      //! @name Packing and Unpacking

      /*!
        \brief Return unique ParObject id

        every class implementing ParObject needs a unique id defined at the
        top of drt_parobject.H (this file) and should return it in this method.
      */
      int UniqueParObjectId() const override
      {
        return CubicConstitutiveLawType::Instance().UniqueParObjectId();
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
       *  \todo  Fully implement this class
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
       *  \todo Fully implement this class
       */
      void Unpack(const std::vector<char>& data) override;

      //@}

      //! @name Access methods

      /// contact constitutive law type
      INPAR::CONTACT::ConstitutiveLawType GetConstitutiveLawType() const override
      {
        return INPAR::CONTACT::ConstitutiveLawType::colaw_cubic;
      }

      /// return copy of this ConstitutiveLaw object
      Teuchos::RCP<ConstitutiveLaw> Clone() const override
      {
        return Teuchos::rcp(new CubicConstitutiveLaw(*this));
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
