/*----------------------------------------------------------------------*/
/*! \file

\brief implements a simple power law as contact constitutive law

\level 1

\maintainer Nora Hagmeyer
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

    /** \brief DRT::ParObjectType for the PowerConstitutiveLaw \f$ Ax^B \f$ relating the gap to the
     * contact pressure
     *
     */
    class PowerConstitutiveLawType : public ConstitutiveLawType
    {
     public:
      /// return name of the contact constitutive law
      std::string Name() const { return "PowerConstitutiveLaw"; }

      static PowerConstitutiveLawType& Instance() { return instance_; };

     private:
      static PowerConstitutiveLawType instance_;
    };

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

      /**
       *  \brief Return unique ParObject id
       *
       * every class implementing ParObject needs a unique id defined at the
       *  top of drt_parobject.H (this file) and should return it in this method.
       */
      int UniqueParObjectId() const override
      {
        return PowerConstitutiveLawType::Instance().UniqueParObjectId();
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
        return INPAR::CONTACT::ConstitutiveLawType::colaw_power;
      }

      /// return copy of this contact constitutive law object
      Teuchos::RCP<ConstitutiveLaw> Clone() const override
      {
        return Teuchos::rcp(new PowerConstitutiveLaw(*this));
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
