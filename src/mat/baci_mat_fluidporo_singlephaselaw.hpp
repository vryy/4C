/*----------------------------------------------------------------------*/
/*! \file
 \brief a material defining the pressure-saturation relationship of
        fluid phase within a multiphase porous fluid

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_FLUIDPORO_SINGLEPHASELAW_HPP
#define FOUR_C_MAT_FLUIDPORO_SINGLEPHASELAW_HPP

#include "baci_config.hpp"

#include "baci_mat_par_parameter.hpp"
#include "baci_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  namespace UTILS
  {
    class FunctionOfAnything;
  }  // namespace UTILS
}  // namespace DRT

namespace MAT
{
  namespace PAR
  {
    //! interface class for generic phase (pressure-saturation) law
    class FluidPoroPhaseLaw : public Parameter
    {
     public:
      /// standard constructor
      explicit FluidPoroPhaseLaw(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! build the phase law
      static FluidPoroPhaseLaw* CreatePhaseLaw(int phaselawId);

      /// create material instance of matching type with my parameters
      Teuchos::RCP<MAT::Material> CreateMaterial() override { return Teuchos::null; };

      /// initialize
      virtual void Initialize() = 0;

      /// pressure ids involved in phase law
      virtual const std::vector<int>* PresIds() = 0;

      //! evaluate saturation of phase
      virtual double EvaluateSaturation(const std::vector<double>& pressure) = 0;

      //! evaluate derivative of saturation w.r.t. pressure
      virtual double EvaluateDerivOfSaturationWrtPressure(
          int doftoderive, const std::vector<double>& pressure) = 0;

      //! evaluate second derivative of saturation w.r.t. pressure
      virtual double EvaluateSecondDerivOfSaturationWrtPressure(
          int firstdoftoderive, int seconddoftoderive, const std::vector<double>& pressure) = 0;

      //! evaluate generalized pressure related to phase law
      virtual double EvaluateGenPressure(double saturation) = 0;

      //! evaluate derivative of pressure w.r.t. saturation
      virtual double EvaluateDerivOfPressureWrtSaturation(int doftoderive, double saturation) = 0;
    };

    /*----------------------------------------------------------------------*/
    //! a linear pressure-saturation relation. Only for testing, no physical meaning.
    class FluidPoroPhaseLawLinear : public FluidPoroPhaseLaw
    {
     public:
      /// standard constructor
      explicit FluidPoroPhaseLawLinear(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// initialize
      void Initialize() override { return; };

      /// @name material parameters
      //@{
      //! number of degrees of freedom
      const int numdof_;
      //! coefficients defining the pressures involved int the pressure-saturation law
      const std::vector<int>* presids_;
      //! 'relative tension' parameter (proportionality coefficient)
      const double reltensions_;
      //! saturation value at zero pressure
      const double sat0_;
      //@}

      /// pressure ids involved in phase law
      const std::vector<int>* PresIds() override { return presids_; };

      //! evaluate saturation of phase
      double EvaluateSaturation(const std::vector<double>& pressure) override;

      //! evaluate generalized pressure related to phase law
      double EvaluateGenPressure(double saturation) override;

      //! evaluate derivative of saturation w.r.t. pressure
      double EvaluateDerivOfSaturationWrtPressure(
          int doftoderive, const std::vector<double>& pressure) override;

      //! evaluate second derivative of saturation w.r.t. pressure
      double EvaluateSecondDerivOfSaturationWrtPressure(int firstdoftoderive, int seconddoftoderive,
          const std::vector<double>& pressure) override;

      //! evaluate derivative of pressure w.r.t. saturtion
      double EvaluateDerivOfPressureWrtSaturation(int doftoderive, double saturation) override;
    };
    // class FluidPoroPhaseLawLinear

    /*----------------------------------------------------------------------*/
    //! tangent pressure-saturation relationship
    class FluidPoroPhaseLawTangent : public FluidPoroPhaseLaw
    {
     public:
      /// standard constructor
      explicit FluidPoroPhaseLawTangent(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// initialize
      void Initialize() override { return; };

      /// @name material parameters
      //@{
      //! number of degrees of freedom
      const int numdof_;
      //! coefficients defining the pressures involved int the pressure-saturation law
      const std::vector<int>* presids_;
      //! relative tension coefficient
      const double reltensions_;
      //! exponent within phase law
      const double exp_;
      //! saturation value at zero pressure
      const double sat0_;
      //@}

      /// pressure ids involved in phase law
      const std::vector<int>* PresIds() override { return presids_; };

      //! evaluate saturation of phase
      double EvaluateSaturation(const std::vector<double>& pressure) override;

      //! evaluate derivative of saturation w.r.t. pressure
      double EvaluateDerivOfSaturationWrtPressure(
          int doftoderive, const std::vector<double>& pressure) override;

      //! evaluate second derivative of saturation w.r.t. pressure
      double EvaluateSecondDerivOfSaturationWrtPressure(int firstdoftoderive, int seconddoftoderive,
          const std::vector<double>& pressure) override;

      //! evaluate generalized pressure related to phase law
      double EvaluateGenPressure(double saturation) override;

      //! evaluate derivative of pressure w.r.t. saturation
      double EvaluateDerivOfPressureWrtSaturation(int doftoderive, double saturation) override;
    };
    // class FluidPoroPhaseLawTangent

    /*----------------------------------------------------------------------*/
    //! a phase law indicating that the saturation is calculated from
    //! saturation constraint (i.e. the sum of all saturations is equal to one)
    class FluidPoroPhaseLawConstraint : public FluidPoroPhaseLaw
    {
     public:
      /// standard constructor
      explicit FluidPoroPhaseLawConstraint(Teuchos::RCP<MAT::PAR::Material> matdata)
          : FluidPoroPhaseLaw(matdata){};

      /// initialize
      void Initialize() override { return; };

      /// pressure ids involved in phase law
      const std::vector<int>* PresIds() override
      {
        FOUR_C_THROW(
            "The constraint phase law does not have pressure coupling! \n "
            "Combining Saturation DOF and constraint phase law is invalid!");
        return nullptr;
      };


      //! evaluate saturation of phase
      double EvaluateSaturation(const std::vector<double>& pressure) override
      {
        FOUR_C_THROW("The constraint phase law does not implement evaluation routines!");
        return 0.0;
      };

      //! evaluate derivative of saturation w.r.t. pressure
      double EvaluateDerivOfSaturationWrtPressure(
          int doftoderive, const std::vector<double>& pressure) override
      {
        FOUR_C_THROW("The constraint phase law does not implement evaluation routines!");
        return 0.0;
      };

      //! evaluate second derivative of saturation w.r.t. pressure
      double EvaluateSecondDerivOfSaturationWrtPressure(
          int firstdoftoderive, int seconddoftoderive, const std::vector<double>& pressure) override
      {
        FOUR_C_THROW("The constraint phase law does not implement evaluation routines!");
        return 0.0;
      };

      //! evaluate generalized pressure related to phase law
      double EvaluateGenPressure(double saturation) override
      {
        FOUR_C_THROW("The constraint phase law does not implement evaluation routines!");
        return 0.0;
      };

      //! evaluate derivative of pressure w.r.t. saturation
      double EvaluateDerivOfPressureWrtSaturation(int doftoderive, double saturation) override
      {
        FOUR_C_THROW("The constraint phase law does not implement evaluation routines!");
        return 0.0;
      };
    };
    // class FluidPoroPhaseLawTangent

    /*----------------------------------------------------------------------*/
    //! pressure-saturation relationship defined by function
    class FluidPoroPhaseLawByFunction : public FluidPoroPhaseLaw
    {
     public:
      /// standard constructor
      explicit FluidPoroPhaseLawByFunction(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// initialize
      void Initialize() override;

      /// @name material parameters
      //@{
      //! number of degrees of freedom
      const int numdof_;
      //! coefficients defining the pressures involved int the pressure-saturation law
      const std::vector<int>* presids_;
      //! function ID for evaluation of saturation
      const int functionID_saturation_;
      //! function ID for evaluation of pressure
      const int functionID_pressure_;
      //@}

      /// pressure ids involved in phase law
      const std::vector<int>* PresIds() override { return presids_; };

      //! evaluate saturation of phase
      double EvaluateSaturation(const std::vector<double>& pressure) override;

      //! evaluate derivative of saturation w.r.t. pressure
      double EvaluateDerivOfSaturationWrtPressure(
          int doftoderive, const std::vector<double>& pressure) override;

      //! evaluate second derivative of saturation w.r.t. pressure
      double EvaluateSecondDerivOfSaturationWrtPressure(int firstdoftoderive, int seconddoftoderive,
          const std::vector<double>& pressure) override;

      //! evaluate generalized pressure related to phase law
      double EvaluateGenPressure(double saturation) override;

      //! evaluate derivative of pressure w.r.t. saturation
      double EvaluateDerivOfPressureWrtSaturation(int doftoderive, double saturation) override;

     private:
      //! templated internal Initialize implementation
      template <int dim>
      void InitializeInternal();

      //! templated internal EvaluateSaturation implementation
      template <int dim>
      double EvaluateSaturationInternal(const std::vector<double>& pressure);

      //! templated internal EvaluateDerivOfSaturationWrtPressure implementation
      template <int dim>
      double EvaluateDerivOfSaturationWrtPressureInternal(
          int doftoderive, const std::vector<double>& pressure);

      //! templated internal EvaluateDerivOfPressureWrtSaturation implementation
      template <int dim>
      double EvaluateDerivOfPressureWrtSaturationInternal(int doftoderive, double saturation);

      //! templated internal EvaluateGenPressure implementation
      template <int dim>
      double EvaluateGenPressureInternal(double saturation);


      //! vector for input of differential pressure to function
      std::vector<std::pair<std::string, double>> dp_;
      //! vector for input of saturation to function
      std::vector<std::pair<std::string, double>> s_;
    };
    // class FluidPoroPhaseLawTangent

  }  // namespace PAR
}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
