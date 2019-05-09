/*----------------------------------------------------------------------*/
/*!
 \brief a material defining the pressure-saturation relationship of
        fluid phase within a multiphase porous fluid

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_MAT_FLUIDPORO_SINGLEPHASELAW_H_
#define SRC_DRT_MAT_FLUIDPORO_SINGLEPHASELAW_H_

#include "matpar_parameter.H"

namespace DRT
{
  namespace UTILS
  {
    class VariableExprFunction;
  }
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

      /// destructor
      virtual ~FluidPoroPhaseLaw() { ; }

      //! build the phase law
      static FluidPoroPhaseLaw* CreatePhaseLaw(int phaselawId);

      /// create material instance of matching type with my parameters
      virtual Teuchos::RCP<MAT::Material> CreateMaterial() { return Teuchos::null; };

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

      /// destructor
      virtual ~FluidPoroPhaseLawLinear() { ; }

      /// initialize
      virtual void Initialize() { return; };

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
      virtual const std::vector<int>* PresIds() { return presids_; };

      //! evaluate saturation of phase
      virtual double EvaluateSaturation(const std::vector<double>& pressure);

      //! evaluate generalized pressure related to phase law
      virtual double EvaluateGenPressure(double saturation);

      //! evaluate derivative of saturation w.r.t. pressure
      virtual double EvaluateDerivOfSaturationWrtPressure(
          int doftoderive, const std::vector<double>& pressure);

      //! evaluate second derivative of saturation w.r.t. pressure
      virtual double EvaluateSecondDerivOfSaturationWrtPressure(
          int firstdoftoderive, int seconddoftoderive, const std::vector<double>& pressure);

      //! evaluate derivative of pressure w.r.t. saturtion
      virtual double EvaluateDerivOfPressureWrtSaturation(int doftoderive, double saturation);
    };
    // class FluidPoroPhaseLawLinear

    /*----------------------------------------------------------------------*/
    //! tangent pressure-saturation relationship
    class FluidPoroPhaseLawTangent : public FluidPoroPhaseLaw
    {
     public:
      /// standard constructor
      explicit FluidPoroPhaseLawTangent(Teuchos::RCP<MAT::PAR::Material> matdata);

      /// destructor
      virtual ~FluidPoroPhaseLawTangent() { ; }

      /// initialize
      virtual void Initialize() { return; };

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
      virtual const std::vector<int>* PresIds() { return presids_; };

      //! evaluate saturation of phase
      virtual double EvaluateSaturation(const std::vector<double>& pressure);

      //! evaluate derivative of saturation w.r.t. pressure
      virtual double EvaluateDerivOfSaturationWrtPressure(
          int doftoderive, const std::vector<double>& pressure);

      //! evaluate second derivative of saturation w.r.t. pressure
      virtual double EvaluateSecondDerivOfSaturationWrtPressure(
          int firstdoftoderive, int seconddoftoderive, const std::vector<double>& pressure);

      //! evaluate generalized pressure related to phase law
      virtual double EvaluateGenPressure(double saturation);

      //! evaluate derivative of pressure w.r.t. saturation
      virtual double EvaluateDerivOfPressureWrtSaturation(int doftoderive, double saturation);
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

      /// destructor
      virtual ~FluidPoroPhaseLawConstraint() { ; }

      /// initialize
      virtual void Initialize() { return; };

      /// pressure ids involved in phase law
      virtual const std::vector<int>* PresIds()
      {
        dserror(
            "The constraint phase law does not have pressure coupling! \n "
            "Combining Saturation DOF and constraint phase law is invalid!");
        return NULL;
      };


      //! evaluate saturation of phase
      virtual double EvaluateSaturation(const std::vector<double>& pressure)
      {
        dserror("The constraint phase law does not implement evaluation routines!");
        return 0.0;
      };

      //! evaluate derivative of saturation w.r.t. pressure
      virtual double EvaluateDerivOfSaturationWrtPressure(
          int doftoderive, const std::vector<double>& pressure)
      {
        dserror("The constraint phase law does not implement evaluation routines!");
        return 0.0;
      };

      //! evaluate second derivative of saturation w.r.t. pressure
      virtual double EvaluateSecondDerivOfSaturationWrtPressure(
          int firstdoftoderive, int seconddoftoderive, const std::vector<double>& pressure)
      {
        dserror("The constraint phase law does not implement evaluation routines!");
        return 0.0;
      };

      //! evaluate generalized pressure related to phase law
      virtual double EvaluateGenPressure(double saturation)
      {
        dserror("The constraint phase law does not implement evaluation routines!");
        return 0.0;
      };

      //! evaluate derivative of pressure w.r.t. saturation
      virtual double EvaluateDerivOfPressureWrtSaturation(int doftoderive, double saturation)
      {
        dserror("The constraint phase law does not implement evaluation routines!");
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

      /// destructor
      virtual ~FluidPoroPhaseLawByFunction() { ; }

      /// initialize
      virtual void Initialize();

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
      virtual const std::vector<int>* PresIds() { return presids_; };

      //! evaluate saturation of phase
      virtual double EvaluateSaturation(const std::vector<double>& pressure);

      //! evaluate derivative of saturation w.r.t. pressure
      virtual double EvaluateDerivOfSaturationWrtPressure(
          int doftoderive, const std::vector<double>& pressure);

      //! evaluate second derivative of saturation w.r.t. pressure
      virtual double EvaluateSecondDerivOfSaturationWrtPressure(
          int firstdoftoderive, int seconddoftoderive, const std::vector<double>& pressure);

      //! evaluate generalized pressure related to phase law
      virtual double EvaluateGenPressure(double saturation);

      //! evaluate derivative of pressure w.r.t. saturation
      virtual double EvaluateDerivOfPressureWrtSaturation(int doftoderive, double saturation);

     private:
      //! helper method for function access
      inline DRT::UTILS::VariableExprFunction& Function(int functnum) const;

      //! vector for input of differential pressure to function
      std::vector<std::pair<std::string, double>> dp_;
      //! vector for input of saturation to function
      std::vector<std::pair<std::string, double>> S_;
    };
    // class FluidPoroPhaseLawTangent

  }  // namespace PAR
}  // namespace MAT

#endif /* SRC_DRT_MAT_FLUIDPORO_SINGLEPHASELAW_H_ */
