/*----------------------------------------------------------------------*/
/*! \file
\brief abstract interface for electrode and electrolyte materials carrying concentration and
electric potential as degrees of freedom


\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_MAT_ELCHSINGLEMAT_HPP
#define FOUR_C_MAT_ELCHSINGLEMAT_HPP

#include "4C_config.hpp"

#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MAT
{
  namespace PAR
  {
    //! parameters for abstract battery material
    class ElchSingleMat : public CORE::MAT::PAR::Parameter
    {
     public:
      //! @name parameters for abstract battery material
      //@{
      //! function number to describe concentration dependence of diffusion coefficient
      const int diffusion_coefficient_concentration_dependence_funct_num_;

      //! function number defining the temperature scaling of the diffusion coefficient
      const int diffusion_coefficient_temperature_scaling_funct_num_;

      //! number of parameters for diffusion coefficient
      const int number_diffusion_coefficent_params_;

      //! parameters for diffusion coefficient
      const std::vector<double> diffusion_coefficent_params_;

      //! number of parameters for scaling function describing temperature dependence of diffusion
      //! coefficient
      const int number_diffusion_temp_scale_funct_params_;

      //! parameters for scaling function describing temperature dependence of diffusion coefficient
      const std::vector<double> diffusion_temp_scale_funct_params_;

      //! function number to describe concentration dependence of conductivity
      const int conductivity_concentration_dependence_funct_num_;

      //! function number defining the temperature scaling of conductivity
      const int conductivity_temperature_scaling_funct_num_;

      //! number of parameters for conductivity
      const int number_conductivity_params_;

      //! parameters for conductivity
      const std::vector<double> conductivity_params_;

      //! number of parameters for scaling function describing temperature dependence of
      //! conductivity
      const int number_conductivity_temp_scale_funct_params_;

      //! parameters for scaling function describing temperature dependence conductivity
      const std::vector<double> conductivity_temp_scale_funct_params_;

      //! universal gas constant for evaluation of diffusion coefficient by means of
      //! Arrhenius-ansatz
      const double R_;
      //@}

     protected:
      //! constructor
      explicit ElchSingleMat(Teuchos::RCP<CORE::MAT::PAR::Material> matdata);

      //! check whether number of parameters is consistent with curve number
      void CheckProvidedParams(int functnr, const std::vector<double>& functparams);
    };  // class MAT::PAR::ElchSingleMat
  }     // namespace PAR


  /*----------------------------------------------------------------------*/
  //! wrapper for abstract battery material
  class ElchSingleMat : public CORE::MAT::Material
  {
   public:
    //! @name packing and unpacking
    /*!
      \brief Return unique ParObject id

      Every class implementing ParObject needs a unique id defined at the
      top of parobject.H (this file) and should return it in this method.
    */
    int UniqueParObjectId() const override = 0;

    /*!
      \brief Pack this class so it can be communicated

      Resizes the vector data and stores all information of a class in it.
      The first information to be stored in data has to be the
      unique ParObject ID delivered by UniqueParObjectId() which will then
      identify the exact class on the receiving processor.

      \param data (in/out): char vector to store class information
    */
    void Pack(CORE::COMM::PackBuffer& data) const override = 0;

    /*!
      \brief Unpack data from a char vector into this class

      The vector data contains all information to rebuild the
      exact copy of an instance of a class on a different processor.
      The first entry in data has to be an integer which is the unique
      parobject id defined at the top of this file and delivered by
      UniqueParObjectId().

      \param data (in) : vector storing all data to be unpacked into this instance.
    */
    void Unpack(const std::vector<char>& data) override = 0;
    //@}

    //! compute diffusion coefficient accounting for concentration and temperature dependence
    virtual double ComputeDiffusionCoefficient(double concentration, double temperature) const;

    //! compute concentration dependent diffusion coefficient according to function number
    double ComputeDiffusionCoefficientConcentrationDependent(double concentration) const;

    //! compute first derivative of diffusion coefficient w.r.t. concentration
    virtual double ComputeConcentrationDerivativeOfDiffusionCoefficient(
        double concentration, double temperature) const;

    //! compute first derivative of diffusion coefficient w.r.t. temperature
    double ComputeTemperatureDerivativeOfDiffusionCoefficient(
        double concentration, double temperature) const;

    //! compute conductivity accounting for concentration and temperature dependence
    double ComputeConductivity(double concentration, double temperature) const;

    //! compute concentration dependent conductivity according to function number
    double ComputeConductivityConcentrationDependent(double concentration) const;

    //! compute first derivative of conductivity w.r.t. concentration
    double ComputeConcentrationDerivativeOfConductivity(
        double concentration, double temperature) const;

    //! compute first derivative of conductivity w.r.t. temperature
    double ComputeTemperatureDerivativeOfConductivity(
        double concentration, double temperature) const;

    //! abbreviations for pre-defined functions
    //@{
    static constexpr int CONSTANT_FUNCTION = -1;
    static constexpr int LINEAR_FUNCTION = -2;
    static constexpr int QUADRATIC_FUNCTION = -3;
    static constexpr int POWER_FUNCTION = -4;
    static constexpr int CONDUCT = -5;
    static constexpr int MOD_CUBIC_FUNCTION = -6;
    static constexpr int CUBIC_FUNCTION = -7;
    static constexpr int NYMAN = -8;
    static constexpr int DEBYE_HUECKEL = -9;
    static constexpr int KOHLRAUSCH_SQUAREROOT = -10;
    static constexpr int GOLDIN = -11;
    static constexpr int STEWART_NEWMAN = -12;
    static constexpr int TDF = -13;
    static constexpr int ARRHENIUS = -14;
    static constexpr int INVERSE_LINEAR = -15;
    //@}

   protected:
    //! compute temperature dependent scale factor
    double ComputeTemperatureDependentScaleFactor(
        double temperature, int functionNumber, const std::vector<double>& functionParams) const;

    //! compute derivative of temperature dependent scale factor w.r.t. temperature
    double ComputeTemperatureDependentScaleFactorDeriv(
        double temperature, int functionNumber, const std::vector<double>& functionParams) const;

    //! return function number describing concentration dependence of the diffusion coefficient
    int DiffusionCoefficientConcentrationDependenceFunctNum() const
    {
      return dynamic_cast<MAT::PAR::ElchSingleMat*>(Parameter())
          ->diffusion_coefficient_concentration_dependence_funct_num_;
    };

    //! return the function number describing the temperature scaling of the diffusion coefficient
    int DiffusionCoefficientTemperatureScalingFunctNum() const
    {
      return dynamic_cast<MAT::PAR::ElchSingleMat*>(Parameter())
          ->diffusion_coefficient_temperature_scaling_funct_num_;
    };

    //! return function number describing concentration dependence of the conductivity
    int ConductivityConcentrationDependenceFunctNum() const
    {
      return dynamic_cast<MAT::PAR::ElchSingleMat*>(Parameter())
          ->conductivity_concentration_dependence_funct_num_;
    };

    //! return the function number describing the temperature scaling of the conductivity
    int ConductivityTemperatureScalingFunctNum() const
    {
      return dynamic_cast<MAT::PAR::ElchSingleMat*>(Parameter())
          ->conductivity_temperature_scaling_funct_num_;
    };

    //! return parameters for diffusion coefficient
    const std::vector<double>& DiffusionCoefficientParams() const
    {
      return dynamic_cast<MAT::PAR::ElchSingleMat*>(Parameter())->diffusion_coefficent_params_;
    };

    //! return parameters for temperature scaling function for diffusion coefficient
    const std::vector<double>& TempScaleFunctionParamsDiff() const
    {
      return dynamic_cast<MAT::PAR::ElchSingleMat*>(Parameter())
          ->diffusion_temp_scale_funct_params_;
    };

    //! return parameters for conductivity
    const std::vector<double>& ConductivityParams() const
    {
      return dynamic_cast<MAT::PAR::ElchSingleMat*>(Parameter())->conductivity_params_;
    };

    //! return parameters for temperature scaling function for conductivity
    const std::vector<double>& TempScaleFunctionParamsCond() const
    {
      return dynamic_cast<MAT::PAR::ElchSingleMat*>(Parameter())
          ->conductivity_temp_scale_funct_params_;
    };

    //! evaluate value as predefined function of any scalar (e.g. concentration, temperature)
    //!
    //! \param functnr      negative function number to be evaluated
    //! \param scalar       scalar value to insert into function
    //! \param functparams  constants that define the functions
    //! \return             function evaluated at value of scalar
    double EvalPreDefinedFunct(
        int functnr, double scalar, const std::vector<double>& functparams) const;

    //! evaluate first derivative of predefined function of any scalar (e.g. concentration,
    //! temperature)
    double EvalFirstDerivPreDefinedFunct(
        int functnr, double scalar, const std::vector<double>& functparams) const;
  };
}  // namespace MAT
FOUR_C_NAMESPACE_CLOSE

#endif
