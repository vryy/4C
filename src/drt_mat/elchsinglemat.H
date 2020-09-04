/*----------------------------------------------------------------------*/
/*! \file
\brief abstract interface for electrode and electrolyte materials carrying concentration and
electric potential as degrees of freedom


\level 2
*/
/*----------------------------------------------------------------------*/
#ifndef ELCHSINGLEMAT_H
#define ELCHSINGLEMAT_H

#include "material.H"
#include "matpar_parameter.H"

namespace MAT
{
  namespace PAR
  {
    //! parameters for abstract battery material
    class ElchSingleMat : public Parameter
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

      //! function number to describe concentration dependence of conductivity
      const int conductivity_concentration_dependence_funct_num_;

      //! function number defining the temperature scaling of the conductivity
      const int conductivity_temperature_scaling_funct_num_;

      //! number of parameters for conductivity
      const int number_conductivity_params_;

      //! parameters for conductivity
      const std::vector<double> conductivity_params_;
      //@}

     protected:
      //! constructor
      explicit ElchSingleMat(Teuchos::RCP<MAT::PAR::Material> matdata);

      //! check whether number of parameters is consistent with curve number
      void CheckProvidedParams(int functnr, const std::vector<double>& functparams);
    };  // class MAT::PAR::ElchSingleMat
  }     // namespace PAR


  /*----------------------------------------------------------------------*/
  //! wrapper for abstract battery material
  class ElchSingleMat : public Material
  {
   public:
    //! @name packing and unpacking
    /*!
      \brief Return unique ParObject id

      Every class implementing ParObject needs a unique id defined at the
      top of drt_parobject.H (this file) and should return it in this method.
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
    void Pack(DRT::PackBuffer& data) const override = 0;

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
    double ComputeDiffusionCoefficient(double concentration, double temperature) const;

    //! compute concentration dependent diffusion coefficient according to function number
    double ComputeDiffusionCoefficientConcentrationDependent(double concentration) const;

    //! compute temperature dependent scale factor based on provided function
    double ComputeTemperatureDependentScaleFactor(double temperature, int functionNumber) const;

    //! compute first derivative of diffusion coefficient according to function number
    double ComputeFirstDerivDiffCoeff(double concentration) const;

    //! compute conductivity accounting for concentration and temperature dependence
    double ComputeConductivity(double concentration, double temperature) const;

    //! compute concentration dependent conductivity according to function number
    double ComputeConductivityConcentrationDependent(double concentration) const;

    //! compute first derivative of conductivity according to function number
    double ComputeFirstDerivCond(double concentration) const;

   protected:
    //! return function number describing concentration dependence of the diffusion coefficient
    int DiffusionCoefficientConcentrationDependenceFunctNum() const
    {
      return static_cast<MAT::PAR::ElchSingleMat*>(Parameter())
          ->diffusion_coefficient_concentration_dependence_funct_num_;
    };

    //! return the function number describing the temperature scaling of the diffusion coefficient
    int DiffusionCoefficientTemperatureScalingFunctNum() const
    {
      return static_cast<MAT::PAR::ElchSingleMat*>(Parameter())
          ->diffusion_coefficient_temperature_scaling_funct_num_;
    };

    //! return function number describing concentration dependence of the conductivity
    int ConductivityConcentrationDependenceFunctNum() const
    {
      return static_cast<MAT::PAR::ElchSingleMat*>(Parameter())
          ->conductivity_concentration_dependence_funct_num_;
    };

    //! return the function number describing the temperature scaling of the conductivity
    int ConductivityTemperatureScalingFunctNum() const
    {
      return static_cast<MAT::PAR::ElchSingleMat*>(Parameter())
          ->conductivity_temperature_scaling_funct_num_;
    };

    //! return parameters for diffusion coefficient
    const std::vector<double>& DiffusionCoefficientParams() const
    {
      return static_cast<MAT::PAR::ElchSingleMat*>(Parameter())->diffusion_coefficent_params_;
    };

    //! return parameters for conductivity
    const std::vector<double>& ConductivityParams() const
    {
      return static_cast<MAT::PAR::ElchSingleMat*>(Parameter())->conductivity_params_;
    };

    //! evaluate value of predefined function
    double EvalPreDefinedFunctValue(
        int functnr, double concentration, const std::vector<double>& functparams) const;

    //! evaluate first derivative of predefined function
    double EvalFirstDerivPreDefinedFunctValue(
        int functnr, double concentration, const std::vector<double>& functparams) const;
  };
}  // namespace MAT
#endif
