// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_ELCHSINGLEMAT_HPP
#define FOUR_C_MAT_ELCHSINGLEMAT_HPP

#include "4C_config.hpp"

#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  namespace PAR
  {
    //! parameters for abstract battery material
    class ElchSingleMat : public Core::Mat::PAR::Parameter
    {
     public:
      //! @name parameters for abstract battery material
      ///@{
      //! function number to describe concentration dependence of diffusion coefficient
      const int diffusion_coefficient_concentration_dependence_funct_num_;

      //! function number defining the temperature scaling of the diffusion coefficient
      const int diffusion_coefficient_temperature_scaling_funct_num_;

      //! number of parameters for diffusion coefficient
      const int number_diffusion_coefficient_params_;

      //! parameters for diffusion coefficient
      const std::vector<double> diffusion_coefficient_params_;

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
      ///@}

     protected:
      //! constructor
      explicit ElchSingleMat(const Core::Mat::PAR::Parameter::Data& matdata);

      //! check whether the number of @p functparams is consistent with the function chosen by
      //! @p functnr
      void check_provided_params(int functnr, const std::vector<double>& functparams);
    };  // class Mat::PAR::ElchSingleMat
  }  // namespace PAR


  /*----------------------------------------------------------------------*/
  //! wrapper for abstract battery material
  class ElchSingleMat : public Core::Mat::Material
  {
   public:
    //! @name packing and unpacking
    ///@{
    [[nodiscard]] int unique_par_object_id() const override = 0;

    void pack(Core::Communication::PackBuffer& data) const override = 0;

    void unpack(Core::Communication::UnpackBuffer& buffer) override = 0;
    ///@}

    //! compute diffusion coefficient accounting for concentration and temperature dependence
    [[nodiscard]] virtual double compute_diffusion_coefficient(
        double concentration, double temperature) const;

    //! compute concentration-dependent diffusion coefficient according to function number
    [[nodiscard]] double compute_diffusion_coefficient_concentration_dependent(
        double concentration) const;

    //! compute first derivative of diffusion coefficient w.r.t. concentration
    [[nodiscard]] virtual double compute_concentration_derivative_of_diffusion_coefficient(
        double concentration, double temperature) const;

    //! compute first derivative of diffusion coefficient w.r.t. temperature
    [[nodiscard]] double compute_temperature_derivative_of_diffusion_coefficient(
        double concentration, double temperature) const;

    //! compute conductivity accounting for concentration and temperature dependence
    [[nodiscard]] double compute_conductivity(double concentration, double temperature) const;

    //! compute concentration-dependent conductivity according to function number
    [[nodiscard]] double compute_conductivity_concentration_dependent(double concentration) const;

    //! compute the first derivative of conductivity w.r.t. concentration
    [[nodiscard]] double compute_concentration_derivative_of_conductivity(
        double concentration, double temperature) const;

    //! compute the first derivative of conductivity w.r.t. temperature
    [[nodiscard]] double compute_temperature_derivative_of_conductivity(
        double concentration, double temperature) const;

    //! abbreviations for pre-defined functions
    ///@{
    static constexpr int CONSTANT_FUNCTION = -1;
    static constexpr int GOLDIN = -11;
    static constexpr int ARRHENIUS = -14;
    ///@}

   protected:
    //! compute temperature dependent scale factor
    [[nodiscard]] double compute_temperature_dependent_scale_factor(
        double temperature, int functionNumber, const std::vector<double>& functionParams) const;

    //! compute derivative of temperature dependent scale factor w.r.t. temperature
    [[nodiscard]] double compute_temperature_dependent_scale_factor_deriv(
        double temperature, int functionNumber, const std::vector<double>& functionParams) const;

    //! return function number describing concentration dependence of the diffusion coefficient
    [[nodiscard]] int diffusion_coefficient_concentration_dependence_funct_num() const
    {
      return dynamic_cast<Mat::PAR::ElchSingleMat*>(parameter())
          ->diffusion_coefficient_concentration_dependence_funct_num_;
    }

    //! return the function number describing the temperature scaling of the diffusion coefficient
    [[nodiscard]] int diffusion_coefficient_temperature_scaling_funct_num() const
    {
      return dynamic_cast<Mat::PAR::ElchSingleMat*>(parameter())
          ->diffusion_coefficient_temperature_scaling_funct_num_;
    }

    //! return function number describing concentration dependence of the conductivity
    [[nodiscard]] int conductivity_concentration_dependence_funct_num() const
    {
      return dynamic_cast<Mat::PAR::ElchSingleMat*>(parameter())
          ->conductivity_concentration_dependence_funct_num_;
    }

    //! return the function number describing the temperature scaling of the conductivity
    [[nodiscard]] int conductivity_temperature_scaling_funct_num() const
    {
      return dynamic_cast<Mat::PAR::ElchSingleMat*>(parameter())
          ->conductivity_temperature_scaling_funct_num_;
    }

    //! return parameters for diffusion coefficient
    [[nodiscard]] const std::vector<double>& diffusion_coefficient_params() const
    {
      return dynamic_cast<Mat::PAR::ElchSingleMat*>(parameter())->diffusion_coefficient_params_;
    }

    //! return parameters for temperature scaling function for diffusion coefficient
    [[nodiscard]] const std::vector<double>& temp_scale_function_params_diff() const
    {
      return dynamic_cast<Mat::PAR::ElchSingleMat*>(parameter())
          ->diffusion_temp_scale_funct_params_;
    }

    //! return parameters for conductivity
    [[nodiscard]] const std::vector<double>& conductivity_params() const
    {
      return dynamic_cast<Mat::PAR::ElchSingleMat*>(parameter())->conductivity_params_;
    }

    //! return parameters for temperature scaling function for conductivity
    [[nodiscard]] const std::vector<double>& temp_scale_function_params_cond() const
    {
      return dynamic_cast<Mat::PAR::ElchSingleMat*>(parameter())
          ->conductivity_temp_scale_funct_params_;
    }

    //! evaluate value as predefined function of any scalar (e.g. concentration, temperature)
    //!
    //! \param functnr      negative function number to be evaluated
    //! \param scalar       scalar value to insert into function
    //! \param functparams  constants that define the functions
    //! \return             function evaluated at value of scalar
    [[nodiscard]] double eval_pre_defined_funct(
        int functnr, double scalar, const std::vector<double>& functparams) const;

    //! evaluate the first derivative of a predefined function of any scalar (e.g. concentration,
    //! temperature)
    [[nodiscard]] double eval_first_deriv_pre_defined_funct(
        int functnr, double scalar, const std::vector<double>& functparams) const;
  };
}  // namespace Mat
FOUR_C_NAMESPACE_CLOSE

#endif
