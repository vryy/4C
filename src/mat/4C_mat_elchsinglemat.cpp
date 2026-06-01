// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_elchsinglemat.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_utils_enum.hpp"
#include "4C_utils_function_of_time.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Mat::PAR::ElchSingleMat::ElchSingleMat(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      diffusion_coefficient_concentration_dependence_funct_num_(
          matdata.parameters.get<int>("DIFF_COEF_CONC_DEP_FUNCT")),
      diffusion_coefficient_temperature_scaling_funct_num_(
          matdata.parameters.get<int>("DIFF_COEF_TEMP_SCALE_FUNCT")),
      number_diffusion_coefficient_params_(matdata.parameters.get<int>("DIFF_PARA_NUM")),
      diffusion_coefficient_params_(matdata.parameters.get<std::vector<double>>("DIFF_PARA")),
      number_diffusion_temp_scale_funct_params_(
          matdata.parameters.get<int>("DIFF_COEF_TEMP_SCALE_FUNCT_PARA_NUM")),
      diffusion_temp_scale_funct_params_(
          matdata.parameters.get<std::vector<double>>("DIFF_COEF_TEMP_SCALE_FUNCT_PARA")),
      conductivity_concentration_dependence_funct_num_(
          matdata.parameters.get<int>("COND_CONC_DEP_FUNCT")),
      conductivity_temperature_scaling_funct_num_(
          matdata.parameters.get<int>("COND_TEMP_SCALE_FUNCT")),
      number_conductivity_params_(matdata.parameters.get<int>("COND_PARA_NUM")),
      conductivity_params_(matdata.parameters.get<std::vector<double>>("COND_PARA")),
      number_conductivity_temp_scale_funct_params_(
          matdata.parameters.get<int>("COND_TEMP_SCALE_FUNCT_PARA_NUM")),
      conductivity_temp_scale_funct_params_(
          matdata.parameters.get<std::vector<double>>("COND_TEMP_SCALE_FUNCT_PARA"))
{
  // safety checks
  if (number_diffusion_coefficient_params_ !=
      static_cast<int>(diffusion_coefficient_params_.size()))
    FOUR_C_THROW("Mismatch in number of parameters for diffusion coefficient!");
  if (number_conductivity_params_ != static_cast<int>(conductivity_params_.size()))
    FOUR_C_THROW("Mismatch in number of parameters for conductivity!");
  if (number_diffusion_temp_scale_funct_params_ !=
      static_cast<int>(diffusion_temp_scale_funct_params_.size()))
    FOUR_C_THROW(
        "Mismatch in number of parameters for temp scale function for diffusion coefficient!");
  if (number_conductivity_temp_scale_funct_params_ !=
      static_cast<int>(conductivity_temp_scale_funct_params_.size()))
    FOUR_C_THROW("Mismatch in number of parameters for temp scale function for conductivity!");
  check_provided_params(
      diffusion_coefficient_concentration_dependence_funct_num_, diffusion_coefficient_params_);
  check_provided_params(conductivity_concentration_dependence_funct_num_, conductivity_params_);
  check_provided_params(
      diffusion_coefficient_temperature_scaling_funct_num_, diffusion_temp_scale_funct_params_);
  check_provided_params(
      conductivity_temperature_scaling_funct_num_, conductivity_temp_scale_funct_params_);
}


/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void Mat::PAR::ElchSingleMat::check_provided_params(
    const int functnr, const std::vector<double>& functparams)
{
  // name of specified curve
  std::string functionname;

  // check set of implemented functions with negative curve number
  if (functnr < 0)
  {
    // expected number of parameters for specified curve
    unsigned int nfunctparams = 0;

    switch (functnr)
    {
      case Mat::ElchSingleMat::CONSTANT_FUNCTION:
      {
        // constant value: functval=functparams[0];
        functionname = "'constant function'";
        nfunctparams = 1;
        break;
      }
      default:
      {
        FOUR_C_THROW("Curve number {} is not implemented", functnr);
      }
    }

    // safety check
    if (functparams.size() != nfunctparams)
    {
      FOUR_C_THROW(
          "Number of provided parameters does not match number of expected parameters for function "
          "with curve number {} ({})!",
          functnr, functionname.c_str());
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_diffusion_coefficient(
    const double concentration, const double temperature) const
{
  // D(c)
  double diffusionCoefficient =
      compute_diffusion_coefficient_concentration_dependent(concentration);

  // D(c)*D(T)
  diffusionCoefficient *= compute_temperature_dependent_scale_factor(temperature,
      diffusion_coefficient_temperature_scaling_funct_num(), temp_scale_function_params_diff());

  return diffusionCoefficient;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_diffusion_coefficient_concentration_dependent(
    const double concentration) const
{
  double diffusionCoefficient(0.0);

  // evaluate pre implemented concentration dependent diffusion coefficient
  if (diffusion_coefficient_concentration_dependence_funct_num() < 0)
  {
    diffusionCoefficient =
        eval_pre_defined_funct(diffusion_coefficient_concentration_dependence_funct_num(),
            concentration, diffusion_coefficient_params());
  }
  else if (diffusion_coefficient_concentration_dependence_funct_num() == 0)
  {
    FOUR_C_THROW(
        "'DIFF_COEF_CONC_DEP_FUNCT' must not be 0! Either set it to a negative value to use one of "
        "the implemented models, or set a positive value and make use of the function framework!");
  }
  // diffusion coefficient is a function of the concentration as defined in the input file
  else
  {
    diffusionCoefficient = Global::Problem::instance()
                               ->function_by_id<Core::Utils::FunctionOfTime>(
                                   diffusion_coefficient_concentration_dependence_funct_num())
                               .evaluate(concentration);
  }

  return diffusionCoefficient;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_temperature_dependent_scale_factor(const double temperature,
    const int functionNumber, const std::vector<double>& functionParams) const
{
  double temperatureDependentScaleFactor(1.0);

  if (functionNumber == 0)
  {
    // do nothing
  }
  else if (functionNumber < 0)
  {
    temperatureDependentScaleFactor =
        eval_pre_defined_funct(functionNumber, temperature, functionParams);
  }
  else if (functionNumber > 0)
  {
    temperatureDependentScaleFactor =
        Global::Problem::instance()
            ->function_by_id<Core::Utils::FunctionOfTime>(functionNumber)
            .evaluate(temperature);
  }
  else
  {
    FOUR_C_THROW(
        "You have to set a reasonable function number for the temperature dependence.\n This can "
        "be either 0 if no temperature dependence is desired, the number of the function in which "
        "you defined the temperature dependence, or a negative number representing a predefined "
        "function");
  }

  return temperatureDependentScaleFactor;
}

/*------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_concentration_derivative_of_diffusion_coefficient(
    const double concentration, const double temperature) const
{
  // derivative of diffusion coefficient w.r.t. concentration
  double diffusion_coeff_conc_deriv(0.0);

  // evaluate derivative w.r.t. concentration of pre implemented concentration dependent diffusion
  // coefficient
  if (diffusion_coefficient_concentration_dependence_funct_num() < 0)
  {
    diffusion_coeff_conc_deriv = eval_first_deriv_pre_defined_funct(
        diffusion_coefficient_concentration_dependence_funct_num(), concentration,
        diffusion_coefficient_params());
  }
  else if (diffusion_coefficient_concentration_dependence_funct_num() == 0)
  {
    FOUR_C_THROW(
        "'DIFF_COEF_CONC_DEP_FUNCT' must not be 0! Either set it to a negative value to use one of "
        "the implemented models, or set a positive value and make use of the function framework!");
  }
  // evaluate derivative w.r.t. concentration of diffusion coefficient that is a function of the
  // concentration as defined in the input file
  else
  {
    diffusion_coeff_conc_deriv = (Global::Problem::instance()
            ->function_by_id<Core::Utils::FunctionOfTime>(
                diffusion_coefficient_concentration_dependence_funct_num())
            .evaluate_derivative(concentration));
  }

  // do the temperature dependent scaling
  diffusion_coeff_conc_deriv *= compute_temperature_dependent_scale_factor(temperature,
      diffusion_coefficient_temperature_scaling_funct_num(), temp_scale_function_params_diff());

  return diffusion_coeff_conc_deriv;
}


/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_temperature_derivative_of_diffusion_coefficient(
    const double concentration, const double temperature) const
{
  // Computation D(c)
  double diffusion_coeff_temp_deriv =
      compute_diffusion_coefficient_concentration_dependent(concentration);

  // Computation D(c)*D'(T)
  diffusion_coeff_temp_deriv *= compute_temperature_dependent_scale_factor_deriv(temperature,
      diffusion_coefficient_temperature_scaling_funct_num(), temp_scale_function_params_diff());

  return diffusion_coeff_temp_deriv;
}

/*---------------------------------------------------------------------*

 *---------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_temperature_dependent_scale_factor_deriv(
    const double temperature, const int functionNumber,
    const std::vector<double>& functionParams) const
{
  double temperatureDependentScaleFactorDeriv(0.0);

  if (functionNumber == 0)
  {
    // do nothing
  }
  else if (functionNumber < 0)
  {
    temperatureDependentScaleFactorDeriv =
        eval_first_deriv_pre_defined_funct(functionNumber, temperature, functionParams);
  }
  else if (functionNumber > 0)
  {
    temperatureDependentScaleFactorDeriv =
        Global::Problem::instance()
            ->function_by_id<Core::Utils::FunctionOfTime>(functionNumber)
            .evaluate_derivative(temperature);
  }
  else
  {
    FOUR_C_THROW(
        "You have to set a reasonable function number for the temperature dependence.\n This can "
        "be either 0 if no temperature dependence is desired or the number of the function in "
        "which you defined the temperature dependence.");
  }

  return temperatureDependentScaleFactorDeriv;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_conductivity(
    const double concentration, const double temperature) const
{
  double conductivity = compute_conductivity_concentration_dependent(concentration);

  // do the temperature dependent scaling
  conductivity *= compute_temperature_dependent_scale_factor(
      temperature, conductivity_temperature_scaling_funct_num(), temp_scale_function_params_cond());

  return conductivity;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_conductivity_concentration_dependent(
    const double concentration) const
{
  double conductivity(0.0);

  // evaluate pre implemented concentration dependent conductivity
  if (conductivity_concentration_dependence_funct_num() < 0)
  {
    conductivity = eval_pre_defined_funct(
        conductivity_concentration_dependence_funct_num(), concentration, conductivity_params());
  }
  else if (conductivity_concentration_dependence_funct_num() == 0)
  {
    FOUR_C_THROW(
        "'COND_CONC_DEP_FUNCT' must not be 0! Either set it to a negative value to use one of the "
        "implemented models, or set a positive value and make use of the function framework!");
  }
  // conductivity is a function of the concentration as defined in the input file
  else
  {
    conductivity = Global::Problem::instance()
                       ->function_by_id<Core::Utils::FunctionOfTime>(
                           conductivity_concentration_dependence_funct_num())
                       .evaluate(concentration);
  }

  return conductivity;
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_concentration_derivative_of_conductivity(
    const double concentration, const double temperature) const
{
  // derivative of conductivity w.r.t. concentration
  double conductivity_conc_deriv(0.0);

  // evaluate derivative w.r.t. concentration of pre implemented concentration dependent
  // conductivity
  if (conductivity_concentration_dependence_funct_num() < 0)
  {
    conductivity_conc_deriv = eval_first_deriv_pre_defined_funct(
        conductivity_concentration_dependence_funct_num(), concentration, conductivity_params());
  }
  else if (conductivity_concentration_dependence_funct_num() == 0)
  {
    FOUR_C_THROW(
        "'COND_CONC_DEP_FUNCT' must not be 0! Either set it to a negative value to use one of the "
        "implemented models, or set a positive value and make use of the function framework!");
  }
  // evaluate derivative w.r.t. concentration of conductivity that is a function of the
  // concentration as defined in the input file
  else
  {
    conductivity_conc_deriv = (Global::Problem::instance()
            ->function_by_id<Core::Utils::FunctionOfTime>(
                conductivity_concentration_dependence_funct_num())
            .evaluate_derivative(concentration));
  }

  // do the temperature dependent scaling
  conductivity_conc_deriv *= compute_temperature_dependent_scale_factor(
      temperature, conductivity_temperature_scaling_funct_num(), temp_scale_function_params_cond());

  return conductivity_conc_deriv;
}


/*-----------------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------------*/
double Mat::ElchSingleMat::compute_temperature_derivative_of_conductivity(
    const double concentration, const double temperature) const
{
  // Cond(c)
  double conductivity_temp_deriv = compute_conductivity_concentration_dependent(concentration);

  // Cond(c)*Cond'(T)
  // derivative of temp dependent scaling factor
  conductivity_temp_deriv *= compute_temperature_dependent_scale_factor_deriv(
      temperature, conductivity_temperature_scaling_funct_num(), temp_scale_function_params_cond());

  return conductivity_temp_deriv;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::ElchSingleMat::eval_pre_defined_funct(
    const int functnr, const double scalar, const std::vector<double>& functparams) const
{
  double functval(0.0);

  switch (functnr)
  {
    // a0
    case CONSTANT_FUNCTION:
      functval = functparams[0];
      break;

    default:
    {
      FOUR_C_THROW("Curve number {} is not implemented!", functnr);
    }
  }

  return functval;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Mat::ElchSingleMat::eval_first_deriv_pre_defined_funct(
    const int functnr, const double scalar, const std::vector<double>& functparams) const
{
  double firstderivfunctval(0.0);

  switch (functnr)
  {
    // d/dc: a0
    case CONSTANT_FUNCTION:
      firstderivfunctval = 0.0;
      break;

    default:
    {
      FOUR_C_THROW("Curve number {} is not implemented!", functnr);
    }
  }

  return firstderivfunctval;
}

FOUR_C_NAMESPACE_CLOSE
