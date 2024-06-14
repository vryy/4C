/*--------------------------------------------------------------------------*/
/*! \file
\brief calculation classes for evaluation of constitutive relation for lubrication

\level 3

*/
/*--------------------------------------------------------------------------*/

#include "4C_mat_lubrication_law.hpp"

#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Mat::PAR::LubricationLaw::LubricationLaw(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Mat::PAR::LubricationLawConstant::LubricationLawConstant(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : LubricationLaw(matdata), viscosity_(matdata.parameters.Get<double>("VISCOSITY"))
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::LubricationLawConstant::create_material()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::PAR::LubricationLawConstant::ComputeViscosity(const double& press, double& viscosity)
{
  viscosity = viscosity_;
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::PAR::LubricationLawConstant::constitutive_derivatives(
    const double& press, const double& viscosity, double& dviscosity_dp)
{
  dviscosity_dp = 0.0;

  return;
}

/*---------------------------------------------------------------------*
 *  Method definitions for Barus viscosity
 *---------------------------------------------------------------------*/

// Standard Constructor
Mat::PAR::LubricationLawBarus::LubricationLawBarus(const Core::Mat::PAR::Parameter::Data& matdata)
    : LubricationLaw(matdata),
      ABSViscosity_(matdata.parameters.Get<double>("ABSViscosity")),
      PreVisCoeff_(matdata.parameters.Get<double>("PreVisCoeff"))
{
  return;
}

// Create material instance of matching type with my parameters
Teuchos::RCP<Core::Mat::Material> Mat::PAR::LubricationLawBarus::create_material()
{
  return Teuchos::null;
}

// Calculate the current viscosity
void Mat::PAR::LubricationLawBarus::ComputeViscosity(const double& press, double& viscosity)
{
  viscosity = ABSViscosity_ * (std::exp(PreVisCoeff_ * press));

  return;
}

// Evaluate constitutive relation for viscosity and compute derivatives
void Mat::PAR::LubricationLawBarus::constitutive_derivatives(
    const double& press, const double& viscosity, double& dviscosity_dp)
{
  dviscosity_dp = viscosity * PreVisCoeff_;

  return;
}

/*---------------------------------------------------------------------*
 *  Method definitions for Roeland viscosity
 *---------------------------------------------------------------------*/

// Standard Constructor
Mat::PAR::LubricationLawRoeland::LubricationLawRoeland(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : LubricationLaw(matdata),
      ABSViscosity_(matdata.parameters.Get<double>("ABSViscosity")),
      PreVisCoeff_(matdata.parameters.Get<double>("PreVisCoeff")),
      RefVisc_(matdata.parameters.Get<double>("RefVisc")),
      RefPress_(matdata.parameters.Get<double>("RefPress"))
{
  z_ = (PreVisCoeff_ * RefPress_) / (log(ABSViscosity_ / RefVisc_));
  return;
}

// Create material instance of matching type with my parameters
Teuchos::RCP<Core::Mat::Material> Mat::PAR::LubricationLawRoeland::create_material()
{
  return Teuchos::null;
}

// Calculate the current viscosity
void Mat::PAR::LubricationLawRoeland::ComputeViscosity(const double& press, double& viscosity)
{
  // double z = (PreVisCoeff_ * RefPress_) / (log ( ABSViscosity_ / RefVisc_ ));

  viscosity =
      ABSViscosity_ * exp(log(ABSViscosity_ / RefVisc_) * (pow((1 + press / RefPress_), z_) - 1));

  return;
}

// Evaluate constitutive relation for viscosity and compute derivatives
void Mat::PAR::LubricationLawRoeland::constitutive_derivatives(
    const double& press, const double& viscosity, double& dviscosity_dp)
{
  // double z = (PreVisCoeff_ * RefPress_ ) / (log ( ABSViscosity_ / RefVisc_ ));

  dviscosity_dp = viscosity * log(ABSViscosity_ / RefVisc_) * z_ *
                  pow((1 + press / RefPress_), (z_ - 1)) * (1 / RefPress_);

  return;
}

FOUR_C_NAMESPACE_CLOSE
