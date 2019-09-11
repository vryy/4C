/*--------------------------------------------------------------------------*/
/*! \file
\brief calculation classes for evaluation of constitutive relation for lubrication

\level 3

\maintainer Mostafa Faraji
*/
/*--------------------------------------------------------------------------*/

#include "matpar_bundle.H"
#include "lubrication_law.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::LubricationLaw::LubricationLaw(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::LubricationLawConstant::LubricationLawConstant(Teuchos::RCP<MAT::PAR::Material> matdata)
    : LubricationLaw(matdata), viscosity_(matdata->GetDouble("VISCOSITY"))
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::LubricationLawConstant::CreateMaterial()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::PAR::LubricationLawConstant::ComputeViscosity(const double& press, double& viscosity)
{
  viscosity = viscosity_;
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::PAR::LubricationLawConstant::ConstitutiveDerivatives(
    const double& press, const double& viscosity, double& dviscosity_dp)
{
  dviscosity_dp = 0.0;

  return;
}

/*---------------------------------------------------------------------*
 *  Method definitions for Barus viscosity
 *---------------------------------------------------------------------*/

// Standard Constructor
MAT::PAR::LubricationLawBarus::LubricationLawBarus(Teuchos::RCP<MAT::PAR::Material> matdata)
    : LubricationLaw(matdata),
      ABSViscosity_(matdata->GetDouble("ABSViscosity")),
      PreVisCoeff_(matdata->GetDouble("PreVisCoeff"))
{
  return;
}

// Create material instance of matching type with my parameters
Teuchos::RCP<MAT::Material> MAT::PAR::LubricationLawBarus::CreateMaterial()
{
  return Teuchos::null;
}

// Calculate the current viscosity
void MAT::PAR::LubricationLawBarus::ComputeViscosity(const double& press, double& viscosity)
{
  viscosity = ABSViscosity_ * (std::exp(PreVisCoeff_ * press));

  return;
}

// Evaluate constitutive relation for viscosity and compute derivatives
void MAT::PAR::LubricationLawBarus::ConstitutiveDerivatives(
    const double& press, const double& viscosity, double& dviscosity_dp)
{
  dviscosity_dp = viscosity * PreVisCoeff_;

  return;
}

/*---------------------------------------------------------------------*
 *  Method definitions for Roeland viscosity
 *---------------------------------------------------------------------*/

// Standard Constructor
MAT::PAR::LubricationLawRoeland::LubricationLawRoeland(Teuchos::RCP<MAT::PAR::Material> matdata)
    : LubricationLaw(matdata),
      ABSViscosity_(matdata->GetDouble("ABSViscosity")),
      PreVisCoeff_(matdata->GetDouble("PreVisCoeff")),
      RefVisc_(matdata->GetDouble("RefVisc")),
      RefPress_(matdata->GetDouble("RefPress"))
{
  z_ = (PreVisCoeff_ * RefPress_) / (log(ABSViscosity_ / RefVisc_));
  return;
}

// Create material instance of matching type with my parameters
Teuchos::RCP<MAT::Material> MAT::PAR::LubricationLawRoeland::CreateMaterial()
{
  return Teuchos::null;
}

// Calculate the current viscosity
void MAT::PAR::LubricationLawRoeland::ComputeViscosity(const double& press, double& viscosity)
{
  // double z = (PreVisCoeff_ * RefPress_) / (log ( ABSViscosity_ / RefVisc_ ));

  viscosity =
      ABSViscosity_ * exp(log(ABSViscosity_ / RefVisc_) * (pow((1 + press / RefPress_), z_) - 1));

  return;
}

// Evaluate constitutive relation for viscosity and compute derivatives
void MAT::PAR::LubricationLawRoeland::ConstitutiveDerivatives(
    const double& press, const double& viscosity, double& dviscosity_dp)
{
  // double z = (PreVisCoeff_ * RefPress_ ) / (log ( ABSViscosity_ / RefVisc_ ));

  dviscosity_dp = viscosity * log(ABSViscosity_ / RefVisc_) * z_ *
                  pow((1 + press / RefPress_), (z_ - 1)) * (1 / RefPress_);

  return;
}
