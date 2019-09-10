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
    double press, double viscosity, double& dviscosity_dp)
{
  dviscosity_dp = 0.0;
  return;
}
