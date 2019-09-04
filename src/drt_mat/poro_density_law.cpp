/*----------------------------------------------------------------------*/
/*! \file
 \brief calculation classes for evaluation of constitutive relation for (microscopic) density in
 porous media

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/



#include "poro_density_law.H"

#include "matpar_bundle.H"

#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::PoroDensityLaw* MAT::PAR::PoroDensityLaw::CreateDensityLaw(int matID)
{
  // initialize null pointer
  MAT::PAR::PoroDensityLaw* densitylaw = NULL;

  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (DRT::Problem::Instance(probinst)->Materials() == Teuchos::null)
    dserror("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (DRT::Problem::Instance(probinst)->Materials()->Num() == 0)
    dserror("List of materials in the global problem instance is empty.");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(matID);

  switch (curmat->Type())
  {
    case INPAR::MAT::m_poro_densitylaw_constant:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::PoroDensityLawConstant(curmat));
      densitylaw = static_cast<MAT::PAR::PoroDensityLawConstant*>(curmat->Parameter());
      break;
    }
    case INPAR::MAT::m_poro_densitylaw_exp:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::PoroDensityLawExp(curmat));
      densitylaw = static_cast<MAT::PAR::PoroDensityLawExp*>(curmat->Parameter());
      break;
    }
    default:
      dserror("invalid material for density law %d", curmat->Type());
      break;
  }

  return densitylaw;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::PoroDensityLawExp::PoroDensityLawExp(Teuchos::RCP<MAT::PAR::Material> matdata)
    : PoroDensityLaw(matdata), bulkmodulus_(matdata->GetDouble("BULKMODULUS"))
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::PoroDensityLawExp::CreateMaterial() { return Teuchos::null; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::PoroDensityLawExp::ComputeCurDensity(const double& refdensity, const double& press)
{
  return refdensity * exp(press / bulkmodulus_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::PoroDensityLawExp::ComputeRefDensityToCurDensity(const double& press)
{
  return exp(-1.0 * press / bulkmodulus_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::PoroDensityLawExp::ComputeRefDensityToCurDensityDerivative(const double& press)
{
  return -1.0 / bulkmodulus_ * exp(-1.0 * press / bulkmodulus_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::PoroDensityLawExp::ComputeRefDensityToCurDensitySecondDerivative(
    const double& press)
{
  return 1.0 / (bulkmodulus_ * bulkmodulus_) * exp(-1.0 * press / bulkmodulus_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::PoroDensityLawExp::ComputeCurDensityDerivative(
    const double& refdensity, const double& press)
{
  return refdensity / bulkmodulus_ * exp(press / bulkmodulus_);
}
