/*----------------------------------------------------------------------*/
/*! \file
 \brief calculation classes for evaluation of constitutive relation for (microscopic) density in
 porous media

   \level 3

 *----------------------------------------------------------------------*/



#include "4C_mat_poro_density_law.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::PoroDensityLaw* MAT::PAR::PoroDensityLaw::CreateDensityLaw(int matID)
{
  // initialize null pointer
  MAT::PAR::PoroDensityLaw* densitylaw = nullptr;

  // retrieve problem instance to read from
  const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (GLOBAL::Problem::Instance(probinst)->Materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (GLOBAL::Problem::Instance(probinst)->Materials()->Num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  // retrieve validated input line of material ID in question
  auto* curmat = GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matID);

  switch (curmat->Type())
  {
    case CORE::Materials::m_poro_densitylaw_constant:
    {
      densitylaw = static_cast<MAT::PAR::PoroDensityLawConstant*>(curmat);
      break;
    }
    case CORE::Materials::m_poro_densitylaw_exp:
    {
      densitylaw = static_cast<MAT::PAR::PoroDensityLawExp*>(curmat);
      break;
    }
    default:
      FOUR_C_THROW("invalid material for density law %d", curmat->Type());
      break;
  }

  return densitylaw;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::PoroDensityLawExp::PoroDensityLawExp(Teuchos::RCP<CORE::MAT::PAR::Material> matdata)
    : PoroDensityLaw(matdata), bulkmodulus_(matdata->Get<double>("BULKMODULUS"))
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::MAT::Material> MAT::PAR::PoroDensityLawExp::create_material()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::PoroDensityLawExp::ComputeCurDensity(const double& refdensity, const double& press)
{
  return refdensity * exp(press / bulkmodulus_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::PoroDensityLawExp::compute_ref_density_to_cur_density(const double& press)
{
  return exp(-1.0 * press / bulkmodulus_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::PoroDensityLawExp::compute_ref_density_to_cur_density_derivative(
    const double& press)
{
  return -1.0 / bulkmodulus_ * exp(-1.0 * press / bulkmodulus_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::PoroDensityLawExp::compute_ref_density_to_cur_density_second_derivative(
    const double& press)
{
  return 1.0 / (bulkmodulus_ * bulkmodulus_) * exp(-1.0 * press / bulkmodulus_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::PoroDensityLawExp::compute_cur_density_derivative(
    const double& refdensity, const double& press)
{
  return refdensity / bulkmodulus_ * exp(press / bulkmodulus_);
}

FOUR_C_NAMESPACE_CLOSE
