/*----------------------------------------------------------------------*/
/*! \file
 \brief calculation classes for evaluation of constitutive relation for
        relative permeability for multiphase porous flow

   \level 3

 *----------------------------------------------------------------------*/



#include "4C_mat_fluidporo_relpermeability_law.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Mat::PAR::FluidPoroRelPermeabilityLaw*
Mat::PAR::FluidPoroRelPermeabilityLaw::create_rel_permeability_law(int matID)
{
  // initialize null pointer
  Mat::PAR::FluidPoroRelPermeabilityLaw* relpermeabilitylaw = nullptr;

  // retrieve problem instance to read from
  const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (Global::Problem::Instance(probinst)->Materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (Global::Problem::Instance(probinst)->Materials()->Num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  // retrieve validated input line of material ID in question
  auto* curmat = Global::Problem::Instance(probinst)->Materials()->ParameterById(matID);

  switch (curmat->Type())
  {
    case Core::Materials::m_fluidporo_relpermeabilitylaw_constant:
    {
      relpermeabilitylaw = static_cast<Mat::PAR::FluidPoroRelPermeabilityLawConstant*>(curmat);
      break;
    }
    case Core::Materials::m_fluidporo_relpermeabilitylaw_exp:
    {
      relpermeabilitylaw = static_cast<Mat::PAR::FluidPoroRelPermeabilityLawExponent*>(curmat);
      break;
    }
    default:
      FOUR_C_THROW("invalid material for permeability law %d", curmat->Type());
      break;
  }

  return relpermeabilitylaw;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Mat::PAR::FluidPoroRelPermeabilityLawConstant::FluidPoroRelPermeabilityLawConstant(
    Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : FluidPoroRelPermeabilityLaw(matdata, true), relpermeability_(matdata->Get<double>("VALUE"))
{
  if (relpermeability_ > 1.0)
    FOUR_C_THROW(
        "relative permeability (actually the sum of the relative permeabilites) of phase cannot be "
        "greater than 1.0");
  return;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Mat::PAR::FluidPoroRelPermeabilityLawExponent::FluidPoroRelPermeabilityLawExponent(
    Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : FluidPoroRelPermeabilityLaw(matdata, false),
      exp_(matdata->Get<double>("EXP")),
      minsat_(matdata->Get<double>("MIN_SAT"))
{
  if (exp_ <= 1.0)
    FOUR_C_THROW("exponent in relative permeability phase law has to be bigger than 1.0");
  // if(minsat_ < 0.0 or minsat_ > 1.0)
  //  FOUR_C_THROW("minimal saturation has to be between 0 and 1");
  return;
}

FOUR_C_NAMESPACE_CLOSE
