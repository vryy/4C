/*----------------------------------------------------------------------*/
/*! \file
 \brief calculation classes for evaluation of constitutive relation for
        relative permeability for multiphase porous flow

   \level 3

 *----------------------------------------------------------------------*/



#include "baci_mat_fluidporo_relpermeability_law.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_par_bundle.hpp"

BACI_NAMESPACE_OPEN
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroRelPermeabilityLaw*
MAT::PAR::FluidPoroRelPermeabilityLaw::CreateRelPermeabilityLaw(int matID)
{
  // initialize null pointer
  MAT::PAR::FluidPoroRelPermeabilityLaw* relpermeabilitylaw = nullptr;

  // retrieve problem instance to read from
  const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (GLOBAL::Problem::Instance(probinst)->Materials() == Teuchos::null)
    dserror("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (GLOBAL::Problem::Instance(probinst)->Materials()->Num() == 0)
    dserror("List of materials in the global problem instance is empty.");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      GLOBAL::Problem::Instance(probinst)->Materials()->ById(matID);

  switch (curmat->Type())
  {
    case INPAR::MAT::m_fluidporo_relpermeabilitylaw_constant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroRelPermeabilityLawConstant(curmat));
      relpermeabilitylaw =
          static_cast<MAT::PAR::FluidPoroRelPermeabilityLawConstant*>(curmat->Parameter());
      break;
    }
    case INPAR::MAT::m_fluidporo_relpermeabilitylaw_exp:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroRelPermeabilityLawExponent(curmat));
      relpermeabilitylaw =
          static_cast<MAT::PAR::FluidPoroRelPermeabilityLawExponent*>(curmat->Parameter());
      break;
    }
    default:
      dserror("invalid material for permeability law %d", curmat->Type());
      break;
  }

  return relpermeabilitylaw;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroRelPermeabilityLawConstant::FluidPoroRelPermeabilityLawConstant(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : FluidPoroRelPermeabilityLaw(matdata, true), relpermeability_(*matdata->Get<double>("VALUE"))
{
  if (relpermeability_ > 1.0)
    dserror(
        "relative permeability (actually the sum of the relative permeabilites) of phase cannot be "
        "greater than 1.0");
  return;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroRelPermeabilityLawExponent::FluidPoroRelPermeabilityLawExponent(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : FluidPoroRelPermeabilityLaw(matdata, false),
      exp_(*matdata->Get<double>("EXP")),
      minsat_(*matdata->Get<double>("MIN_SAT"))
{
  if (exp_ <= 1.0) dserror("exponent in relative permeability phase law has to be bigger than 1.0");
  // if(minsat_ < 0.0 or minsat_ > 1.0)
  //  dserror("minimal saturation has to be between 0 and 1");
  return;
}

BACI_NAMESPACE_CLOSE
