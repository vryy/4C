/*----------------------------------------------------------------------*/
/*!
 \brief calculation classes for evaluation of constitutive relation for
        relative permeability for multiphase porous flow

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/



#include "fluidporo_relpermeability_law.H"

#include "matpar_bundle.H"

#include "../drt_lib/drt_globalproblem.H"
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroRelPermeabilityLaw*
MAT::PAR::FluidPoroRelPermeabilityLaw::CreateRelPermeabilityLaw(int matID)
{
  // initialize null pointer
  MAT::PAR::FluidPoroRelPermeabilityLaw* relpermeabilitylaw = NULL;

  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (DRT::Problem::Instance(probinst)->Materials() == Teuchos::null)
    dserror("Sorry dude, cannot work out problem instance.");
  // yet another safety check
  if (DRT::Problem::Instance(probinst)->Materials()->Num() == 0)
    dserror("Sorry dude, no materials defined.");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(matID);

  switch (curmat->Type())
  {
    case INPAR::MAT::m_fluidporo_relpermeabilitylaw_constant:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::PAR::FluidPoroRelPermeabilityLawConstant(curmat));
      relpermeabilitylaw =
          static_cast<MAT::PAR::FluidPoroRelPermeabilityLawConstant*>(curmat->Parameter());
      break;
    }
    case INPAR::MAT::m_fluidporo_relpermeabilitylaw_exp:
    {
      if (curmat->Parameter() == NULL)
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
    : FluidPoroRelPermeabilityLaw(matdata, true), relpermeability_(matdata->GetDouble("VALUE"))
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
      exp_(matdata->GetDouble("EXP")),
      minsat_(matdata->GetDouble("MIN_SAT"))
{
  if (exp_ <= 1.0) dserror("exponent in relative permeability phase law has to be bigger than 1.0");
  // if(minsat_ < 0.0 or minsat_ > 1.0)
  //  dserror("minimal saturation has to be between 0 and 1");
  return;
}
