/*----------------------------------------------------------------------*/
/*!
 \file fluidporo_viscosity_law.cpp

 \brief calculation classes for evaluation of constitutive relation for viscosity for multiphase porous flow

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/



#include "fluidporo_viscosity_law.H"

#include "matpar_bundle.H"

#include "../drt_lib/drt_globalproblem.H"
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroViscosityLaw* MAT::PAR::FluidPoroViscosityLaw::CreateViscosityLaw(int matID)
{
  // initialize null pointer
  MAT::PAR::FluidPoroViscosityLaw* viscositylaw = NULL;

  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (DRT::Problem::Instance(probinst)->Materials() == Teuchos::null)
    dserror("Sorry dude, cannot work out problem instance.");
  // yet another safety check
  if (DRT::Problem::Instance(probinst)->Materials()->Num() == 0)
    dserror("Sorry dude, no materials defined.");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat = DRT::Problem::Instance(probinst)->Materials()->ById(matID);

  switch (curmat->Type())
  {
  case INPAR::MAT::m_fluidporo_viscositylaw_constant:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::FluidPoroViscosityLawConstant(curmat));
    viscositylaw = static_cast<MAT::PAR::FluidPoroViscosityLawConstant*>(curmat->Parameter());
    break;
  }
  case INPAR::MAT::m_fluidporo_viscositylaw_celladh:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::PAR::FluidPoroViscosityLawCellAdherence(curmat));
    viscositylaw = static_cast<MAT::PAR::FluidPoroViscosityLawCellAdherence*>(curmat->Parameter());
    break;
  }
  default:
    dserror("invalid material for viscosity law %d", curmat->Type());
    break;
  }

  return viscositylaw;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroViscosityLawConstant::FluidPoroViscosityLawConstant(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: FluidPoroViscosityLaw(matdata,true),
  viscosity_(matdata->GetDouble("VALUE"))
{

  return;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroViscosityLawCellAdherence::FluidPoroViscosityLawCellAdherence(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: FluidPoroViscosityLaw(matdata,false),
  visc0_(matdata->GetDouble("VISC_0")),
  xi_(matdata->GetDouble("XI")),
  psi_(matdata->GetDouble("PSI"))

{
  if(visc0_ <= 0.0)
    dserror("VISC_0 cannot be smaller or equal to zero!");
  if(xi_ <= 0.0)
    dserror("XI cannot be smaller or equal to zero!");
  if(psi_ <= 0.0)
    dserror("PSI cannot be smaller or equal to zero!");

  return;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroViscosityLawCellAdherence::GetViscosity(
    const double abspressgrad
    ) const
{

  // visc = visc0 / ((1 - xi)*(1 - psi / | grad(pressure) |) * heaviside(1 - psi / | grad(pressure) |) + xi)

  if(abspressgrad <= psi_) // case heaviside(1 - psi / | grad(pressure) |) = 0
    return visc0_/xi_;
  else // case heaviside(1 - psi / | grad(pressure) |) = 1
    return visc0_/((1.0-xi_)*(1.0-psi_/abspressgrad)+xi_);

}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroViscosityLawCellAdherence::GetDerivOfViscosityWrtAbsPressGrad(
    const double abspressgrad
    ) const
{

  // visc = visc0 / ((1 - xi)*(1 - psi / | grad(pressure) |) * heaviside(1 - psi / | grad(pressure) |) + xi)
  // case heaviside = 1:
  // dvisc / d|grad(pressure) = visc0 * (xi - 1.0) * psi / ((xi - 1.0) * psi + | grad(pressure) |) ^ 2

  if(abspressgrad <= psi_) // case heaviside(1 - psi / | grad(pressure) |) = 0
    return 0.0;
  else // case heaviside(1 - psi / | grad(pressure) |) = 1
    return visc0_*(xi_-1.0)*psi_/((xi_-1.0)*psi_+abspressgrad)/((xi_-1.0)*psi_+abspressgrad);

}
