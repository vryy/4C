/*----------------------------------------------------------------------*/
/*! \file
 \brief calculation classes for evaluation of constitutive relation for viscosity for multiphase
 porous flow

   \level 3

 *----------------------------------------------------------------------*/



#include "4C_mat_fluidporo_viscosity_law.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroViscosityLaw* MAT::PAR::FluidPoroViscosityLaw::CreateViscosityLaw(int matID)
{
  // initialize null pointer
  MAT::PAR::FluidPoroViscosityLaw* viscositylaw = nullptr;

  // retrieve problem instance to read from
  const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (GLOBAL::Problem::Instance(probinst)->Materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (GLOBAL::Problem::Instance(probinst)->Materials()->Num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<CORE::MAT::PAR::Material> curmat =
      GLOBAL::Problem::Instance(probinst)->Materials()->ById(matID);

  switch (curmat->Type())
  {
    case CORE::Materials::m_fluidporo_viscositylaw_constant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroViscosityLawConstant(curmat));
      viscositylaw = static_cast<MAT::PAR::FluidPoroViscosityLawConstant*>(curmat->Parameter());
      break;
    }
    case CORE::Materials::m_fluidporo_viscositylaw_celladh:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::FluidPoroViscosityLawCellAdherence(curmat));
      viscositylaw =
          static_cast<MAT::PAR::FluidPoroViscosityLawCellAdherence*>(curmat->Parameter());
      break;
    }
    default:
      FOUR_C_THROW("invalid material for viscosity law %d", curmat->Type());
      break;
  }

  return viscositylaw;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroViscosityLawConstant::FluidPoroViscosityLawConstant(
    Teuchos::RCP<CORE::MAT::PAR::Material> matdata)
    : FluidPoroViscosityLaw(matdata, true), viscosity_(matdata->Get<double>("VALUE"))
{
  return;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::PAR::FluidPoroViscosityLawCellAdherence::FluidPoroViscosityLawCellAdherence(
    Teuchos::RCP<CORE::MAT::PAR::Material> matdata)
    : FluidPoroViscosityLaw(matdata, false),
      visc0_(matdata->Get<double>("VISC_0")),
      xi_(matdata->Get<double>("XI")),
      psi_(matdata->Get<double>("PSI"))

{
  if (visc0_ <= 0.0) FOUR_C_THROW("VISC_0 cannot be smaller or equal to zero!");
  if (xi_ <= 0.0) FOUR_C_THROW("XI cannot be smaller or equal to zero!");
  if (psi_ <= 0.0) FOUR_C_THROW("PSI cannot be smaller or equal to zero!");

  return;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroViscosityLawCellAdherence::GetViscosity(const double abspressgrad) const
{
  // visc = visc0 / ((1 - xi)*(1 - psi / | grad(pressure) |) * heaviside(1 - psi / | grad(pressure)
  // |) + xi)

  if (abspressgrad <= psi_)  // case heaviside(1 - psi / | grad(pressure) |) = 0
    return visc0_ / xi_;
  else  // case heaviside(1 - psi / | grad(pressure) |) = 1
    return visc0_ / ((1.0 - xi_) * (1.0 - psi_ / abspressgrad) + xi_);
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double MAT::PAR::FluidPoroViscosityLawCellAdherence::get_deriv_of_viscosity_wrt_abs_press_grad(
    const double abspressgrad) const
{
  // visc = visc0 / ((1 - xi)*(1 - psi / | grad(pressure) |) * heaviside(1 - psi / | grad(pressure)
  // |) + xi) case heaviside = 1: dvisc / d|grad(pressure) = visc0 * (xi - 1.0) * psi / ((xi - 1.0)
  // * psi + | grad(pressure) |) ^ 2

  if (abspressgrad <= psi_)  // case heaviside(1 - psi / | grad(pressure) |) = 0
    return 0.0;
  else  // case heaviside(1 - psi / | grad(pressure) |) = 1
    return visc0_ * (xi_ - 1.0) * psi_ / ((xi_ - 1.0) * psi_ + abspressgrad) /
           ((xi_ - 1.0) * psi_ + abspressgrad);
}

FOUR_C_NAMESPACE_CLOSE
