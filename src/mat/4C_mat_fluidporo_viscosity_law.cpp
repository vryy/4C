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
Mat::PAR::FluidPoroViscosityLaw* Mat::PAR::FluidPoroViscosityLaw::CreateViscosityLaw(int matID)
{
  // initialize null pointer
  Mat::PAR::FluidPoroViscosityLaw* viscositylaw = nullptr;

  // retrieve problem instance to read from
  const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (Global::Problem::Instance(probinst)->Materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (Global::Problem::Instance(probinst)->Materials()->Num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  auto* curmat = Global::Problem::Instance(probinst)->Materials()->ParameterById(matID);

  switch (curmat->Type())
  {
    case Core::Materials::m_fluidporo_viscositylaw_constant:
    {
      viscositylaw = static_cast<Mat::PAR::FluidPoroViscosityLawConstant*>(curmat);
      break;
    }
    case Core::Materials::m_fluidporo_viscositylaw_celladh:
    {
      viscositylaw = static_cast<Mat::PAR::FluidPoroViscosityLawCellAdherence*>(curmat);
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
Mat::PAR::FluidPoroViscosityLawConstant::FluidPoroViscosityLawConstant(
    Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : FluidPoroViscosityLaw(matdata, true), viscosity_(matdata->Get<double>("VALUE"))
{
  return;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Mat::PAR::FluidPoroViscosityLawCellAdherence::FluidPoroViscosityLawCellAdherence(
    Teuchos::RCP<Core::Mat::PAR::Material> matdata)
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
double Mat::PAR::FluidPoroViscosityLawCellAdherence::GetViscosity(const double abspressgrad) const
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
double Mat::PAR::FluidPoroViscosityLawCellAdherence::get_deriv_of_viscosity_wrt_abs_press_grad(
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
