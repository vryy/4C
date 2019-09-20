/*---------------------------------------------------------------------------*/
/*! \file
\brief equation of state handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_equationofstate.H"

#include "particle_interaction_utils.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHEquationOfStateBase::SPHEquationOfStateBase()
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHEquationOfStateBase::Init()
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHEquationOfStateBase::Setup()
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHEquationOfStateBase::WriteRestart(
    const int step, const double time) const
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHEquationOfStateBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

PARTICLEINTERACTION::SPHEquationOfStateGenTait::SPHEquationOfStateGenTait(
    const double& speedofsound, const double& refdensfac, const double& exponent)
    : PARTICLEINTERACTION::SPHEquationOfStateBase(),
      speedofsound_(speedofsound),
      refdensfac_(refdensfac),
      exponent_(exponent)
{
  // empty constructor
}

double PARTICLEINTERACTION::SPHEquationOfStateGenTait::DensityToPressure(
    const double& density, const double& density0) const
{
  if (exponent_ == 1)
    return UTILS::pow<2>(speedofsound_) * (density - refdensfac_ * density0);
  else
  {
    double initPressure = UTILS::pow<2>(speedofsound_) * density0 / exponent_;
    return initPressure * (std::pow((density / density0), exponent_) - refdensfac_);
  }
}

double PARTICLEINTERACTION::SPHEquationOfStateGenTait::PressureToDensity(
    const double& pressure, const double& density0) const
{
  if (exponent_ == 1)
    return pressure / UTILS::pow<2>(speedofsound_) + refdensfac_ * density0;
  else
  {
    double initPressure = UTILS::pow<2>(speedofsound_) * density0 / exponent_;
    return density0 * std::pow(((pressure / initPressure) + refdensfac_), (1.0 / exponent_));
  }
}

double PARTICLEINTERACTION::SPHEquationOfStateGenTait::DensityToEnergy(
    const double& density, const double& mass, const double& density0) const
{
  // thermodynamic energy E with p=-dE/dV, T=dE/dS (see Espanol2003, Eq.(5))
  // Attention: currently only the pressure-dependent contribution of the thermodynamic energy is
  // implemented! Thus, it is only valid for isentrop problems, i.e. dE/dS=0! Remark: integration of
  // the pressure law with the relation V=mass/density and integration constant from initial
  // condition E(V=mass/initDensity)
  if (exponent_ == 1)
    return -UTILS::pow<2>(speedofsound_) * mass *
           (std::log(UTILS::pow<2>(mass) / (density0 * density)) -
               refdensfac_ * (1 + (density0 / density)));
  else
  {
    double initPressure = UTILS::pow<2>(speedofsound_) * density0 / exponent_;
    return -initPressure *
           ((1.0 / (1 - exponent_)) *
                   (mass / (std::pow(density0, exponent_) * std::pow(density, (1 - exponent_))) +
                       mass / density0) -
               refdensfac_ * (mass / density0 + mass / density));
  }
}

PARTICLEINTERACTION::SPHEquationOfStateIdealGas::SPHEquationOfStateIdealGas(
    const double& speedofsound)
    : PARTICLEINTERACTION::SPHEquationOfStateBase(), speedofsound_(speedofsound)
{
  // empty constructor
}

double PARTICLEINTERACTION::SPHEquationOfStateIdealGas::DensityToPressure(
    const double& density, const double& density0) const
{
  return UTILS::pow<2>(speedofsound_) * density;
}

double PARTICLEINTERACTION::SPHEquationOfStateIdealGas::PressureToDensity(
    const double& pressure, const double& density0) const
{
  return pressure / UTILS::pow<2>(speedofsound_);
}

double PARTICLEINTERACTION::SPHEquationOfStateIdealGas::DensityToEnergy(
    const double& density, const double& mass, const double& density0) const
{
  // thermodynamic energy E with p=-dE/dV, T=dE/dS (see Espanol2003, Eq.(5))
  // Attention: currently only the pressure-dependent contribution of the thermodynamic energy is
  // implemented! Thus, it is only valid for isentrop problems, i.e. dE/dS=0! Remark: integration of
  // the pressure law with the relation V=mass/density and integration constant from initial
  // condition E(V=mass/initDensity)
  return -UTILS::pow<2>(speedofsound_) * mass *
         std::log(UTILS::pow<2>(mass) / (density0 * density));
}
