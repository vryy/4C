/*----------------------------------------------------------------------*/
/*!
\file particle_sph_eos.cpp

\brief equation of state handler for SPH methods

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                 sfuchs 11/17 |
 *----------------------------------------------------------------------*/
#include "particle_sph_eos.H"
#include <cmath>


/*----------------------------------------------------------------------*
 | constructor                                             sfuchs 11/17 |
 *----------------------------------------------------------------------*/
PARTICLE::EquationOfState_GenTait::EquationOfState_GenTait(
    const double& speedOfSound,
    const double& initDensity,
    const double& refDensFac,
    const double& exponent) : EquationOfState_Base(),
    speedOfSound_(speedOfSound),
    initDensity_(initDensity),
    refDensFac_(refDensFac),
    exponent_(exponent)
{
  initPressure_ = std::pow(speedOfSound_,2) * initDensity_ / exponent_;

  return;
}

/*----------------------------------------------------------------------*
 | determine the pressure                                  sfuchs 11/17 |
 *----------------------------------------------------------------------*/
double PARTICLE::EquationOfState_GenTait::DensityToPressure(
    const double& density
    )
{
  if (exponent_ == 1)
    return std::pow(speedOfSound_,2)*(density-refDensFac_*initDensity_);
  else
    return initPressure_*(std::pow((density/initDensity_),exponent_)-refDensFac_);
}

/*----------------------------------------------------------------------*
 | determine the density                                   sfuchs 11/17 |
 *----------------------------------------------------------------------*/
double PARTICLE::EquationOfState_GenTait::PressureToDensity(
    const double& pressure
    )
{
  if (exponent_ == 1)
    return pressure*std::pow(speedOfSound_,-2)+refDensFac_*initDensity_;
  else
    return initDensity_*std::pow(((pressure/initPressure_)+refDensFac_),(1.0/exponent_));
}

/*----------------------------------------------------------------------*
 | determine the energy                                    sfuchs 11/17 |
 *----------------------------------------------------------------------*/
double PARTICLE::EquationOfState_GenTait::DensityToEnergy(
    const double& density,
    const double& mass
    )
{
  // thermodynamic energy E with p=-dE/dV, T=dE/dS (see Espanol2003, Eq.(5))
  // Attention: currently only the pressure-dependent contribution of the thermodynamic energy is implemented! Thus, it is only valid for isentrop problems, i.e. dE/dS=0!
  // Remark: integration of the pressure law with the relation V=mass/density and integration constant from initial condition E(V=mass/initDensity)
  if (exponent_ == 1)
    return -std::pow(speedOfSound_,2)*mass*( std::log(std::pow(mass,2)/(initDensity_*density)) - refDensFac_*(1+(initDensity_/density)) );
  else
    return -initPressure_*( (1.0/(1-exponent_))*( mass/(std::pow(initDensity_,exponent_)*std::pow(density,(1-exponent_))) + mass/initDensity_ ) - refDensFac_*(mass/initDensity_+mass/density) );
}


/*----------------------------------------------------------------------*
 | constructor                                             sfuchs 11/17 |
 *----------------------------------------------------------------------*/
PARTICLE::EquationOfState_IdealGas::EquationOfState_IdealGas(
    const double& speedOfSound,
    const double& initDensity) : EquationOfState_Base(),
        speedOfSound_(speedOfSound),
        initDensity_(initDensity)
{
  return;
}

/*----------------------------------------------------------------------*
 | determine the pressure                                  sfuchs 11/17 |
 *----------------------------------------------------------------------*/
double PARTICLE::EquationOfState_IdealGas::DensityToPressure(
    const double& density
    )
{
  return std::pow(speedOfSound_,2)*density;
}

/*----------------------------------------------------------------------*
 | determine the density                                   sfuchs 11/17 |
 *----------------------------------------------------------------------*/
double PARTICLE::EquationOfState_IdealGas::PressureToDensity(
    const double& pressure
    )
{
  return std::pow(speedOfSound_,-2)*pressure;
}

/*----------------------------------------------------------------------*
 | determine the energy                                    sfuchs 11/17 |
 *----------------------------------------------------------------------*/
double PARTICLE::EquationOfState_IdealGas::DensityToEnergy(
    const double& density,
    const double& mass
    )
{
  // thermodynamic energy E with p=-dE/dV, T=dE/dS (see Espanol2003, Eq.(5))
  // Attention: currently only the pressure-dependent contribution of the thermodynamic energy is implemented! Thus, it is only valid for isentrop problems, i.e. dE/dS=0!
  // Remark: integration of the pressure law with the relation V=mass/density and integration constant from initial condition E(V=mass/initDensity)
  return -std::pow(speedOfSound_,2)*mass*std::log(std::pow(mass,2)/(initDensity_*density));
}
