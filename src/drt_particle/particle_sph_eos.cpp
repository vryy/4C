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
    const double& refDensFac,
    const double& exponent) : EquationOfState_Base(),
    speedOfSound_(speedOfSound),
    refDensFac_(refDensFac),
    exponent_(exponent)
{
  return;
}

/*----------------------------------------------------------------------*
 | determine the pressure                                  sfuchs 11/17 |
 *----------------------------------------------------------------------*/
double PARTICLE::EquationOfState_GenTait::DensityToPressure(
    const double& density,
    const double& density0
    )
{

  if (exponent_ == 1)
    return std::pow(speedOfSound_,2)*(density-refDensFac_*density0);
  else
  {
    double initPressure = std::pow(speedOfSound_,2) * density0 / exponent_;
    return initPressure*(std::pow((density/density0),exponent_)-refDensFac_);
  }
}

/*----------------------------------------------------------------------*
 | determine the density                                   sfuchs 11/17 |
 *----------------------------------------------------------------------*/
double PARTICLE::EquationOfState_GenTait::PressureToDensity(
    const double& pressure,
    const double& density0
    )
{
  if (exponent_ == 1)
    return pressure*std::pow(speedOfSound_,-2)+refDensFac_*density0;
  else
  {
    double initPressure = std::pow(speedOfSound_,2) * density0 / exponent_;
    return density0*std::pow(((pressure/initPressure)+refDensFac_),(1.0/exponent_));
  }
}

/*----------------------------------------------------------------------*
 | determine the energy                                    sfuchs 11/17 |
 *----------------------------------------------------------------------*/
double PARTICLE::EquationOfState_GenTait::DensityToEnergy(
    const double& density,
    const double& mass,
    const double& density0
    )
{
  // thermodynamic energy E with p=-dE/dV, T=dE/dS (see Espanol2003, Eq.(5))
  // Attention: currently only the pressure-dependent contribution of the thermodynamic energy is implemented! Thus, it is only valid for isentrop problems, i.e. dE/dS=0!
  // Remark: integration of the pressure law with the relation V=mass/density and integration constant from initial condition E(V=mass/initDensity)
  if (exponent_ == 1)
    return -std::pow(speedOfSound_,2)*mass*( std::log(std::pow(mass,2)/(density0*density)) - refDensFac_*(1+(density0/density)) );
  else
  {
    double initPressure = std::pow(speedOfSound_,2) * density0 / exponent_;
    return -initPressure*( (1.0/(1-exponent_))*( mass/(std::pow(density0,exponent_)*std::pow(density,(1-exponent_))) + mass/density0 ) - refDensFac_*(mass/density0+mass/density) );
  }
}


/*----------------------------------------------------------------------*
 | constructor                                             sfuchs 11/17 |
 *----------------------------------------------------------------------*/
PARTICLE::EquationOfState_IdealGas::EquationOfState_IdealGas(
    const double& speedOfSound) : EquationOfState_Base(),
        speedOfSound_(speedOfSound)
{
  return;
}

/*----------------------------------------------------------------------*
 | determine the pressure                                  sfuchs 11/17 |
 *----------------------------------------------------------------------*/
double PARTICLE::EquationOfState_IdealGas::DensityToPressure(
    const double& density,
    const double& density0
    )
{
  return std::pow(speedOfSound_,2)*density;
}

/*----------------------------------------------------------------------*
 | determine the density                                   sfuchs 11/17 |
 *----------------------------------------------------------------------*/
double PARTICLE::EquationOfState_IdealGas::PressureToDensity(
    const double& pressure,
    const double& density0
    )
{
  return std::pow(speedOfSound_,-2)*pressure;
}

/*----------------------------------------------------------------------*
 | determine the energy                                    sfuchs 11/17 |
 *----------------------------------------------------------------------*/
double PARTICLE::EquationOfState_IdealGas::DensityToEnergy(
    const double& density,
    const double& mass,
    const double& density0
    )
{
  // thermodynamic energy E with p=-dE/dV, T=dE/dS (see Espanol2003, Eq.(5))
  // Attention: currently only the pressure-dependent contribution of the thermodynamic energy is implemented! Thus, it is only valid for isentrop problems, i.e. dE/dS=0!
  // Remark: integration of the pressure law with the relation V=mass/density and integration constant from initial condition E(V=mass/initDensity)
  return -std::pow(speedOfSound_,2)*mass*std::log(std::pow(mass,2)/(density0*density));
}
