/*---------------------------------------------------------------------------*/
/*! \file
\brief adhesion surface energy handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_dem_adhesion_surface_energy.H"

#include "../drt_inpar/inpar_particle.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyBase::DEMAdhesionSurfaceEnergyBase(
    const Teuchos::ParameterList& params)
    : params_dem_(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyBase::Init()
{
  // nothing to do
}

void PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyBase::Setup()
{
  // safety check
  if (not(params_dem_.get<double>("ADHESION_SURFACE_ENERGY") > 0.0))
    dserror("non-positive adhesion surface energy!");
}

PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyConstant::DEMAdhesionSurfaceEnergyConstant(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyBase(params)
{
  // empty constructor
}

PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionBase::
    DEMAdhesionSurfaceEnergyDistributionBase(const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyBase(params),
      variance_(params_dem_.get<double>("ADHESION_SURFACE_ENERGY_DISTRIBUTION_VAR")),
      cutofffactor_(params_dem_.get<double>("ADHESION_SURFACE_ENERGY_DISTRIBUTION_CUTOFF_FACTOR"))
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionBase::Setup()
{
  // call base class setup
  DEMAdhesionSurfaceEnergyBase::Setup();

  // safety checks
  if (variance_ < 0.0) dserror("negative variance for adhesion surface energy distribution!");
  if (cutofffactor_ < 0.0) dserror("negative cutoff factor of adhesion surface energy!");
}

void PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionBase::
    AdjustSurfaceEnergyToAllowedBounds(
        const double& mean_surface_energy, double& surface_energy) const
{
  const double adhesion_surface_energy_min =
      std::min(0.0, mean_surface_energy - cutofffactor_ * variance_);
  const double adhesion_surface_energy_max = mean_surface_energy + cutofffactor_ * variance_;

  if (surface_energy > adhesion_surface_energy_max)
    surface_energy = adhesion_surface_energy_max;
  else if (surface_energy < adhesion_surface_energy_min)
    surface_energy = adhesion_surface_energy_min;
}

PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionNormal::
    DEMAdhesionSurfaceEnergyDistributionNormal(const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionBase(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionNormal::AdhesionSurfaceEnergy(
    const double& mean_surface_energy, double& surface_energy) const
{
  // initialize random number generator
  DRT::Problem::Instance()->Random()->SetMeanVariance(mean_surface_energy, variance_);

  // set normal distributed random value for surface energy
  surface_energy = DRT::Problem::Instance()->Random()->Normal();

  // adjust surface energy to allowed bounds
  AdjustSurfaceEnergyToAllowedBounds(mean_surface_energy, surface_energy);
}

PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionLogNormal::
    DEMAdhesionSurfaceEnergyDistributionLogNormal(const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionBase(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionLogNormal::AdhesionSurfaceEnergy(
    const double& mean_surface_energy, double& surface_energy) const
{
  // initialize random number generator
  DRT::Problem::Instance()->Random()->SetMeanVariance(std::log(mean_surface_energy), variance_);

  // set log-normal distributed random value for surface energy
  surface_energy = std::exp(DRT::Problem::Instance()->Random()->Normal());

  // adjust surface energy to allowed bounds
  AdjustSurfaceEnergyToAllowedBounds(mean_surface_energy, surface_energy);
}
