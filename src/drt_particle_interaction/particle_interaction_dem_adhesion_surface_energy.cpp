/*---------------------------------------------------------------------------*/
/*!
\brief adhesion surface energy handler for discrete element method (DEM) interactions

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_dem_adhesion_surface_energy.H"

#include "../drt_inpar/inpar_particle.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyBase::DEMAdhesionSurfaceEnergyBase(
    const Teuchos::ParameterList& params)
    : params_dem_(params),
      adhesion_surface_energy_(params_dem_.get<double>("ADHESION_SURFACE_ENERGY"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init adhesion surface energy handler                       sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyBase::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup adhesion surface energy handler                      sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyBase::Setup()
{
  // safety check
  if (not(adhesion_surface_energy_ > 0.0)) dserror("non-positive adhesion surface energy!");
}

/*---------------------------------------------------------------------------*
 | write restart of adhesion surface energy handler           sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyBase::WriteRestart(
    const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of adhesion surface energy handler            sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyConstant::DEMAdhesionSurfaceEnergyConstant(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyBase(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionBase::
    DEMAdhesionSurfaceEnergyDistributionBase(const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyBase(params),
      variance_(params_dem_.get<double>("ADHESION_SURFACE_ENERGY_DISTRIBUTION_VAR")),
      adhesion_surface_energy_min_(
          adhesion_surface_energy_ -
          params_dem_.get<double>("ADHESION_SURFACE_ENERGY_DISTRIBUTION_CUTOFF_FACTOR") *
              variance_),
      adhesion_surface_energy_max_(
          adhesion_surface_energy_ +
          params_dem_.get<double>("ADHESION_SURFACE_ENERGY_DISTRIBUTION_CUTOFF_FACTOR") * variance_)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | setup adhesion surface energy handler                      sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionBase::Setup()
{
  // call base class setup
  DEMAdhesionSurfaceEnergyBase::Setup();

  // safety checks
  if (variance_ < 0.0) dserror("negative variance for adhesion surface energy distribution!");

  if (adhesion_surface_energy_min_ < 0.0) dserror("negative minimum adhesion surface energy!");
  if (adhesion_surface_energy_max_ < 0.0) dserror("negative maximum adhesion surface energy!");

  if (adhesion_surface_energy_min_ > adhesion_surface_energy_max_)
    dserror("minimum adhesion surface energy larger than maximum adhesion surface energy!");
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionNormal::
    DEMAdhesionSurfaceEnergyDistributionNormal(const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionBase(params),
      mean_(adhesion_surface_energy_)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | get adhesion surface energy                                sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionNormal::AdhesionSurfaceEnergy(
    double& adhesion_surface_energy) const
{
  // initialize random number generator
  DRT::Problem::Instance()->Random()->SetMeanVariance(mean_, variance_);

  // set normal distributed random value for surface energy
  adhesion_surface_energy = DRT::Problem::Instance()->Random()->Normal();

  // adjust surface energy to allowed bounds
  if (adhesion_surface_energy > adhesion_surface_energy_max_)
    adhesion_surface_energy = adhesion_surface_energy_max_;
  else if (adhesion_surface_energy < adhesion_surface_energy_min_)
    adhesion_surface_energy = adhesion_surface_energy_min_;
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionLogNormal::
    DEMAdhesionSurfaceEnergyDistributionLogNormal(const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionBase(params),
      mean_(std::log(adhesion_surface_energy_))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | get adhesion surface energy                                sfuchs 07/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyDistributionLogNormal::AdhesionSurfaceEnergy(
    double& adhesion_surface_energy) const
{
  // initialize random number generator
  DRT::Problem::Instance()->Random()->SetMeanVariance(mean_, variance_);

  // set log-normal distributed random value for surface energy
  adhesion_surface_energy = std::exp(DRT::Problem::Instance()->Random()->Normal());

  // adjust surface energy to allowed bounds
  if (adhesion_surface_energy > adhesion_surface_energy_max_)
    adhesion_surface_energy = adhesion_surface_energy_max_;
  else if (adhesion_surface_energy < adhesion_surface_energy_min_)
    adhesion_surface_energy = adhesion_surface_energy_min_;
}
