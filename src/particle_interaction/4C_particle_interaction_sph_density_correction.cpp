/*---------------------------------------------------------------------------*/
/*! \file
\brief density correction handler in smoothed particle hydrodynamics (SPH)
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_interaction_sph_density_correction.hpp"

#include "4C_inpar_particle.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
ParticleInteraction::SPHDensityCorrectionBase::SPHDensityCorrectionBase()
{
  // empty constructor
}

void ParticleInteraction::SPHDensityCorrectionBase::Init()
{
  // nothing to do
}

void ParticleInteraction::SPHDensityCorrectionBase::Setup()
{
  // nothing to do
}

void ParticleInteraction::SPHDensityCorrectionBase::corrected_density_interior(
    const double* denssum, double* dens) const
{
  dens[0] = denssum[0];
}

ParticleInteraction::SPHDensityCorrectionInterior::SPHDensityCorrectionInterior()
    : ParticleInteraction::SPHDensityCorrectionBase()
{
  // empty constructor
}

bool ParticleInteraction::SPHDensityCorrectionInterior::ComputeDensityBC() const { return false; }

void ParticleInteraction::SPHDensityCorrectionInterior::corrected_density_free_surface(
    const double* denssum, const double* colorfield, const double* dens_bc, double* dens) const
{
  // density of free surface particles is not corrected
}

ParticleInteraction::SPHDensityCorrectionNormalized::SPHDensityCorrectionNormalized()
    : ParticleInteraction::SPHDensityCorrectionBase()
{
  // empty constructor
}

bool ParticleInteraction::SPHDensityCorrectionNormalized::ComputeDensityBC() const { return false; }

void ParticleInteraction::SPHDensityCorrectionNormalized::corrected_density_free_surface(
    const double* denssum, const double* colorfield, const double* dens_bc, double* dens) const
{
  dens[0] = denssum[0] / colorfield[0];
}

ParticleInteraction::SPHDensityCorrectionRandles::SPHDensityCorrectionRandles()
    : ParticleInteraction::SPHDensityCorrectionBase()
{
  // empty constructor
}

bool ParticleInteraction::SPHDensityCorrectionRandles::ComputeDensityBC() const { return true; }

void ParticleInteraction::SPHDensityCorrectionRandles::corrected_density_free_surface(
    const double* denssum, const double* colorfield, const double* dens_bc, double* dens) const
{
  dens[0] = denssum[0] + dens_bc[0] * (1.0 - colorfield[0]);
}

FOUR_C_NAMESPACE_CLOSE
