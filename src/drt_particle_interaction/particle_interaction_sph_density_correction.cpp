/*---------------------------------------------------------------------------*/
/*! \file
\brief density correction handler in smoothed particle hydrodynamics (SPH)
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_density_correction.H"

#include "../drt_inpar/inpar_particle.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHDensityCorrectionBase::SPHDensityCorrectionBase()
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHDensityCorrectionBase::Init()
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHDensityCorrectionBase::Setup()
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHDensityCorrectionBase::WriteRestart(
    const int step, const double time) const
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHDensityCorrectionBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
  | set corrected density of interior particles               sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityCorrectionBase::CorrectedDensityInterior(
    const double* denssum, double* dens) const
{
  dens[0] = denssum[0];
}

PARTICLEINTERACTION::SPHDensityCorrectionInterior::SPHDensityCorrectionInterior()
    : PARTICLEINTERACTION::SPHDensityCorrectionBase()
{
  // empty constructor
}

bool PARTICLEINTERACTION::SPHDensityCorrectionInterior::ComputeDensityBC() const { return false; }

/*---------------------------------------------------------------------------*
  | set corrected density of free surface particles           sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityCorrectionInterior::CorrectedDensityFreeSurface(
    const double* denssum, const double* colorfield, const double* dens_bc, double* dens) const
{
  // density of free surface particles is not corrected
}

PARTICLEINTERACTION::SPHDensityCorrectionNormalized::SPHDensityCorrectionNormalized()
    : PARTICLEINTERACTION::SPHDensityCorrectionBase()
{
  // empty constructor
}

bool PARTICLEINTERACTION::SPHDensityCorrectionNormalized::ComputeDensityBC() const { return false; }

/*---------------------------------------------------------------------------*
  | set corrected density of free surface particles           sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityCorrectionNormalized::CorrectedDensityFreeSurface(
    const double* denssum, const double* colorfield, const double* dens_bc, double* dens) const
{
  dens[0] = denssum[0] / colorfield[0];
}

PARTICLEINTERACTION::SPHDensityCorrectionRandles::SPHDensityCorrectionRandles()
    : PARTICLEINTERACTION::SPHDensityCorrectionBase()
{
  // empty constructor
}

bool PARTICLEINTERACTION::SPHDensityCorrectionRandles::ComputeDensityBC() const { return true; }

void PARTICLEINTERACTION::SPHDensityCorrectionRandles::CorrectedDensityFreeSurface(
    const double* denssum, const double* colorfield, const double* dens_bc, double* dens) const
{
  dens[0] = denssum[0] + dens_bc[0] * (1.0 - colorfield[0]);
}
