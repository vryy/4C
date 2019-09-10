/*---------------------------------------------------------------------------*/
/*! \file
\brief density correction handler in smoothed particle hydrodynamics (SPH)

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_density_correction.H"

#include "../drt_inpar/inpar_particle.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHDensityCorrectionBase::SPHDensityCorrectionBase()
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init density correction handler                            sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityCorrectionBase::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup density correction handler                           sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityCorrectionBase::Setup()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | write restart of density correction handler                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityCorrectionBase::WriteRestart(
    const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of density correction handler                 sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
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

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHDensityCorrectionInterior::SPHDensityCorrectionInterior()
    : PARTICLEINTERACTION::SPHDensityCorrectionBase()
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | density boundary condition is needed                       sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
bool PARTICLEINTERACTION::SPHDensityCorrectionInterior::ComputeDensityBC() const { return false; }

/*---------------------------------------------------------------------------*
  | set corrected density of free surface particles           sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityCorrectionInterior::CorrectedDensityFreeSurface(
    const double* denssum, const double* colorfield, const double* dens_bc, double* dens) const
{
  // density of free surface particles is not corrected
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHDensityCorrectionNormalized::SPHDensityCorrectionNormalized()
    : PARTICLEINTERACTION::SPHDensityCorrectionBase()
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | density boundary condition is needed                       sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
bool PARTICLEINTERACTION::SPHDensityCorrectionNormalized::ComputeDensityBC() const { return false; }

/*---------------------------------------------------------------------------*
  | set corrected density of free surface particles           sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityCorrectionNormalized::CorrectedDensityFreeSurface(
    const double* denssum, const double* colorfield, const double* dens_bc, double* dens) const
{
  dens[0] = denssum[0] / colorfield[0];
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHDensityCorrectionRandles::SPHDensityCorrectionRandles()
    : PARTICLEINTERACTION::SPHDensityCorrectionBase()
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | density boundary condition is needed                       sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
bool PARTICLEINTERACTION::SPHDensityCorrectionRandles::ComputeDensityBC() const { return true; }

/*---------------------------------------------------------------------------*
 | set corrected density of free surface particles            sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHDensityCorrectionRandles::CorrectedDensityFreeSurface(
    const double* denssum, const double* colorfield, const double* dens_bc, double* dens) const
{
  dens[0] = denssum[0] + dens_bc[0] * (1.0 - colorfield[0]);
}
