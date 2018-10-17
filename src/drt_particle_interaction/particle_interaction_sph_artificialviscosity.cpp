/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_artificialviscosity.cpp

\brief artificial viscosity handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_artificialviscosity.H"

#include "../drt_lib/drt_dserror.H"

#include <cmath>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHArtificialViscosity::SPHArtificialViscosity()
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init artificial viscosity handler                          sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHArtificialViscosity::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup artificial viscosity handler                         sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHArtificialViscosity::Setup()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | write restart of artificial viscosity handler              sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHArtificialViscosity::WriteRestart(
    const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of artificial viscosity handler               sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHArtificialViscosity::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | evaluate artificial viscosity                              sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHArtificialViscosity::ArtificialViscosity(const double* dens_i,
    const double* dens_j, const double* vel_i, const double* vel_j, const double* mass_j,
    const double& artificialviscosity, const double& dWdrij, const double& h_i, const double& h_j,
    const double& c_i, const double& c_j, const double& abs_rij, const double* e_ij,
    double* acc_i) const
{
  // particle averaged smoothing length
  const double h_ij = 0.5 * (h_i + h_j);

  // particle averaged speed of sound
  const double c_ij = 0.5 * (c_i + c_j);

  // particle averaged density
  const double dens_ij = 0.5 * (dens_i[0] + dens_j[0]);

  const double e_ij_vrel_ij = ((vel_i[0] - vel_j[0]) * e_ij[0] + (vel_i[1] - vel_j[1]) * e_ij[1] +
                               (vel_i[2] - vel_j[2]) * e_ij[2]);

  // avoid division by zero for close particles
  const double epsilon = 0.01;

  const double fac = h_ij * c_ij * e_ij_vrel_ij * abs_rij /
                     (dens_ij * (std::pow(abs_rij, 2) + epsilon * std::pow(h_ij, 2)));

  for (int i = 0; i < 3; ++i) acc_i[i] += artificialviscosity * mass_j[0] * dWdrij * fac * e_ij[i];
}
