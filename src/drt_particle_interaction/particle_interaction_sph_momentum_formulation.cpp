/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_momentum_formulation.cpp

\brief momentum formulation handler for smoothed particle hydrodynamics (SPH) interactions

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
#include "particle_interaction_sph_momentum_formulation.H"

#include "../drt_lib/drt_dserror.H"

#include <cmath>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHMomentumFormulationBase::SPHMomentumFormulationBase()
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init momentum formulation handler                          sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationBase::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup momentum formulation handler                         sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationBase::Setup()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | write restart of momentum formulation handler              sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationBase::WriteRestart(
    const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of momentum formulation handler               sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHMomentumFormulationMonaghan::SPHMomentumFormulationMonaghan()
    : PARTICLEINTERACTION::SPHMomentumFormulationBase()
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | evaluate specific coefficient                              sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationMonaghan::SpecificCoefficient(const double* dens_i,
    const double* dens_j, const double* mass_i, const double* mass_j, const double& dWdrij,
    double& specificcoefficient_ij) const
{
  specificcoefficient_ij = (dWdrij * mass_j[0]);
}

/*---------------------------------------------------------------------------*
 | evaluate pressure gradient                                 sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationMonaghan::PressureGradient(const double* dens_i,
    const double* dens_j, const double* press_i, const double* press_j,
    const double& specificcoefficient_ij, const double* e_ij, double* acc_i) const
{
  const double pressuregradient =
      -specificcoefficient_ij *
      (press_i[0] / std::pow(dens_i[0], 2) + press_j[0] / std::pow(dens_j[0], 2));

  for (int i = 0; i < 3; ++i) acc_i[i] += pressuregradient * e_ij[i];
}

/*---------------------------------------------------------------------------*
 | evaluate shear forces                                      sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationMonaghan::ShearForces(const double* dens_i,
    const double* dens_j, const double* vel_i, const double* vel_j, const double& kernelfac,
    const double& visc_i, const double& visc_j, const double& bulk_visc_i,
    const double& bulk_visc_j, const double& abs_rij, const double& specificcoefficient_ij,
    const double* e_ij, double* acc_i) const
{
  double viscosity = 0.0;
  if (visc_i > 0.0 and visc_j > 0.0) viscosity = (2.0 * visc_i * visc_j / (visc_i + visc_j));

  double bulkviscosity = 0.0;
  if (bulk_visc_i > 0.0 and bulk_visc_j > 0.0)
    bulkviscosity = (2.0 * bulk_visc_i * bulk_visc_j / (bulk_visc_i + bulk_visc_j));

  const double convectioncoefficient = kernelfac * (bulkviscosity + viscosity / 3.0);
  const double diffusioncoefficient = 5.0 * viscosity / 3.0 - bulkviscosity;

  // safety check
  if (diffusioncoefficient < 0.0) dserror("diffusion coefficient is negative!");

  const double inv_densi_densj_absdist = 1.0 / (dens_i[0] * dens_j[0] * abs_rij);
  const double e_ij_vrel_ij = ((vel_i[0] - vel_j[0]) * e_ij[0] + (vel_i[1] - vel_j[1]) * e_ij[1] +
                               (vel_i[2] - vel_j[2]) * e_ij[2]);

  for (int i = 0; i < 3; ++i)
  {
    // diffusion
    acc_i[i] += specificcoefficient_ij * diffusioncoefficient * inv_densi_densj_absdist *
                (vel_i[i] - vel_j[i]);

    // convection
    acc_i[i] += specificcoefficient_ij * convectioncoefficient * e_ij_vrel_ij *
                inv_densi_densj_absdist * e_ij[i];
  }
}

/*---------------------------------------------------------------------------*
 | evaluate background pressure (standard formulation)        sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationMonaghan::StandardBackgroundPressure(
    const double* dens_i, const double* dens_j, const double& backgroundpressure,
    const double& specificcoefficient_ij, const double* e_ij, double* acc_i) const
{
  const double fac = -specificcoefficient_ij * backgroundpressure *
                     (1.0 / std::pow(dens_i[0], 2) + 1.0 / std::pow(dens_j[0], 2));

  for (int i = 0; i < 3; ++i) acc_i[i] += fac * e_ij[i];
}

/*---------------------------------------------------------------------------*
 | evaluate background pressure (generalized formulation)     sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationMonaghan::GeneralizedBackgroundPressure(
    const double* dens_i, const double* mass_i, const double* mass_j,
    const double& backgroundpressure_tilde, const double& dWdrij_tilde, const double* e_ij,
    double* acc_i) const
{
  const double fac =
      -backgroundpressure_tilde * (mass_j[0] / std::pow(dens_i[0], 2)) * dWdrij_tilde;

  for (int i = 0; i < 3; ++i) acc_i[i] += fac * e_ij[i];
}

/*---------------------------------------------------------------------------*
 | evaluate modified velocity contribution                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationMonaghan::ModifiedVelocityContribution(
    const double* dens_i, const double* dens_j, const double* vel_i, const double* vel_j,
    const double* mod_vel_i, const double* mod_vel_j, const double& specificcoefficient_ij,
    const double* e_ij, double* acc_i) const
{
  double A_i[3][3];
  double A_j[3][3];

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
      A_i[i][j] = dens_i[0] * vel_i[i] * (mod_vel_i[j] - vel_i[j]);

      if (vel_j and mod_vel_j) A_j[i][j] = dens_j[0] * vel_j[i] * (mod_vel_j[j] - vel_j[j]);
    }

  double A_ij_e_ij[3];
  for (int i = 0; i < 3; ++i)
  {
    A_ij_e_ij[i] = (1.0 / std::pow(dens_i[0], 2)) *
                   (A_i[i][0] * e_ij[0] + A_i[i][1] * e_ij[1] + A_i[i][2] * e_ij[2]);

    if (vel_j and mod_vel_j)
      A_ij_e_ij[i] += (1.0 / std::pow(dens_j[0], 2)) *
                      (A_j[i][0] * e_ij[0] + A_j[i][1] * e_ij[1] + A_j[i][2] * e_ij[2]);
  }


  for (int i = 0; i < 3; ++i) acc_i[i] += specificcoefficient_ij * A_ij_e_ij[i];
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHMomentumFormulationAdami::SPHMomentumFormulationAdami()
    : PARTICLEINTERACTION::SPHMomentumFormulationBase()
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | evaluate specific coefficient                              sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationAdami::SpecificCoefficient(const double* dens_i,
    const double* dens_j, const double* mass_i, const double* mass_j, const double& dWdrij,
    double& specificcoefficient_ij) const
{
  specificcoefficient_ij =
      (std::pow((mass_i[0] / dens_i[0]), 2) + std::pow((mass_j[0] / dens_j[0]), 2)) *
      (dWdrij / mass_i[0]);
}

/*---------------------------------------------------------------------------*
 | evaluate pressure gradient                                 sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationAdami::PressureGradient(const double* dens_i,
    const double* dens_j, const double* press_i, const double* press_j,
    const double& specificcoefficient_ij, const double* e_ij, double* acc_i) const
{
  const double pressuregradient = -specificcoefficient_ij *
                                  (dens_i[0] * press_j[0] + dens_j[0] * press_i[0]) /
                                  (dens_i[0] + dens_j[0]);

  for (int i = 0; i < 3; ++i) acc_i[i] += pressuregradient * e_ij[i];
}

/*---------------------------------------------------------------------------*
 | evaluate shear forces                                      sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationAdami::ShearForces(const double* dens_i,
    const double* dens_j, const double* vel_i, const double* vel_j, const double& kernelfac,
    const double& visc_i, const double& visc_j, const double& bulk_visc_i,
    const double& bulk_visc_j, const double& abs_rij, const double& specificcoefficient_ij,
    const double* e_ij, double* acc_i) const
{
  double viscosity = 0.0;
  if (visc_i > 0.0 and visc_j > 0.0)
    viscosity = (2.0 * visc_i * visc_j / (visc_i + visc_j));
  else
    return;

  const double laminarviscosity = specificcoefficient_ij * viscosity / abs_rij;

  for (int i = 0; i < 3; ++i) acc_i[i] += laminarviscosity * (vel_i[i] - vel_j[i]);
}

/*---------------------------------------------------------------------------*
 | evaluate background pressure (standard formulation)        sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationAdami::StandardBackgroundPressure(
    const double* dens_i, const double* dens_j, const double& backgroundpressure,
    const double& specificcoefficient_ij, const double* e_ij, double* acc_i) const
{
  for (int i = 0; i < 3; ++i) acc_i[i] -= specificcoefficient_ij * backgroundpressure * e_ij[i];
}

/*---------------------------------------------------------------------------*
 | evaluate background pressure (generalized formulation)     sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationAdami::GeneralizedBackgroundPressure(
    const double* dens_i, const double* mass_i, const double* mass_j,
    const double& backgroundpressure_tilde, const double& dWdrij_tilde, const double* e_ij,
    double* acc_i) const
{
  const double fac =
      -(backgroundpressure_tilde / mass_i[0]) * std::pow((mass_i[0] / dens_i[0]), 2) * dWdrij_tilde;

  for (int i = 0; i < 3; ++i) acc_i[i] += fac * e_ij[i];
}

/*---------------------------------------------------------------------------*
 | evaluate modified velocity contribution                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationAdami::ModifiedVelocityContribution(
    const double* dens_i, const double* dens_j, const double* vel_i, const double* vel_j,
    const double* mod_vel_i, const double* mod_vel_j, const double& specificcoefficient_ij,
    const double* e_ij, double* acc_i) const
{
  double A_i[3][3];
  double A_j[3][3];

  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {
      A_i[i][j] = dens_i[0] * vel_i[i] * (mod_vel_i[j] - vel_i[j]);

      if (vel_j and mod_vel_j) A_j[i][j] = dens_j[0] * vel_j[i] * (mod_vel_j[j] - vel_j[j]);
    }

  double A_ij_e_ij[3];
  for (int i = 0; i < 3; ++i)
  {
    A_ij_e_ij[i] = 0.5 * (A_i[i][0] * e_ij[0] + A_i[i][1] * e_ij[1] + A_i[i][2] * e_ij[2]);

    if (vel_j and mod_vel_j)
      A_ij_e_ij[i] += 0.5 * (A_j[i][0] * e_ij[0] + A_j[i][1] * e_ij[1] + A_j[i][2] * e_ij[2]);
  }

  for (int i = 0; i < 3; ++i) acc_i[i] += specificcoefficient_ij * A_ij_e_ij[i];
}
