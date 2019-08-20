/*---------------------------------------------------------------------------*/
/*! \file
\brief momentum formulation handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_momentum_formulation.H"

#include "particle_interaction_utils.H"

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
    const double& dWdrji, double& speccoeff_ij, double& speccoeff_ji) const
{
  speccoeff_ij = (dWdrij * mass_j[0]);
  speccoeff_ji = (dWdrji * mass_i[0]);
}

/*---------------------------------------------------------------------------*
 | evaluate pressure gradient                                 sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationMonaghan::PressureGradient(const double* dens_i,
    const double* dens_j, const double* press_i, const double* press_j, const double& speccoeff_ij,
    const double& speccoeff_ji, const double* e_ij, double* acc_i, double* acc_j) const
{
  const double fac =
      (press_i[0] / UTILS::pow<2>(dens_i[0]) + press_j[0] / UTILS::pow<2>(dens_j[0]));

  if (acc_i) UTILS::vec_addscale(acc_i, -speccoeff_ij * fac, e_ij);
  if (acc_j) UTILS::vec_addscale(acc_j, speccoeff_ji * fac, e_ij);
}

/*---------------------------------------------------------------------------*
 | evaluate shear forces                                      sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationMonaghan::ShearForces(const double* dens_i,
    const double* dens_j, const double* vel_i, const double* vel_j, const double& kernelfac,
    const double& visc_i, const double& visc_j, const double& bulk_visc_i,
    const double& bulk_visc_j, const double& abs_rij, const double& speccoeff_ij,
    const double& speccoeff_ji, const double* e_ij, double* acc_i, double* acc_j) const
{
  double scaled_viscosity = 0.0;
  if (visc_i > 0.0 and visc_j > 0.0)
    scaled_viscosity = (2.0 * visc_i * visc_j / (3.0 * (visc_i + visc_j)));

  double bulk_viscosity = 0.0;
  if (bulk_visc_i > 0.0 and bulk_visc_j > 0.0)
    bulk_viscosity = (2.0 * bulk_visc_i * bulk_visc_j / (bulk_visc_i + bulk_visc_j));

  const double convection_coeff = kernelfac * (bulk_viscosity + scaled_viscosity);
  const double diffusion_coeff = 5.0 * scaled_viscosity - bulk_viscosity;

  // safety check
  if (diffusion_coeff < 0.0) dserror("diffusion coefficient is negative!");

  double vel_ij[3];
  UTILS::vec_set(vel_ij, vel_i);
  UTILS::vec_sub(vel_ij, vel_j);

  const double inv_densi_densj_absdist = 1.0 / (dens_i[0] * dens_j[0] * abs_rij);

  // diffusion
  const double fac_diff = diffusion_coeff * inv_densi_densj_absdist;
  if (acc_i) UTILS::vec_addscale(acc_i, speccoeff_ij * fac_diff, vel_ij);
  if (acc_j) UTILS::vec_addscale(acc_j, -speccoeff_ji * fac_diff, vel_ij);

  // convection
  const double fac_conv = convection_coeff * UTILS::vec_dot(vel_ij, e_ij) * inv_densi_densj_absdist;
  if (acc_i) UTILS::vec_addscale(acc_i, speccoeff_ij * fac_conv, e_ij);
  if (acc_j) UTILS::vec_addscale(acc_j, -speccoeff_ji * fac_conv, e_ij);
}

/*---------------------------------------------------------------------------*
 | evaluate background pressure (standard formulation)        sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationMonaghan::StandardBackgroundPressure(
    const double* dens_i, const double* dens_j, const double& bg_press_i, const double& bg_press_j,
    const double& speccoeff_ij, const double& speccoeff_ji, const double* e_ij, double* mod_acc_i,
    double* mod_acc_j) const
{
  const double fac = (1.0 / UTILS::pow<2>(dens_i[0]) + 1.0 / UTILS::pow<2>(dens_j[0]));

  if (mod_acc_i) UTILS::vec_addscale(mod_acc_i, -speccoeff_ij * bg_press_i * fac, e_ij);
  if (mod_acc_j) UTILS::vec_addscale(mod_acc_j, speccoeff_ji * bg_press_j * fac, e_ij);
}

/*---------------------------------------------------------------------------*
 | evaluate background pressure (generalized formulation)     sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationMonaghan::GeneralizedBackgroundPressure(
    const double* dens_i, const double* dens_j, const double* mass_i, const double* mass_j,
    const double& mod_bg_press_i, const double& mod_bg_press_j, const double& mod_dWdrij,
    const double& mod_dWdrji, const double* e_ij, double* mod_acc_i, double* mod_acc_j) const
{
  if (mod_acc_i)
    UTILS::vec_addscale(
        mod_acc_i, -mod_bg_press_i * (mass_j[0] / UTILS::pow<2>(dens_i[0])) * mod_dWdrij, e_ij);

  if (mod_acc_j)
    UTILS::vec_addscale(
        mod_acc_j, mod_bg_press_j * (mass_i[0] / UTILS::pow<2>(dens_j[0])) * mod_dWdrji, e_ij);
}

/*---------------------------------------------------------------------------*
 | evaluate modified velocity contribution                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationMonaghan::ModifiedVelocityContribution(
    const double* dens_i, const double* dens_j, const double* vel_i, const double* vel_j,
    const double* mod_vel_i, const double* mod_vel_j, const double& speccoeff_ij,
    const double& speccoeff_ji, const double* e_ij, double* acc_i, double* acc_j) const
{
  double A_ij_e_ij[3] = {0.0};

  if (mod_vel_i)
  {
    double modvel_ii[3];
    UTILS::vec_set(modvel_ii, mod_vel_i);
    UTILS::vec_sub(modvel_ii, vel_i);
    UTILS::vec_addscale(A_ij_e_ij, (UTILS::vec_dot(modvel_ii, e_ij) / dens_i[0]), vel_i);
  }

  if (mod_vel_j)
  {
    double modvel_jj[3];
    UTILS::vec_set(modvel_jj, mod_vel_j);
    UTILS::vec_sub(modvel_jj, vel_j);
    UTILS::vec_addscale(A_ij_e_ij, (UTILS::vec_dot(modvel_jj, e_ij) / dens_j[0]), vel_j);
  }

  if (acc_i) UTILS::vec_addscale(acc_i, speccoeff_ij, A_ij_e_ij);
  if (acc_j) UTILS::vec_addscale(acc_j, -speccoeff_ji, A_ij_e_ij);
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
    const double& dWdrji, double& speccoeff_ij, double& speccoeff_ji) const
{
  const double fac = (UTILS::pow<2>(mass_i[0] / dens_i[0]) + UTILS::pow<2>(mass_j[0] / dens_j[0]));

  speccoeff_ij = fac * (dWdrij / mass_i[0]);
  speccoeff_ji = fac * (dWdrji / mass_j[0]);
}

/*---------------------------------------------------------------------------*
 | evaluate pressure gradient                                 sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationAdami::PressureGradient(const double* dens_i,
    const double* dens_j, const double* press_i, const double* press_j, const double& speccoeff_ij,
    const double& speccoeff_ji, const double* e_ij, double* acc_i, double* acc_j) const
{
  const double fac = (dens_i[0] * press_j[0] + dens_j[0] * press_i[0]) / (dens_i[0] + dens_j[0]);

  if (acc_i) UTILS::vec_addscale(acc_i, -speccoeff_ij * fac, e_ij);
  if (acc_j) UTILS::vec_addscale(acc_j, speccoeff_ji * fac, e_ij);
}

/*---------------------------------------------------------------------------*
 | evaluate shear forces                                      sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationAdami::ShearForces(const double* dens_i,
    const double* dens_j, const double* vel_i, const double* vel_j, const double& kernelfac,
    const double& visc_i, const double& visc_j, const double& bulk_visc_i,
    const double& bulk_visc_j, const double& abs_rij, const double& speccoeff_ij,
    const double& speccoeff_ji, const double* e_ij, double* acc_i, double* acc_j) const
{
  double viscosity = 0.0;
  if (visc_i > 0.0 and visc_j > 0.0)
    viscosity = (2.0 * visc_i * visc_j / (visc_i + visc_j));
  else
    return;

  double vel_ij[3];
  UTILS::vec_set(vel_ij, vel_i);
  UTILS::vec_sub(vel_ij, vel_j);

  const double fac = viscosity / abs_rij;

  if (acc_i) UTILS::vec_addscale(acc_i, speccoeff_ij * fac, vel_ij);
  if (acc_j) UTILS::vec_addscale(acc_j, -speccoeff_ji * fac, vel_ij);
}

/*---------------------------------------------------------------------------*
 | evaluate background pressure (standard formulation)        sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationAdami::StandardBackgroundPressure(
    const double* dens_i, const double* dens_j, const double& bg_press_i, const double& bg_press_j,
    const double& speccoeff_ij, const double& speccoeff_ji, const double* e_ij, double* mod_acc_i,
    double* mod_acc_j) const
{
  if (mod_acc_i) UTILS::vec_addscale(mod_acc_i, -speccoeff_ij * bg_press_i, e_ij);
  if (mod_acc_j) UTILS::vec_addscale(mod_acc_j, speccoeff_ji * bg_press_j, e_ij);
}

/*---------------------------------------------------------------------------*
 | evaluate background pressure (generalized formulation)     sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationAdami::GeneralizedBackgroundPressure(
    const double* dens_i, const double* dens_j, const double* mass_i, const double* mass_j,
    const double& mod_bg_press_i, const double& mod_bg_press_j, const double& mod_dWdrij,
    const double& mod_dWdrji, const double* e_ij, double* mod_acc_i, double* mod_acc_j) const
{
  if (mod_acc_i)
    UTILS::vec_addscale(
        mod_acc_i, -(mod_bg_press_i * mass_i[0] * mod_dWdrij) / UTILS::pow<2>(dens_i[0]), e_ij);

  if (mod_acc_j)
    UTILS::vec_addscale(
        mod_acc_j, (mod_bg_press_j * mass_j[0] * mod_dWdrji) / UTILS::pow<2>(dens_j[0]), e_ij);
}

/*---------------------------------------------------------------------------*
 | evaluate modified velocity contribution                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHMomentumFormulationAdami::ModifiedVelocityContribution(
    const double* dens_i, const double* dens_j, const double* vel_i, const double* vel_j,
    const double* mod_vel_i, const double* mod_vel_j, const double& speccoeff_ij,
    const double& speccoeff_ji, const double* e_ij, double* acc_i, double* acc_j) const
{
  double A_ij_e_ij[3] = {0.0};

  if (mod_vel_i)
  {
    double modvel_ii[3];
    UTILS::vec_set(modvel_ii, mod_vel_i);
    UTILS::vec_sub(modvel_ii, vel_i);
    UTILS::vec_addscale(A_ij_e_ij, 0.5 * dens_i[0] * UTILS::vec_dot(modvel_ii, e_ij), vel_i);
  }

  if (mod_vel_j)
  {
    double modvel_jj[3];
    UTILS::vec_set(modvel_jj, mod_vel_j);
    UTILS::vec_sub(modvel_jj, vel_j);
    UTILS::vec_addscale(A_ij_e_ij, 0.5 * dens_j[0] * UTILS::vec_dot(modvel_jj, e_ij), vel_j);
  }

  if (acc_i) UTILS::vec_addscale(acc_i, speccoeff_ij, A_ij_e_ij);
  if (acc_j) UTILS::vec_addscale(acc_j, -speccoeff_ji, A_ij_e_ij);
}
