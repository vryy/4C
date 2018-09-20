/*---------------------------------------------------------------------------*/
/*!
\file unit_particle_interaction_sph_momentum_formulation.H

\brief unittests for momentum formulation handler in smoothed particle hydrodynamics (SPH)

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#ifndef UNIT_PARTICLE_INTERACTION_SPH_MOMENTUM_FORMULATION_H
#define UNIT_PARTICLE_INTERACTION_SPH_MOMENTUM_FORMULATION_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#include <cxxtest/TestSuite.h>
#include "src/drt_particle_interaction/particle_interaction_sph_momentum_formulation.H"

/*---------------------------------------------------------------------------*
 | forward declaration                                        sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHMomentumFormulationMonaghan_TestSuite;
  class SPHMomentumFormulationAdami_TestSuite;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | momentum formulation handler test suite                    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
class PARTICLEINTERACTION::SPHMomentumFormulationMonaghan_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::SPHMomentumFormulationMonaghan> momentumformulation_;

 public:
  void setUp()
  {
    // create momentum formulation handler
    momentumformulation_ = std::unique_ptr<PARTICLEINTERACTION::SPHMomentumFormulationMonaghan>(
        new PARTICLEINTERACTION::SPHMomentumFormulationMonaghan());

    // init momentum formulation handler
    momentumformulation_->Init();

    // setup momentum formulation handler
    momentumformulation_->Setup();
  }

  void tearDown() { momentumformulation_ = nullptr; }

  // note: the public functions Init() and Setup() of class SPHMomentumFormulationMonaghan are
  // called in setUp()
  //       and thus implicitly tested by all following unittests

  void test_SpecificCoefficient()
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double mass_i = 0.25;
    const double mass_j = 0.23;
    const double dWdrij = 0.3;
    double specificcoefficient_ij = 0.0;

    momentumformulation_->SpecificCoefficient(
        &dens_i, &dens_j, &mass_i, &mass_j, dWdrij, specificcoefficient_ij);

    // compute reference solution
    double specificcoefficient_ij_ref = dWdrij * mass_j;

    // compare results
    TS_ASSERT_EQUALS(specificcoefficient_ij, specificcoefficient_ij_ref);
  }

  void test_PressureGradient()
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double press_i = 0.45;
    const double press_j = 0.43;
    const double specificcoefficient_ij = 12.42;

    double e_ij[3];
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);

    double acc_i[3];
    acc_i[0] = 0.0;
    acc_i[1] = 0.0;
    acc_i[2] = 0.0;

    momentumformulation_->PressureGradient(
        &dens_i, &dens_j, &press_i, &press_j, specificcoefficient_ij, e_ij, acc_i);

    // compute reference solution
    const double pressuregradient =
        -specificcoefficient_ij * (press_i / std::pow(dens_i, 2) + press_j / std::pow(dens_j, 2));

    double acc_i_ref[3];
    acc_i_ref[0] = 0.0;
    acc_i_ref[1] = 0.0;
    acc_i_ref[2] = 0.0;
    for (int i = 0; i < 3; ++i) acc_i_ref[i] = pressuregradient * e_ij[i];

    // compare results
    for (int i = 0; i < 3; ++i) TS_ASSERT_EQUALS(acc_i[i], acc_i_ref[i]);
  }

  void test_ShearForces()
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;

    double vel_i[3];
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3];
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;

    const double kernelfac = 3;

    const double visc_i = 0.03;
    const double visc_j = 0.05;
    const double bulk_visc_i = 0.015;
    const double bulk_visc_j = 0.006;

    const double abs_rij = 0.3;
    const double specificcoefficient_ij = 12.42;

    double e_ij[3];
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);

    double acc_i[3];
    acc_i[0] = 0.0;
    acc_i[1] = 0.0;
    acc_i[2] = 0.0;

    momentumformulation_->ShearForces(&dens_i, &dens_j, vel_i, vel_j, kernelfac, visc_i, visc_j,
        bulk_visc_i, bulk_visc_j, abs_rij, specificcoefficient_ij, e_ij, acc_i);

    // compute reference solution
    double viscosity = 0.0;
    if (visc_i > 0.0 and visc_j > 0.0) viscosity = (2.0 * visc_i * visc_j / (visc_i + visc_j));

    double bulkviscosity = 0.0;
    if (bulk_visc_i > 0.0 and bulk_visc_j > 0.0)
      bulkviscosity = (2.0 * bulk_visc_i * bulk_visc_j / (bulk_visc_i + bulk_visc_j));

    const double convectioncoefficient = kernelfac * (bulkviscosity + viscosity / 3.0);
    const double diffusioncoefficient = 5.0 * viscosity / 3.0 - bulkviscosity;

    const double inv_densi_densj_absdist = 1.0 / (dens_i * dens_j * abs_rij);
    const double e_ij_vrel_ij = ((vel_i[0] - vel_j[0]) * e_ij[0] + (vel_i[1] - vel_j[1]) * e_ij[1] +
                                 (vel_i[2] - vel_j[2]) * e_ij[2]);

    double acc_i_ref[3];
    acc_i_ref[0] = 0.0;
    acc_i_ref[1] = 0.0;
    acc_i_ref[2] = 0.0;
    for (int i = 0; i < 3; ++i)
    {
      // diffusion
      acc_i_ref[i] = specificcoefficient_ij * diffusioncoefficient * inv_densi_densj_absdist *
                     (vel_i[i] - vel_j[i]);

      // convection
      acc_i_ref[i] += specificcoefficient_ij * convectioncoefficient * e_ij_vrel_ij *
                      inv_densi_densj_absdist * e_ij[i];
    }

    // compare results
    for (int i = 0; i < 3; ++i) TS_ASSERT_EQUALS(acc_i[i], acc_i_ref[i]);
  }

  void test_StandardBackgroundPressure()
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;

    const double backgroundpressure = 10.0;
    const double specificcoefficient_ij = 12.42;

    double e_ij[3];
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);

    double acc_i[3];
    acc_i[0] = 0.0;
    acc_i[1] = 0.0;
    acc_i[2] = 0.0;

    momentumformulation_->StandardBackgroundPressure(
        &dens_i, &dens_j, backgroundpressure, specificcoefficient_ij, e_ij, acc_i);

    // compute reference solution
    const double fac = -specificcoefficient_ij * backgroundpressure *
                       (1.0 / std::pow(dens_i, 2) + 1.0 / std::pow(dens_j, 2));

    double acc_i_ref[3];
    acc_i_ref[0] = 0.0;
    acc_i_ref[1] = 0.0;
    acc_i_ref[2] = 0.0;
    for (int i = 0; i < 3; ++i) acc_i_ref[i] = fac * e_ij[i];

    // compare results
    for (int i = 0; i < 3; ++i) TS_ASSERT_EQUALS(acc_i[i], acc_i_ref[i]);
  }

  void test_GeneralizedBackgroundPressure()
  {
    const double dens_i = 1.02;
    const double mass_i = 0.25;
    const double mass_j = 0.23;

    const double backgroundpressure_tilde = 8.4;
    const double dWdrij_tilde = 0.32;

    double e_ij[3];
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);

    double acc_i[3];
    acc_i[0] = 0.0;
    acc_i[1] = 0.0;
    acc_i[2] = 0.0;

    momentumformulation_->GeneralizedBackgroundPressure(
        &dens_i, &mass_i, &mass_j, backgroundpressure_tilde, dWdrij_tilde, e_ij, acc_i);

    // compute reference solution
    const double fac = -backgroundpressure_tilde * (mass_j / std::pow(dens_i, 2)) * dWdrij_tilde;

    double acc_i_ref[3];
    acc_i_ref[0] = 0.0;
    acc_i_ref[1] = 0.0;
    acc_i_ref[2] = 0.0;
    for (int i = 0; i < 3; ++i) acc_i_ref[i] = fac * e_ij[i];

    // compare results
    for (int i = 0; i < 3; ++i) TS_ASSERT_EQUALS(acc_i[i], acc_i_ref[i]);
  }

  void test_ModifiedVelocityContribution()
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;

    double vel_i[3];
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3];
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;

    double mod_vel_i[3];
    mod_vel_i[0] = 0.1;
    mod_vel_i[1] = 2.0;
    mod_vel_i[2] = -0.4;
    double mod_vel_j[3];
    mod_vel_j[0] = 0.27;
    mod_vel_j[1] = -1.3;
    mod_vel_j[2] = 0.5;

    const double specificcoefficient_ij = 12.42;

    double e_ij[3];
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);

    double acc_i[3];
    acc_i[0] = 0.0;
    acc_i[1] = 0.0;
    acc_i[2] = 0.0;

    momentumformulation_->ModifiedVelocityContribution(
        &dens_i, &dens_j, vel_i, vel_j, mod_vel_i, mod_vel_j, specificcoefficient_ij, e_ij, acc_i);

    // compute reference solution
    double A_i[3][3];
    double A_j[3][3];

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
        A_i[i][j] = dens_i * vel_i[i] * (mod_vel_i[j] - vel_i[j]);

        if (vel_j and mod_vel_j) A_j[i][j] = dens_j * vel_j[i] * (mod_vel_j[j] - vel_j[j]);
      }

    double A_ij_e_ij[3];
    for (int i = 0; i < 3; ++i)
    {
      A_ij_e_ij[i] = (1.0 / std::pow(dens_i, 2)) *
                     (A_i[i][0] * e_ij[0] + A_i[i][1] * e_ij[1] + A_i[i][2] * e_ij[2]);

      if (vel_j and mod_vel_j)
        A_ij_e_ij[i] += (1.0 / std::pow(dens_j, 2)) *
                        (A_j[i][0] * e_ij[0] + A_j[i][1] * e_ij[1] + A_j[i][2] * e_ij[2]);
    }

    double acc_i_ref[3];
    acc_i_ref[0] = 0.0;
    acc_i_ref[1] = 0.0;
    acc_i_ref[2] = 0.0;
    for (int i = 0; i < 3; ++i) acc_i_ref[i] = specificcoefficient_ij * A_ij_e_ij[i];

    // compare results
    for (int i = 0; i < 3; ++i) TS_ASSERT_EQUALS(acc_i[i], acc_i_ref[i]);
  }
};

class PARTICLEINTERACTION::SPHMomentumFormulationAdami_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::SPHMomentumFormulationAdami> momentumformulation_;

 public:
  void setUp()
  {
    // create momentum formulation handler
    momentumformulation_ = std::unique_ptr<PARTICLEINTERACTION::SPHMomentumFormulationAdami>(
        new PARTICLEINTERACTION::SPHMomentumFormulationAdami());

    // init momentum formulation handler
    momentumformulation_->Init();

    // setup momentum formulation handler
    momentumformulation_->Setup();
  }

  void tearDown() { momentumformulation_ = nullptr; }

  // note: the public functions Init() and Setup() of class SPHMomentumFormulationAdami are called
  // in setUp()
  //       and thus implicitly tested by all following unittests

  void test_SpecificCoefficient()
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double mass_i = 0.25;
    const double mass_j = 0.23;
    const double dWdrij = 0.3;
    double specificcoefficient_ij = 0.0;

    momentumformulation_->SpecificCoefficient(
        &dens_i, &dens_j, &mass_i, &mass_j, dWdrij, specificcoefficient_ij);

    // compute reference solution
    double specificcoefficient_ij_ref =
        (std::pow((mass_i / dens_i), 2) + std::pow((mass_j / dens_j), 2)) * (dWdrij / mass_i);

    // compare results
    TS_ASSERT_EQUALS(specificcoefficient_ij, specificcoefficient_ij_ref);
  }

  void test_PressureGradient()
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double press_i = 0.45;
    const double press_j = 0.43;
    const double specificcoefficient_ij = 12.42;

    double e_ij[3];
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);

    double acc_i[3];
    acc_i[0] = 0.0;
    acc_i[1] = 0.0;
    acc_i[2] = 0.0;

    momentumformulation_->PressureGradient(
        &dens_i, &dens_j, &press_i, &press_j, specificcoefficient_ij, e_ij, acc_i);

    // compute reference solution
    const double pressuregradient =
        -specificcoefficient_ij * (dens_i * press_j + dens_j * press_i) / (dens_i + dens_j);

    double acc_i_ref[3];
    acc_i_ref[0] = 0.0;
    acc_i_ref[1] = 0.0;
    acc_i_ref[2] = 0.0;
    for (int i = 0; i < 3; ++i) acc_i_ref[i] = pressuregradient * e_ij[i];

    // compare results
    for (int i = 0; i < 3; ++i) TS_ASSERT_EQUALS(acc_i[i], acc_i_ref[i]);
  }

  void test_ShearForces()
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;

    double vel_i[3];
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3];
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;

    const double kernelfac = 3;

    const double visc_i = 0.03;
    const double visc_j = 0.05;
    const double bulk_visc_i = 0.015;
    const double bulk_visc_j = 0.006;

    const double abs_rij = 0.3;
    const double specificcoefficient_ij = 12.42;

    double e_ij[3];
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);

    double acc_i[3];
    acc_i[0] = 0.0;
    acc_i[1] = 0.0;
    acc_i[2] = 0.0;

    momentumformulation_->ShearForces(&dens_i, &dens_j, vel_i, vel_j, kernelfac, visc_i, visc_j,
        bulk_visc_i, bulk_visc_j, abs_rij, specificcoefficient_ij, e_ij, acc_i);

    // compute reference solution
    double viscosity = 0.0;
    if (visc_i > 0.0 and visc_j > 0.0) viscosity = (2.0 * visc_i * visc_j / (visc_i + visc_j));

    const double laminarviscosity = specificcoefficient_ij * viscosity / abs_rij;

    double acc_i_ref[3];
    acc_i_ref[0] = 0.0;
    acc_i_ref[1] = 0.0;
    acc_i_ref[2] = 0.0;
    for (int i = 0; i < 3; ++i) acc_i_ref[i] = laminarviscosity * (vel_i[i] - vel_j[i]);

    // compare results
    for (int i = 0; i < 3; ++i) TS_ASSERT_EQUALS(acc_i[i], acc_i_ref[i]);
  }

  void test_StandardBackgroundPressure()
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;

    const double backgroundpressure = 10.0;
    const double specificcoefficient_ij = 12.42;

    double e_ij[3];
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);

    double acc_i[3];
    acc_i[0] = 0.0;
    acc_i[1] = 0.0;
    acc_i[2] = 0.0;

    momentumformulation_->StandardBackgroundPressure(
        &dens_i, &dens_j, backgroundpressure, specificcoefficient_ij, e_ij, acc_i);

    // compute reference solution
    const double fac = -specificcoefficient_ij * backgroundpressure;

    double acc_i_ref[3];
    acc_i_ref[0] = 0.0;
    acc_i_ref[1] = 0.0;
    acc_i_ref[2] = 0.0;
    for (int i = 0; i < 3; ++i) acc_i_ref[i] = fac * e_ij[i];

    // compare results
    for (int i = 0; i < 3; ++i) TS_ASSERT_EQUALS(acc_i[i], acc_i_ref[i]);
  }

  void test_GeneralizedBackgroundPressure()
  {
    const double dens_i = 1.02;
    const double mass_i = 0.25;
    const double mass_j = 0.23;

    const double backgroundpressure_tilde = 8.4;
    const double dWdrij_tilde = 0.32;

    double e_ij[3];
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);

    double acc_i[3];
    acc_i[0] = 0.0;
    acc_i[1] = 0.0;
    acc_i[2] = 0.0;

    momentumformulation_->GeneralizedBackgroundPressure(
        &dens_i, &mass_i, &mass_j, backgroundpressure_tilde, dWdrij_tilde, e_ij, acc_i);

    // compute reference solution
    const double fac =
        -(backgroundpressure_tilde / mass_i) * std::pow((mass_i / dens_i), 2) * dWdrij_tilde;

    double acc_i_ref[3];
    acc_i_ref[0] = 0.0;
    acc_i_ref[1] = 0.0;
    acc_i_ref[2] = 0.0;
    for (int i = 0; i < 3; ++i) acc_i_ref[i] = fac * e_ij[i];

    // compare results
    for (int i = 0; i < 3; ++i) TS_ASSERT_EQUALS(acc_i[i], acc_i_ref[i]);
  }

  void test_ModifiedVelocityContribution()
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;

    double vel_i[3];
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3];
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;

    double mod_vel_i[3];
    mod_vel_i[0] = 0.1;
    mod_vel_i[1] = 2.0;
    mod_vel_i[2] = -0.4;
    double mod_vel_j[3];
    mod_vel_j[0] = 0.27;
    mod_vel_j[1] = -1.3;
    mod_vel_j[2] = 0.5;

    const double specificcoefficient_ij = 12.42;

    double e_ij[3];
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);

    double acc_i[3];
    acc_i[0] = 0.0;
    acc_i[1] = 0.0;
    acc_i[2] = 0.0;

    momentumformulation_->ModifiedVelocityContribution(
        &dens_i, &dens_j, vel_i, vel_j, mod_vel_i, mod_vel_j, specificcoefficient_ij, e_ij, acc_i);

    // compute reference solution
    double A_i[3][3];
    double A_j[3][3];

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
        A_i[i][j] = dens_i * vel_i[i] * (mod_vel_i[j] - vel_i[j]);

        if (vel_j and mod_vel_j) A_j[i][j] = dens_j * vel_j[i] * (mod_vel_j[j] - vel_j[j]);
      }

    double A_ij_e_ij[3];
    for (int i = 0; i < 3; ++i)
    {
      A_ij_e_ij[i] = 0.5 * (A_i[i][0] * e_ij[0] + A_i[i][1] * e_ij[1] + A_i[i][2] * e_ij[2]);

      if (vel_j and mod_vel_j)
        A_ij_e_ij[i] += 0.5 * (A_j[i][0] * e_ij[0] + A_j[i][1] * e_ij[1] + A_j[i][2] * e_ij[2]);
    }

    double acc_i_ref[3];
    acc_i_ref[0] = 0.0;
    acc_i_ref[1] = 0.0;
    acc_i_ref[2] = 0.0;
    for (int i = 0; i < 3; ++i) acc_i_ref[i] = specificcoefficient_ij * A_ij_e_ij[i];

    // compare results
    for (int i = 0; i < 3; ++i) TS_ASSERT_EQUALS(acc_i[i], acc_i_ref[i]);
  }
};

/*---------------------------------------------------------------------------*/
#endif
