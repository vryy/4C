/*---------------------------------------------------------------------------*/
/*! \file
\brief unittests for momentum formulation handler in smoothed particle hydrodynamics (SPH)
\level 3
*/
/*---------------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "4C_particle_interaction_sph_momentum_formulation.hpp"

#include "4C_unittest_utils_assertions_test.hpp"

#include <cmath>


namespace
{
  using namespace FourC;

  class SPHMomentumFormulationMonaghanTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<ParticleInteraction::SPHMomentumFormulationMonaghan> momentumformulation_;

    SPHMomentumFormulationMonaghanTest()
    {
      // create momentum formulation handler
      momentumformulation_ =
          std::make_unique<ParticleInteraction::SPHMomentumFormulationMonaghan>();

      // init momentum formulation handler
      momentumformulation_->init();

      // setup momentum formulation handler
      momentumformulation_->setup();
    }
    // note: the public functions init() and setup() of class SPHMomentumFormulationMonaghan are
    // called in the constructor and thus implicitly tested by all following unittests
  };
  TEST_F(SPHMomentumFormulationMonaghanTest, SpecificCoefficient)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double mass_i = 0.25;
    const double mass_j = 0.23;
    const double dWdrij = 0.3;
    const double dWdrji = 0.85;

    double speccoeff_ij = 0.0;
    double speccoeff_ji = 0.0;

    momentumformulation_->SpecificCoefficient(
        &dens_i, &dens_j, &mass_i, &mass_j, dWdrij, dWdrji, &speccoeff_ij, &speccoeff_ji);

    // compute reference solution
    const double speccoeff_ij_ref = dWdrij * mass_j;
    const double speccoeff_ji_ref = dWdrji * mass_i;

    // compare results
    EXPECT_NEAR(speccoeff_ij, speccoeff_ij_ref, 1.0e-14);
    EXPECT_NEAR(speccoeff_ji, speccoeff_ji_ref, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationMonaghanTest, PressureGradient)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double press_i = 0.45;
    const double press_j = 0.43;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;

    double acc_i[3] = {};
    double acc_j[3] = {};

    momentumformulation_->PressureGradient(
        &dens_i, &dens_j, &press_i, &press_j, speccoeff_ij, speccoeff_ji, e_ij, acc_i, acc_j);

    // compute reference solution
    double acc_i_ref[3] = {};
    double acc_j_ref[3] = {};

    const double fac = (press_i / std::pow(dens_i, 2) + press_j / std::pow(dens_j, 2));
    for (int i = 0; i < 3; ++i) acc_i_ref[i] = -speccoeff_ij * fac * e_ij[i];
    for (int i = 0; i < 3; ++i) acc_j_ref[i] = speccoeff_ji * fac * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_i, acc_i_ref, 3, 1.0e-14);
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_j, acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationMonaghanTest, PressureGradientNullptrAccI)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double press_i = 0.45;
    const double press_j = 0.43;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;

    double acc_j[3] = {};

    momentumformulation_->PressureGradient(
        &dens_i, &dens_j, &press_i, &press_j, speccoeff_ij, speccoeff_ji, e_ij, nullptr, acc_j);

    // compute reference solution
    double acc_j_ref[3] = {};

    const double fac = (press_i / std::pow(dens_i, 2) + press_j / std::pow(dens_j, 2));
    for (int i = 0; i < 3; ++i) acc_j_ref[i] = speccoeff_ji * fac * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_j, acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationMonaghanTest, PressureGradientNullptrAccJ)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double press_i = 0.45;
    const double press_j = 0.43;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;

    double acc_i[3] = {};

    momentumformulation_->PressureGradient(
        &dens_i, &dens_j, &press_i, &press_j, speccoeff_ij, speccoeff_ji, e_ij, acc_i, nullptr);

    // compute reference solution
    double acc_i_ref[3] = {};

    const double fac = (press_i / std::pow(dens_i, 2) + press_j / std::pow(dens_j, 2));
    for (int i = 0; i < 3; ++i) acc_i_ref[i] = -speccoeff_ij * fac * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_i, acc_i_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationMonaghanTest, ShearForces)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double vel_i[3] = {};
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3] = {};
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;
    const double abs_rij = 0.3;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;
    const double visc_i = 0.03;
    const double visc_j = 0.05;
    const double bulk_visc_i = 0.015;
    const double bulk_visc_j = 0.006;

    const double kernelfac = 3;

    double acc_i[3] = {};
    double acc_j[3] = {};

    momentumformulation_->ShearForces(&dens_i, &dens_j, vel_i, vel_j, kernelfac, visc_i, visc_j,
        bulk_visc_i, bulk_visc_j, abs_rij, speccoeff_ij, speccoeff_ji, e_ij, acc_i, acc_j);

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

    double acc_i_ref[3] = {};
    double acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i)
    {
      // diffusion
      acc_i_ref[i] =
          speccoeff_ij * diffusioncoefficient * inv_densi_densj_absdist * (vel_i[i] - vel_j[i]);
      acc_j_ref[i] =
          -speccoeff_ji * diffusioncoefficient * inv_densi_densj_absdist * (vel_i[i] - vel_j[i]);

      // convection
      acc_i_ref[i] +=
          speccoeff_ij * convectioncoefficient * e_ij_vrel_ij * inv_densi_densj_absdist * e_ij[i];
      acc_j_ref[i] +=
          -speccoeff_ji * convectioncoefficient * e_ij_vrel_ij * inv_densi_densj_absdist * e_ij[i];
    }

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_i, acc_i_ref, 3, 1.0e-14);
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_j, acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationMonaghanTest, ShearForcesNullptrAccI)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double vel_i[3] = {};
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3] = {};
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;
    const double abs_rij = 0.3;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;
    const double visc_i = 0.03;
    const double visc_j = 0.05;
    const double bulk_visc_i = 0.015;
    const double bulk_visc_j = 0.006;

    const double kernelfac = 3;

    double acc_j[3] = {};

    momentumformulation_->ShearForces(&dens_i, &dens_j, vel_i, vel_j, kernelfac, visc_i, visc_j,
        bulk_visc_i, bulk_visc_j, abs_rij, speccoeff_ij, speccoeff_ji, e_ij, nullptr, acc_j);

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

    double acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i)
    {
      // diffusion
      acc_j_ref[i] =
          -speccoeff_ji * diffusioncoefficient * inv_densi_densj_absdist * (vel_i[i] - vel_j[i]);

      // convection
      acc_j_ref[i] +=
          -speccoeff_ji * convectioncoefficient * e_ij_vrel_ij * inv_densi_densj_absdist * e_ij[i];
    }

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_j, acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationMonaghanTest, ShearForcesNullptrAccJ)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double vel_i[3] = {};
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3] = {};
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;
    const double abs_rij = 0.3;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;
    const double visc_i = 0.03;
    const double visc_j = 0.05;
    const double bulk_visc_i = 0.015;
    const double bulk_visc_j = 0.006;

    const double kernelfac = 3;

    double acc_i[3] = {};

    momentumformulation_->ShearForces(&dens_i, &dens_j, vel_i, vel_j, kernelfac, visc_i, visc_j,
        bulk_visc_i, bulk_visc_j, abs_rij, speccoeff_ij, speccoeff_ji, e_ij, acc_i, nullptr);

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

    double acc_i_ref[3] = {};

    for (int i = 0; i < 3; ++i)
    {
      // diffusion
      acc_i_ref[i] =
          speccoeff_ij * diffusioncoefficient * inv_densi_densj_absdist * (vel_i[i] - vel_j[i]);

      // convection
      acc_i_ref[i] +=
          speccoeff_ij * convectioncoefficient * e_ij_vrel_ij * inv_densi_densj_absdist * e_ij[i];
    }

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_i, acc_i_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationMonaghanTest, standard_background_pressure)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;
    const double bg_press_i = 10.0;
    const double bg_press_j = 4.5;

    double mod_acc_i[3] = {};
    double mod_acc_j[3] = {};

    momentumformulation_->standard_background_pressure(&dens_i, &dens_j, bg_press_i, bg_press_j,
        speccoeff_ij, speccoeff_ji, e_ij, mod_acc_i, mod_acc_j);

    // compute reference solution
    const double fac = (1.0 / std::pow(dens_i, 2) + 1.0 / std::pow(dens_j, 2));

    double mod_acc_i_ref[3] = {};
    double mod_acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i) mod_acc_i_ref[i] = -speccoeff_ij * bg_press_i * fac * e_ij[i];
    for (int i = 0; i < 3; ++i) mod_acc_j_ref[i] = speccoeff_ji * bg_press_j * fac * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(mod_acc_i, mod_acc_i_ref, 3, 1.0e-14);
    FOUR_C_EXPECT_ITERABLE_NEAR(mod_acc_j, mod_acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationMonaghanTest, StandardBackgroundPressureNullptrModAccI)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;
    const double bg_press_i = 10.0;
    const double bg_press_j = 4.5;

    double mod_acc_j[3] = {};

    momentumformulation_->standard_background_pressure(&dens_i, &dens_j, bg_press_i, bg_press_j,
        speccoeff_ij, speccoeff_ji, e_ij, nullptr, mod_acc_j);

    // compute reference solution
    const double fac = (1.0 / std::pow(dens_i, 2) + 1.0 / std::pow(dens_j, 2));

    double mod_acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i) mod_acc_j_ref[i] = speccoeff_ji * bg_press_j * fac * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(mod_acc_j, mod_acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationMonaghanTest, StandardBackgroundPressureNullptrModAccJ)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;
    const double bg_press_i = 10.0;
    const double bg_press_j = 4.5;

    double mod_acc_i[3] = {};

    momentumformulation_->standard_background_pressure(&dens_i, &dens_j, bg_press_i, bg_press_j,
        speccoeff_ij, speccoeff_ji, e_ij, mod_acc_i, nullptr);

    // compute reference solution
    const double fac = (1.0 / std::pow(dens_i, 2) + 1.0 / std::pow(dens_j, 2));

    double mod_acc_i_ref[3] = {};

    for (int i = 0; i < 3; ++i) mod_acc_i_ref[i] = -speccoeff_ij * bg_press_i * fac * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(mod_acc_i, mod_acc_i_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationMonaghanTest, generalized_background_pressure)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double mass_i = 0.25;
    const double mass_j = 0.23;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double mod_bg_press_i = 8.4;
    const double mod_bg_press_j = 4.5;
    const double mod_dWdrij = 0.32;
    const double mod_dWdrji = 0.97;

    double mod_acc_i[3] = {};
    double mod_acc_j[3] = {};

    momentumformulation_->generalized_background_pressure(&dens_i, &dens_j, &mass_i, &mass_j,
        mod_bg_press_i, mod_bg_press_j, mod_dWdrij, mod_dWdrji, e_ij, mod_acc_i, mod_acc_j);

    // compute reference solution
    double mod_acc_i_ref[3] = {};
    double mod_acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i)
      mod_acc_i_ref[i] = -mod_bg_press_i * (mass_j / std::pow(dens_i, 2)) * mod_dWdrij * e_ij[i];

    for (int i = 0; i < 3; ++i)
      mod_acc_j_ref[i] = mod_bg_press_j * (mass_i / std::pow(dens_j, 2)) * mod_dWdrji * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(mod_acc_i, mod_acc_i_ref, 3, 1.0e-14);
    FOUR_C_EXPECT_ITERABLE_NEAR(mod_acc_j, mod_acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationMonaghanTest, GeneralizedBackgroundPressureNullptrModAccI)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double mass_i = 0.25;
    const double mass_j = 0.23;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double mod_bg_press_i = 8.4;
    const double mod_bg_press_j = 4.5;
    const double mod_dWdrij = 0.32;
    const double mod_dWdrji = 0.97;

    double mod_acc_j[3] = {};

    momentumformulation_->generalized_background_pressure(&dens_i, &dens_j, &mass_i, &mass_j,
        mod_bg_press_i, mod_bg_press_j, mod_dWdrij, mod_dWdrji, e_ij, nullptr, mod_acc_j);

    // compute reference solution
    double mod_acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i)
      mod_acc_j_ref[i] = mod_bg_press_j * (mass_i / std::pow(dens_j, 2)) * mod_dWdrji * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(mod_acc_j, mod_acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationMonaghanTest, GeneralizedBackgroundPressureNullptrModAccJ)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double mass_i = 0.25;
    const double mass_j = 0.23;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double mod_bg_press_i = 8.4;
    const double mod_bg_press_j = 4.5;
    const double mod_dWdrij = 0.32;
    const double mod_dWdrji = 0.97;

    double mod_acc_i[3] = {};

    momentumformulation_->generalized_background_pressure(&dens_i, &dens_j, &mass_i, &mass_j,
        mod_bg_press_i, mod_bg_press_j, mod_dWdrij, mod_dWdrji, e_ij, mod_acc_i, nullptr);

    // compute reference solution
    double mod_acc_i_ref[3] = {};

    for (int i = 0; i < 3; ++i)
      mod_acc_i_ref[i] = -mod_bg_press_i * (mass_j / std::pow(dens_i, 2)) * mod_dWdrij * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(mod_acc_i, mod_acc_i_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationMonaghanTest, modified_velocity_contribution)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double vel_i[3] = {};
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3] = {};
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;
    double mod_vel_i[3] = {};
    mod_vel_i[0] = 0.12;
    mod_vel_i[1] = 1.9;
    mod_vel_i[2] = -0.42;
    double mod_vel_j[3] = {};
    mod_vel_j[0] = 0.3;
    mod_vel_j[1] = -1.2;
    mod_vel_j[2] = 0.77;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;

    double acc_i[3] = {};
    double acc_j[3] = {};

    momentumformulation_->modified_velocity_contribution(&dens_i, &dens_j, vel_i, vel_j, mod_vel_i,
        mod_vel_j, speccoeff_ij, speccoeff_ji, e_ij, acc_i, acc_j);

    // compute reference solution
    double A_i[3][3] = {};
    double A_j[3][3] = {};

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
        A_i[i][j] = dens_i * vel_i[i] * (mod_vel_i[j] - vel_i[j]);

        A_j[i][j] = dens_j * vel_j[i] * (mod_vel_j[j] - vel_j[j]);
      }

    double A_ij_e_ij[3] = {};
    for (int i = 0; i < 3; ++i)
    {
      A_ij_e_ij[i] = (1.0 / std::pow(dens_i, 2)) *
                     (A_i[i][0] * e_ij[0] + A_i[i][1] * e_ij[1] + A_i[i][2] * e_ij[2]);

      A_ij_e_ij[i] += (1.0 / std::pow(dens_j, 2)) *
                      (A_j[i][0] * e_ij[0] + A_j[i][1] * e_ij[1] + A_j[i][2] * e_ij[2]);
    }

    double acc_i_ref[3] = {};
    double acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i) acc_i_ref[i] = speccoeff_ij * A_ij_e_ij[i];
    for (int i = 0; i < 3; ++i) acc_j_ref[i] = -speccoeff_ji * A_ij_e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_i, acc_i_ref, 3, 1.0e-14);
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_j, acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationMonaghanTest, ModifiedVelocityContributionNullptrAccI)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double vel_i[3] = {};
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3] = {};
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;
    double mod_vel_i[3] = {};
    mod_vel_i[0] = 0.12;
    mod_vel_i[1] = 1.9;
    mod_vel_i[2] = -0.42;
    double mod_vel_j[3] = {};
    mod_vel_j[0] = 0.3;
    mod_vel_j[1] = -1.2;
    mod_vel_j[2] = 0.77;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;

    double acc_j[3] = {};

    momentumformulation_->modified_velocity_contribution(&dens_i, &dens_j, vel_i, vel_j, mod_vel_i,
        mod_vel_j, speccoeff_ij, speccoeff_ji, e_ij, nullptr, acc_j);

    // compute reference solution
    double A_i[3][3] = {};
    double A_j[3][3] = {};

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
        A_i[i][j] = dens_i * vel_i[i] * (mod_vel_i[j] - vel_i[j]);

        A_j[i][j] = dens_j * vel_j[i] * (mod_vel_j[j] - vel_j[j]);
      }

    double A_ij_e_ij[3] = {};
    for (int i = 0; i < 3; ++i)
    {
      A_ij_e_ij[i] = (1.0 / std::pow(dens_i, 2)) *
                     (A_i[i][0] * e_ij[0] + A_i[i][1] * e_ij[1] + A_i[i][2] * e_ij[2]);

      A_ij_e_ij[i] += (1.0 / std::pow(dens_j, 2)) *
                      (A_j[i][0] * e_ij[0] + A_j[i][1] * e_ij[1] + A_j[i][2] * e_ij[2]);
    }

    double acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i) acc_j_ref[i] = -speccoeff_ji * A_ij_e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_j, acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationMonaghanTest, ModifiedVelocityContributionNullptrAccJ)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double vel_i[3] = {};
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3] = {};
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;
    double mod_vel_i[3] = {};
    mod_vel_i[0] = 0.12;
    mod_vel_i[1] = 1.9;
    mod_vel_i[2] = -0.42;
    double mod_vel_j[3] = {};
    mod_vel_j[0] = 0.3;
    mod_vel_j[1] = -1.2;
    mod_vel_j[2] = 0.77;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;

    double acc_i[3] = {};

    momentumformulation_->modified_velocity_contribution(&dens_i, &dens_j, vel_i, vel_j, mod_vel_i,
        mod_vel_j, speccoeff_ij, speccoeff_ji, e_ij, acc_i, nullptr);

    // compute reference solution
    double A_i[3][3] = {};
    double A_j[3][3] = {};

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
        A_i[i][j] = dens_i * vel_i[i] * (mod_vel_i[j] - vel_i[j]);

        A_j[i][j] = dens_j * vel_j[i] * (mod_vel_j[j] - vel_j[j]);
      }

    double A_ij_e_ij[3] = {};
    for (int i = 0; i < 3; ++i)
    {
      A_ij_e_ij[i] = (1.0 / std::pow(dens_i, 2)) *
                     (A_i[i][0] * e_ij[0] + A_i[i][1] * e_ij[1] + A_i[i][2] * e_ij[2]);

      A_ij_e_ij[i] += (1.0 / std::pow(dens_j, 2)) *
                      (A_j[i][0] * e_ij[0] + A_j[i][1] * e_ij[1] + A_j[i][2] * e_ij[2]);
    }

    double acc_i_ref[3] = {};

    for (int i = 0; i < 3; ++i) acc_i_ref[i] = speccoeff_ij * A_ij_e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_i, acc_i_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationMonaghanTest, ModifiedVelocityContributionNullptrAccIModVelI)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double vel_i[3] = {};
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3] = {};
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;
    double mod_vel_j[3] = {};
    mod_vel_j[0] = 0.3;
    mod_vel_j[1] = -1.2;
    mod_vel_j[2] = 0.77;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;

    double acc_j[3] = {};

    momentumformulation_->modified_velocity_contribution(&dens_i, &dens_j, vel_i, vel_j, nullptr,
        mod_vel_j, speccoeff_ij, speccoeff_ji, e_ij, nullptr, acc_j);

    // compute reference solution
    double A_j[3][3] = {};

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) A_j[i][j] = dens_j * vel_j[i] * (mod_vel_j[j] - vel_j[j]);

    double A_ij_e_ij[3] = {};
    for (int i = 0; i < 3; ++i)
      A_ij_e_ij[i] += (1.0 / std::pow(dens_j, 2)) *
                      (A_j[i][0] * e_ij[0] + A_j[i][1] * e_ij[1] + A_j[i][2] * e_ij[2]);

    double acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i) acc_j_ref[i] = -speccoeff_ji * A_ij_e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_j, acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationMonaghanTest, ModifiedVelocityContributionNullptrAccJModVelJ)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double vel_i[3] = {};
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3] = {};
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;
    double mod_vel_i[3] = {};
    mod_vel_i[0] = 0.12;
    mod_vel_i[1] = 1.9;
    mod_vel_i[2] = -0.42;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;

    double acc_i[3] = {};

    momentumformulation_->modified_velocity_contribution(&dens_i, &dens_j, vel_i, vel_j, mod_vel_i,
        nullptr, speccoeff_ij, speccoeff_ji, e_ij, acc_i, nullptr);

    // compute reference solution
    double A_i[3][3] = {};

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) A_i[i][j] = dens_i * vel_i[i] * (mod_vel_i[j] - vel_i[j]);

    double A_ij_e_ij[3] = {};
    for (int i = 0; i < 3; ++i)
      A_ij_e_ij[i] = (1.0 / std::pow(dens_i, 2)) *
                     (A_i[i][0] * e_ij[0] + A_i[i][1] * e_ij[1] + A_i[i][2] * e_ij[2]);

    double acc_i_ref[3] = {};

    for (int i = 0; i < 3; ++i) acc_i_ref[i] = speccoeff_ij * A_ij_e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_i, acc_i_ref, 3, 1.0e-14);
  }


  class SPHMomentumFormulationAdamiTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<ParticleInteraction::SPHMomentumFormulationAdami> momentumformulation_;

    SPHMomentumFormulationAdamiTest()
    {
      // create momentum formulation handler
      momentumformulation_ = std::make_unique<ParticleInteraction::SPHMomentumFormulationAdami>();

      // init momentum formulation handler
      momentumformulation_->init();

      // setup momentum formulation handler
      momentumformulation_->setup();
    }
    // note: the public functions init() and setup() of class SPHMomentumFormulationAdami are called
    // in setup() and thus implicitly tested by all following unittests
  };

  TEST_F(SPHMomentumFormulationAdamiTest, SpecificCoefficient)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double mass_i = 0.25;
    const double mass_j = 0.23;
    const double dWdrij = 0.3;
    const double dWdrji = 0.85;

    double speccoeff_ij = 0.0;
    double speccoeff_ji = 0.0;

    momentumformulation_->SpecificCoefficient(
        &dens_i, &dens_j, &mass_i, &mass_j, dWdrij, dWdrji, &speccoeff_ij, &speccoeff_ji);

    // compute reference solution
    const double fac = (std::pow((mass_i / dens_i), 2) + std::pow((mass_j / dens_j), 2));
    const double speccoeff_ij_ref = fac * (dWdrij / mass_i);
    const double speccoeff_ji_ref = fac * (dWdrji / mass_j);

    // compare results
    EXPECT_NEAR(speccoeff_ij, speccoeff_ij_ref, 1.0e-14);
    EXPECT_NEAR(speccoeff_ji, speccoeff_ji_ref, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationAdamiTest, PressureGradient)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double press_i = 0.45;
    const double press_j = 0.43;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;

    double acc_i[3] = {};
    double acc_j[3] = {};

    momentumformulation_->PressureGradient(
        &dens_i, &dens_j, &press_i, &press_j, speccoeff_ij, speccoeff_ji, e_ij, acc_i, acc_j);

    // compute reference solution
    double acc_i_ref[3] = {};
    double acc_j_ref[3] = {};

    const double fac = (dens_i * press_j + dens_j * press_i) / (dens_i + dens_j);
    for (int i = 0; i < 3; ++i) acc_i_ref[i] = -speccoeff_ij * fac * e_ij[i];
    for (int i = 0; i < 3; ++i) acc_j_ref[i] = speccoeff_ji * fac * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_i, acc_i_ref, 3, 1.0e-14);
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_j, acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationAdamiTest, PressureGradientNullptrAccI)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double press_i = 0.45;
    const double press_j = 0.43;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;

    double acc_j[3] = {};

    momentumformulation_->PressureGradient(
        &dens_i, &dens_j, &press_i, &press_j, speccoeff_ij, speccoeff_ji, e_ij, nullptr, acc_j);

    // compute reference solution
    double acc_j_ref[3] = {};

    const double fac = (dens_i * press_j + dens_j * press_i) / (dens_i + dens_j);
    for (int i = 0; i < 3; ++i) acc_j_ref[i] = speccoeff_ji * fac * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_j, acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationAdamiTest, PressureGradientNullptrAccJ)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double press_i = 0.45;
    const double press_j = 0.43;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;

    double acc_i[3] = {};

    momentumformulation_->PressureGradient(
        &dens_i, &dens_j, &press_i, &press_j, speccoeff_ij, speccoeff_ji, e_ij, acc_i, nullptr);

    // compute reference solution
    double acc_i_ref[3] = {};

    const double fac = (dens_i * press_j + dens_j * press_i) / (dens_i + dens_j);
    for (int i = 0; i < 3; ++i) acc_i_ref[i] = -speccoeff_ij * fac * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_i, acc_i_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationAdamiTest, ShearForces)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double vel_i[3] = {};
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3] = {};
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;
    const double abs_rij = 0.3;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;
    const double visc_i = 0.03;
    const double visc_j = 0.05;
    const double bulk_visc_i = 0.015;
    const double bulk_visc_j = 0.006;

    const double kernelfac = 3;

    double acc_i[3] = {};
    double acc_j[3] = {};

    momentumformulation_->ShearForces(&dens_i, &dens_j, vel_i, vel_j, kernelfac, visc_i, visc_j,
        bulk_visc_i, bulk_visc_j, abs_rij, speccoeff_ij, speccoeff_ji, e_ij, acc_i, acc_j);

    // compute reference solution
    double viscosity = 0.0;
    if (visc_i > 0.0 and visc_j > 0.0) viscosity = (2.0 * visc_i * visc_j / (visc_i + visc_j));

    const double fac = viscosity / abs_rij;

    double acc_i_ref[3] = {};
    double acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i) acc_i_ref[i] = speccoeff_ij * fac * (vel_i[i] - vel_j[i]);
    for (int i = 0; i < 3; ++i) acc_j_ref[i] = -speccoeff_ji * fac * (vel_i[i] - vel_j[i]);

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_i, acc_i_ref, 3, 1.0e-14);
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_j, acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationAdamiTest, ShearForcesNullptrAccI)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double vel_i[3] = {};
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3] = {};
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;
    const double abs_rij = 0.3;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;
    const double visc_i = 0.03;
    const double visc_j = 0.05;
    const double bulk_visc_i = 0.015;
    const double bulk_visc_j = 0.006;

    const double kernelfac = 3;

    double acc_j[3] = {};

    momentumformulation_->ShearForces(&dens_i, &dens_j, vel_i, vel_j, kernelfac, visc_i, visc_j,
        bulk_visc_i, bulk_visc_j, abs_rij, speccoeff_ij, speccoeff_ji, e_ij, nullptr, acc_j);

    // compute reference solution
    double viscosity = 0.0;
    if (visc_i > 0.0 and visc_j > 0.0) viscosity = (2.0 * visc_i * visc_j / (visc_i + visc_j));

    const double fac = viscosity / abs_rij;

    double acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i) acc_j_ref[i] = -speccoeff_ji * fac * (vel_i[i] - vel_j[i]);

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_j, acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationAdamiTest, ShearForcesNullptrAccJ)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double vel_i[3] = {};
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3] = {};
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;
    const double abs_rij = 0.3;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;
    const double visc_i = 0.03;
    const double visc_j = 0.05;
    const double bulk_visc_i = 0.015;
    const double bulk_visc_j = 0.006;

    const double kernelfac = 3;

    double acc_i[3] = {};

    momentumformulation_->ShearForces(&dens_i, &dens_j, vel_i, vel_j, kernelfac, visc_i, visc_j,
        bulk_visc_i, bulk_visc_j, abs_rij, speccoeff_ij, speccoeff_ji, e_ij, acc_i, nullptr);

    // compute reference solution
    double viscosity = 0.0;
    if (visc_i > 0.0 and visc_j > 0.0) viscosity = (2.0 * visc_i * visc_j / (visc_i + visc_j));

    const double fac = viscosity / abs_rij;

    double acc_i_ref[3] = {};

    for (int i = 0; i < 3; ++i) acc_i_ref[i] = speccoeff_ij * fac * (vel_i[i] - vel_j[i]);

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_i, acc_i_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationAdamiTest, standard_background_pressure)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;
    const double bg_press_i = 10.0;
    const double bg_press_j = 4.5;

    double mod_acc_i[3] = {};
    double mod_acc_j[3] = {};

    momentumformulation_->standard_background_pressure(&dens_i, &dens_j, bg_press_i, bg_press_j,
        speccoeff_ij, speccoeff_ji, e_ij, mod_acc_i, mod_acc_j);

    // compute reference solution
    double mod_acc_i_ref[3] = {};
    double mod_acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i) mod_acc_i_ref[i] = -speccoeff_ij * bg_press_i * e_ij[i];
    for (int i = 0; i < 3; ++i) mod_acc_j_ref[i] = speccoeff_ji * bg_press_j * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(mod_acc_i, mod_acc_i_ref, 3, 1.0e-14);
    FOUR_C_EXPECT_ITERABLE_NEAR(mod_acc_j, mod_acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationAdamiTest, StandardBackgroundPressureNullptrModAccI)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;
    const double bg_press_i = 10.0;
    const double bg_press_j = 4.5;

    double mod_acc_j[3] = {};

    momentumformulation_->standard_background_pressure(&dens_i, &dens_j, bg_press_i, bg_press_j,
        speccoeff_ij, speccoeff_ji, e_ij, nullptr, mod_acc_j);

    // compute reference solution
    double mod_acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i) mod_acc_j_ref[i] = speccoeff_ji * bg_press_j * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(mod_acc_j, mod_acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationAdamiTest, StandardBackgroundPressureNullptrModAccJ)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;
    const double bg_press_i = 10.0;
    const double bg_press_j = 4.5;

    double mod_acc_i[3] = {};

    momentumformulation_->standard_background_pressure(&dens_i, &dens_j, bg_press_i, bg_press_j,
        speccoeff_ij, speccoeff_ji, e_ij, mod_acc_i, nullptr);

    // compute reference solution
    double mod_acc_i_ref[3] = {};

    for (int i = 0; i < 3; ++i) mod_acc_i_ref[i] = -speccoeff_ij * bg_press_i * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(mod_acc_i, mod_acc_i_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationAdamiTest, generalized_background_pressure)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double mass_i = 0.25;
    const double mass_j = 0.23;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double mod_bg_press_i = 8.4;
    const double mod_bg_press_j = 4.5;
    const double mod_dWdrij = 0.32;
    const double mod_dWdrji = 0.97;

    double mod_acc_i[3] = {};
    double mod_acc_j[3] = {};

    momentumformulation_->generalized_background_pressure(&dens_i, &dens_j, &mass_i, &mass_j,
        mod_bg_press_i, mod_bg_press_j, mod_dWdrij, mod_dWdrji, e_ij, mod_acc_i, mod_acc_j);

    // compute reference solution
    double mod_acc_i_ref[3] = {};
    double mod_acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i)
      mod_acc_i_ref[i] =
          -(mod_bg_press_i / mass_i) * std::pow((mass_i / dens_i), 2) * mod_dWdrij * e_ij[i];

    for (int i = 0; i < 3; ++i)
      mod_acc_j_ref[i] =
          (mod_bg_press_j / mass_j) * std::pow((mass_j / dens_j), 2) * mod_dWdrji * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(mod_acc_i, mod_acc_i_ref, 3, 1.0e-14);
    FOUR_C_EXPECT_ITERABLE_NEAR(mod_acc_j, mod_acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationAdamiTest, GeneralizedBackgroundPressureNullptrModAccI)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double mass_i = 0.25;
    const double mass_j = 0.23;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double mod_bg_press_i = 8.4;
    const double mod_bg_press_j = 4.5;
    const double mod_dWdrij = 0.32;
    const double mod_dWdrji = 0.97;

    double mod_acc_j[3] = {};

    momentumformulation_->generalized_background_pressure(&dens_i, &dens_j, &mass_i, &mass_j,
        mod_bg_press_i, mod_bg_press_j, mod_dWdrij, mod_dWdrji, e_ij, nullptr, mod_acc_j);

    // compute reference solution
    double mod_acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i)
      mod_acc_j_ref[i] =
          (mod_bg_press_j / mass_j) * std::pow((mass_j / dens_j), 2) * mod_dWdrji * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(mod_acc_j, mod_acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationAdamiTest, GeneralizedBackgroundPressureNullptrModAccJ)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    const double mass_i = 0.25;
    const double mass_j = 0.23;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double mod_bg_press_i = 8.4;
    const double mod_bg_press_j = 4.5;
    const double mod_dWdrij = 0.32;
    const double mod_dWdrji = 0.97;

    double mod_acc_i[3] = {};

    momentumformulation_->generalized_background_pressure(&dens_i, &dens_j, &mass_i, &mass_j,
        mod_bg_press_i, mod_bg_press_j, mod_dWdrij, mod_dWdrji, e_ij, mod_acc_i, nullptr);

    // compute reference solution
    double mod_acc_i_ref[3] = {};

    for (int i = 0; i < 3; ++i)
      mod_acc_i_ref[i] =
          -(mod_bg_press_i / mass_i) * std::pow((mass_i / dens_i), 2) * mod_dWdrij * e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(mod_acc_i, mod_acc_i_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationAdamiTest, modified_velocity_contribution)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double vel_i[3] = {};
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3] = {};
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;
    double mod_vel_i[3] = {};
    mod_vel_i[0] = 0.12;
    mod_vel_i[1] = 1.9;
    mod_vel_i[2] = -0.42;
    double mod_vel_j[3] = {};
    mod_vel_j[0] = 0.3;
    mod_vel_j[1] = -1.2;
    mod_vel_j[2] = 0.77;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;

    double acc_i[3] = {};
    double acc_j[3] = {};

    momentumformulation_->modified_velocity_contribution(&dens_i, &dens_j, vel_i, vel_j, mod_vel_i,
        mod_vel_j, speccoeff_ij, speccoeff_ji, e_ij, acc_i, acc_j);

    // compute reference solution
    double A_i[3][3] = {};
    double A_j[3][3] = {};

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
        A_i[i][j] = dens_i * vel_i[i] * (mod_vel_i[j] - vel_i[j]);

        A_j[i][j] = dens_j * vel_j[i] * (mod_vel_j[j] - vel_j[j]);
      }

    double A_ij_e_ij[3] = {};
    for (int i = 0; i < 3; ++i)
    {
      A_ij_e_ij[i] = 0.5 * (A_i[i][0] * e_ij[0] + A_i[i][1] * e_ij[1] + A_i[i][2] * e_ij[2]);

      A_ij_e_ij[i] += 0.5 * (A_j[i][0] * e_ij[0] + A_j[i][1] * e_ij[1] + A_j[i][2] * e_ij[2]);
    }

    double acc_i_ref[3] = {};
    double acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i) acc_i_ref[i] = speccoeff_ij * A_ij_e_ij[i];
    for (int i = 0; i < 3; ++i) acc_j_ref[i] = -speccoeff_ji * A_ij_e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_i, acc_i_ref, 3, 1.0e-14);
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_j, acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationAdamiTest, ModifiedVelocityContributionNullptrAccI)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double vel_i[3] = {};
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3] = {};
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;
    double mod_vel_i[3] = {};
    mod_vel_i[0] = 0.12;
    mod_vel_i[1] = 1.9;
    mod_vel_i[2] = -0.42;
    double mod_vel_j[3] = {};
    mod_vel_j[0] = 0.3;
    mod_vel_j[1] = -1.2;
    mod_vel_j[2] = 0.77;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;

    double acc_j[3] = {};

    momentumformulation_->modified_velocity_contribution(&dens_i, &dens_j, vel_i, vel_j, mod_vel_i,
        mod_vel_j, speccoeff_ij, speccoeff_ji, e_ij, nullptr, acc_j);

    // compute reference solution
    double A_i[3][3] = {};
    double A_j[3][3] = {};

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
        A_i[i][j] = dens_i * vel_i[i] * (mod_vel_i[j] - vel_i[j]);

        A_j[i][j] = dens_j * vel_j[i] * (mod_vel_j[j] - vel_j[j]);
      }

    double A_ij_e_ij[3] = {};
    for (int i = 0; i < 3; ++i)
    {
      A_ij_e_ij[i] = 0.5 * (A_i[i][0] * e_ij[0] + A_i[i][1] * e_ij[1] + A_i[i][2] * e_ij[2]);

      A_ij_e_ij[i] += 0.5 * (A_j[i][0] * e_ij[0] + A_j[i][1] * e_ij[1] + A_j[i][2] * e_ij[2]);
    }

    double acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i) acc_j_ref[i] = -speccoeff_ji * A_ij_e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_j, acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationAdamiTest, ModifiedVelocityContributionNullptrAccJ)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double vel_i[3] = {};
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3] = {};
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;
    double mod_vel_i[3] = {};
    mod_vel_i[0] = 0.12;
    mod_vel_i[1] = 1.9;
    mod_vel_i[2] = -0.42;
    double mod_vel_j[3] = {};
    mod_vel_j[0] = 0.3;
    mod_vel_j[1] = -1.2;
    mod_vel_j[2] = 0.77;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;

    double acc_i[3] = {};

    momentumformulation_->modified_velocity_contribution(&dens_i, &dens_j, vel_i, vel_j, mod_vel_i,
        mod_vel_j, speccoeff_ij, speccoeff_ji, e_ij, acc_i, nullptr);

    // compute reference solution
    double A_i[3][3] = {};
    double A_j[3][3] = {};

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
      {
        A_i[i][j] = dens_i * vel_i[i] * (mod_vel_i[j] - vel_i[j]);

        A_j[i][j] = dens_j * vel_j[i] * (mod_vel_j[j] - vel_j[j]);
      }

    double A_ij_e_ij[3] = {};
    for (int i = 0; i < 3; ++i)
    {
      A_ij_e_ij[i] = 0.5 * (A_i[i][0] * e_ij[0] + A_i[i][1] * e_ij[1] + A_i[i][2] * e_ij[2]);

      A_ij_e_ij[i] += 0.5 * (A_j[i][0] * e_ij[0] + A_j[i][1] * e_ij[1] + A_j[i][2] * e_ij[2]);
    }

    double acc_i_ref[3] = {};

    for (int i = 0; i < 3; ++i) acc_i_ref[i] = speccoeff_ij * A_ij_e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_i, acc_i_ref, 3, 1.0e-14);
  }


  TEST_F(SPHMomentumFormulationAdamiTest, ModifiedVelocityContributionNullptrAccIModVelI)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double vel_i[3] = {};
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3] = {};
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;
    double mod_vel_j[3] = {};
    mod_vel_j[0] = 0.3;
    mod_vel_j[1] = -1.2;
    mod_vel_j[2] = 0.77;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;

    double acc_j[3] = {};

    momentumformulation_->modified_velocity_contribution(&dens_i, &dens_j, vel_i, vel_j, nullptr,
        mod_vel_j, speccoeff_ij, speccoeff_ji, e_ij, nullptr, acc_j);

    // compute reference solution
    double A_j[3][3] = {};

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) A_j[i][j] = dens_j * vel_j[i] * (mod_vel_j[j] - vel_j[j]);

    double A_ij_e_ij[3] = {};
    for (int i = 0; i < 3; ++i)
      A_ij_e_ij[i] += 0.5 * (A_j[i][0] * e_ij[0] + A_j[i][1] * e_ij[1] + A_j[i][2] * e_ij[2]);

    double acc_j_ref[3] = {};

    for (int i = 0; i < 3; ++i) acc_j_ref[i] = -speccoeff_ji * A_ij_e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_j, acc_j_ref, 3, 1.0e-14);
  }

  TEST_F(SPHMomentumFormulationAdamiTest, ModifiedVelocityContributionNullptrAccIModVelJ)
  {
    const double dens_i = 1.02;
    const double dens_j = 0.97;
    double vel_i[3] = {};
    vel_i[0] = 0.1;
    vel_i[1] = 2.0;
    vel_i[2] = -0.4;
    double vel_j[3] = {};
    vel_j[0] = 0.27;
    vel_j[1] = -1.3;
    vel_j[2] = 0.5;
    double mod_vel_i[3] = {};
    mod_vel_i[0] = 0.12;
    mod_vel_i[1] = 1.9;
    mod_vel_i[2] = -0.42;
    double e_ij[3] = {};
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);
    const double speccoeff_ij = 12.42;
    const double speccoeff_ji = 0.35;

    double acc_i[3] = {};

    momentumformulation_->modified_velocity_contribution(&dens_i, &dens_j, vel_i, vel_j, mod_vel_i,
        nullptr, speccoeff_ij, speccoeff_ji, e_ij, acc_i, nullptr);

    // compute reference solution
    double A_i[3][3] = {};

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) A_i[i][j] = dens_i * vel_i[i] * (mod_vel_i[j] - vel_i[j]);

    double A_ij_e_ij[3] = {};
    for (int i = 0; i < 3; ++i)
      A_ij_e_ij[i] = 0.5 * (A_i[i][0] * e_ij[0] + A_i[i][1] * e_ij[1] + A_i[i][2] * e_ij[2]);

    double acc_i_ref[3] = {};

    for (int i = 0; i < 3; ++i) acc_i_ref[i] = speccoeff_ij * A_ij_e_ij[i];

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_i, acc_i_ref, 3, 1.0e-14);
  }
}  // namespace
