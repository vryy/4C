/*---------------------------------------------------------------------------*/
/*! \file
\brief unittests for equation of state handler for smoothed particle hydrodynamics (SPH)
interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "src/drt_particle_interaction/particle_interaction_sph_equationofstate.H"

namespace
{

  class SPHEquationOfStateGenTaitTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<PARTICLEINTERACTION::SPHEquationOfStateGenTait> equationofstate_;
    std::unique_ptr<PARTICLEINTERACTION::SPHEquationOfStateGenTait> equationofstate_special_;

    SPHEquationOfStateGenTaitTest()
    {
      const double speedofsound = 3.5;
      const double refdensfac = 0.9;
      const double exponent = 7;

      // create equation of state handler
      equationofstate_ = std::make_unique<PARTICLEINTERACTION::SPHEquationOfStateGenTait>(
          speedofsound, refdensfac, exponent);
      equationofstate_special_ = std::make_unique<PARTICLEINTERACTION::SPHEquationOfStateGenTait>(
          speedofsound, refdensfac, 1.0);

      // init equation of state handler
      equationofstate_->Init();
      equationofstate_special_->Init();

      // setup equation of state handler
      equationofstate_->Setup();
      equationofstate_special_->Setup();
    }
    // note: the public functions Init() and Setup() of class SPHEquationOfStateGenTait are called
    // in SetUp() and thus implicitly tested by all following unittests
  };

  TEST_F(SPHEquationOfStateGenTaitTest, DensityToPressure)
  {
    const double density = 0.75;
    const double density0 = 0.4;

    EXPECT_NEAR(equationofstate_->DensityToPressure(density, density0), 0.5640046918e2, 1.0e-08);
    EXPECT_NEAR(equationofstate_special_->DensityToPressure(density, density0), 0.47775e1, 1.0e-08);
  }
  TEST_F(SPHEquationOfStateGenTaitTest, PressureToDensity)
  {
    const double pressure = 1.8;
    const double density0 = 0.4;

    EXPECT_NEAR(equationofstate_->PressureToDensity(pressure, density0), 0.477832244, 1.0e-08);
    EXPECT_NEAR(
        equationofstate_special_->PressureToDensity(pressure, density0), 0.506938775, 1.0e-08);
  }

  TEST_F(SPHEquationOfStateGenTaitTest, DensityToEnergy)
  {
    const double density = 0.75;
    const double mass = 1.3;
    const double density0 = 0.4;

    EXPECT_NEAR(
        equationofstate_->DensityToEnergy(density, mass, density0), 0.1999413554e2, 1.0e-08);
    EXPECT_NEAR(equationofstate_special_->DensityToEnergy(density, mass, density0), -0.5553068732e1,
        1.0e-08);
  }


  class SPHEquationOfStateIdealGasTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<PARTICLEINTERACTION::SPHEquationOfStateIdealGas> equationofstate_;

    SPHEquationOfStateIdealGasTest()
    {
      const double speedofsound = 3.5;

      // create equation of state handler
      equationofstate_ =
          std::make_unique<PARTICLEINTERACTION::SPHEquationOfStateIdealGas>(speedofsound);

      // init equation of state handler
      equationofstate_->Init();

      // setup equation of state handler
      equationofstate_->Setup();
    }
    // note: the public functions Init() and Setup() of class SPHEquationOfStateIdealGas are called
    // in SetUp() and thus implicitly tested by all following unittests
  };

  TEST_F(SPHEquationOfStateIdealGasTest, DensityToPressure)
  {
    const double density = 0.75;
    const double density0 = 0.4;

    EXPECT_NEAR(equationofstate_->DensityToPressure(density, density0), 9.1875, 1.0e-08);
  }

  TEST_F(SPHEquationOfStateIdealGasTest, PressureToDensity)
  {
    const double pressure = 1.8;
    const double density0 = 0.4;

    EXPECT_NEAR(equationofstate_->PressureToDensity(pressure, density0), 0.146938775, 1.0e-08);
  }

  TEST_F(SPHEquationOfStateIdealGasTest, DensityToEnergy)
  {
    const double density = 0.75;
    const double mass = 1.3;
    const double density0 = 0.4;

    EXPECT_NEAR(
        equationofstate_->DensityToEnergy(density, mass, density0), -0.2752956873e2, 1.0e-08);
  }
}  // namespace
