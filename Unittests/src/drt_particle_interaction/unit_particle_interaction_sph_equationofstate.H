/*---------------------------------------------------------------------------*/
/*! \file
\brief unittests for equation of state handler for smoothed particle hydrodynamics (SPH)
interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#ifndef UNIT_PARTICLE_INTERACTION_SPH_EQUATIONOFSTATE_H
#define UNIT_PARTICLE_INTERACTION_SPH_EQUATIONOFSTATE_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#include <cxxtest/TestSuite.h>
#include "src/drt_particle_interaction/particle_interaction_sph_equationofstate.H"

/*---------------------------------------------------------------------------*
 | forward declaration                                        sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHEquationOfStateGenTait_TestSuite;
  class SPHEquationOfStateIdealGas_TestSuite;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | kernel handler test suite                                  sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
class PARTICLEINTERACTION::SPHEquationOfStateGenTait_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::SPHEquationOfStateGenTait> equationofstate_;
  std::unique_ptr<PARTICLEINTERACTION::SPHEquationOfStateGenTait> equationofstate_special_;

 public:
  void setUp()
  {
    const double speedofsound = 3.5;
    const double refdensfac = 0.9;
    const double exponent = 7;

    // create equation of state handler
    equationofstate_ = std::unique_ptr<PARTICLEINTERACTION::SPHEquationOfStateGenTait>(
        new PARTICLEINTERACTION::SPHEquationOfStateGenTait(speedofsound, refdensfac, exponent));
    equationofstate_special_ = std::unique_ptr<PARTICLEINTERACTION::SPHEquationOfStateGenTait>(
        new PARTICLEINTERACTION::SPHEquationOfStateGenTait(speedofsound, refdensfac, 1.0));

    // init equation of state handler
    equationofstate_->Init();
    equationofstate_special_->Init();

    // setup equation of state handler
    equationofstate_->Setup();
    equationofstate_special_->Setup();
  }

  void tearDown()
  {
    equationofstate_ = nullptr;
    equationofstate_special_ = nullptr;
  }

  // note: the public functions Init() and Setup() of class SPHEquationOfStateGenTait are called in
  // setUp()
  //       and thus implicitly tested by all following unittests

  void test_DensityToPressure()
  {
    const double density = 0.75;
    const double density0 = 0.4;

    TS_ASSERT_DELTA(
        equationofstate_->DensityToPressure(density, density0), 0.5640046918e2, 1.0e-08);
    TS_ASSERT_DELTA(
        equationofstate_special_->DensityToPressure(density, density0), 0.47775e1, 1.0e-08);
  }

  void test_PressureToDensity()
  {
    const double pressure = 1.8;
    const double density0 = 0.4;

    TS_ASSERT_DELTA(equationofstate_->PressureToDensity(pressure, density0), 0.477832244, 1.0e-08);
    TS_ASSERT_DELTA(
        equationofstate_special_->PressureToDensity(pressure, density0), 0.506938775, 1.0e-08);
  }

  void test_DensityToEnergy()
  {
    const double density = 0.75;
    const double mass = 1.3;
    const double density0 = 0.4;

    TS_ASSERT_DELTA(
        equationofstate_->DensityToEnergy(density, mass, density0), 0.1999413554e2, 1.0e-08);
    TS_ASSERT_DELTA(equationofstate_special_->DensityToEnergy(density, mass, density0),
        -0.5553068732e1, 1.0e-08);
  }
};

class PARTICLEINTERACTION::SPHEquationOfStateIdealGas_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::SPHEquationOfStateIdealGas> equationofstate_;

 public:
  void setUp()
  {
    const double speedofsound = 3.5;

    // create equation of state handler
    equationofstate_ = std::unique_ptr<PARTICLEINTERACTION::SPHEquationOfStateIdealGas>(
        new PARTICLEINTERACTION::SPHEquationOfStateIdealGas(speedofsound));

    // init equation of state handler
    equationofstate_->Init();

    // setup equation of state handler
    equationofstate_->Setup();
  }

  void tearDown() { equationofstate_ = nullptr; }

  // note: the public functions Init() and Setup() of class SPHEquationOfStateIdealGas are called in
  // setUp()
  //       and thus implicitly tested by all following unittests

  void test_DensityToPressure()
  {
    const double density = 0.75;
    const double density0 = 0.4;

    TS_ASSERT_DELTA(equationofstate_->DensityToPressure(density, density0), 9.1875, 1.0e-08);
  }

  void test_PressureToDensity()
  {
    const double pressure = 1.8;
    const double density0 = 0.4;

    TS_ASSERT_DELTA(equationofstate_->PressureToDensity(pressure, density0), 0.146938775, 1.0e-08);
  }

  void test_DensityToEnergy()
  {
    const double density = 0.75;
    const double mass = 1.3;
    const double density0 = 0.4;

    TS_ASSERT_DELTA(
        equationofstate_->DensityToEnergy(density, mass, density0), -0.2752956873e2, 1.0e-08);
  }
};

/*---------------------------------------------------------------------------*/
#endif
