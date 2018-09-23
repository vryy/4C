/*---------------------------------------------------------------------------*/
/*!
\file unit_particle_interaction_sph_density_correction.H

\brief unittests for density correction handler in smoothed particle hydrodynamics (SPH)

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#ifndef UNIT_PARTICLE_INTERACTION_SPH_DENSITYCORRECTION_H
#define UNIT_PARTICLE_INTERACTION_SPH_DENSITYCORRECTION_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include <cxxtest/TestSuite.h>
#include "src/drt_particle_interaction/particle_interaction_sph_density_correction.H"

/*---------------------------------------------------------------------------*
 | forward declaration                                        sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHDensityCorrectionInterior_TestSuite;
  class SPHDensityCorrectionNormalized_TestSuite;
  class SPHDensityCorrectionRandles_TestSuite;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | density correction handler test suite                      sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
class PARTICLEINTERACTION::SPHDensityCorrectionInterior_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::SPHDensityCorrectionInterior> densitycorrection_;

 public:
  void setUp()
  {
    // create density correction handler
    densitycorrection_ = std::unique_ptr<PARTICLEINTERACTION::SPHDensityCorrectionInterior>(
        new PARTICLEINTERACTION::SPHDensityCorrectionInterior());

    // init density correction handler
    densitycorrection_->Init();

    // setup density correction handler
    densitycorrection_->Setup();
  }

  void tearDown() { densitycorrection_ = nullptr; }

  // note: the public functions Init() and Setup() of class SPHEquationOfStateGenTait are called in
  // setUp()
  //       and thus implicitly tested by all following unittests

  void test_ComputeDensityBC() { TS_ASSERT_EQUALS(densitycorrection_->ComputeDensityBC(), false); }

  void test_CorrectedDensityInterior()
  {
    const double denssum = 1.07;
    double dens = 0.98;

    densitycorrection_->CorrectedDensityInterior(&denssum, &dens);

    TS_ASSERT_EQUALS(dens, denssum);
  }

  void test_CorrectedDensityFreeSurface()
  {
    const double denssum = 1.07;
    const double colorfield = 0.82;
    const double dens_bc = 1.05;
    double dens = 0.78;

    const double dens_ref = dens;

    densitycorrection_->CorrectedDensityFreeSurface(&denssum, &colorfield, &dens_bc, &dens);

    TS_ASSERT_EQUALS(dens, dens_ref);
  }
};

class PARTICLEINTERACTION::SPHDensityCorrectionNormalized_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::SPHDensityCorrectionNormalized> densitycorrection_;

 public:
  void setUp()
  {
    // create density correction handler
    densitycorrection_ = std::unique_ptr<PARTICLEINTERACTION::SPHDensityCorrectionNormalized>(
        new PARTICLEINTERACTION::SPHDensityCorrectionNormalized());

    // init density correction handler
    densitycorrection_->Init();

    // setup density correction handler
    densitycorrection_->Setup();
  }

  void tearDown() { densitycorrection_ = nullptr; }

  // note: the public functions Init() and Setup() of class SPHEquationOfStateGenTait are called in
  // setUp()
  //       and thus implicitly tested by all following unittests

  void test_ComputeDensityBC() { TS_ASSERT_EQUALS(densitycorrection_->ComputeDensityBC(), false); }

  void test_CorrectedDensityInterior()
  {
    const double denssum = 1.07;
    double dens = 0.98;

    densitycorrection_->CorrectedDensityInterior(&denssum, &dens);

    TS_ASSERT_EQUALS(dens, denssum);
  }

  void test_CorrectedDensityFreeSurface()
  {
    const double denssum = 1.07;
    const double colorfield = 0.82;
    const double dens_bc = 1.05;
    double dens = 0.78;

    const double dens_ref = denssum / colorfield;

    densitycorrection_->CorrectedDensityFreeSurface(&denssum, &colorfield, &dens_bc, &dens);

    TS_ASSERT_EQUALS(dens, dens_ref);
  }
};

class PARTICLEINTERACTION::SPHDensityCorrectionRandles_TestSuite : public CxxTest::TestSuite
{
 private:
  std::unique_ptr<PARTICLEINTERACTION::SPHDensityCorrectionRandles> densitycorrection_;

 public:
  void setUp()
  {
    // create density correction handler
    densitycorrection_ = std::unique_ptr<PARTICLEINTERACTION::SPHDensityCorrectionRandles>(
        new PARTICLEINTERACTION::SPHDensityCorrectionRandles());

    // init density correction handler
    densitycorrection_->Init();

    // setup density correction handler
    densitycorrection_->Setup();
  }

  void tearDown() { densitycorrection_ = nullptr; }

  // note: the public functions Init() and Setup() of class SPHEquationOfStateGenTait are called in
  // setUp()
  //       and thus implicitly tested by all following unittests

  void test_ComputeDensityBC() { TS_ASSERT_EQUALS(densitycorrection_->ComputeDensityBC(), true); }

  void test_CorrectedDensityInterior()
  {
    const double denssum = 1.07;
    double dens = 0.98;

    densitycorrection_->CorrectedDensityInterior(&denssum, &dens);

    TS_ASSERT_EQUALS(dens, denssum);
  }

  void test_CorrectedDensityFreeSurface()
  {
    const double denssum = 1.07;
    const double colorfield = 0.82;
    const double dens_bc = 1.05;
    double dens = 0.78;

    const double dens_ref = denssum + dens_bc * (1.0 - colorfield);

    densitycorrection_->CorrectedDensityFreeSurface(&denssum, &colorfield, &dens_bc, &dens);

    TS_ASSERT_EQUALS(dens, dens_ref);
  }
};

/*---------------------------------------------------------------------------*/
#endif
