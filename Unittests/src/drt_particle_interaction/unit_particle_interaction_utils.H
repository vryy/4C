/*---------------------------------------------------------------------------*/
/*!
\file unit_particle_interaction_utils.H

\brief unittests for utils for particle interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
#ifndef UNIT_PARTICLE_INTERACTION_UTILS_H
#define UNIT_PARTICLE_INTERACTION_UTILS_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
#include <cxxtest/TestSuite.h>
#include "src/drt_particle_interaction/particle_interaction_utils.H"

/*---------------------------------------------------------------------------*
 | forward declaration                                        sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  namespace UTILS
  {
    class PowHelper_TestSuite;
  }
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | pow helper test suite                                      sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
class PARTICLEINTERACTION::UTILS::PowHelper_TestSuite : public CxxTest::TestSuite
{
 public:
  void setUp() {}

  void tearDown() {}

  void test_pow()
  {
    TS_ASSERT_DELTA(UTILS::pow<2>(1.34), 1.7956, 1e-14);
    TS_ASSERT_DELTA(UTILS::pow<5>(0.8), 0.32768, 1e-14);
    TS_ASSERT_DELTA(UTILS::pow<4>(3.5), 150.0625, 1e-14);
  }
};

/*---------------------------------------------------------------------------*/
#endif
