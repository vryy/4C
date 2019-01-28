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

  void test_vec_set()
  {
    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    UTILS::vec_set(c, a);

    TS_ASSERT_DELTA(c[0], 1.0, 1e-14);
    TS_ASSERT_DELTA(c[1], -2.0, 1e-14);
    TS_ASSERT_DELTA(c[2], 4.25, 1e-14);
  }

  void test_vec_add()
  {
    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    UTILS::vec_add(c, a);

    TS_ASSERT_DELTA(c[0], 3.5, 1e-14);
    TS_ASSERT_DELTA(c[1], 5.5, 1e-14);
    TS_ASSERT_DELTA(c[2], 2.45, 1e-14);
  }

  void test_vec_sub()
  {
    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    UTILS::vec_sub(c, a);

    TS_ASSERT_DELTA(c[0], 1.5, 1e-14);
    TS_ASSERT_DELTA(c[1], 9.5, 1e-14);
    TS_ASSERT_DELTA(c[2], -6.05, 1e-14);
  }

  void test_vec_setscale()
  {
    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    UTILS::vec_setscale(c, 1.8, a);

    TS_ASSERT_DELTA(c[0], 1.8, 1e-14);
    TS_ASSERT_DELTA(c[1], -3.6, 1e-14);
    TS_ASSERT_DELTA(c[2], 7.65, 1e-14);
  }

  void test_vec_addscale()
  {
    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    UTILS::vec_addscale(c, 1.8, a);

    TS_ASSERT_DELTA(c[0], 4.3, 1e-14);
    TS_ASSERT_DELTA(c[1], 3.9, 1e-14);
    TS_ASSERT_DELTA(c[2], 5.85, 1e-14);
  }

  void test_vec_addcross()
  {
    const double a[3] = {1.0, -2.0, 4.25};
    const double b[3] = {-0.5, -3.5, 0.2};
    double c[3] = {2.5, 7.5, -1.8};
    UTILS::vec_addcross(c, a, b);

    TS_ASSERT_DELTA(c[0], 16.9750, 1e-14);
    TS_ASSERT_DELTA(c[1], 5.175, 1e-14);
    TS_ASSERT_DELTA(c[2], -6.3, 1e-14);
  }

  void test_vec_dot()
  {
    const double a[3] = {1.0, -2.0, 4.25};
    const double b[3] = {-0.5, -3.5, 0.2};

    const double a_dot_b = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

    TS_ASSERT_DELTA(UTILS::vec_dot(a, b), a_dot_b, 1e-14);
  }

  void test_vec_norm2()
  {
    const double a[3] = {1.0, -2.0, 4.25};

    const double a_norm2 = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

    TS_ASSERT_DELTA(UTILS::vec_norm2(a), a_norm2, 1e-14);
  }
};

/*---------------------------------------------------------------------------*/
#endif
