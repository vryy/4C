/*---------------------------------------------------------------------------*/
/*! \file
\brief unittests for utils for particle interactions
\level 3
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
#include "src/common/unit_cxx_test_wrapper.H"
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
class PARTICLEINTERACTION::UTILS::PowHelper_TestSuite : public BACICxxTestWrapper
{
 public:
  void Setup() override {}

  void TearDown() override {}

  void TestPow()
  {
    TS_ASSERT_DELTA(UTILS::Pow<2>(1.34), 1.7956, 1.0e-14);
    TS_ASSERT_DELTA(UTILS::Pow<5>(0.8), 0.32768, 1.0e-14);
    TS_ASSERT_DELTA(UTILS::Pow<4>(3.5), 150.0625, 1.0e-14);
  }

  void TestVecClear()
  {
    double c[3] = {2.5, 7.5, -1.8};
    UTILS::VecClear(c);

    TS_ASSERT_DELTA(c[0], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(c[1], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(c[2], 0.0, 1.0e-14);
  }

  void TestVecSet()
  {
    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    UTILS::VecSet(c, a);

    TS_ASSERT_DELTA(c[0], 1.0, 1.0e-14);
    TS_ASSERT_DELTA(c[1], -2.0, 1.0e-14);
    TS_ASSERT_DELTA(c[2], 4.25, 1.0e-14);
  }

  void TestVecAdd()
  {
    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    UTILS::VecAdd(c, a);

    TS_ASSERT_DELTA(c[0], 3.5, 1.0e-14);
    TS_ASSERT_DELTA(c[1], 5.5, 1.0e-14);
    TS_ASSERT_DELTA(c[2], 2.45, 1.0e-14);
  }

  void TestVecSub()
  {
    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    UTILS::VecSub(c, a);

    TS_ASSERT_DELTA(c[0], 1.5, 1.0e-14);
    TS_ASSERT_DELTA(c[1], 9.5, 1.0e-14);
    TS_ASSERT_DELTA(c[2], -6.05, 1.0e-14);
  }

  void TestVecScale()
  {
    double c[3] = {2.5, 7.5, -1.8};
    UTILS::VecScale(c, 1.8);

    TS_ASSERT_DELTA(c[0], 4.5, 1.0e-14);
    TS_ASSERT_DELTA(c[1], 13.5, 1.0e-14);
    TS_ASSERT_DELTA(c[2], -3.24, 1.0e-14);
  }

  void TestVecSetScale()
  {
    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    UTILS::VecSetScale(c, 1.8, a);

    TS_ASSERT_DELTA(c[0], 1.8, 1.0e-14);
    TS_ASSERT_DELTA(c[1], -3.6, 1.0e-14);
    TS_ASSERT_DELTA(c[2], 7.65, 1.0e-14);
  }

  void TestVecAddScale()
  {
    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    UTILS::VecAddScale(c, 1.8, a);

    TS_ASSERT_DELTA(c[0], 4.3, 1.0e-14);
    TS_ASSERT_DELTA(c[1], 3.9, 1.0e-14);
    TS_ASSERT_DELTA(c[2], 5.85, 1.0e-14);
  }

  void TestVecSetCross()
  {
    const double a[3] = {1.0, -2.0, 4.25};
    const double b[3] = {-0.5, -3.5, 0.2};
    double c[3] = {2.5, 7.5, -1.8};
    UTILS::VecSetCross(c, a, b);

    TS_ASSERT_DELTA(c[0], 14.4750, 1.0e-14);
    TS_ASSERT_DELTA(c[1], -2.325, 1.0e-14);
    TS_ASSERT_DELTA(c[2], -4.5, 1.0e-14);
  }

  void TestVecAddCross()
  {
    const double a[3] = {1.0, -2.0, 4.25};
    const double b[3] = {-0.5, -3.5, 0.2};
    double c[3] = {2.5, 7.5, -1.8};
    UTILS::VecAddCross(c, a, b);

    TS_ASSERT_DELTA(c[0], 16.9750, 1.0e-14);
    TS_ASSERT_DELTA(c[1], 5.175, 1.0e-14);
    TS_ASSERT_DELTA(c[2], -6.3, 1.0e-14);
  }

  void TestVecDot()
  {
    const double a[3] = {1.0, -2.0, 4.25};
    const double b[3] = {-0.5, -3.5, 0.2};

    const double a_dot_b = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

    TS_ASSERT_DELTA(UTILS::VecDot(a, b), a_dot_b, 1.0e-14);
  }

  void TestVecNormTwo()
  {
    const double a[3] = {1.0, -2.0, 4.25};

    const double a_norm2 = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

    TS_ASSERT_DELTA(UTILS::VecNormTwo(a), a_norm2, 1.0e-14);
  }

  void TestUnitSurfaceTangents()
  {
    double n[3] = {1.0, -2.0, 4.25};
    UTILS::VecScale(n, 1.0 / UTILS::VecNormTwo(n));

    double t1[3] = {0.0};
    double t2[3] = {0.0};

    UTILS::UnitSurfaceTangents(n, t1, t2);

    TS_ASSERT_DELTA(UTILS::VecNormTwo(t1), 1.0, 1.0e-14);
    TS_ASSERT_DELTA(UTILS::VecNormTwo(t2), 1.0, 1.0e-14);

    TS_ASSERT_DELTA(UTILS::VecDot(n, t1), 0.0, 1.0e-14);
    TS_ASSERT_DELTA(UTILS::VecDot(n, t2), 0.0, 1.0e-14);
    TS_ASSERT_DELTA(UTILS::VecDot(t1, t2), 0.0, 1.0e-14);

    double n_ref[3] = {0.0};
    VecSetCross(n_ref, t1, t2);

    TS_ASSERT_DELTA(n_ref[0], n[0], 1.0e-14);
    TS_ASSERT_DELTA(n_ref[1], n[1], 1.0e-14);
    TS_ASSERT_DELTA(n_ref[2], n[2], 1.0e-14);
  }

  void TestLinTransLower() { TS_ASSERT_DELTA(UTILS::LinTrans(0.5, 1.2, 3.8), 0.0, 1.0e-14); }

  void TestLinTransIn() { TS_ASSERT_DELTA(UTILS::LinTrans(2.24, 1.2, 3.8), 0.4, 1.0e-14); }

  void TestLinTransUpper() { TS_ASSERT_DELTA(UTILS::LinTrans(4.0, 1.2, 3.8), 1.0, 1.0e-14); }
};

/*---------------------------------------------------------------------------*/
#endif
