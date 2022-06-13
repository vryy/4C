/*---------------------------------------------------------------------------*/
/*! \file
\brief unittests for utils for particle interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "src/drt_particle_interaction/particle_interaction_utils.H"

namespace
{
  TEST(PowHelperTest, Pow)
  {
    EXPECT_NEAR(PARTICLEINTERACTION::UTILS::Pow<2>(1.34), 1.7956, 1.0e-14);
    EXPECT_NEAR(PARTICLEINTERACTION::UTILS::Pow<5>(0.8), 0.32768, 1.0e-14);
    EXPECT_NEAR(PARTICLEINTERACTION::UTILS::Pow<4>(3.5), 150.0625, 1.0e-14);
  }

  TEST(PowHelperTest, VecClear)
  {
    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecClear(c);

    EXPECT_NEAR(c[0], 0.0, 1.0e-14);
    EXPECT_NEAR(c[1], 0.0, 1.0e-14);
    EXPECT_NEAR(c[2], 0.0, 1.0e-14);
  }

  TEST(PowHelperTest, VecSet)
  {
    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecSet(c, a);

    EXPECT_NEAR(c[0], 1.0, 1.0e-14);
    EXPECT_NEAR(c[1], -2.0, 1.0e-14);
    EXPECT_NEAR(c[2], 4.25, 1.0e-14);
  }

  TEST(PowHelperTest, VecAdd)
  {
    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecAdd(c, a);

    EXPECT_NEAR(c[0], 3.5, 1.0e-14);
    EXPECT_NEAR(c[1], 5.5, 1.0e-14);
    EXPECT_NEAR(c[2], 2.45, 1.0e-14);
  }

  TEST(PowHelperTest, VecSub)
  {
    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecSub(c, a);

    EXPECT_NEAR(c[0], 1.5, 1.0e-14);
    EXPECT_NEAR(c[1], 9.5, 1.0e-14);
    EXPECT_NEAR(c[2], -6.05, 1.0e-14);
  }

  TEST(PowHelperTest, VecScale)
  {
    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecScale(c, 1.8);

    EXPECT_NEAR(c[0], 4.5, 1.0e-14);
    EXPECT_NEAR(c[1], 13.5, 1.0e-14);
    EXPECT_NEAR(c[2], -3.24, 1.0e-14);
  }

  TEST(PowHelperTest, VecSetScale)
  {
    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecSetScale(c, 1.8, a);

    EXPECT_NEAR(c[0], 1.8, 1.0e-14);
    EXPECT_NEAR(c[1], -3.6, 1.0e-14);
    EXPECT_NEAR(c[2], 7.65, 1.0e-14);
  }

  TEST(PowHelperTest, VecAddScale)
  {
    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecAddScale(c, 1.8, a);

    EXPECT_NEAR(c[0], 4.3, 1.0e-14);
    EXPECT_NEAR(c[1], 3.9, 1.0e-14);
    EXPECT_NEAR(c[2], 5.85, 1.0e-14);
  }

  TEST(PowHelperTest, VecSetCross)
  {
    const double a[3] = {1.0, -2.0, 4.25};
    const double b[3] = {-0.5, -3.5, 0.2};
    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecSetCross(c, a, b);

    EXPECT_NEAR(c[0], 14.4750, 1.0e-14);
    EXPECT_NEAR(c[1], -2.325, 1.0e-14);
    EXPECT_NEAR(c[2], -4.5, 1.0e-14);
  }

  TEST(PowHelperTest, VecAddCross)
  {
    const double a[3] = {1.0, -2.0, 4.25};
    const double b[3] = {-0.5, -3.5, 0.2};
    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecAddCross(c, a, b);

    EXPECT_NEAR(c[0], 16.9750, 1.0e-14);
    EXPECT_NEAR(c[1], 5.175, 1.0e-14);
    EXPECT_NEAR(c[2], -6.3, 1.0e-14);
  }

  TEST(PowHelperTest, VecDot)
  {
    const double a[3] = {1.0, -2.0, 4.25};
    const double b[3] = {-0.5, -3.5, 0.2};

    const double a_dot_b = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

    EXPECT_NEAR(PARTICLEINTERACTION::UTILS::VecDot(a, b), a_dot_b, 1.0e-14);
  }

  TEST(PowHelperTest, VecNormTwo)
  {
    const double a[3] = {1.0, -2.0, 4.25};

    const double a_norm2 = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

    EXPECT_NEAR(PARTICLEINTERACTION::UTILS::VecNormTwo(a), a_norm2, 1.0e-14);
  }

  TEST(PowHelperTest, UnitSurfaceTangents)
  {
    double n[3] = {1.0, -2.0, 4.25};
    PARTICLEINTERACTION::UTILS::VecScale(n, 1.0 / PARTICLEINTERACTION::UTILS::VecNormTwo(n));

    double t1[3] = {0.0};
    double t2[3] = {0.0};

    PARTICLEINTERACTION::UTILS::UnitSurfaceTangents(n, t1, t2);

    EXPECT_NEAR(PARTICLEINTERACTION::UTILS::VecNormTwo(t1), 1.0, 1.0e-14);
    EXPECT_NEAR(PARTICLEINTERACTION::UTILS::VecNormTwo(t2), 1.0, 1.0e-14);

    EXPECT_NEAR(PARTICLEINTERACTION::UTILS::VecDot(n, t1), 0.0, 1.0e-14);
    EXPECT_NEAR(PARTICLEINTERACTION::UTILS::VecDot(n, t2), 0.0, 1.0e-14);
    EXPECT_NEAR(PARTICLEINTERACTION::UTILS::VecDot(t1, t2), 0.0, 1.0e-14);

    double n_ref[3] = {0.0};
    PARTICLEINTERACTION::UTILS::VecSetCross(n_ref, t1, t2);

    EXPECT_NEAR(n_ref[0], n[0], 1.0e-14);
    EXPECT_NEAR(n_ref[1], n[1], 1.0e-14);
    EXPECT_NEAR(n_ref[2], n[2], 1.0e-14);
  }

  TEST(PowHelperTest, LinTransLower)
  {
    EXPECT_NEAR(PARTICLEINTERACTION::UTILS::LinTrans(0.5, 1.2, 3.8), 0.0, 1.0e-14);
  }

  TEST(PowHelperTest, LinTransIn)
  {
    EXPECT_NEAR(PARTICLEINTERACTION::UTILS::LinTrans(2.24, 1.2, 3.8), 0.4, 1.0e-14);
  }

  TEST(PowHelperTest, LinTransUpper)
  {
    EXPECT_NEAR(PARTICLEINTERACTION::UTILS::LinTrans(4.0, 1.2, 3.8), 1.0, 1.0e-14);
  }

}  // namespace
