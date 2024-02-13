/*---------------------------------------------------------------------------*/
/*! \file
\brief unittests for utils for particle interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "baci_particle_interaction_utils.hpp"

#include "baci_unittest_utils_assertions_test.hpp"

namespace
{
  using namespace BACI;

  TEST(PowHelperTest, Pow)
  {
    EXPECT_NEAR(PARTICLEINTERACTION::UTILS::Pow<2>(1.34), 1.7956, 1.0e-14);
    EXPECT_NEAR(PARTICLEINTERACTION::UTILS::Pow<5>(0.8), 0.32768, 1.0e-14);
    EXPECT_NEAR(PARTICLEINTERACTION::UTILS::Pow<4>(3.5), 150.0625, 1.0e-14);
  }

  TEST(PowHelperTest, VecClear)
  {
    const double c_ref[3] = {0.0, 0.0, 0.0};

    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecClear(c);

    BACI_EXPECT_ITERABLE_NEAR(c, c_ref, 3, 1.0e-14);
  }

  TEST(PowHelperTest, VecSet)
  {
    const double c_ref[3] = {1.0, -2.0, 4.25};

    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecSet(c, a);

    BACI_EXPECT_ITERABLE_NEAR(c, c_ref, 3, 1.0e-14);
  }

  TEST(PowHelperTest, VecAdd)
  {
    const double c_ref[3] = {3.5, 5.5, 2.45};

    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecAdd(c, a);

    BACI_EXPECT_ITERABLE_NEAR(c, c_ref, 3, 1.0e-14);
  }

  TEST(PowHelperTest, VecSub)
  {
    const double c_ref[3] = {1.5, 9.5, -6.05};

    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecSub(c, a);

    BACI_EXPECT_ITERABLE_NEAR(c, c_ref, 3, 1.0e-14);
  }

  TEST(PowHelperTest, VecScale)
  {
    const double c_ref[3] = {4.5, 13.5, -3.24};

    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecScale(c, 1.8);

    BACI_EXPECT_ITERABLE_NEAR(c, c_ref, 3, 1.0e-14);
  }

  TEST(PowHelperTest, VecSetScale)
  {
    const double c_ref[3] = {1.8, -3.6, 7.65};

    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecSetScale(c, 1.8, a);

    BACI_EXPECT_ITERABLE_NEAR(c, c_ref, 3, 1.0e-14);
  }

  TEST(PowHelperTest, VecAddScale)
  {
    const double c_ref[3] = {4.3, 3.9, 5.85};

    const double a[3] = {1.0, -2.0, 4.25};
    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecAddScale(c, 1.8, a);

    BACI_EXPECT_ITERABLE_NEAR(c, c_ref, 3, 1.0e-14);
  }

  TEST(PowHelperTest, VecSetCross)
  {
    const double c_ref[3] = {14.4750, -2.325, -4.5};

    const double a[3] = {1.0, -2.0, 4.25};
    const double b[3] = {-0.5, -3.5, 0.2};
    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecSetCross(c, a, b);

    BACI_EXPECT_ITERABLE_NEAR(c, c_ref, 3, 1.0e-14);
  }

  TEST(PowHelperTest, VecAddCross)
  {
    const double c_ref[3] = {16.9750, 5.175, -6.3};

    const double a[3] = {1.0, -2.0, 4.25};
    const double b[3] = {-0.5, -3.5, 0.2};
    double c[3] = {2.5, 7.5, -1.8};
    PARTICLEINTERACTION::UTILS::VecAddCross(c, a, b);

    BACI_EXPECT_ITERABLE_NEAR(c, c_ref, 3, 1.0e-14);
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

    BACI_EXPECT_ITERABLE_NEAR(n, n_ref, 3, 1.0e-14);
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
