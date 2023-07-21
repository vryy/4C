/*---------------------------------------------------------------------------*/
/*! \file
\brief unittests for utils for rigid bodies
\level 3
*/
/*---------------------------------------------------------------------------*/

#include <gtest/gtest.h>
#include "baci_unittest_utils_assertions.h"
#include "baci_particle_rigidbody_utils.H"

namespace
{
  TEST(QuaternionTest, Clear)
  {
    const double q_ref[4] = {0.0, 0.0, 0.0, 1.0};
    double q[4] = {1.0, 2.0, 3.0, 4.0};
    PARTICLERIGIDBODY::UTILS::QuaternionClear(q);

    BACI_EXPECT_ITERABLE_NEAR(q, q_ref, 4, 1.0e-14);
  }

  TEST(QuaternionTest, Set)
  {
    double q1[4] = {0.0, 0.0, 0.0, 1.0};
    const double q2[4] = {1.0, 2.0, 3.0, 4.0};

    PARTICLERIGIDBODY::UTILS::QuaternionSet(q1, q2);

    BACI_EXPECT_ITERABLE_NEAR(q1, q2, 4, 1.0e-14);
  }

  TEST(QuaternionTest, Invert)
  {
    double q1[4] = {0.0, 0.0, 0.0, 1.0};
    const double q2[4] = {1.0, 2.0, 3.0, 4.0};
    const double q_ref[4] = {-q2[0], -q2[1], -q2[2], q2[3]};

    PARTICLERIGIDBODY::UTILS::QuaternionInvert(q1, q2);

    BACI_EXPECT_ITERABLE_NEAR(q1, q_ref, 4, 1.0e-14);
  }

  TEST(QuaternionTest, Product)
  {
    const double q12_ref[4] = {
        -0.01121419126499877, 0.9977058985744629, -0.05089858263600289, 0.04320319605818204};
    const double phi1[3] = {0.1, -2.0, 0.3};
    const double phi2[3] = {-0.8, 5.0, 0.0};

    double q1[4];
    PARTICLERIGIDBODY::UTILS::QuaternionFromAngle(q1, phi1);

    double q2[4];
    PARTICLERIGIDBODY::UTILS::QuaternionFromAngle(q2, phi2);

    double q12[4];
    PARTICLERIGIDBODY::UTILS::QuaternionProduct(q12, q2, q1);

    BACI_EXPECT_ITERABLE_NEAR(q12, q12_ref, 4, 1.0e-14);
  }

  TEST(QuaternionTest, FromAngleZero)
  {
    const double q_ref[4] = {0.0, 0.0, 0.0, 1.0};
    const double phi[3] = {0.0, 0.0, 0.0};

    double q[4];
    PARTICLERIGIDBODY::UTILS::QuaternionFromAngle(q, phi);

    BACI_EXPECT_ITERABLE_NEAR(q, q_ref, 4, 1.0e-14);
  }

  TEST(QuaternionTest, FromAngleXAxis)
  {
    const double q_ref[4] = {std::sin(M_PI / 4), 0.0, 0.0, std::cos(M_PI / 4)};
    const double phi[3] = {M_PI / 2, 0.0, 0.0};

    double q[4];
    PARTICLERIGIDBODY::UTILS::QuaternionFromAngle(q, phi);

    BACI_EXPECT_ITERABLE_NEAR(q, q_ref, 4, 1.0e-14);
  }

  TEST(QuaternionTest, FromAngleYAxis)
  {
    const double q_ref[4] = {0.0, std::sin(M_PI / 4), 0.0, std::cos(M_PI / 4)};
    const double phi[3] = {0.0, M_PI / 2, 0.0};

    double q[4];
    PARTICLERIGIDBODY::UTILS::QuaternionFromAngle(q, phi);

    BACI_EXPECT_ITERABLE_NEAR(q, q_ref, 4, 1.0e-14);
  }

  TEST(QuaternionTest, FromAngleZAxis)
  {
    const double q_ref[4] = {0.0, 0.0, std::sin(M_PI / 4), std::cos(M_PI / 4)};
    const double phi[3] = {0.0, 0.0, M_PI / 2};

    double q[4];
    PARTICLERIGIDBODY::UTILS::QuaternionFromAngle(q, phi);

    BACI_EXPECT_ITERABLE_NEAR(q, q_ref, 4, 1.0e-14);
  }

  TEST(QuaternionTest, FromAngleGeneral)
  {
    const double q_ref[4] = {
        -0.2759788075111623, 0.8279364225334871, 0.4139682112667435, -0.2588190451025209};
    const double phi[3] = {-M_PI / 3, M_PI, M_PI / 2};

    double q[4];
    PARTICLERIGIDBODY::UTILS::QuaternionFromAngle(q, phi);

    BACI_EXPECT_ITERABLE_NEAR(q, q_ref, 4, 1.0e-14);
  }

  TEST(QuaternionTest, RotateVectorXUnitAroundZAxis)
  {
    const double w_ref[3] = {0.0, 1.0, 0.0};
    double v[3] = {1.0, 0.0, 0.0};
    double w[3] = {0.0, 0.0, 0.0};

    const double q[4] = {0.0, 0.0, std::sin(M_PI / 4), std::cos(M_PI / 4)};

    PARTICLERIGIDBODY::UTILS::QuaternionRotateVector(w, q, v);

    BACI_EXPECT_ITERABLE_NEAR(w, w_ref, 3, 1.0e-14);
  }

  TEST(QuaternionTest, RotateVectorZUnitAroundYAxis)
  {
    const double w_ref[3] = {1.0, 0.0, 0.0};
    double v[3] = {0.0, 0.0, 1.0};
    double w[3] = {0.0, 0.0, 0.0};

    const double q[4] = {0.0, std::sin(M_PI / 4), 0.0, std::cos(M_PI / 4)};

    PARTICLERIGIDBODY::UTILS::QuaternionRotateVector(w, q, v);

    BACI_EXPECT_ITERABLE_NEAR(w, w_ref, 3, 1.0e-14);
  }

  TEST(QuaternionTest, RotateVectorYUnitAroundXAxis)
  {
    const double w_ref[3] = {0.0, 0.0, 1.0};
    double v[3] = {0.0, 1.0, 0.0};
    double w[3] = {0.0, 0.0, 0.0};

    const double q[4] = {std::sin(M_PI / 4), 0.0, 0.0, std::cos(M_PI / 4)};

    PARTICLERIGIDBODY::UTILS::QuaternionRotateVector(w, q, v);

    BACI_EXPECT_ITERABLE_NEAR(w, w_ref, 3, 1.0e-14);
  }

  TEST(QuaternionTest, RotateVectorGeneral)
  {
    const double w_ref[3] = {0.7145801717316358, -0.9159468817988596, 1.97494721141881};
    double v[3] = {0.5, 1.0, -2.0};
    double w[3] = {0.0, 0.0, 0.0};

    const double phi[3] = {-M_PI / 3, M_PI, M_PI / 2};

    double q[4];
    PARTICLERIGIDBODY::UTILS::QuaternionFromAngle(q, phi);

    PARTICLERIGIDBODY::UTILS::QuaternionRotateVector(w, q, v);

    BACI_EXPECT_ITERABLE_NEAR(w, w_ref, 3, 1.0e-14);
  }
}  // namespace
