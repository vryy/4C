/*---------------------------------------------------------------------------*/
/*! \file
\brief unittests for utils for rigid bodies
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef UNIT_PARTICLE_RIGIDBODY_UTILS_H
#define UNIT_PARTICLE_RIGIDBODY_UTILS_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include <cxxtest/TestSuite.h>
#include "src/drt_particle_rigidbody/particle_rigidbody_utils.H"

/*---------------------------------------------------------------------------*
 | forward declaration                                                       |
 *---------------------------------------------------------------------------*/
namespace PARTICLERIGIDBODY
{
  namespace UTILS
  {
    class Quaternion_TestSuite;
  }
}  // namespace PARTICLERIGIDBODY

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
class PARTICLERIGIDBODY::UTILS::Quaternion_TestSuite : public CxxTest::TestSuite
{
 public:
  void setUp() {}

  void tearDown() {}

  void test_quaternion_clear()
  {
    double q[4] = {1.0, 2.0, 3.0, 4.0};
    UTILS::quaternion_clear(q);

    TS_ASSERT_DELTA(q[0], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(q[1], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(q[2], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(q[3], 1.0, 1.0e-14);
  }

  void test_quaternion_set()
  {
    double q1[4] = {0.0, 0.0, 0.0, 1.0};
    const double q2[4] = {1.0, 2.0, 3.0, 4.0};

    UTILS::quaternion_set(q1, q2);

    TS_ASSERT_DELTA(q1[0], q2[0], 1.0e-14);
    TS_ASSERT_DELTA(q1[1], q2[1], 1.0e-14);
    TS_ASSERT_DELTA(q1[2], q2[2], 1.0e-14);
    TS_ASSERT_DELTA(q1[3], q2[3], 1.0e-14);
  }

  void test_quaternion_invert()
  {
    double q1[4] = {0.0, 0.0, 0.0, 1.0};
    const double q2[4] = {1.0, 2.0, 3.0, 4.0};

    UTILS::quaternion_invert(q1, q2);

    TS_ASSERT_DELTA(q1[0], -q2[0], 1.0e-14);
    TS_ASSERT_DELTA(q1[1], -q2[1], 1.0e-14);
    TS_ASSERT_DELTA(q1[2], -q2[2], 1.0e-14);
    TS_ASSERT_DELTA(q1[3], q2[3], 1.0e-14);
  }

  void test_quaternion_product()
  {
    const double phi1[3] = {0.1, -2.0, 0.3};
    const double phi2[3] = {-0.8, 5.0, 0.0};

    double q1[4];
    UTILS::quaternion_from_angle(q1, phi1);

    double q2[4];
    UTILS::quaternion_from_angle(q2, phi2);

    double q12[4];
    UTILS::quaternion_product(q12, q2, q1);

    TS_ASSERT_DELTA(q12[0], -0.01121419126499877, 1.0e-14);
    TS_ASSERT_DELTA(q12[1], 0.9977058985744629, 1.0e-14);
    TS_ASSERT_DELTA(q12[2], -0.05089858263600289, 1.0e-14);
    TS_ASSERT_DELTA(q12[3], 0.04320319605818204, 1.0e-14);
  }

  void test_quaternion_from_angle_zero()
  {
    const double phi[3] = {0.0, 0.0, 0.0};

    double q[4];
    UTILS::quaternion_from_angle(q, phi);

    TS_ASSERT_DELTA(q[0], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(q[1], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(q[2], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(q[3], 1.0, 1.0e-14);
  }

  void test_quaternion_from_angle_x()
  {
    const double phi[3] = {M_PI / 2, 0.0, 0.0};

    double q[4];
    UTILS::quaternion_from_angle(q, phi);

    TS_ASSERT_DELTA(q[0], std::sin(M_PI / 4), 1.0e-14);
    TS_ASSERT_DELTA(q[1], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(q[2], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(q[3], std::cos(M_PI / 4), 1.0e-14);
  }

  void test_quaternion_from_angle_y()
  {
    const double phi[3] = {0.0, M_PI / 2, 0.0};

    double q[4];
    UTILS::quaternion_from_angle(q, phi);

    TS_ASSERT_DELTA(q[0], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(q[1], std::sin(M_PI / 4), 1.0e-14);
    TS_ASSERT_DELTA(q[2], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(q[3], std::cos(M_PI / 4), 1.0e-14);
  }

  void test_quaternion_from_angle_z()
  {
    const double phi[3] = {0.0, 0.0, M_PI / 2};

    double q[4];
    UTILS::quaternion_from_angle(q, phi);

    TS_ASSERT_DELTA(q[0], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(q[1], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(q[2], std::sin(M_PI / 4), 1.0e-14);
    TS_ASSERT_DELTA(q[3], std::cos(M_PI / 4), 1.0e-14);
  }

  void test_quaternion_from_angle_general()
  {
    const double phi[3] = {-M_PI / 3, M_PI, M_PI / 2};

    double q[4];
    UTILS::quaternion_from_angle(q, phi);

    TS_ASSERT_DELTA(q[0], -0.2759788075111623, 1.0e-14);
    TS_ASSERT_DELTA(q[1], 0.8279364225334871, 1.0e-14);
    TS_ASSERT_DELTA(q[2], 0.4139682112667435, 1.0e-14);
    TS_ASSERT_DELTA(q[3], -0.2588190451025209, 1.0e-14);
  }

  void test_quaternion_rotate_vector_x_around_z()
  {
    double v[3] = {1.0, 0.0, 0.0};
    double w[3] = {0.0, 0.0, 0.0};

    const double q[4] = {0.0, 0.0, std::sin(M_PI / 4), std::cos(M_PI / 4)};

    UTILS::quaternion_rotate_vector(w, q, v);

    TS_ASSERT_DELTA(w[0], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(w[1], 1.0, 1.0e-14);
    TS_ASSERT_DELTA(w[2], 0.0, 1.0e-14);
  }

  void test_quaternion_rotate_vector_z_around_y()
  {
    double v[3] = {0.0, 0.0, 1.0};
    double w[3] = {0.0, 0.0, 0.0};

    const double q[4] = {0.0, std::sin(M_PI / 4), 0.0, std::cos(M_PI / 4)};

    UTILS::quaternion_rotate_vector(w, q, v);

    TS_ASSERT_DELTA(w[0], 1.0, 1.0e-14);
    TS_ASSERT_DELTA(w[1], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(w[2], 0.0, 1.0e-14);
  }

  void test_quaternion_rotate_vector_y_around_x()
  {
    double v[3] = {0.0, 1.0, 0.0};
    double w[3] = {0.0, 0.0, 0.0};

    const double q[4] = {std::sin(M_PI / 4), 0.0, 0.0, std::cos(M_PI / 4)};

    UTILS::quaternion_rotate_vector(w, q, v);

    TS_ASSERT_DELTA(w[0], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(w[1], 0.0, 1.0e-14);
    TS_ASSERT_DELTA(w[2], 1.0, 1.0e-14);
  }

  void test_quaternion_rotate_vector_general()
  {
    double v[3] = {0.5, 1.0, -2.0};
    double w[3] = {0.0, 0.0, 0.0};

    const double phi[3] = {-M_PI / 3, M_PI, M_PI / 2};

    double q[4];
    UTILS::quaternion_from_angle(q, phi);

    UTILS::quaternion_rotate_vector(w, q, v);

    TS_ASSERT_DELTA(w[0], 0.7145801717316358, 1.0e-14);
    TS_ASSERT_DELTA(w[1], -0.9159468817988596, 1.0e-14);
    TS_ASSERT_DELTA(w[2], 1.97494721141881, 1.0e-14);
  }
};

/*---------------------------------------------------------------------------*/
#endif
