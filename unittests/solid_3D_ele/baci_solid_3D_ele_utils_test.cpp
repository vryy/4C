/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests for the solid element utilities (e.g. stress/strain conversions)

\level 1
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "baci_solid_3D_ele_utils.hpp"

#include "baci_unittest_utils_assertions_test.hpp"

namespace
{
  using namespace FourC;

  CORE::LINALG::Matrix<3, 3> GetF()
  {
    CORE::LINALG::Matrix<3, 3> F(true);
    F(0, 0) = 1.1;
    F(0, 1) = 0.2;
    F(0, 2) = 0.5;
    F(1, 0) = 0.14;
    F(1, 1) = 1.2;
    F(1, 2) = 0.3;
    F(2, 0) = 0.05;
    F(2, 1) = 0.2;
    F(2, 2) = 1.3;

    return F;
  }

  TEST(TestStressStrainMeasures, GreenLagrangeToEulerAlmansi)
  {
    CORE::LINALG::Matrix<6, 1> green_lagrange_strain(
        std::array<double, 6>{0.11605, 0.26, 0.515, 0.398, 0.72, 0.657}.data());

    CORE::LINALG::Matrix<6, 1> euler_almansi_strain =
        STR::UTILS::GreenLagrangeToEulerAlmansi(green_lagrange_strain, GetF());

    CORE::LINALG::Matrix<6, 1> euler_almansi_strain_ref(
        std::array<double, 6>{0.055233442151184, 0.101134166403205, 0.104112596224498,
            0.182642289473823, 0.214768580862521, 0.315358749090858}
            .data());

    BACI_EXPECT_NEAR(euler_almansi_strain, euler_almansi_strain_ref, 1e-13);
  }

  TEST(TestStressStrainMeasures, GreenLagrangeToLogStrain)
  {
    CORE::LINALG::Matrix<6, 1> green_lagrange_strain(
        std::array<double, 6>{0.11605, 0.26, 0.515, 0.398, 0.72, 0.657}.data());

    CORE::LINALG::Matrix<6, 1> log_strain =
        STR::UTILS::GreenLagrangeToLogStrain(green_lagrange_strain);

    CORE::LINALG::Matrix<6, 1> log_strain_ref(
        std::array<double, 6>{0.039139830823291, 0.150129540734586, 0.281109187392933,
            0.218832208837098, 0.400808067245772, 0.400940161591198}
            .data());

    BACI_EXPECT_NEAR(log_strain, log_strain_ref, 1e-13);
  }

  TEST(TestStressStrainMeasures, SecondPiolaKirchhoffToCauchy)
  {
    CORE::LINALG::Matrix<6, 1> pk2(std::array<double, 6>{283.6946919505318, 195.86721709838096,
        202.01904686970775, 142.72731871521245, 182.86374040756576, 278.020938548381}
                                       .data());

    CORE::LINALG::Matrix<6, 1> cauchy(true);
    STR::UTILS::Pk2ToCauchy(pk2, GetF(), cauchy);

    CORE::LINALG::Matrix<6, 1> cauchy_ref(
        std::array<double, 6>{504.0646185061422, 317.85764952017706, 302.4131750725638,
            340.6815203116966, 306.97914008976466, 411.0514636046741}
            .data());

    BACI_EXPECT_NEAR(cauchy, cauchy_ref, 1e-12);
  }
}  // namespace