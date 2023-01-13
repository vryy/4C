/*---------------------------------------------------------------------------*/
/*! \file

\brief Unittests for the cubic spline interpolation

\level 1
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include <Teuchos_RCP.hpp>
#include "lib_cubic_spline_interpolation.H"

namespace
{
  class CubicSplineInterpolationTest : public ::testing::Test
  {
   protected:
    CubicSplineInterpolationTest()
    {
      const std::vector<double> x = {0.30, 0.35, 0.40, 0.45};
      const std::vector<double> y = {4.40, 4.30, 4.25, 4.10};

      cubic_spline_ = Teuchos::rcp(new DRT::UTILS::CubicSplineInterpolation(x, y));
    }

    Teuchos::RCP<DRT::UTILS::CubicSplineInterpolation> cubic_spline_;
  };

  TEST_F(CubicSplineInterpolationTest, InputArgumentSortedAscending)
  {
    const std::vector<double> x = {0.3, 0.6, 0.5};
    const std::vector<double> y(x.size(), 0.0);

    EXPECT_THROW(DRT::UTILS::CubicSplineInterpolation(x, y), std::runtime_error);
  }

  TEST_F(CubicSplineInterpolationTest, InputArgumentsDifferentLength)
  {
    const std::vector<double> x = {0.3, 0.5, 0.7};
    const std::vector<double> y(x.size() - 1, 0.0);

    EXPECT_THROW(DRT::UTILS::CubicSplineInterpolation(x, y), std::runtime_error);
  }

  TEST_F(CubicSplineInterpolationTest, EvaluateOutsideValidityBounds)
  {
    const std::vector<double> x_test = {0.2, 0.5};

    for (std::size_t i = 0; i < x_test.size(); ++i)
    {
      EXPECT_THROW(cubic_spline_->EvaluateScalar(x_test[i]), std::runtime_error);
      EXPECT_THROW(cubic_spline_->EvaluateScalarFirstDerivative(x_test[i]), std::runtime_error);
      EXPECT_THROW(cubic_spline_->EvaluateScalarSecondDerivative(x_test[i]), std::runtime_error);
    }
  }

  TEST_F(CubicSplineInterpolationTest, EvaluateScalar)
  {
    const std::vector<double> x_test = {0.33, 0.36, 0.4, 0.42};
    const std::vector<double> reference_solution = {4.33232, 4.29, 4.25, 4.20152};

    for (std::size_t i = 0; i < x_test.size(); ++i)
      EXPECT_NEAR(cubic_spline_->EvaluateScalar(x_test[i]), reference_solution[i], 1.0e-12);
  }

  TEST_F(CubicSplineInterpolationTest, EvaluateScalarFirstDerivative)
  {
    const std::vector<double> x_test = {0.33, 0.36, 0.4, 0.42};
    const std::vector<double> reference_solution = {-1.968, -0.84, -1.8, -2.952};

    for (std::size_t i = 0; i < x_test.size(); ++i)
      EXPECT_NEAR(
          cubic_spline_->EvaluateScalarFirstDerivative(x_test[i]), reference_solution[i], 1.0e-12);
  }

  TEST_F(CubicSplineInterpolationTest, EvaluateScalarSecondDerivative)
  {
    const std::vector<double> x_test = {0.33, 0.36, 0.4, 0.42};
    const std::vector<double> reference_solution = {28.8, 24.0, -72.0, -43.2};

    for (std::size_t i = 0; i < x_test.size(); ++i)
      EXPECT_NEAR(
          cubic_spline_->EvaluateScalarSecondDerivative(x_test[i]), reference_solution[i], 1.0e-12);
  }
}  // namespace