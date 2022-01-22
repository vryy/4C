/*---------------------------------------------------------------------------*/
/*! \file

\brief unittests for the cubic spline interpolation

\level 1

*/


#ifndef UNIT_DRT_CUBIC_SPLINE_INTERPOLATION_H
#define UNIT_DRT_CUBIC_SPLINE_INTERPOLATION_H

#include "Unittests/src/common/unit_cxx_test_wrapper.H"
#include "../../../src/drt_lib/drt_cubic_spline_interpolation.H"

namespace DRT
{
  namespace UTILS
  {
    class CubicSplineInterpolation_TestSuite;
  }
}  // namespace DRT

class DRT::UTILS::CubicSplineInterpolation_TestSuite : public BACICxxTestWrapper
{
 public:
  void Setup() override
  {
    const std::vector<double> x = {0.30, 0.35, 0.40, 0.45};
    const std::vector<double> y = {4.40, 4.30, 4.25, 4.10};

    cubic_spline_ = Teuchos::rcp(new DRT::UTILS::CubicSplineInterpolation(x, y));
  }

  void TearDown() override { cubic_spline_ = Teuchos::null; }

  void TestInputArgumentSortedAscending()
  {
    const std::vector<double> x = {0.3, 0.6, 0.5};
    const std::vector<double> y(x.size(), 0.0);

    TS_ASSERT_THROWS(DRT::UTILS::CubicSplineInterpolation(x, y), const std::runtime_error &)
  }

  void TestInputArgumentsDifferentLength()
  {
    const std::vector<double> x = {0.3, 0.5, 0.7};
    const std::vector<double> y(x.size() - 1, 0.0);

    TS_ASSERT_THROWS(DRT::UTILS::CubicSplineInterpolation(x, y), const std::runtime_error &)
  }

  void TestEvaluateOutsideValidityBounds()
  {
    const std::vector<double> x_test = {0.2, 0.5};

    for (std::size_t i = 0; i < x_test.size(); ++i)
    {
      TS_ASSERT_THROWS(cubic_spline_->EvaluateScalar(x_test[i]), const std::runtime_error &)
      TS_ASSERT_THROWS(
          cubic_spline_->EvaluateScalarFirstDerivative(x_test[i]), const std::runtime_error &)
      TS_ASSERT_THROWS(
          cubic_spline_->EvaluateScalarSecondDerivative(x_test[i]), const std::runtime_error &)
    }
  }

  void TestEvaluateScalar()
  {
    const std::vector<double> x_test = {0.33, 0.36, 0.4, 0.42};
    const std::vector<double> reference_solution = {4.33232, 4.29, 4.25, 4.20152};

    for (std::size_t i = 0; i < x_test.size(); ++i)
      TS_ASSERT_DELTA(cubic_spline_->EvaluateScalar(x_test[i]), reference_solution[i], 1.0e-12)
  }

  void TestEvaluateScalarFirstDerivative()
  {
    const std::vector<double> x_test = {0.33, 0.36, 0.4, 0.42};
    const std::vector<double> reference_solution = {-1.968, -0.84, -1.8, -2.952};

    for (std::size_t i = 0; i < x_test.size(); ++i)
      TS_ASSERT_DELTA(
          cubic_spline_->EvaluateScalarFirstDerivative(x_test[i]), reference_solution[i], 1.0e-12)
  }

  void TestEvaluateScalarSecondDerivative()
  {
    const std::vector<double> x_test = {0.33, 0.36, 0.4, 0.42};
    const std::vector<double> reference_solution = {28.8, 24.0, -72.0, -43.2};

    for (std::size_t i = 0; i < x_test.size(); ++i)
      TS_ASSERT_DELTA(
          cubic_spline_->EvaluateScalarSecondDerivative(x_test[i]), reference_solution[i], 1.0e-12)
  }

 private:
  Teuchos::RCP<CubicSplineInterpolation> cubic_spline_;
};

#endif
