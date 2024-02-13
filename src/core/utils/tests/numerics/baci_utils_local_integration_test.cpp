/*----------------------------------------------------------------------*/
/*! \file
 *
\brief Testcases for the core numerics integration library

\level 3
*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "baci_utils_local_integration.hpp"

#include <Sacado.hpp>
#include <Sacado_Fad_DFad.hpp>

BACI_NAMESPACE_OPEN

namespace
{
  using FADdouble = Sacado::Fad::DFad<double>;

  template <int coeff0, int coeff1, int coeff2>
  constexpr auto QUADRATIC_FUNCTION = [](auto x)
  { return std::make_tuple(x, coeff2 * std::pow(x, 2) + coeff1 * x + coeff0); };

  template <int coeff0, int coeff1>
  constexpr auto LINEAR_FUNCTION = [](auto x) { return std::make_tuple(x, coeff1 * x + coeff0); };

  TEST(CoreUtilsLocalIntegrationTest, IntegrateSimpsonStep)
  {
    constexpr auto f1 = QUADRATIC_FUNCTION<1, 2, 1>;
    auto value1 = CORE::UTILS::IntegrateSimpsonStep(
        0.5, std::get<1>(f1(0.0)), std::get<1>(f1(0.5)), std::get<1>(f1(1.0)));
    EXPECT_NEAR(value1, 7.0 / 3.0, 1e-8);

    constexpr auto f2 = QUADRATIC_FUNCTION<2, -2, 5>;
    auto value2 = CORE::UTILS::IntegrateSimpsonStep(
        0.5, std::get<1>(f2(0.0)), std::get<1>(f2(0.5)), std::get<1>(f2(1.0)));
    EXPECT_NEAR(value2, 8.0 / 3.0, 1e-8);
  }

  TEST(CoreUtilsLocalIntegrationTest, IntegrateSimpsonEquidistant)
  {
    constexpr auto f1 = QUADRATIC_FUNCTION<1, 2, 1>;
    auto value1 = CORE::UTILS::IntegrateSimpsonStep(f1(0.0), f1(0.5), f1(1.0));
    EXPECT_NEAR(value1, 7.0 / 3.0, 1e-8);

    constexpr auto f2 = QUADRATIC_FUNCTION<2, -2, 5>;
    auto value2 = CORE::UTILS::IntegrateSimpsonStep(f2(0.0), f2(0.5), f2(1.0));
    EXPECT_NEAR(value2, 8.0 / 3.0, 1e-8);
  }

  TEST(CoreUtilsLocalIntegrationTest, IntegrateSimpsonNonEquidistant)
  {
    constexpr auto f1 = QUADRATIC_FUNCTION<1, 2, 1>;
    auto value1 = CORE::UTILS::IntegrateSimpsonStep(f1(0.0), f1(0.9), f1(1.0));
    EXPECT_NEAR(value1, 7.0 / 3.0, 1e-8);

    constexpr auto f2 = QUADRATIC_FUNCTION<2, -2, 5>;
    auto value2 = CORE::UTILS::IntegrateSimpsonStep(f2(0.0), f2(0.9), f2(1.0));
    EXPECT_NEAR(value2, 8.0 / 3.0, 1e-8);
  }

  TEST(CoreUtilsLocalIntegrationTest, IntegrateSimpsonBCEquidistant)
  {
    constexpr auto f1 = QUADRATIC_FUNCTION<1, 2, 1>;
    auto value1 = CORE::UTILS::IntegrateSimpsonStepBC(f1(0.0), f1(0.5), f1(1.0));
    EXPECT_NEAR(value1, 37.0 / 24.0, 1e-8);

    constexpr auto f2 = QUADRATIC_FUNCTION<2, -2, 5>;
    auto value2 = CORE::UTILS::IntegrateSimpsonStepBC(f2(0.0), f2(0.5), f2(1.0));
    EXPECT_NEAR(value2, 41.0 / 24.0, 1e-8);
  }

  TEST(CoreUtilsLocalIntegrationTest, IntegrateSimpsonBCNonEquidistant)
  {
    constexpr auto f1 = QUADRATIC_FUNCTION<1, 2, 1>;
    auto value1 = CORE::UTILS::IntegrateSimpsonStepBC(f1(0.0), f1(0.9), f1(1.0));
    EXPECT_NEAR(value1, 1141.0 / 3000.0, 1e-8);

    constexpr auto f2 = QUADRATIC_FUNCTION<2, -2, 5>;
    auto value2 = CORE::UTILS::IntegrateSimpsonStepBC(f2(0.0), f2(0.9), f2(1.0));
    EXPECT_NEAR(value2, 277.0 / 600.0, 1e-8);
  }

  TEST(CoreUtilsLocalIntegrationTest, IntegrateSimpsonBCEquidistantWithDerivative)
  {
    constexpr auto f1 = QUADRATIC_FUNCTION<2, -2, 5>;
    auto [value, derivative] =
        CORE::UTILS::IntegrateSimpsonStepBCAndReturnDerivativeC(f1(0.0), f1(0.5), f1(1.0));
    EXPECT_NEAR(0.20833333333333333, derivative, 1e-8);
  }

  TEST(CoreUtilsLocalIntegrationTest, IntegrateSimpsonBCNonEquidistantWithDerivative)
  {
    constexpr auto f1 = QUADRATIC_FUNCTION<2, -2, 5>;
    auto [value, derivative] =
        CORE::UTILS::IntegrateSimpsonStepBCAndReturnDerivativeC(f1(0.0), f1(0.9), f1(1.0));
    EXPECT_NEAR(0.048333333333333318, derivative, 1e-8);
  }

  TEST(CoreUtilsLocalIntegrationTest, IntegrateSimpsonBCNonEquidistantWithDerivativeWithFad)
  {
    auto function_compute_derivative = [](double x) -> std::tuple<double, FADdouble>
    {
      if (x == 1.2)
      {
        return {x, FADdouble(1, 0, x)};
      }

      return {x, x};
    };
    FADdouble fad_derivative =
        CORE::UTILS::IntegrateSimpsonStepBC(function_compute_derivative(-0.1),
            function_compute_derivative(0.2), function_compute_derivative(1.2));


    constexpr auto f1 = QUADRATIC_FUNCTION<1, 2, 1>;

    auto [value, derivative] =
        CORE::UTILS::IntegrateSimpsonStepBCAndReturnDerivativeC(f1(-0.1), f1(0.2), f1(1.2));

    EXPECT_NEAR(value, 223.0 / 75.0, 1e-8);
    EXPECT_NEAR(fad_derivative.dx(0), derivative, 1e-8);
  }

  TEST(CoreUtilsLocalIntegrationTest, IntegrateTrapezoidal)
  {
    constexpr auto f = LINEAR_FUNCTION<1, 2>;
    auto value = CORE::UTILS::IntegrateTrapezoidalStep(f(0.0), f(1.0));
    EXPECT_NEAR(2, value, 1e-8);
  }

  TEST(CoreUtilsLocalIntegrationTest, IntegrateTrapezoidalDerivative)
  {
    constexpr auto f = LINEAR_FUNCTION<1, 2>;
    auto [value, derivative] =
        CORE::UTILS::IntegrateTrapezoidalStepAndReturnDerivativeB(f(0.0), f(1.0));
    EXPECT_NEAR(2, value, 1e-8);
    EXPECT_NEAR(0.5, derivative, 1e-8);
  }

  TEST(CoreUtilsLocalIntegrationTest, IntegrateTrapezoidalDerivativeWithFad)
  {
    auto function_compute_derivative = [](double x) -> std::tuple<double, FADdouble>
    {
      if (x == 1.2)
      {
        return {x, FADdouble(1, 0, x)};
      }

      return {x, x};
    };
    FADdouble fad_derivative = CORE::UTILS::IntegrateTrapezoidalStep(
        function_compute_derivative(-0.1), function_compute_derivative(1.2));
    constexpr auto f = LINEAR_FUNCTION<1, 2>;
    auto [value, derivative] =
        CORE::UTILS::IntegrateTrapezoidalStepAndReturnDerivativeB(f(-0.1), f(1.2));

    EXPECT_NEAR(value, 2.73, 1e-8);
    EXPECT_NEAR(fad_derivative.dx(0), derivative, 1e-8);
  }

  TEST(CoreUtilsLocalIntegrationTest, IntegrateIntegrateSimpsonTrapezoidal2Items)
  {
    std::array<double, 2> times = {0.0, 1.0};
    auto value = CORE::UTILS::IntegrateSimpsonTrapezoidal(times, LINEAR_FUNCTION<1, 2>);
    EXPECT_NEAR(2, value, 1e-8);
  }

  TEST(CoreUtilsLocalIntegrationTest, IntegrateIntegrateSimpsonTrapezoidalEvenItems)
  {
    std::array times = {0.0, 0.1, 0.3, 0.5, 0.88, 1.0};
    auto value = CORE::UTILS::IntegrateSimpsonTrapezoidal(times, QUADRATIC_FUNCTION<1, 2, 1>);
    EXPECT_NEAR(7.0 / 3.0, value, 1e-8);
  }

  TEST(CoreUtilsLocalIntegrationTest, IntegrateIntegrateSimpsonTrapezoidalOddItems)
  {
    std::array times = {0.0, 0.1, 0.3, 0.5, 0.88, 0.89, 1.0};
    auto value = CORE::UTILS::IntegrateSimpsonTrapezoidal(times, QUADRATIC_FUNCTION<1, 2, 1>);
    EXPECT_NEAR(7.0 / 3.0, value, 1e-8);
  }

  TEST(CoreUtilsLocalIntegrationTest, IntegrateIntegrateSimpsonTrapezoidalOneItem)
  {
    std::array times = {0.0};
    auto value = CORE::UTILS::IntegrateSimpsonTrapezoidal(times, QUADRATIC_FUNCTION<1, 2, 1>);
    EXPECT_NEAR(0, value, 1e-8);
  }

  TEST(CoreUtilsLocalIntegrationTest, IntegrateIntegrateSimpsonTrapezoidalZeroItem)
  {
    auto function = [](double x) -> std::tuple<double, double> {
      return {x, std::pow(x, 2) + 2 * x + 1};
    };
    std::array<double, 0> times = {};
    auto value = CORE::UTILS::IntegrateSimpsonTrapezoidal(times, function);
    EXPECT_NEAR(0, value, 1e-8);
  }

  TEST(CoreUtilsLocalIntegrationTest, IntegrateIntegrateSimpsonTrapezoidalMultipleIntervals)
  {
    auto function = [](double x) -> std::tuple<double, double>
    {
      if (x <= 1.0) return QUADRATIC_FUNCTION<1, 2, 1>(x);
      return QUADRATIC_FUNCTION<1, -2, 5>(x);
    };
    std::array times = {0.0, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0};
    auto value = CORE::UTILS::IntegrateSimpsonTrapezoidal(times, function);
    EXPECT_NEAR((7.0 + 29.0) / 3.0, value, 1e-8);
  }
}  // namespace
BACI_NAMESPACE_CLOSE
