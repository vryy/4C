/*----------------------------------------------------------------------*/
/*! \file
 *
\brief Testcases for the core utils

\level 3
*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "4C_utils_local_numeric_methods.hpp"

#include "4C_unittest_utils_assertions_test.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  TEST(CoreUtilsTest, Bisection)
  {
    std::function<double(double)> function = [](double x) { return std::pow(x, 2) - 1; };

    double root = CORE::UTILS::Bisection(function, 0.0, 5.0, 1e-12, 200);

    EXPECT_NEAR(root, 1.0, 1e-12);
  }

  TEST(CoreUtilsTest, FirstDerivativeCentralDifferences)
  {
    const double x = 1.5;
    const double h = 0.001;

    // lets say f(x) = x^2 + 1
    const double f_xplus = std::pow(x + h, 2) + 1;
    const double f_xminus = std::pow(x - h, 2) + 1;

    // correct dfdx would be dfdx = 2x
    const double ref_dfdx = 2 * x;

    double dfdx = CORE::UTILS::FirstDerivativeCentralDifferences(f_xminus, f_xplus, h);

    EXPECT_NEAR(dfdx, ref_dfdx, std::pow(h, 2.0));
  }

  TEST(CoreUtilsTest, SecondDerivativeCentralDifferences)
  {
    const double x = 1.5;
    const double h = 0.001;

    // lets say f(x) = x^2 + 1
    const double f_x = std::pow(x, 2) + 1;
    const double f_xplus = std::pow(x + h, 2) + 1;
    const double f_xminus = std::pow(x - h, 2) + 1;

    // correct dfdx would be ddfddx = 2
    const double ref_ddfddx = 2;

    double ddfddx = CORE::UTILS::SecondDerivativeCentralDifferences(f_xminus, f_x, f_xplus, h);

    EXPECT_NEAR(ddfddx, ref_ddfddx, std::pow(h, 2.0));
  }
}  // namespace
FOUR_C_NAMESPACE_CLOSE
