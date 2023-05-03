/*----------------------------------------------------------------------*/
/*! \file
 *
\brief Testcases for the core utils

\level 3
*/
/*----------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "utils_local_numeric_methods.H"
#include "unittests_assertions.h"

namespace
{
  TEST(CoreUtilsTest, NewtonScalar)
  {
    std::function<CORE::UTILS::ValuesFunctAndFunctDeriv(double)> function = [](double x)
    {
      return CORE::UTILS::ValuesFunctAndFunctDeriv{
          .val_funct = std::pow(x, 2) - 1, .val_deriv_funct = 2 * x};
    };

    double root = CORE::UTILS::NewtonScalar(function, 5.0, 1e-12, 200);

    EXPECT_NEAR(root, 1.0, 1e-12);
  }

  TEST(CoreUtilsTest, Bisection)
  {
    std::function<CORE::UTILS::ValuesFunctAndFunctDeriv(double)> function = [](double x)
    {
      return CORE::UTILS::ValuesFunctAndFunctDeriv{
          .val_funct = std::pow(x, 2) - 1, .val_deriv_funct = 0.0};
    };

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