/*----------------------------------------------------------------------*/
/*! \file

\brief Collection of small local numerical methods.

\level 2
*/
/*----------------------------------------------------------------------*/

#include "baci_utils_local_numeric_methods.hpp"

#include "baci_utils_exceptions.hpp"

#include <cmath>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double CORE::UTILS::Bisection(const std::function<double(double)> &funct, const double a_init,
    const double b_init, const double tol, const int maxiter)
{
  double a = a_init;
  double b = b_init;
  double c;
  int numiter = 0;

  while (numiter <= maxiter)
  {
    numiter++;
    c = (a + b) / 2;

    auto f_c = funct(c);
    auto f_a = funct(a);

    if (std::fabs(f_c) < tol)
    {
      return c;
    }
    // sgn(f(a)) == sgn(f(c))
    else if (f_c * f_a > 0)
    {
      a = c;
    }
    else
    {
      b = c;
    }
  }

  dserror(
      "Maximal number of iterations reached for Bisection method (%d iterations). Error is "
      "still at %14.14f",
      numiter, (b - a) / 2);

  return c;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::UTILS::ValuesFunctAndFunctDerivs
CORE::UTILS::EvaluateFunctionAndDerivativesCentralDifferences(
    const std::function<double(double)> &func, const double x, const double delta_x)
{
  ValuesFunctAndFunctDerivs f_df_ddf;

  f_df_ddf.val_funct = func(x);

  double f_i_minus_1 = func(x - delta_x);
  double f_i_plus_1 = func(x + delta_x);
  f_df_ddf.val_deriv_funct = FirstDerivativeCentralDifferences(f_i_minus_1, f_i_plus_1, delta_x);
  f_df_ddf.val_deriv_deriv_funct =
      SecondDerivativeCentralDifferences(f_i_minus_1, f_df_ddf.val_funct, f_i_plus_1, delta_x);

  return f_df_ddf;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double CORE::UTILS::FirstDerivativeCentralDifferences(
    const double f_i_minus_1, const double f_i_plus_1, const double delta_x)
{
  double dfdx = (f_i_plus_1 - f_i_minus_1) / (2 * delta_x);
  return dfdx;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double CORE::UTILS::SecondDerivativeCentralDifferences(
    const double f_i_minus_1, const double f_i, const double f_i_plus_1, const double delta_x)
{
  double ddfddx = (f_i_plus_1 - 2 * f_i + f_i_minus_1) / (delta_x * delta_x);
  return ddfddx;
}

FOUR_C_NAMESPACE_CLOSE
