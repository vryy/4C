/*----------------------------------------------------------------------*/
/*! \file

\brief Cubic spline interpolation

\level 1

*/
/*----------------------------------------------------------------------*/

#include <Epetra_SerialDenseSolver.h>
#include <utility>

#include "lib_cubic_spline_interpolation.H"
#include "lib_dserror.H"

/*----------------------------------------------------------------------*/
DRT::UTILS::CubicSplineInterpolation::CubicSplineInterpolation(
    std::vector<double> x, std::vector<double> y)
    : a_(std::move(y)), x_(std::move(x))
// zeroth-order coefficients a_ equal the function values y
{
  // safety checks
  if (x_.size() != a_.size())
    dserror("Length of vectors provided to cubic spline interpolation must match!");
  if (not std::is_sorted(x_.begin(), x_.end()))
    dserror("Data points must be sorted in ascending order!");

  // temp variables for solution of the linear system
  const int N = static_cast<int>(x_.size());
  Epetra_SerialDenseMatrix A(N, N);
  Epetra_SerialDenseVector c(N), b(N);

  BuildMatrixAndRhs(N, A, b);
  SolveLinearSystem(A, c, b);
  SetupInternalVectors(c);
}

/*----------------------------------------------------------------------*/
void DRT::UTILS::CubicSplineInterpolation::BuildMatrixAndRhs(
    const int N, Epetra_SerialDenseMatrix &A, Epetra_SerialDenseVector &b) const
{
  // fill everything except the boundary condition lines
  for (int i = 1; i < N - 1; ++i)
  {
    // a_ = y, see initialization list in constructor
    const double delta_y_p = a_[i + 1] - a_[i];
    const double delta_y_m = a_[i] - a_[i - 1];
    const double delta_x_p = x_[i + 1] - x_[i];
    const double delta_x_m = x_[i] - x_[i - 1];
    const double delta_x_pm = x_[i + 1] - x_[i - 1];

    const double alpha_i = delta_x_m / delta_x_pm;
    const double beta_i = delta_x_p / delta_x_pm;
    const double gamma_i = 3.0 / delta_x_pm * ((delta_y_p / delta_x_p) - (delta_y_m / delta_x_m));

    A(i, i - 1) = alpha_i;
    A(i, i) = 2.0;
    A(i, i + 1) = beta_i;
    b(i) = gamma_i;
  }

  // set values due to natural boundary conditions
  A(0, 0) = A(N - 1, N - 1) = 2.0;
  b(0) = b(N - 1) = 0.0;
}

/*----------------------------------------------------------------------*/
double DRT::UTILS::CubicSplineInterpolation::EvaluateScalar(const double x) const
{
  // safety check
  if (x < x_.front() or x > x_.back())
    dserror("Sampling point x = %lf lies outside sampling point range!", x);

  double value(0.0);

  // evaluate cubic spline interpolation
  for (std::size_t i = 0; i < x_.size() - 1; ++i)
  {
    if (x <= x_[i + 1])
    {
      const double delta_x = x - x_[i];
      value =
          a_[i] + b_[i] * delta_x + c_[i] * std::pow(delta_x, 2.0) + d_[i] * std::pow(delta_x, 3.0);
      break;
    }
  }

  return value;
}

/*----------------------------------------------------------------------*/
double DRT::UTILS::CubicSplineInterpolation::EvaluateScalarFirstDerivative(const double x) const
{
  // safety check
  if (x < x_.front() or x > x_.back())
    dserror("Sampling point x = %lf lies outside sampling point range!", x);

  double first_derivative(0.0);

  // evaluate cubic spline interpolation
  for (std::size_t i = 0; i < x_.size() - 1; ++i)
  {
    if (x <= x_[i + 1])
    {
      const double delta_x = x - x_[i];
      first_derivative = b_[i] + 2.0 * c_[i] * delta_x + 3.0 * d_[i] * std::pow(delta_x, 2.0);
      break;
    }
  }

  return first_derivative;
}

/*----------------------------------------------------------------------*/
double DRT::UTILS::CubicSplineInterpolation::EvaluateScalarSecondDerivative(const double x) const
{
  // safety check
  if (x < x_.front() or x > x_.back())
    dserror("Sampling point x = %lf lies outside sampling point range!", x);

  double second_derivative(0.0);

  // evaluate cubic spline interpolation
  for (std::size_t i = 0; i < x_.size() - 1; ++i)
  {
    if (x <= x_[i + 1])
    {
      const double delta_x = x - x_[i];
      second_derivative = 2.0 * c_[i] + 6.0 * d_[i] * delta_x;
      break;
    }
  }

  return second_derivative;
}

/*----------------------------------------------------------------------*/
void DRT::UTILS::CubicSplineInterpolation::SolveLinearSystem(
    Epetra_SerialDenseMatrix &A, Epetra_SerialDenseVector &c, Epetra_SerialDenseVector &b) const
{
  // solve for third-order coefficients for cubic spline interpolation
  Epetra_SerialDenseSolver solver;
  solver.SetMatrix(A);
  solver.SetVectors(c, b);
  solver.FactorWithEquilibration(true);
  solver.SolveToRefinedSolution(true);
  if (solver.Factor() or solver.Solve()) dserror("Solution of linear system of equations failed!");
}

/*----------------------------------------------------------------------*/
void DRT::UTILS::CubicSplineInterpolation::SetupInternalVectors(const Epetra_SerialDenseVector &c)
{
  const std::size_t system_size = x_.size();

  // resize of the internal member vectors
  c_.resize(system_size, 0.0);
  b_.resize(system_size - 1, 0.0), d_.resize(system_size - 1, 0.0);

  for (auto i = 0; i < static_cast<int>(system_size); ++i) c_[i] = c(i);

  for (std::size_t i = 0; i < system_size - 1; ++i)
  {
    // a_ = y, see initialization list in constructor
    const double delta_y_p = a_[i + 1] - a_[i];
    const double delta_x_p = x_[i + 1] - x_[i];

    b_[i] = delta_y_p / delta_x_p - 1.0 / 3.0 * delta_x_p * (2.0 * c_[i] + c_[i + 1]);
    d_[i] = (c_[i + 1] - c_[i]) / (3.0 * delta_x_p);
  }
}