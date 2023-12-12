/*----------------------------------------------------------------------*/
/*! \file

\brief Cubic spline interpolation

\level 1

*/
/*----------------------------------------------------------------------*/

#include "baci_utils_cubic_spline_interpolation.H"

#include "baci_utils_exceptions.H"

#include <Teuchos_SerialDenseSolver.hpp>

#include <utility>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
CORE::UTILS::CubicSplineInterpolation::CubicSplineInterpolation(
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
  CORE::LINALG::SerialDenseMatrix A(N, N);
  CORE::LINALG::SerialDenseVector c(N), b(N);

  BuildMatrixAndRhs(N, A, b);
  SolveLinearSystem(A, c, b);
  SetupInternalVectors(c);
}

/*----------------------------------------------------------------------*/
void CORE::UTILS::CubicSplineInterpolation::BuildMatrixAndRhs(
    const int N, CORE::LINALG::SerialDenseMatrix &A, CORE::LINALG::SerialDenseVector &b) const
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
double CORE::UTILS::CubicSplineInterpolation::Evaluate(const double x) const
{
  // safety check
  if (x < x_.front() or x > x_.back())
    dserror("Sampling point x = %lf lies outside sampling point range!", x);

  auto greater_equal_x = [x](const double val) { return val >= x; };
  auto right_position = std::find_if(x_.begin(), x_.end(), greater_equal_x);
  dsassert(right_position != x_.begin(), "Internal error.");
  // the left side of the sought interval is found by deleting 1
  const auto left_position = std::distance(x_.begin(), right_position) - 1;

  const double delta_x = x - x_[left_position];
  return a_[left_position] + b_[left_position] * delta_x + c_[left_position] * delta_x * delta_x +
         d_[left_position] * delta_x * delta_x * delta_x;
}

/*----------------------------------------------------------------------*/
double CORE::UTILS::CubicSplineInterpolation::EvaluateDerivative(
    const double x, const int deriv_order) const
{
  // safety check
  if (x < x_.front() or x > x_.back())
    dserror("Sampling point x = %lf lies outside sampling point range!", x);

  auto greater_equal_x = [x](const double val) { return val >= x; };
  auto right_position = std::find_if(x_.begin(), x_.end(), greater_equal_x);
  dsassert(right_position != x_.begin(), "Internal error.");
  // the left side of the sought interval is found by deleting 1
  const auto left_position = std::distance(x_.begin(), right_position) - 1;

  switch (deriv_order)
  {
    case 1:
    {
      const double delta_x = x - x_[left_position];
      return b_[left_position] + 2.0 * c_[left_position] * delta_x +
             3.0 * d_[left_position] * delta_x * delta_x;
    }
    case 2:
    {
      const double delta_x = x - x_[left_position];
      return 2.0 * c_[left_position] + 6.0 * d_[left_position] * delta_x;
    }
  }

  dserror(
      "Evaluation of %i derivative is not implemented for CubicSplineInterpolation", deriv_order);
}

/*----------------------------------------------------------------------*/
void CORE::UTILS::CubicSplineInterpolation::SolveLinearSystem(CORE::LINALG::SerialDenseMatrix &A,
    CORE::LINALG::SerialDenseVector &c, CORE::LINALG::SerialDenseVector &b) const
{
  // solve for third-order coefficients for cubic spline interpolation
  using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
  using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
  Teuchos::SerialDenseSolver<ordinalType, scalarType> solver;
  solver.setMatrix(Teuchos::rcpFromRef(A));
  solver.setVectors(Teuchos::rcpFromRef(c), Teuchos::rcpFromRef(b));
  solver.factorWithEquilibration(true);
  solver.solveToRefinedSolution(true);
  if (solver.factor() or solver.solve()) dserror("Solution of linear system of equations failed!");
}

/*----------------------------------------------------------------------*/
void CORE::UTILS::CubicSplineInterpolation::SetupInternalVectors(
    const CORE::LINALG::SerialDenseVector &c)
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
BACI_NAMESPACE_CLOSE
