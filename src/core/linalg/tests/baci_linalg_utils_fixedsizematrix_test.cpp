/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for the fixed size matrix.

\level 0
*/


#include <gtest/gtest.h>

#include "baci_linalg_fixedsizematrix.hpp"
#include "baci_utils_fad.hpp"

BACI_NAMESPACE_OPEN

namespace
{
  auto GetTestValuesAssignmentOperators()
  {
    static const unsigned int n_dim = 2;
    using fad_type = Sacado::Fad::DFad<double>;

    CORE::LINALG::Matrix<n_dim, 1, double> X;
    X(0) = 3.0;
    X(1) = 4.1;

    CORE::LINALG::Matrix<n_dim, 1, fad_type> u;
    u(0) = CORE::FADUTILS::HigherOrderFadValue<fad_type>::apply(n_dim, 0, 0.4);
    u(1) = CORE::FADUTILS::HigherOrderFadValue<fad_type>::apply(n_dim, 1, 0.3);

    CORE::LINALG::Matrix<n_dim, 1, fad_type> x_ref;
    x_ref(0) = CORE::FADUTILS::HigherOrderFadValue<fad_type>::apply(n_dim, 0, 3.4);
    x_ref(1) = CORE::FADUTILS::HigherOrderFadValue<fad_type>::apply(n_dim, 1, 4.4);

    return std::tuple{X, u, x_ref};
  };

  template <typename T, typename V>
  void CheckTestResults(const T& x, const V& x_ref)
  {
    const double eps = 1e-12;
    const unsigned int n_dim = x.numRows();

    // Check the values of the result as well as the first derivatives
    for (unsigned int i = 0; i < n_dim; i++)
    {
      EXPECT_NEAR(CORE::FADUTILS::CastToDouble(x(i)), CORE::FADUTILS::CastToDouble(x_ref(i)), eps);
      for (unsigned int j = 0; j < (unsigned int)x(i).length(); j++)
      {
        EXPECT_NEAR(CORE::FADUTILS::CastToDouble(x(i).dx(j)),
            CORE::FADUTILS::CastToDouble(x_ref(i).dx(j)), eps);
      }
    }
  }

  TEST(FixedSizeMatrixTest, AssignmentOperatorPlusEqualDifferentTypes)
  {
    const auto [X, u, x_ref] = GetTestValuesAssignmentOperators();
    using result_type = typename std::decay<decltype(x_ref)>::type;
    result_type x;

    x = u;
    x += X;

    CheckTestResults(x, x_ref);
  }

  TEST(FixedSizeMatrixTest, AssignmentOperatorMinusEqualDifferentTypes)
  {
    const auto [X, u, x_ref] = GetTestValuesAssignmentOperators();
    using result_type = typename std::decay<decltype(x_ref)>::type;
    result_type x;

    x = u;
    x.Scale(-1.0);
    x -= X;
    x.Scale(-1.0);

    CheckTestResults(x, x_ref);
  }

  TEST(FixedSizeMatrixTest, UpdateDifferentTypes)
  {
    const auto [X, u, x_ref] = GetTestValuesAssignmentOperators();
    using result_type = typename std::decay<decltype(x_ref)>::type;
    result_type x;

    x.Update(X);
    x += u;

    CheckTestResults(x, x_ref);
  }

  TEST(FixedSizeMatrixTest, MultiplyDifferentTypes)
  {
    using fad_type = Sacado::Fad::DFad<double>;

    CORE::LINALG::Matrix<2, 4> shape_function_matrix(true);
    shape_function_matrix(0, 0) = 0.75;
    shape_function_matrix(1, 1) = 0.75;
    shape_function_matrix(0, 2) = 0.25;
    shape_function_matrix(1, 3) = 0.25;

    CORE::LINALG::Matrix<4, 1, fad_type> nodal_dof;
    nodal_dof(0) = CORE::FADUTILS::HigherOrderFadValue<fad_type>::apply(4, 0, 0.4);
    nodal_dof(1) = CORE::FADUTILS::HigherOrderFadValue<fad_type>::apply(4, 1, 1.4);
    nodal_dof(2) = CORE::FADUTILS::HigherOrderFadValue<fad_type>::apply(4, 2, 2.4);
    nodal_dof(3) = CORE::FADUTILS::HigherOrderFadValue<fad_type>::apply(4, 3, 3.4);

    CORE::LINALG::Matrix<2, 1, fad_type> u;
    u.Multiply(shape_function_matrix, nodal_dof);

    CORE::LINALG::Matrix<2, 1, fad_type> u_ref;
    u_ref(0) = CORE::FADUTILS::HigherOrderFadValue<fad_type>::apply(4, 0, 0.9);
    u_ref(1) = CORE::FADUTILS::HigherOrderFadValue<fad_type>::apply(4, 1, 1.9);
    u_ref(0).fastAccessDx(0) = 0.75;
    u_ref(0).fastAccessDx(2) = 0.25;
    u_ref(1).fastAccessDx(1) = 0.75;
    u_ref(1).fastAccessDx(3) = 0.25;

    CheckTestResults(u, u_ref);
  }
}  // namespace

BACI_NAMESPACE_CLOSE
