/*----------------------------------------------------------------------*/
/*! \file
 *
\brief Testcases for the core utils

\level 3
*/
/*----------------------------------------------------------------------*/

#include "gtest/gtest.h"
#include "utils_local_newton.H"

namespace
{
  TEST(CoreUtilsLocalNewtonTest, NewtonScalar)
  {
    auto residuum_and_jacobian = [](double x) -> std::tuple<double, double>
    {
      double residuum = 3 * std::pow(x, 3) + 2 * std ::pow(x, 2) + 4 * x + 1;
      double jacobian = 9 * std::pow(x, 2) + 4 * x + 4;
      return {residuum, jacobian};
    };

    auto [x_with_jac, jacobian] =
        CORE::UTILS::SolveLocalNewtonAndReturnJacobian(residuum_and_jacobian, 0.0, 1e-9);

    EXPECT_NEAR(x_with_jac, -0.271887376775884, 1e-8);
    EXPECT_NEAR(jacobian, 3.577755203747107557, 1e-8);

    auto x = CORE::UTILS::SolveLocalNewton(residuum_and_jacobian, 0.0, 1e-9);

    EXPECT_NEAR(x, -0.271887376775884, 1e-8);
  }


  TEST(CoreUtilsLocalNewtonTest, NewtonVector)
  {
    auto residuum_and_jacobian =
        [](const LINALG::Matrix<2, 1>& x) -> std::tuple<LINALG::Matrix<2, 1>, LINALG::Matrix<2, 2>>
    {
      LINALG::Matrix<2, 1> residuum(false);
      LINALG::Matrix<2, 2> jacobian(false);

      residuum(0) = std::exp(x(0) * x(1)) - 1;
      residuum(1) = 2 * std::exp(x(0) * x(1)) - x(0);

      jacobian(0, 0) = x(1) * std::exp(x(0) * x(1));
      jacobian(0, 1) = x(0) * std::exp(x(0) * x(1));

      jacobian(1, 0) = 2 * x(1) * std::exp(x(0) * x(1)) - 1;
      jacobian(1, 1) = 2 * x(0) * std::exp(x(0) * x(1));

      return {residuum, jacobian};
    };

    LINALG::Matrix<2, 1> x_0(true);
    x_0(0) = 1.1;
    x_0(1) = 0.1;

    auto [x_with_jac, jacobian] =
        CORE::UTILS::SolveLocalNewtonAndReturnJacobian(residuum_and_jacobian, x_0, 1e-9);

    EXPECT_NEAR(x_with_jac(0), 2, 1e-8);
    EXPECT_NEAR(x_with_jac(1), 0, 1e-8);


    auto x = CORE::UTILS::SolveLocalNewton(residuum_and_jacobian, x_0, 1e-9);

    EXPECT_NEAR(x(0), 2, 1e-8);
    EXPECT_NEAR(x(1), 0, 1e-8);
  }

  // Definition of our custom types:
  struct CustomScalarType
  {
    double value_ = 0.0;
    explicit CustomScalarType(double value) : value_(value) {}
  };

  bool operator>(CustomScalarType x, CustomScalarType y) { return x.value_ > y.value_; }

  struct CustomVectorType
  {
    double value_ = 0.0;
    explicit CustomVectorType(double value) : value_(value) {}
  };

  struct CustomJacobianType
  {
    double jacobian_ = 0.0;
    explicit CustomJacobianType(double jacobian) : jacobian_(jacobian) {}
  };

  CustomScalarType L2Norm(CustomVectorType x) { return CustomScalarType{std::abs(x.value_)}; }

  void LocalNewtonIteration(
      CustomVectorType& x, CustomVectorType residuum, CustomJacobianType jacobian)
  {
    x.value_ -= residuum.value_ / jacobian.jacobian_;
  }

  TEST(CoreUtilsLocalNewtonTest, NewtonCustomType)
  {
    auto residuum_and_jacobian =
        [](const CustomVectorType& x) -> std::tuple<CustomVectorType, CustomJacobianType>
    {
      CustomVectorType residuum{
          3 * std::pow(x.value_, 3) + 2 * std ::pow(x.value_, 2) + 4 * x.value_ + 1};
      CustomJacobianType jacobian{9 * std::pow(x.value_, 2) + 4 * x.value_ + 4};
      return {residuum, jacobian};
    };

    CustomVectorType x_0{0.0};
    auto [x_with_jac, jacobian] = CORE::UTILS::SolveLocalNewtonAndReturnJacobian(
        residuum_and_jacobian, x_0, CustomScalarType(1e-9));

    EXPECT_NEAR(x_with_jac.value_, -0.271887376775884, 1e-8);
    EXPECT_NEAR(jacobian.jacobian_, 3.577755203747107557, 1e-8);


    auto x = CORE::UTILS::SolveLocalNewton(residuum_and_jacobian, x_0, CustomScalarType(1e-9));
    EXPECT_NEAR(x.value_, -0.271887376775884, 1e-8);
  }
}  // namespace