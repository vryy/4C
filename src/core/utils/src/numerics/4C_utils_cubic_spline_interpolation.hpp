/*----------------------------------------------------------------------*/
/*! \file

\brief Cubic spline interpolation

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_UTILS_CUBIC_SPLINE_INTERPOLATION_HPP
#define FOUR_C_UTILS_CUBIC_SPLINE_INTERPOLATION_HPP

#include "4C_config.hpp"

#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/

namespace CORE::UTILS
{
  /**
   * @brief Cubic spline interpolation based on two input vectors \f$ \vec{x}, \vec{y} \f$
   *
   * This class defines a cubic spline interpolation connecting \f$ n \f$ points
   * \f$ (x_i,y_i), i \in [0,1,...,n-1] \f$ expressed as splines \f$S_i(x)\f$ of the form:
   * \f[
   * S_i(x) = a_i + b_i(x-x_i) + c_i(x-x_i)^2 + d_i(x-x_i)^3, j \in [0,1,...,n-1],
   * \f]
   * where \f$ S_i(x) \f$ is a third order polynomial on the interval \f$ [x_i, x_{i+1}] \f$.
   *
   * @note The spline is calculated such that the first and second derivatives match at the
   * intersection of two segments. Furthermore, so-called natural boundary conditions, i.e. the
   * second derivatives at the two boundaries are set to zero, are applied.
   */
  class CubicSplineInterpolation
  {
   public:
    /*!
     * @brief Constructor of function defining a cubic spline interpolation created based on the
     * two input vectors \f$ \vec{x}, \vec{y} \f$
     *
     * @param[in] x  argument vector \f$ \vec{x} \f$
     * @param[in] y  value vector \f$ \vec{y} \f$
     */
    CubicSplineInterpolation(std::vector<double> x, std::vector<double> y);

    /*!
     * @brief Evaluate the scalar function
     *
     * @param[in] x  generic scalar
     * @return value of function evaluated at @x
     */
    [[nodiscard]] double Evaluate(double x) const;

    /*!
     * @brief Evaluate the first derivative of the scalar function
     *
     * @param[in] x  generic scalar
     * @return @deriv_order derivative of the function evaluated at @x
     */
    [[nodiscard]] double EvaluateDerivative(double x, int deriv_order) const;

   private:
    /*!
     * @brief Method that creates the matrix and right-hand side of the linear system that has to
     * be solved to calculate the polynomial coefficients of the cubic spline interpolation:
     * \f[
     *     \boldsymbol{A} \vec{c} = \vec{b}
     * \f]
     * with
     * \f[
     *    \boldsymbol{A} = \begin{bmatrix}
     *    2        & 0        &         &              &        & 0           \\
     *    \alpha_1 & 2        & \beta_1 &              &        & \vdots      \\
     *    0        & \alpha_2 &    2    & \beta_2      &        &             \\
     *             &          & \ddots  &  \ddots      & \ddots & 0           \\
     *    \vdots   &          &         & \alpha_{n-2} & 2      & \beta_{n-2} \\
     *    0        &          &  \dots  &              & 0      & 2           \\
     *    \end{bmatrix}
     * \f]
     * and
     * \f[
     *    \vec{b} = \begin{bmatrix} 0 & \gamma_1 & \dots & \gamma_{n-2} & 0 \end{bmatrix}^T
     * \f]
     * with #x_ (\f$ \vec{x} \f$) and #a_ (\f$ \vec{y} \f$):
     *      \f$ \alpha_i = \frac{x_i - x_{i-1}}{x_{i+1} - x_{i-1}} \f$,
     *      \f$ \beta_i  = \frac{x_{i+1} - x_i}{x_{i+1} - x_{i-1}} \f$, and
     *      \f$ \gamma_i = \frac{3}{x_{i+1} - x_{i-1}} \left( \frac{y_{i+1}-y_i}{x_{i+1}-x_i} -
     *      \frac{y_i-y_{i-1}}{x_i-x_{i-1}} \right) \f$
     *
     * @param[in] N      size of linear system
     * @param[out] A     matrix of linear system that defines the polynomial coefficients of the
     *                   cubic spline interpolation
     * @param[out] b     right-hand side vector of linear system that defines the polynomial
     *                   coefficients of the cubic spline interpolation
     *
     * @note The first and last line of the implementation represent the so-called natural
     * boundary conditions, i.e. that the second derivative equals zero at the boundary
     */
    void BuildMatrixAndRhs(
        const int N, CORE::LINALG::SerialDenseMatrix& A, CORE::LINALG::SerialDenseVector& b) const;

    /*!
     * @brief Calculates the coefficient vectors
     *
     * The coefficient vectors #b_ (\f$\vec{b}\f$), #c_ (\f$\vec{c}\f$), and #d_ (\f$\vec{d}\f$)
     * are calculated using #x_ (\f$\vec{x}\f$), #a_ (\f$\vec{y}\f$), and
     * @p c (\f$\vec{c}_{in}\f$) as follows:
     * \f$ b_i = \frac{y_{i+1}-y_i}{x_{i+1}-x_i}-1/3\cdot (x_{i+1}-x_i) \cdot (2 c_i+c_{i+1}) \f$;
     * \f$ \vec{c} = \vec{c}_{in} \f$;
     * \f$ d_i = \frac{c_{i+1}-c_i}{3(x_{i+1}-x_i)} \f$;
     *
     * @param[in] c  solution vector of linear system that defines the polynomial coefficients of
     *               the cubic spline interpolation, i.e. coefficient vector #c_
     */
    void SetupInternalVectors(const CORE::LINALG::SerialDenseVector& c);

    /*!
     * @brief solves linear system \f$ \boldsymbol{A} \vec{c} = \vec{b} \f$
     *
     * @param[in]  A matrix
     * @param[out] c solution vector, equals coefficient vector #c_
     * @param[in]  b right-hand side vector
     */
    void SolveLinearSystem(CORE::LINALG::SerialDenseMatrix& A, CORE::LINALG::SerialDenseVector& c,
        CORE::LINALG::SerialDenseVector& b) const;

    //! zeroth-order coefficients for cubic spline interpolation
    std::vector<double> a_;

    //! first-order coefficients for cubic spline interpolation
    std::vector<double> b_;

    //! second-order coefficients for cubic spline interpolation
    std::vector<double> c_;

    //! third-order coefficients for cubic spline interpolation
    std::vector<double> d_;

    //! sampling points for cubic spline interpolation
    std::vector<double> x_;
  };
}  // namespace CORE::UTILS

FOUR_C_NAMESPACE_CLOSE

#endif