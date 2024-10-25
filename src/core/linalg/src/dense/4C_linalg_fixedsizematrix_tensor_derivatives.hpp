// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later



#ifndef FOUR_C_LINALG_FIXEDSIZEMATRIX_TENSOR_DERIVATIVES_HPP
#define FOUR_C_LINALG_FIXEDSIZEMATRIX_TENSOR_DERIVATIVES_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg::Tensor
{
  /*!
   * @brief Add the derivative of the square of a tensor to a 4th order symmetric material tensor in
   * matrix notation
   *
   * In tensor index notation this method does
   * \f[
   * C_{IJKL} := \text{scalar_this} \cdot C_{IJKL}
   *           + \text{scalar_squared_dx} \cdot \frac{d (X^2)_{IJ}}{d X_{KL}}
   * \f]
   *
   * Compute the derivative of the square of a symmetric 2nd order tensor w.r.t. to the tensor add
   * the result to a 4th order tensor (in Voigt matrix notation!) using the symmetry-conditions
   * inherent to elasticity tensors.
   *
   * @param[in,out] C       Material tangent matrix to be modified
   * @param[in] scalar_squared_dx   Scalar to multiply with dX^2/dX
   * @param[in] X           Dense matrix (3 x 3)
   * @param[in] scalar_this  Scalar to multiply with C before adding dX^2/dX
   */
  void add_derivative_of_squared_tensor(Core::LinAlg::Matrix<6, 6>& C, double scalar_squared_dx,
      Core::LinAlg::Matrix<3, 3> X, double scalar_this);



  /*!
   * @brief Add scaled derivative of invA*B*invA w.r.t. A
   *
   * Add the following contribution to the tensor out(6,6) based on
   * - the inverse of the tensor \f$\mathbf{A}\f$
   * - the term invA*B*invA \f$\mathbf{A}^{-1} \mathbf{B} \mathbf{A}^{-1}\f$ with the tensor B
   *   \f$\mathbf{B}\f$
   *
   * \f[
   *    \text{scalar} \cdot \frac{\partial \mathbf{A}^{-1} \mathbf{B} \mathbf{A}^{-1}}{\partial
   *    \mathbf{A}}
   * \f]
   *
   * wherein the derivative \f$\frac{\partial \mathbf{A}^{-1} \mathbf{B} \mathbf{A}^{-1}}{\partial
   * \mathbf{A}}\f$ is computed to:
   * \f[
   *   - \frac{1}{2} \cdot \left( A^{-1}_{ik} A^{-1}_{jm} B_{mn} A^{-1}_{nl} + A^{-1}_{il}
   *    A^{-1}_{jm} B_{mn} A^{-1}_{nk} + A^{-1}_{jk} A^{-1}_{im} B_{mn} A^{-1}_{nl} + A^{-1}_{jl}
   *    A^{-1}_{im} B_{mn} A^{-1}_{nk} \right)
   * \f]
   *
   * @param[in] fac         Scaling factor
   * @param[in] invA        Inverse of the 2nd order tensor A stress-like Voigt notation
   * @param[in] invABinvA   2nd order tensor product invA*B*invA in stress-like Voigt notation
   * @param[in,out] out     4th order tensor derivative of invA*B*invA w.r.t. A in stress-like Voigt
   *                        notation
   */
  void add_derivative_of_inva_b_inva_product(double const& fac,
      const Core::LinAlg::Matrix<6, 1>& invA, const Core::LinAlg::Matrix<6, 1>& invABinvA,
      Core::LinAlg::Matrix<6, 6>& out);

}  // namespace Core::LinAlg::Tensor


FOUR_C_NAMESPACE_CLOSE

#endif