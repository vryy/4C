/*! \file
\level 1
\brief Util functions for tensor transformations
*/

#ifndef FOUR_C_LINALG_FIXEDSIZEMATRIX_TENSOR_TRANSFORMATION_HPP
#define FOUR_C_LINALG_FIXEDSIZEMATRIX_TENSOR_TRANSFORMATION_HPP

#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINALG::TENSOR
{
  /*!
   * @brief Rotation of all basis vectors of 2-Tensor Tin by the rotation Matrix Q
   *
   * @tparam n dim
   * @param Q (in) : Rotation Matrix for each basis vector.
   * @param Tin (in) : Unrotated tensor
   * @param Tout (out) : Rotated tensor
   */
  template <unsigned int n>
  inline void TensorRotation(const CORE::LINALG::Matrix<n, n>& Q,
      const CORE::LINALG::Matrix<n, n>& Tin, CORE::LINALG::Matrix<n, n>& Tout)
  {
    CORE::LINALG::Matrix<n, n> temp(false);
    temp.MultiplyNN(1.0, Q, Tin);
    Tout.MultiplyNT(1.0, temp, Q);
  }

  /*!
   * @brief Inverse rotation of all basis vectors of 2-Tensor Tin by the rotation Matrix Q
   *
   * @tparam n dim
   * @param Q (in) : Rotation Matrix for each basis vector.
   * @param Tin (in) : Rotated tensor
   * @param Tout (out) : Unrotated tensor
   */
  template <unsigned int n>
  inline void InverseTensorRotation(const CORE::LINALG::Matrix<n, n>& Q,
      const CORE::LINALG::Matrix<n, n>& Tin, CORE::LINALG::Matrix<n, n>& Tout)
  {
    CORE::LINALG::Matrix<n, n> temp(false);
    temp.MultiplyTN(1.0, Q, Tin);
    Tout.MultiplyNN(1.0, temp, Q);
  }

  /*!
   * @brief Evaluate transformation matrix for fourth order tensors in Voigt notation
   *
   * This transformation could for example be used to rotate the linearization tensor of 2.
   * Piola-Kirchhoff stresses w.r.t. right Cauchy Green strains by the rotation matrix Q.
   *
   * @param Q  (in) : Rotation matrix for each basis vector
   * @param Qfourth (out) : Transformation matrix for fourth order tensors in Voigt notation
   */
  inline void FourthTensorRotationMatrix(
      const CORE::LINALG::Matrix<3, 3>& Q, CORE::LINALG::Matrix<6, 6>& Qfourth)
  {
    Qfourth(0, 0) = Q(0, 0) * Q(0, 0);
    Qfourth(0, 1) = Q(0, 1) * Q(0, 1);
    Qfourth(0, 2) = Q(0, 2) * Q(0, 2);
    Qfourth(0, 3) = Q(0, 0) * Q(0, 1);
    Qfourth(0, 4) = Q(0, 1) * Q(0, 2);
    Qfourth(0, 5) = Q(0, 0) * Q(0, 2);

    Qfourth(1, 0) = Q(1, 0) * Q(1, 0);
    Qfourth(1, 1) = Q(1, 1) * Q(1, 1);
    Qfourth(1, 2) = Q(1, 2) * Q(1, 2);
    Qfourth(1, 3) = Q(1, 0) * Q(1, 1);
    Qfourth(1, 4) = Q(1, 1) * Q(1, 2);
    Qfourth(1, 5) = Q(1, 0) * Q(1, 2);

    Qfourth(2, 0) = Q(2, 0) * Q(2, 0);
    Qfourth(2, 1) = Q(2, 1) * Q(2, 1);
    Qfourth(2, 2) = Q(2, 2) * Q(2, 2);
    Qfourth(2, 3) = Q(2, 0) * Q(2, 1);
    Qfourth(2, 4) = Q(2, 1) * Q(2, 2);
    Qfourth(2, 5) = Q(2, 0) * Q(2, 2);

    Qfourth(3, 0) = 2.0 * Q(0, 0) * Q(1, 0);
    Qfourth(3, 1) = 2.0 * Q(0, 1) * Q(1, 1);
    Qfourth(3, 2) = 2.0 * Q(0, 2) * Q(1, 2);
    Qfourth(3, 3) = Q(0, 0) * Q(1, 1) + Q(1, 0) * Q(0, 1);
    Qfourth(3, 4) = Q(0, 1) * Q(1, 2) + Q(1, 1) * Q(0, 2);
    Qfourth(3, 5) = Q(0, 0) * Q(1, 2) + Q(1, 0) * Q(0, 2);

    Qfourth(4, 0) = 2.0 * Q(1, 0) * Q(2, 0);
    Qfourth(4, 1) = 2.0 * Q(1, 1) * Q(2, 1);
    Qfourth(4, 2) = 2.0 * Q(1, 2) * Q(2, 2);
    Qfourth(4, 3) = Q(1, 0) * Q(2, 1) + Q(2, 0) * Q(1, 1);
    Qfourth(4, 4) = Q(1, 1) * Q(2, 2) + Q(2, 1) * Q(1, 2);
    Qfourth(4, 5) = Q(1, 0) * Q(2, 2) + Q(2, 0) * Q(1, 2);

    Qfourth(5, 0) = 2.0 * Q(0, 0) * Q(2, 0);
    Qfourth(5, 1) = 2.0 * Q(0, 1) * Q(2, 1);
    Qfourth(5, 2) = 2.0 * Q(0, 2) * Q(2, 2);
    Qfourth(5, 3) = Q(0, 0) * Q(2, 1) + Q(2, 0) * Q(0, 1);
    Qfourth(5, 4) = Q(0, 1) * Q(2, 2) + Q(2, 1) * Q(0, 2);
    Qfourth(5, 5) = Q(0, 0) * Q(2, 2) + Q(2, 0) * Q(0, 2);
  }
  /*!
   * @brief Inverse rotation of all basis vectors of Tin by the rotation matrix Q
   *
   * @param Q (in) : Rotation matrix
   * @param fourthTensorGlobal (in) : Rotated fourth order tensor
   * @param fourthTensorLocal (out) : Unrotated fourth order tensor
   */
  inline void InverseFourthTensorRotation(const CORE::LINALG::Matrix<3, 3>& Q,
      const CORE::LINALG::Matrix<6, 6>& Tin, CORE::LINALG::Matrix<6, 6>& Tout)
  {
    CORE::LINALG::Matrix<6, 6> Qfourth;
    FourthTensorRotationMatrix(Q, Qfourth);

    CORE::LINALG::Matrix<6, 6> temp(false);
    temp.MultiplyTN(1.0, Qfourth, Tin);
    Tout.MultiplyNN(1.0, temp, Qfourth);
  }

}  // namespace CORE::LINALG::TENSOR


FOUR_C_NAMESPACE_CLOSE

#endif