/*----------------------------------------------------------------------*/
/*! \file
\brief Various service routines related to materials

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_SERVICE_HPP
#define FOUR_C_MAT_SERVICE_HPP


#include "4C_config.hpp"

#include "4C_linalg_fixedsizematrix.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  template <int dim>
  class FourTensor;
}

namespace CORE::MAT::PAR
{
  class Parameter;
}

namespace MAT
{

  /*!
   * @brief Multiply two 2nd order tensors A x B and add the result to a 4th order symmetric
   * material tensor in matrix notation
   *
   * In tensor index notation this method does
   * \f[
   * C_{IJKL} := \text{ScalarThis} \cdot C_{IJKL}
   *          + \frac{1}{2} \cdot \text{ScalarAB} \cdot \left( A_{IJ} \cdot B_{KL} \right)
   * \f]
   *
   * Compute the "elasticity tensor product" A x B of two 2nd order tensors (in matrix notation) and
   * add the result to a 4th order tensor (in Voigt matrix notation!) using the symmetry-conditions
   * inherent to elasticity tensors. This is another version of this function using the fixed size
   * matrix.
   *
   * @param[in,out] C       Material tangent matrix to be modified
   * @param[in] ScalarAB    Scalar to multiply with A x B
   * @param[in] A           Dense matrix (3 x 3) as 2nd order tensor A
   * @param[in] B           Dense matrix (3 x 3) as 2nd order tensor B
   * @param[in] ScalarThis  Scalar to multiply with C before adding A x B
   */
  void ElastSymTensorMultiply(CORE::LINALG::Matrix<6, 6>& C, const double ScalarAB,
      const CORE::LINALG::Matrix<3, 3>& A, const CORE::LINALG::Matrix<3, 3>& B,
      const double ScalarThis);

  /*!
   * @brief Multiply two 2nd order tensors (A x B + B x A) and add the result to a 4th order
   * symmetric material tensor in matrix notation
   *
   * In tensor index notation this method does
   * \f[
   * C_{IJKL} := \text{ScalarThis} \cdot C_{IJKL}
   *          + \text{ScalarAB} \cdot \left( A_{IJ} \cdot B_{KL} + B_{IJ} \cdot A_{KL} \right)
   * \f]
   *
   * Compute the "elasticity tensor product" (A x B + B x A) of two 2nd order tensors (in matrix
   * notation) and add the result to a 4th order tensor (in Voigt matrix notation!) using the
   * symmetry-conditions inherent to elasticity tensors
   *
   * @param[in,out] C       Material tangent matrix to be modified
   * @param[in] ScalarAB    Scalar to multiply with (A x B + B x A)
   * @param[in] A           Dense matrix (3 x 3) as 2nd order tensor A
   * @param[in] B           Dense matrix (3 x 3) as 2nd order tensor B
   * @param[in] ScalarThis  Scalar to multiply with C before adding (A x B + B x A)
   */
  void ElastSymTensorMultiplyAddSym(CORE::LINALG::Matrix<6, 6>& C, const double ScalarAB,
      const CORE::LINALG::Matrix<3, 3>& A, const CORE::LINALG::Matrix<3, 3>& B,
      const double ScalarThis);

  /*!
   * @brief Multiply two 2nd order tensors A o B and add the result to a 4th order symmetric
   * material tensor in matrix notation.
   *
   * In tensor index notation this method does
   * \f[
   * C_{IJKL} := \text{ScalarThis} \cdot C_{IJKL} + \frac{1}{2} \cdot \text{ScalarAB} \cdot \left(
   *             A_{IK} \cdot B_{JL} + A_{IL} \cdot B_{JK} \right) \f]
   *
   * Compute the "material tensor product" A o B (also known as kronecker-tensor-product) of two 2nd
   * order tensors (in matrix notation) and add the result to a 4th order tensor (also in matrix
   * notation) using the symmetry-conditions inherent to material tensors, or tangent matrices,
   * respectively AND the Voigt notation of E,S, and C with the famous factor 2! The implementation
   * uses the fixed size matrix.ss
   *
   * @param[in,out] C       Material tangent matrix to be modified
   * @param[in] ScalarAB    Scalar to multiply with A o B
   * @param[in] A           Dense matrix (3 x 3) as 2nd order tensor A
   * @param[in] B           Dense matrix (3 x 3) as 2nd order tensor B
   * @param[in] ScalarThis  Scalar to multiply with C before adding A o B
   */
  void ElastSymTensor_o_Multiply(CORE::LINALG::Matrix<6, 6>& C, const double ScalarAB,
      const CORE::LINALG::Matrix<3, 3>& A, const CORE::LINALG::Matrix<3, 3>& B,
      const double ScalarThis);

  /*!
   * @brief Add 'Holzapfel product' contribution to constitutive tensor using Voigt notation
   *
   * This function adds the following contribution to the given constitutive
   * matrix \f$\mathbb{C}\f$ (#cmat) based on the inverse of the right Cauchy-Green vector
   * \f$\boldsymbol{C}^{-1}\f$ (#invc):
   *
   * In compact tensor notation this method performs
   * \f[
   *    \mathbb{C} := \mathbb{C} + \text{scalar} \cdot \boldsymbol{C}^{-1} \odot
   *                           \boldsymbol{C}^{-1} \f]
   *
   * which is related to the derivative of the inverse right Cauchy--Green tensor
   * with respect to its non-inverted sibling, i.e.
   * \f[
   *    \frac{\partial \boldsymbol{C}^{-1}}{\partial \boldsymbol{C}}
   *    = - \boldsymbol{C}^{-1} \odot \boldsymbol{C}^{-1}
   * \f]
   * as written by Holzapfel [1], p. 254.
   *
   * In tensor index notation this \f$\odot\f$-product means
   * \f[
   * \mathbb{C}^{ABCD} := \mathbb{C}^{ABCD} + \text{scalar} \left( {C^{-1}}^{AC} {C^{-1}}^{BD} +
   * {C^{-1}}^{AD} {C^{-1}}^{BC} \right) \f]
   *
   * References:
   * [1] G.A. Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
   *
   * @param[in,out] cmat Material tangent matrix(6x6) to be modified, its rows are denoted
   *                     stress-like 6-Voigt, i.e. \f$[ ()_{11}, \; ()_{22}, \; ()_{33}, \; ()_{12},
   *                     \; ()_{23}, \; ()_{31} ]\f$ and its columns strain-like 6-Voigt, i.e. \f$[
   *                     ()_{11}, \; ()_{22}, \; ()_{33}, \; 2()_{12}, \; 2()_{23}, \; 2()_{31} ]\f$
   *
   * @param[in] invC     Inverse right Cauchy-Green vector in stress-like 6-Voigt notation
   *                     \f$[ C^{-1}_{11} \; C^{-1}_{22} \; C^{-1}_{33} \; C^{-1}_{12} \;
   *                     C^{-1}_{23} \; C^{-1}_{31}]\f$
   *
   * @param[in] scalar   Scalar, corresponding to \f$\delta_7\f$ in Holzapfel [1], p.261
   */
  template <typename T>
  void AddtoCmatHolzapfelProduct(CORE::LINALG::Matrix<6, 6, T>& cmat,
      const CORE::LINALG::Matrix<6, 1, T>& invc, const T scalar);

  /*!
   * @brief Put apart the volumetric and isochoric part of the 2nd Piola--Kirchhoff stress and the
   * material constitutive tensor
   *
   * The volumetric 2nd Piola--Kirchhof stress \f$S^{AB}_\text{vol}\f$
   * are obtained by
   * \f[
   *   S^{AB}_\text{vol}
   *   = \frac{1}{3} \left( S^{CD} \, C_{CD} \right) \, {C^{-1}}^{AB}
   * \f]
   *
   * and its corresponding contribution to the constitutive tensor
   * \f$C^{ABCD}_\text{vol}\f$ is
   * \f[
   *   C^{ABCD}_\text{vol}
   *   = \frac{2}{3} {C^{-1}}^{AB} \, S^{CD}_\text{lin}
   *   + \frac{2}{3} \left(S^{EF} C_{EF}\right) \, \left( -\frac{1}{2} \left(
   *     {C^{-1}}^{AC} {C^{-1}}^{BD} + {C^{-1}}^{AD} {C^{-1}}^{BC}
   *   \right) \right)
   * \f]
   *
   * with the 'linearised' 2nd Piola--Kirchhoff stresses
   * \f[
   *   S^{CD}_\text{lin} = S^{CD} + \frac{1}{2} C_{AB} C^{ABCD} .
   * \f]
   *
   * The isochoric 2nd Piola--Kirchhoff stress \f$S^{AB}_\text{iso}\f$
   * and and the isochoric constribution \f$C^{ABCD}_\text{iso}\f$
   * to the constitutive tensor are
   * \f[
   *   S^{AB}_\text{iso} = S^{AB} - S^{AB}_\text{vol}
   * \f]
   *
   * and
   * \f[
   *   C^{ABCD}_\text{iso} = C^{ABCD} - C^{ABCD}_\text{vol}
   * \f].
   *
   * @param[in,out] pk2vol  Volumetric 2nd Piola--Kirchhoff 6-Voigt stress, i.e.
   *                        \f$[ S^{11}\; S^{22}\; S^{33}\; S^{12}\; S^{23}\; S^{31}
   *                        ]_\text{vol}\f$, if != nullptr
   * @param[in] cvol        Volumetric contribution to constitutive tensor, if != nullptr
   * @param[in] pk2iso      Isochoric 2nd Piola--Kirchhoff 6-Voigt stress, i.e.
   *                        \f$[ S^{11}\; S^{22}\; S^{33}\; S^{12}\; S^{23}\; S^{31}
   *                        ]_\text{iso}\f$, if != nullptr
   * @param[in,out] ciso    Isochoric contribution to constitutive tensor, if != nullptr
   * @param[in] gl          Green--Lagrange 6-Voigt strain, i.e. \f$[ E_{11}\; E_{22}\; E_{33}\;
   *                        2E_{12}\; 2E_{23}\; 2E_{31} ]\f$
   * @param[in] pk2         2nd Piola--Kirchhoff 6-Voigt stress, i.e. \f$[ S^{11}\; S^{22}\;
   *                        S^{33}\; S^{12}\; S^{23}\; S^{31} ]\f$
   * @param[in] cmat        Constitutive tensor
   */
  void VolumetrifyAndIsochorify(CORE::LINALG::Matrix<6, 1>* pk2vol,
      CORE::LINALG::Matrix<6, 6>* cvol, CORE::LINALG::Matrix<6, 1>* pk2iso,
      CORE::LINALG::Matrix<6, 6>* ciso, const CORE::LINALG::Matrix<6, 1>& gl,
      const CORE::LINALG::Matrix<6, 1>& pk2, const CORE::LINALG::Matrix<6, 6>& cmat);

  /*!
   * @brief Add the derivative of the square of a tensor to a 4th order symmetric material tensor in
   * matrix notation
   *
   * In tensor index notation this method does
   * \f[
   * C_{IJKL} := \text{ScalarThis} \cdot C_{IJKL}
   *           + \text{ScalarDX2} \cdot \frac{d (X^2)_{IJ}}{d X_{KL}}
   * \f]
   *
   * Compute the derivative of the square of a symmetric 2nd order tensor w.r.t. to the tensor add
   * the result to a 4th order tensor (in Voigt matrix notation!) using the symmetry-conditions
   * inherent to elasticity tensors.
   *
   * @param[in,out] C       Material tangent matrix to be modified
   * @param[in] ScalarDX2   Scalar to multiply with dX^2/dX
   * @param[in] X           Dense matrix (3 x 3)
   * @param[in] ScalarThis  Scalar to multiply with C before adding dX^2/dX
   */
  void AddToCmatDerivTensorSquare(CORE::LINALG::Matrix<6, 6>& C, double ScalarDX2,
      CORE::LINALG::Matrix<3, 3> X, double ScalarThis);

  /*!
   * @brief Add symmetric Holzapfel product to a 4th order tensor in matrix notation
   *
   * In tensor index notation this does:
   * \f[
   *    X_{abcd} += \text{fac} \cdot \left( A_{ca} \cdot B_{db} + A_{da} \cdot B_{cb} + A_{db}
   * \cdot B_{ca} + A_{cb} \cdot B_{da} \right)
   * \f]
   *
   * The result has minor symmetries, but no major symmetry, i.e. symmetric w.r.t. A and B
   *
   * @param[in,out] X 4th order tensor \f[X_{abcd}\f]
   * @param[in] A     2nd order tensor \f[A_{ef}\f]
   * @param[in] B     2nd order tensor \f[B_{gh}\f]
   * @param[in] fac   Scalar factor
   */
  void AddSymmetricHolzapfelProduct(CORE::LINALG::Matrix<6, 6>& X,
      const CORE::LINALG::Matrix<3, 3>& A, const CORE::LINALG::Matrix<3, 3>& B, const double fac);

  /*!
   * @brief Add right non-symmetric Holzapfel product to a 4th order tensor in matrix notation
   *
   * In tensor index notation this does:
   * \f[
   *    \text{out}_{IJKL} += \text{fac} \cdot \left( A_{IK} \cdot B_{JL}
   *                       + A_{JK} \cdot B_{IL} \right)
   * \f]
   *
   * Note that the result is symmetric within the first two indices. Thus the first two indices
   * are stored in stress-like Voigt notation. The 2nd index pair has no symmetries and is
   * therefore stored as 9-vector. Be aware that this is NOT symmetric w.r.t. A and B !
   *
   * @param[in,out] out     4th order tensor \f[\text{out}_{ABCD}\f]
   * @param[in] A           2nd order tensor \f[A_{EF}\f]
   * @param[in] B           2nd order tensor \f[B_{GH}\f]
   * @param[in] fac         Scalar factor
   */
  template <typename T>
  void AddRightNonSymmetricHolzapfelProduct(CORE::LINALG::Matrix<6, 9, T>& out,
      CORE::LINALG::Matrix<3, 3, T> const& A, CORE::LINALG::Matrix<3, 3, T> const& B, T const fac);

  /*!
   * @brief Add right non-symmetric Holzapfel product, where the symmetric part is stored
   * strain-like to a 4th order tensor in matrix notation
   *
   * In tensor index notation this does:
   * \f[
   *    \text{out}_{IJKL} += \text{fac} \cdot \left( A_{IK} \cdot B_{JL}
   *                      + A_{JK} \cdot B_{IL} \right)
   * \f]
   *
   * Note that the result is symmetric within the first two indices. The first two indices are
   * stored in strain-like Voigt notation. The 2nd index pair has no symmetries and is therefore
   * stored as 9-vector. Be aware that this is NOT symmetric w.r.t. A and B !
   *
   * @param[in,out] out     4th order tensor \f[\text{out}_{ABCD}\f]
   * @param[in] A           2nd order tensor \f[A_{EF}\f]
   * @param[in] B           2nd order tensor \f[B_{GH}\f]
   * @param[in] fac         Scalar factor
   */
  template <typename T>
  void AddRightNonSymmetricHolzapfelProductStrainLike(CORE::LINALG::Matrix<6, 9, T>& out,
      CORE::LINALG::Matrix<3, 3, T> const& A, CORE::LINALG::Matrix<3, 3, T> const& B, T const fac);

  /*!
   * @brief Add left non-symmetric Holzapfel product to a 4th order tensor in matrix notation
   *
   * In tensor index notation this does:
   * \f[
   *    \text{out}_{IJKL} += \text{fac} \cdot \left( A_{IK} \cdot B_{JL}
   *                       + A_{IL} \cdot B_{JK} \right)
   * \f]
   *
   * Note that the result is symmetric within the second two indices. Thus the second two indices
   * are stored in stress-like Voigt notation. The 1st index pair has no symmetries and is
   * therefore stored as 9-vector. Be aware that this is NOT symmetric w.r.t. A and B !
   *
   * @param[in,out] out     4th order tensor \f[\text{out}_{ABCD}\f]
   * @param[in] A           2nd order tensor \f[A_{EF}\f]
   * @param[in] B           2nd order tensor \f[B_{GH}\f]
   * @param[in] fac         Scalar factor
   */
  void AddLeftNonSymmetricHolzapfelProduct(CORE::LINALG::Matrix<9, 6>& out,
      CORE::LINALG::Matrix<3, 3> const& A, CORE::LINALG::Matrix<3, 3> const& B, double const fac);

  /*!
   * @brief Add non-symmetric product to a 4th order tensor in matrix notation
   *
   * In tensor index notation this does:
   * \f[
   *    X_{ABCD} = \text{fac} \cdot \left( A_{AC} \cdot B_{DB} \right)
   * \f]
   *
   * There are no symmetries and the matrix is therefore stored in a 9x9 matrix. A and B are not
   * required to be symmetric!!!
   *
   * @param[in] fac         Scalar factor
   * @param[in] A           2nd order tensor \f[A_{AB}\f]
   * @param[in] B           2nd order tensor \f[B_{CD}\f]
   * @param[in,out] out     4th order tensor \f[\text{out}_{EFGH}\f]
   */
  void AddNonSymmetricProduct(double const& fac, CORE::LINALG::Matrix<3, 3> const& A,
      CORE::LINALG::Matrix<3, 3> const& B, CORE::LINALG::Matrix<9, 9>& out);

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
  void AddDerivInvABInvBProduct(double const& fac, const CORE::LINALG::Matrix<6, 1>& invA,
      const CORE::LINALG::Matrix<6, 1>& invABinvA, CORE::LINALG::Matrix<6, 6>& out);

  /*!
   * @brief Evaluates the principal invariants of any symmetric tensor
   *
   * @param[out] prinv  Principal invariants of the tensor
   * @param[in]  tens   Tensor
   */
  void InvariantsPrincipal(
      CORE::LINALG::Matrix<3, 1>& prinv, const CORE::LINALG::Matrix<3, 3>& tens);

  /*!
   * @brief Converts the principal invariants to the modified principal invariants.
   *
   * @param[out] modinv Modified principal invariants of the Right Cauchy-Green strain tensor
   * @param[in] prinv   Principal invariants of the Right Cauchy-Green strain tensor
   */
  void invariants_modified(
      CORE::LINALG::Matrix<3, 1>& modinv, const CORE::LINALG::Matrix<3, 1>& prinv);

  /*!
   * @brief Evaluates principal stretches \f$\lambda_\alpha\f$ and material principal directions
   *
   * @param[out] prstr  Principal stretches
   * @param[out] prdir  Principal directions
   * @param[in] rcg     Right Cauchy-Green strain tensor in strain-like Voigt notation
   */
  void StretchesPrincipal(CORE::LINALG::Matrix<3, 1>& prstr, CORE::LINALG::Matrix<3, 3>& prdir,
      const CORE::LINALG::Matrix<6, 1>& rcg);

  /*!
   * @brief Evaluates the modified principal stretches \f$\bar{\lambda}_\alpha\f$
   *
   * @param[out] modstr Modified principal stretches
   * @param[in] prstr   Principal stretches
   */
  void StretchesModified(
      CORE::LINALG::Matrix<3, 1>& modstr, const CORE::LINALG::Matrix<3, 1>& prstr);

  inline void InvariantsPrincipal(
      CORE::LINALG::Matrix<3, 1>& prinv, const CORE::LINALG::Matrix<3, 3>& tens)
  {
    // 1st invariant, trace tens
    prinv(0) = tens(0, 0) + tens(1, 1) + tens(2, 2);

    // 2nd invariant, 0.5( (trace(tens))^2 - trace(tens^2))
    prinv(1) = tens(0, 0) * tens(1, 1) + tens(1, 1) * tens(2, 2) + tens(0, 0) * tens(2, 2) -
               tens(0, 1) * tens(1, 0) - tens(1, 2) * tens(2, 1) - tens(0, 2) * tens(2, 0);

    // 3rd invariant, determinant tens
    prinv(2) = tens.Determinant();
  }

  template <class T>
  T* CreateMaterialParameterInstance(CORE::MAT::PAR::Parameter* curmat)
  {
    auto* params = dynamic_cast<T*>(curmat);
    return params;
  }

  /// Set every tensor value to zero
  template <int dim>
  void ClearFourTensor(CORE::LINALG::FourTensor<dim>& fourTensor);

  /*!
   * @brief Multiply 4th order tensor with 2nd order matrix
   *
   * @param[in,out] fourTensorResult   Result of 4th order tensor - matrix multiplication
   *                                   \f$C^{ijkl} += A^{ijkm} B_m^l\f$
   * @param[in] fourTensor             4th order tensor \f$A^{ijkm}\f$
   * @param[in] matrix                 2-tensor \f$B_m^l\f$
   * @param[in] clearResultTensor      if true: \f$C^{ijkl}\f$ is cleared when method is called
   */
  template <int dim>
  void MultiplyFourTensorMatrix(CORE::LINALG::FourTensor<dim>& fourTensorResult,
      const CORE::LINALG::FourTensor<dim>& fourTensor, const CORE::LINALG::Matrix<dim, dim>& matrix,
      const bool clearResultTensor = true);

  /*!
   * @brief Multiply 2nd order matrix with 4th order tensor
   *
   * @param[in,out] fourTensorResult   Result of matrix - 4th order tensor multiplication
   *                                   \f$C^{ijkl} += B^i_m A^{mjkl}\f$
   * @param[in] matrix                 2-tensor \f$B^i_m\f$
   * @param[in] fourTensor             4th order tensor \f$A^{mjkl}\f$
   * @param[in] clearResultTensor      if true: \f$C^{ijkl}\f$ is cleared when method is called
   */
  template <int dim>
  void MultiplyMatrixFourTensor(CORE::LINALG::FourTensor<dim>& fourTensorResult,
      const CORE::LINALG::Matrix<dim, dim>& matrix, const CORE::LINALG::FourTensor<dim>& fourTensor,
      const bool clearResultTensor = true);

  /*!
   * @brief Multiply 2nd order matrix with 4th order tensor by the second index
   *
   * @param[in,out] fourTensorResult   Result of matrix - 4th order tensor multiplication
   *                                   \f$C^{ijkl} += B_m^j A^{imkl}\f$
   * @param[in] matrix                 2nd order tensor \f$B_m^j\f$
   * @param[in] fourTensor             4th order tensor \f$A^{imkl}\f$
   * @param[in] clearResultTensor      if true: C_ijkl is cleared when method is called
   */
  template <int dim>
  void MultiplyMatrixFourTensorBySecondIndex(CORE::LINALG::FourTensor<dim>& fourTensorResult,
      const CORE::LINALG::Matrix<dim, dim>& matrix, const CORE::LINALG::FourTensor<dim>& fourTensor,
      const bool clearResultTensor = true);

  /*!
   * @brief Multiply 4th order tensor with 4th order tensor by the first two indices
   *
   * @param[in,out] fourTensorResult   Result of 4th order tensor - 4th order tensor multiplication
   *                                   \f$C^{ijkl} += A^{ij}_{ab} B^{abkl}\f$
   * @param[in] fourTensor1            4th order tensor \f$A^{ij}_{ab}\f$
   * @param[in] fourTensor2            4th order tensor \f$B^{abkl}\f$
   * @param[in] clearResultTensor      if true: C_ijkl is cleared when method is called
   */
  template <int dim>
  void MultiplyFourTensorFourTensor(CORE::LINALG::FourTensor<dim>& fourTensorResult,
      const CORE::LINALG::FourTensor<dim>& fourTensor1,
      const CORE::LINALG::FourTensor<dim>& fourTensor2, const bool clearResultTensor = true);

  /*!
   * @brief Pull back of a symmetric elastic 4th order tensor (in matrix/voigt notation) via the 2nd
   * order deformation gradient (also in matrix notation)
   */
  template <int dim>
  CORE::LINALG::Matrix<6, 6> PullBackFourTensor(
      const CORE::LINALG::Matrix<dim, dim>& defgr, const CORE::LINALG::Matrix<6, 6>& cMatVoigt);

  /*!
   * @brief Pull back the ijkl-th entry of a symmetric elastic 4th order tensor (in matrix/voigt
   * notation) via the 2nd order deformation gradient (also in matrix notation)
   */
  template <int dim>
  double PullBackFourTensorijkl(const CORE::LINALG::Matrix<dim, dim>& defgr,
      const CORE::LINALG::FourTensor<dim>& fourTensor, const int i, const int j, const int k,
      const int l);

  /*!
   * @brief Setup 4th order tensor from 6x6 Voigt notation
   *
   * @note Setup 4th order tensor from 6x6 Voigt matrix (which has to be the representative of a 4
   * tensor with at least minor symmetries)
   *
   * @param[out] fourTensor   4th order tensor that is set up based on matrixVoigt
   * @param[in]  matrixVoigt  4th order tensor in Voigt notation with at least minor symmetries
   */
  template <int dim>
  void SetupFourTensor(
      CORE::LINALG::FourTensor<dim>& fourTensor, const CORE::LINALG::Matrix<6, 6>& matrixVoigt);

  /*!
   * @brief Setup 6x6 Voigt matrix from 4th order tensor with minor symmetries
   *
   * @param[out] matrixVoigt  6x6 Voigt matrix that is set up based on fourTensor
   * @param[in]  fourTensor   4th order tensor with minor symmetries
   */
  template <int dim>
  void Setup6x6VoigtMatrix(
      CORE::LINALG::Matrix<6, 6>& matrixVoigt, const CORE::LINALG::FourTensor<dim>& fourTensor);

  /*!
   * @brief Transpose 4th order tensor w.r.t. to basis vectors 1 and 2
   */
  template <int dim>
  void TransposeFourTensor12(CORE::LINALG::FourTensor<dim>& resultTensor,
      const CORE::LINALG::FourTensor<dim>& inputTensor);

  /*!
   * @brief Adds the dyadic product of two 2nd order tensors (in matrix notation) to a 4th order
   * tensor
   *
   * @param[in,out] fourTensorResult    4th order tensor \f$C^{ijkl} = A^{ij} B^{kl}\f$
   * @param[in]     matrixA             2nd order tensor \f$A^{ij}\f$
   * @param[in]     matrixB             2nd order tensor \f$B^{kl}\f$
   */
  void AddDyadicProductMatrixMatrix(CORE::LINALG::FourTensor<3>& fourTensorResult,
      const CORE::LINALG::Matrix<3, 3>& matrixA, const CORE::LINALG::Matrix<3, 3>& matrixB);

  /*!
   * @brief Adds the contraction of a 2nd order tensor (in matrix notation) and a 4th order tensor
   * to a 2nd order tensor (in matrix notation)
   *
   * @param[in,out] matrixResult    2nd order tensor \f$C^{kl} = A_{ij} B^{ijkl}\f$
   * @param[in]     matrix          2nd order tensor \f$A_{ij}\f$
   * @param[in]     fourTensor      4th order tensor \f$B^{ijkl}\f$
   */
  void AddContractionMatrixFourTensor(CORE::LINALG::Matrix<3, 3>& matrixResult,
      const CORE::LINALG::Matrix<3, 3>& matrix, const CORE::LINALG::FourTensor<3>& fourTensor);

  /*!
   * @brief Returns the double contraction of two 2nd order tensors (in matrix notation)
   *
   * @param[out]    scalarContraction   0th order tensor \f$s = A_{ij} B^{ij}\f$
   * @param[in]     matrixA             2nd order tensor \f$A_{ij}\f$
   * @param[in]     matrixB             2nd order tensor \f$B^{ij}\f$
   */
  double ContractMatrixMatrix(
      const CORE::LINALG::Matrix<3, 3>& matrixA, const CORE::LINALG::Matrix<3, 3>& matrixB);
}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
