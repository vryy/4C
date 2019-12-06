/*----------------------------------------------------------------------*/
/*! \file
\brief Various service routines related to materials

\level 1

\maintainer Amadeus Gebauer
*/
/*----------------------------------------------------------------------*/

#ifndef MATERIAL_SERVICE_H
#define MATERIAL_SERVICE_H



#include <Teuchos_RCP.hpp>
#include <Epetra_SerialDenseVector.h>
#include "../linalg/linalg_fixedsizematrix.H"


/// MAT: materials
namespace MAT
{
  /*! \brief Multiply two 2nd order tensors A x B and add the result to a
   *         4th order symmetric material tensor in matrix notation
   *
   * In tensor index notation this method does
   * \f[
   * C_{IJKL} := ScalarThis \cdot C_{IJKL}
   *           + 1/2 \cdot ScalarAB \cdot ( A_{IJ} \cdot B_{KL} )
   * \f]
   *
   * compute the "elasticity tensor product" A x B of
   * two 2nd order tensors (in matrix notation) and add the result to
   * a 4th order tensor (in Voigt matrix notation!) using the
   * symmetry-conditions inherent to elasticity tensors
   * The implementation is based on the Epetra-Method Matrix.Multiply.
   *
   * \param C (in/out)    : Material tangent matrix to be modified
   * \param ScalarAB (in) : Scalar to multiply with A x B
   * \param A (in)        : Dense matrix (3 by 3) as 2nd order tensor A
   * \param B (in)        : Dense matrix (3 by 3) as 2nd order tensor B
   * \param ScalarThis (in): Scalar to multiply with C before adding A x B
   */
  void ElastSymTensorMultiply(Epetra_SerialDenseMatrix& C, const double ScalarAB,
      const Epetra_SerialDenseMatrix& A, const Epetra_SerialDenseMatrix& B,
      const double ScalarThis);

  /*! \brief Multiply two 2nd order tensors A x B and add the result to a
   *         4th order symmetric material tensor in matrix notation
   *
   * In tensor index notation this method does
   * \f[
   * C_{IJKL} := ScalarThis \cdot C_{IJKL}
   *          + ScalarAB \cdot ( A_{IJ} \cdot B_{KL} )
   * \f]
   *
   * compute the "elasticity tensor product" A x B of
   * two 2nd order tensors (in matrix notation) and add the result to
   * a 4th order tensor (in Voigt matrix notation!) using the
   * symmetry-conditions inherent to elasticity tensors.
   * This is another version of this function using the fixed size matrix.
   *
   * \param C (in/out)    : Material tangent matrix to be modified
   * \param ScalarAB (in) : Scalar to multiply with A x B
   * \param A (in)        : Dense matrix (3 by 3) as 2nd order tensor A
   * \param B (in)        : Dense matrix (3 by 3) as 2nd order tensor B
   * \param ScalarThis (in): Scalar to multiply with C before adding A x B
   */
  void ElastSymTensorMultiply(LINALG::Matrix<6, 6>& C, const double ScalarAB,
      const LINALG::Matrix<3, 3>& A, const LINALG::Matrix<3, 3>& B, const double ScalarThis);

  /*! \brief Multiply two 2nd order tensors (A x B + B x A) and add the result to a
   *         4th order symmetric material tensor in matrix notation
   *
   * In tensor index notation this method does
   * \f[
   * C_{IJKL} := ScalarThis \cdot C_{IJKL}
   *          + ScalarAB \cdot ( A_{IJ} \cdot B_{KL} + B_{IJ} \cdot A_{KL} )
   * \f]
   *
   * compute the "elasticity tensor product" (A x B + B x A) of
   * two 2nd order tensors (in matrix notation) and add the result to
   * a 4th order tensor (in Voigt matrix notation!) using the
   * symmetry-conditions inherent to elasticity tensors
   *
   * \param C (in/out)     : Material tangent matrix to be modified
   * \param ScalarAB (in)  : Scalar to multiply with (A x B + B x A)
   * \param A (in)         : Dense matrix (3 by 3) as 2nd order tensor A
   * \param B (in)         : Dense matrix (3 by 3) as 2nd order tensor B
   * \param ScalarThis (in): Scalar to multiply with C before adding (A x B + B x A)
   */
  void ElastSymTensorMultiplyAddSym(Epetra_SerialDenseMatrix& C, const double ScalarAB,
      const Epetra_SerialDenseMatrix& A, const Epetra_SerialDenseMatrix& B,
      const double ScalarThis);

  /*! \brief Multiply two 2nd order tensors (A x B + B x A) and add the result to a
   *         4th order symmetric material tensor in matrix notation
   *
   * In tensor index notation this method does
   * \f[
   * C_{IJKL} := ScalarThis \cdot C_{IJKL}
   *          + ScalarAB \cdot ( A_{IJ} \cdot B_{KL} + B_{IJ} \cdot A_{KL} )
   * \f]
   *
   * compute the "elasticity tensor product" (A x B + B x A) of
   * two 2nd order tensors (in matrix notation) and add the result to
   * a 4th order tensor (in Voigt matrix notation!) using the
   * symmetry-conditions inherent to elasticity tensors
   * This is another version of this function using the fixed size matrix.
   *
   * \param C (in/out)     : Material tangent matrix to be modified
   * \param ScalarAB (in)  : Scalar to multiply with (A x B + B x A)
   * \param A (in)         : Dense matrix (3 by 3) as 2nd order tensor A
   * \param B (in)         : Dense matrix (3 by 3) as 2nd order tensor B
   * \param ScalarThis (in): Scalar to multiply with C before adding (A x B + B x A)
   */
  void ElastSymTensorMultiplyAddSym(LINALG::Matrix<6, 6>& C, const double ScalarAB,
      const LINALG::Matrix<3, 3>& A, const LINALG::Matrix<3, 3>& B, const double ScalarThis);

  /*! \brief Multiply two 2nd order tensors A o B and add the result to a
   *         4th order symmetric material tensor in matrix notation
   *
   * In tensor index notation this method does
   * \f[
   * C_{IJKL} := ScalarThis \cdot C_{IJKL}
   *          + 1/2 \cdot ScalarAB \cdot ( A_{IK} \cdot B_{JL} + A_{IL} \cdot B_{JK} )
   * \f]
   *
   * \param C (in/out)    : Material tangent matrix to be modified
   * \param ScalarAB (in) : Scalar to multiply with A o B
   * \param A (in)        : Dense matrix (3 by 3) as 2nd order tensor A
   * \param B (in)        : Dense matrix (3 by 3) as 2nd order tensor B
   * \param ScalarThis (in): Scalar to multiply with C before adding A o B
   */
  void ElastSymTensor_o_Multiply(Epetra_SerialDenseMatrix& C, const double ScalarAB,
      const Epetra_SerialDenseMatrix& A, const Epetra_SerialDenseMatrix& B,
      const double ScalarThis);

  /*! \brief Multiply two 2nd order tensors A o B and add the result to a
   *         4th order symmetric material tensor in matrix notation
   *         This is another version of this function using the fixed size matrix.
   *
   * In tensor index notation this method does
   * \f[
   * C_{IJKL} := ScalarThis \cdot C_{IJKL}
   *          + 1/2 \cdot ScalarAB \cdot ( A_{IK} \cdot B_{JL} + A_{IL} \cdot B_{JK} )
   * \f]
   *
   * \param C (in/out)    : Material tangent matrix to be modified
   * \param ScalarAB (in) : Scalar to multiply with A o B
   * \param A (in)        : Dense matrix (3 by 3) as 2nd order tensor A
   * \param B (in)        : Dense matrix (3 by 3) as 2nd order tensor B
   * \param ScalarThis (in): Scalar to multiply with C before adding A o B
   */
  void ElastSymTensor_o_Multiply(LINALG::Matrix<6, 6>& C, const double ScalarAB,
      const LINALG::Matrix<3, 3>& A, const LINALG::Matrix<3, 3>& B, const double ScalarThis);


  /*! \brief Add 'Holzapfel product' contribution to constitutive tensor
   *         using Voigt notation
   *
   * This function adds the following contribution to the given constitutive
   * matrix cmat(6,6) based on the inverse of the right Cauchy-Green vector
   * invc(6):
   *
   * \verbatim
     scalar * ( Cinv boeppel Cinv )
   * \endverbatim
   *
   * For that purpose we need the derivative
   *
   * \verbatim
     \partial tensor(C)^-1
     -----------------------
     \partial tensor(C)
   * \endverbatim
   *
   * which yields the following product
   *
   * \verbatim
     - ( Cinv boeppel Cinv )_{abcd} = 1/2 * ( Cinv_{ac} Cinv_{bd} + Cinv_{ad} Cinv_{bc} )
   * \endverbatim
   *
   * For more details see Holzapfel p. 254
   *
   * \param cmat (in/out): Material tangent matrix(6x6) to be modified
   * \param invc   (in)  : Inverse right Cauchy-Green vector(6)
   * \param scalar (in)  : scalar corresponding to delta7
   */
  void AddtoCmatHolzapfelProduct(
      Epetra_SerialDenseMatrix& cmat, const Epetra_SerialDenseVector& invc, const double scalar);

  /*! \brief Add 'Holzapfel product' contribution to constitutive tensor
   *         using Voigt notation
   *
   * This function adds the following contribution to the given constitutive
   * matrix #cmat based on the inverse of the right Cauchy-Green vector
   * #invc.
   *
   * In compact tensor notation this method performs
   * \f[
   * \boldsymbol{C}_\text{m} := \boldsymbol{C}_\text{m}
   *                + scalar \cdot \boldsymbol{C}^{-1} \odot \boldsymbol{C}^{-1}
   * \f]
   * which is related to the derivative of the inverse right Cauchy--Green tensor
   * with respect to its non-inverted sibling, i.e.
   * \f[
   * \frac{\partial \boldsymbol{C}^{-1}}{\partial \boldsymbol{C}}
   * = - \boldsymbol{C}^{-1} \odot \boldsymbol{C}^{-1}
   * \f]
   * as written by Holzapfel [1], p. 254.
   *
   * In tensor index notation this \f$\odot\f$-product means
   * \f[
   * C^{ABCD} := C^{ABCD}
   *          + scalar ( (C^{-1})^{AC} (C^{-1})^{BD}
   *                     + (C^{-1})^{AD} (C^{-1})^{BC} )
   * \f]
   *
   * References:
   * [1] G.A. Holzapfel, "Nonlinear solid mechanics", Wiley, 2000.
   *
   * \param cmat (in/out): Material tangent matrix(6x6) to be modified,
   *                       its rows are denoted stress-like 6-Voigt, i.e.
   *                       \f$[ ()_{11}, \; ()_{22}, \; ()_{33},
   *                       \; ()_{12}, \; ()_{23}, \; ()_{31} ]\f$
   *                       and its columns strain-like 6-Voigt, i.e.
   *                       \f$[ ()_{11}, \; ()_{22}, \; ()_{33},
   *                       \; 2*()_{12}, \; 2*()_{23}, \; 2*()_{31} ]\f$
   *
   * \param invc   (in)  : inverse right Cauchy-Green vector
   *                       in stress-like 6-Voigt notation
   *                       \f$[ C^{-1}_{11} \; C^{-1}_{22} \; C^{-1}_{33}
   *                       \; C^{-1}_{12} \; C^{-1}_{23} \; C^{-1}_{31}]\f$
   * \param scalar (in)  : scalar, corresponding to \f$\delta_7\f$ in Holzapfel [1], p.261
   */
  template <typename T>
  void AddtoCmatHolzapfelProduct(
      LINALG::Matrix<6, 6, T>& cmat, const LINALG::Matrix<6, 1, T>& invc, const T scalar);


  /*! \brief Put apart the volumetric and isochoric part
   *         of the 2nd Piola--Kirchhoff stress and the
   *         material constitutive tensor
   *
   * The volumetric 2nd Piola--Kirchhof stress \f$S^{AB}_\text{vol}\f$
   * are obtained by
   * \f[
   *   S^{AB}_\text{vol}
   *   = \frac{1}{3} \big( S^{CD} \, C_{CD} \big) \, (C^{-1})^{AB}
   * \f]
   * and its corresponding contribution to the constitutive tensor
   * \f$C^{ABCD}_\text{vol}\f$ is
   * \f[
   *   C^{ABCD}_\text{vol}
   *   = \frac{2}{3} (C^{-1})^{AB} \, S^{CD}_\text{lin}
   *   + \frac{2}{3} \big(S^{EF} C_{EF}\big) \, \Big( -\frac{1}{2} \big(
   *     (C^{-1})^{AC} (C^{-1})^{BD} + (C^{-1})^{AD} (C^{-1})^{BC}
   *   \big) \Big)
   * \f]
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
   * and
   * \f[
   *   C^{ABCD}_\text{iso} = C^{ABCD} - C^{ABCD}_\text{vol}
   * \f]
   *
   * \author bborn
   * \date 05/09
   */
  void VolumetrifyAndIsochorify(
      LINALG::Matrix<6, 1>* pk2vol,  ///< volumetric 2nd Piola--Kirchhoff 6-Voigt stress, i.e.
                                     ///< \f$[ S^{11}\; S^{22}\; S^{33}\; S^{12}\; S^{23}\; S^{31}
                                     ///< ]_\text{vol}\f$, if != NULL
      LINALG::Matrix<6, 6>* cvol,    ///< volumetric contribution to constitutive tensor,
                                     ///< if != NULL
      LINALG::Matrix<6, 1>* pk2iso,  ///< isochoric 2nd Piola--Kirchhoff 6-Voigt stress, i.e.
                                     ///< \f$[ S^{11}\; S^{22}\; S^{33}\; S^{12}\; S^{23}\; S^{31}
                                     ///< ]_\text{iso}\f$, if != NULL
      LINALG::Matrix<6, 6>* ciso,    ///< isochoric contribution to constitutive tensor,
                                     ///< if != NULL
      const LINALG::Matrix<6, 1>&
          gl,  ///< Green--Lagrange 6-Voigt strain, i.e.
               ///< \f$[ E_{11}\; E_{22}\; E_{33}\; 2E_{12}\; 2E_{23}\; 2E_{31} ]\f$
      const LINALG::Matrix<6, 1>&
          pk2,  ///< 2nd Piola--Kirchhoff 6-Voigt stress, i.e.
                ///< \f$[ S^{11}\; S^{22}\; S^{33}\; S^{12}\; S^{23}\; S^{31} ]\f$
      const LINALG::Matrix<6, 6>& cmat  ///< constitutive tensor
  );

  /*! \brief Add the derivative of the square of a tensor to a
   *         4th order symmetric material tensor in matrix notation
   *
   * In tensor index notation this method does
   * \f[
   * C_{IJKL} := ScalarThis \cdot C_{IJKL}
   *          + ScalarDX2 \cdot \frac{d (X^2)_{IJ}}{d X_{KL}}
   * \f]
   *
   * compute the derivative of the square of a symmetric 2nd order
   * tensor w.r.t. to the tensor add the result to
   * a 4th order tensor (in Voigt matrix notation!) using the
   * symmetry-conditions inherent to elasticity tensors.
   *
   * \param C (in/out)    : Material tangent matrix to be modified
   * \param ScalarDX2 (in): Scalar to multiply with dX^2/dX
   * \param X (in)        : Dense matrix (3 by 3)
   * \param ScalarThis (in): Scalar to multiply with C before adding dX^2/dX
   */
  void AddToCmatDerivTensorSquare(
      LINALG::Matrix<6, 6>& C, double ScalarDX2, LINALG::Matrix<3, 3> X, double ScalarThis);

  /**
   * \brief Add symmetric Holzapfel product to a 4th order tensor in matrix notation
   *
   * In tensor index notation this does:
   * \f[ X_{abcd} += fac*( A_{ca}*B_{db} + A_{da}*B_{cb} + A_{db}*B_{ca} + A_{cb}*B_{da} ) \f]
   * The result has minor symmetries, but no major symmetry, i.e. symmetric w.r.t. A and B
   *
   * \param X[in,out] 4th order tensor \f[X_{abcd}\f]
   * \param A[in]     2nd order tensor \f[A_{ef}\f]
   * \param B[in]     2nd order tensor \f[B_{gh}\f]
   * \param fac[in]   scalar factor
   */
  void AddSymmetricHolzapfelProduct(LINALG::Matrix<6, 6>& X, const LINALG::Matrix<3, 3>& A,
      const LINALG::Matrix<3, 3>& B, const double fac);

  /**
   * \brief Add right non-symmetric Holzapfel product to a 4th order tensor in matrix notation
   *
   * In tensor index notation this does:
   * \f[out_{IJKL} += fac*( A_{IK} * B_{JL} + A_{JK} * B_{IL} ) \f]
   * Note that the result is symmetric within the first two indices. Thus the first two indices are
   * stored in stress-like Voigt notation. The 2nd index pair has no symmetries and is therefore
   * stored as 9-vector.
   * Be aware that this is NOT symmetric w.r.t. A and B !
   *
   * \param out[in,out] 4th order tensor \f[out_{ABCD}\f]
   * \param A[in]       2nd order tensor \f[A_{EF}\f]
   * \param B[in]       2nd order tensor \f[B_{GH}\f]
   * \param fac[in]     scalar factor
   */
  template <typename T>
  void AddRightNonSymmetricHolzapfelProduct(LINALG::Matrix<6, 9, T>& out,
      LINALG::Matrix<3, 3, T> const& A, LINALG::Matrix<3, 3, T> const& B, T const fac);

  /**
   * \brief Add right non-symmetric Holzapfel product, where the symmetric part is stored
   * strain-like to a 4th order tensor in matrix notation
   *
   * In tensor index notation this does:
   * \f[ out_{IJKL} += fac*( A_{IK} * B_{JL} + A_{JK} * B_{IL} ) \f]
   * Note that the result is symmetric within the first two indices. The first two indices are
   * stored in strain-like Voigt notation. The 2nd index pair has no symmetries and is therefore
   * stored as 9-vector.
   * Be aware that this is NOT symmetric w.r.t. A and B !
   *
   * \param out[in,out] 4th order tensor \f[out_{ABCD}\f]
   * \param A[in]       2nd order tensor \f[A_{EF}\f]
   * \param B[in]       2nd order tensor \f[B_{GH}\f]
   * \param fac[in]     scalar factor
   */
  template <typename T>
  void AddRightNonSymmetricHolzapfelProductStrainLike(LINALG::Matrix<6, 9, T>& out,
      LINALG::Matrix<3, 3, T> const& A, LINALG::Matrix<3, 3, T> const& B, T const fac);

  /**
   * \brief Add left non-symmetric Holzapfel product to a 4th order tensor in matrix notation
   *
   * In tensor index notation this does:
   * \f[ out_{IJKL} += fac*( A_{IK} * B_{JL} + A_{IL} * B_{JK} ) \f]
   * Note that the result is symmetric within the second two indices. Thus the second two indices
   * are stored in stress-like Voigt notation. The 1st index pair has no symmetries and is therefore
   * stored as 9-vector.
   * Be aware that this is NOT symmetric w.r.t. A and B !
   *
   * \param out[in,out] 4th order tensor \f[out_{ABCD}\f]
   * \param A[in]       2nd order tensor \f[A_{EF}\f]
   * \param B[in]       2nd order tensor \f[B_{GH}\f]
   * \param fac[in]     scalar factor
   */
  void AddLeftNonSymmetricHolzapfelProduct(LINALG::Matrix<9, 6>& out, LINALG::Matrix<3, 3> const& A,
      LINALG::Matrix<3, 3> const& B, double const fac);

  /**
   * \brief Add non-symmetric product to a 4th order tensor in matrix notation
   *
   * In tensor index notation this does:
   * \f[ X_{ABCD} = fac*( A_{AC} * B_{DB} ) \f]
   * There are no symmetries and the matrix is therefore stored in a 9x9 matrix. A and B are not
   * required to be symmetric!!!
   *
   * \param fac[in]     scalar factor
   * \param A[in]       2nd order tensor \f[A_{AB}\f]
   * \param B[in]       2nd order tensor \f[B_{CD}\f]
   * \param out[in,out] 4th order tensor \f[out_{EFGH}\f]
   */
  void AddNonSymmetricProduct(double const& fac,  ///< Scalar factor
      LINALG::Matrix<3, 3> const& A,              ///< Matrix in the above mentioned equation
      LINALG::Matrix<3, 3> const& B,              ///< Matrix in the above mentioned equation
      LINALG::Matrix<9, 9>& out);                 ///< Output

  /*!
   * \brief Evaluates the principal invariants of any symmetric tensor
   *
   * @param prinv (out) : Principal invariants of the tensor
   * @param tens (in) : Tensor
   */
  void InvariantsPrincipal(LINALG::Matrix<3, 1>& prinv,  ///< principal invariants
      const LINALG::Matrix<3, 3>& tens                   ///< symmetric Cartesian 2-tensor
  );

  /*!
   * \brief Converts the principal invariants to the modified principal invariants.
   *
   * @param modinv (out) : Modified principal invariants of the Right Cauchy-Green strain tensor
   * @param prinv (in) : Principal invariants of the Right Cauchy-Green strain tensor
   */
  void InvariantsModified(LINALG::Matrix<3, 1>& modinv, const LINALG::Matrix<3, 1>& prinv);

  /*!
   * \brief Evaluates principal stretches \f$\lambda_\alpha\f$ and material principal directions
   *
   * @param prstr (out) : Principal stretches
   * @param prdir (out) : Principal directions
   * @param rcg (in) : Right Cauchy-Green strain tensor in strain-like Voigt notation
   */
  void StretchesPrincipal(
      LINALG::Matrix<3, 1>& prstr, LINALG::Matrix<3, 3>& prdir, const LINALG::Matrix<6, 1>& rcg);

  /*!
   * \brief Evaluates the modified principal stretches \f$\bar{\lambda}_\alpha\f$
   *
   * @param modstr (out) : Modified principal stretches
   * @param prstr (in) : Principal stretches
   */
  void StretchesModified(LINALG::Matrix<3, 1>& modstr, const LINALG::Matrix<3, 1>& prstr);

  /*!
   * \brief Generates a nxn Identity matrix in tensor notation
   *
   * @param mat Matrix to store the identity tensor
   */
  inline void IdentityMatrix(LINALG::Matrix<3, 3>& mat)
  {
    mat.Clear();
    for (unsigned int i = 0; i < 3; ++i)
    {
      mat(i, i) = 1.0;
    }
  }

  inline void InvariantsPrincipal(LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<3, 3>& tens)
  {
    // 1st invariant, trace tens
    prinv(0) = tens(0, 0) + tens(1, 1) + tens(2, 2);

    // 2nd invariant, 0.5( (trace(tens))^2 - trace(tens^2))
    prinv(1) = tens(0, 0) * tens(1, 1) + tens(1, 1) * tens(2, 2) + tens(0, 0) * tens(2, 2) -
               tens(0, 1) * tens(1, 0) - tens(1, 2) * tens(2, 1) - tens(0, 2) * tens(2, 0);

    // 3rd invariant, determinant tens
    prinv(2) = tens.Determinant();
  }

}  // namespace MAT

#endif
