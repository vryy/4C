/*----------------------------------------------------------------------*/
/*!
\brief Various service routines related to materials

\level 1

\maintainer Fabian Braeu
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
      LINALG::TMatrix<T, 6, 6>& cmat, const LINALG::TMatrix<T, 6, 1>& invc, const T scalar);


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
  void AddRightNonSymmetricHolzapfelProduct(LINALG::TMatrix<T, 6, 9>& out,
      LINALG::TMatrix<T, 3, 3> const& A, LINALG::TMatrix<T, 3, 3> const& B, T const fac);

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
  void AddRightNonSymmetricHolzapfelProductStrainLike(LINALG::TMatrix<T, 6, 9>& out,
      LINALG::TMatrix<T, 3, 3> const& A, LINALG::TMatrix<T, 3, 3> const& B, T const fac);

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

  /// Convert matrix to "stress-like" Voigt notation
  ///          A_00 A_01 A_02
  /// Matrix   A_10 A_11 A_12   is converted to the (stress-like) vector below:
  ///          A_20 A_21 A_22
  ///
  /// Vector   V_0 = A_00; V_1 = A_11; V_2 = A_22; V_3 = 0.5*(A_01 + A_10); V_4 = 0.5*(A_12 + A_21);
  /// V_5 = 0.5*(A_02 + A_20)
  template <typename T>
  void MatrixtoStressLikeVoigtNotation(
      LINALG::TMatrix<T, 3, 3> const& in, LINALG::TMatrix<T, 6, 1>& out);

  /// Convert matrix to "strain-like" Voigt notation
  ///          A_00 A_01 A_02
  /// Matrix   A_10 A_11 A_12   is converted to the (strain-like) vector below:
  ///          A_20 A_21 A_22
  ///
  /// Vector   V_0 = A_00; V_1 = A_11; V_2 = A_22; V_3 = (A_01 + A_10); V_4 = (A_12 + A_21); V_5 =
  /// (A_02 + A_20)
  void MatrixtoStrainLikeVoigtNotation(LINALG::Matrix<3, 3> const& in, LINALG::Matrix<6, 1>& out);

  /// Convert "stress-like" Voigt notation to matrix
  ///
  /// Vector   V_0 V_1 V_2 V_3 V_4 V_5   (stress-like) is converted to the matrix below:
  ///
  /// Matrix   A_00 = V_0; A_11 = V_1; A_22 = V_2; A_01 = A_10 = V_3; A_21 = A_12 = V_4; A_20 = A_02
  /// = V_5
  void StressLikeVoigtNotationtoMatrix(LINALG::Matrix<6, 1> const& in, LINALG::Matrix<3, 3>& out);

  /// convert non-symmetric 2-tensor to 9x1 vector
  ///          A_00 A_01 A_02
  /// Matrix   A_10 A_11 A_12  is converted to the vector below:
  ///          A_20 A_21 A_22
  ///
  /// Vector   V_0 = A_00; V_1 = A_11; V_2 = A_22; V_3 = A_01; V_4 = A_12; V_5 = A_02; V_6 = A_10;
  /// V_7 = A_21; V_8 = A_20
  void Matrix3x3to9x1(LINALG::Matrix<3, 3> const& in, LINALG::Matrix<9, 1>& out);

  /*--------------------------------------------------------------------------*/
  /// Voigt notation types
  enum class Notation
  {
    stress,
    strain
  };

  /*--------------------------------------------------------------------------*/
  /// Index maps
  struct IMap
  {
    IMap() = default;

    /// map second order tensor indices to perturbed Voigt vector notation
    static const unsigned second_[3][3];

    /// map perturbed Voigt matrix notation for the fourth order elasticity
    /// tensor to the fourth order indices
    static const unsigned fourth_[6][6][4];
  };

  /*--------------------------------------------------------------------------*/
  /** \brief Utility routines for the perturbed Voigt tensor notation
   *
   *  \author hiermeier \date 09/17 */
  template <Notation type>
  class VoigtUtils
  {
   public:
    /// constructor
    VoigtUtils() = default;

    /// destructor
    ~VoigtUtils() = default;

    /** \brief compute power of a symmetric 3x3 matrix in perturbed Voigt notation
     *
     *  \f[
     *  [vtensorpow]_{AB} = [vtensor^{pow}]_{AB}
     *  \f]
     *
     *  \param[in] pow          positive integer exponent
     *  \param[in] vtensor      input tensor in Voigt <type> notation
     *  \param[in] vtensor_pow  result, i.e. input tensor to the given power
     *
     *  \author hiermeier \date 09/17 */
    static void PowerOfSymmetricTensor(
        const unsigned pow, const LINALG::Matrix<6, 1>& vtensor, LINALG::Matrix<6, 1>& vtensor_pow);

    /** \brief Compute the inverse tensor in perturbed Voigt notation
     *
     *  \param[in]  vtensor      tensor in Voigt <type> notation
     *  \param[out] vtensor_inv  inverse tensor in Voigt <type> notation
     *
     *  \author hiermeier \date 09/17 */
    static void InverseTensor(
        const LINALG::Matrix<6, 1>& vtensor, LINALG::Matrix<6, 1>& vtensor_inv);

    /** \brief Compute the product of a tensor in perturbed Voigt notation
     *  and a vector
     *
     *  \f$ [vecres]_{A} = vtensor_{AB} vec^{B} \f$
     *
     *  \param[in]  vtensor      tensor in Voigt <type> notation
     *  \param[out] vtensor_inv  inverser tensor in Voigt <type> notation
     *
     *  \author hiermeier \date 09/17 */
    static void MultiplyTensorVector(const LINALG::Matrix<6, 1>& vtensor,
        const LINALG::Matrix<3, 1>& vec, LINALG::Matrix<3, 1>& vec_res);

    /** \brief Compute the symmetric outer product of two vectors
     *
     *  \f$ [abba]_{AB} = [veca]_{A} [vecb]_{B} + [veca]_{B} [vecb]_{A} \f$
     *
     *  \param[in]  vec_a  first vector
     *  \param[in]  vec_b  second vector
     *  \param[out] ab_ba  symmetric outer product of the two input vectors
     *                     in the Voigt <type> notation
     *
     *  \author hiermeier \date 09/17 */
    static void SymmetricOuterProduct(const LINALG::Matrix<3, 1>& vec_a,
        const LINALG::Matrix<3, 1>& vec_b, LINALG::Matrix<6, 1>& ab_ba);

    /** \brief scale off diagonal values
     *
     *  \note This function changes the values only if the strain notation is used.
     *
     *  \param[out] tensor  scale the off-diagonal values of this tensor
     *
     *  \author hiermeier \date 09/17 */
    static void ScaleOffDiagonalVals(LINALG::Matrix<6, 1>& tensor);

    /** \brief unscale off diagonal values
     *
     *  \note This function changes the values only if the strain notation is used.
     *
     *  \param[out] tensor  unscale the off-diagonal values of this tensor
     *
     *  \author hiermeier \date 09/17 */
    static void UnscaleOffDiagonalVals(LINALG::Matrix<6, 1>& strain);

    /** \brief unscale factors for the perturbed Voigt strain notation
     *
     *  These factors are meaningful if the strain convention is followed. */
    static const double unscale_fac_[6];

    /** \brief scale factors for the perturbed Voigt stress notation
     *
     *  These factors are meaningful if the strain convention is followed. */
    static const double scale_fac_[6];
  };

  /// @name type defs for easier use
  /// @{

  /// Utility class for the perturbed Voigt stress notation
  typedef VoigtUtils<Notation::stress> VStressUtils;

  /// Utility class for the perturbed Voigt strain notation
  typedef VoigtUtils<Notation::strain> VStrainUtils;

  /// @}

  /// access scaling factors
  template <Notation type>
  inline double ScaleFactor(unsigned i);
  template <>
  inline double ScaleFactor<Notation::stress>(unsigned i)
  {
    return 1.0;
  };
  template <>
  inline double ScaleFactor<Notation::strain>(unsigned i)
  {
    return VStrainUtils::scale_fac_[i];
  };

  /// access unscaling factors
  template <Notation type>
  inline double UnscaleFactor(unsigned i);
  template <>
  inline double UnscaleFactor<Notation::stress>(unsigned i)
  {
    return 1.0;
  };
  template <>
  inline double UnscaleFactor<Notation::strain>(unsigned i)
  {
    return VStrainUtils::unscale_fac_[i];
  };

  /** \brief convert tensor from perturbed Voigt strain notation to
   *  perturbed Voigt stress notation
   *
   * \author hiermeier \date 09/17 */
  inline void ConvertStrainToStressNotation(LINALG::Matrix<6, 1>& strain)
  {
    VStrainUtils::UnscaleOffDiagonalVals(strain);
  }

  /** \brief convert tensor from perturbed Voigt stress notation to
   *  perturbed Voigt strain notation
   *
   * \author hiermeier \date 09/17 */
  inline void ConvertStressToStrainNotation(LINALG::Matrix<6, 1>& strain)
  {
    VStrainUtils::ScaleOffDiagonalVals(strain);
  }
}  // namespace MAT

#endif
