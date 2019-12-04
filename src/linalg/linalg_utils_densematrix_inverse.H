/*----------------------------------------------------------------------*/
/*! \file

\brief Inverse dense matrices up to 4x4 and other inverse methods.

\level 0
\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/
#ifndef LINALG_UTILS_DENSEMATRIX_INVERSE_H
#define LINALG_UTILS_DENSEMATRIX_INVERSE_H

#include "linalg_utils_densematrix_determinant.H"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "../drt_lib/drt_exporter.H"

namespace LINALG
{
  /**
   * \brief Reorder the matrix entries for the inverse of a nonsymmetric 2x2 matrix.
   *
   * @param A (in/out) Matrix A.
   */
  template <typename T>
  void Inverse2x2ReorderMatrixEntries(LINALG::Matrix<2, 2, T>& A)
  {
    T b00 = A(0, 0);
    T b01 = A(0, 1);
    T b10 = A(1, 0);
    T b11 = A(1, 1);
    A(0, 0) = b11;
    A(1, 0) = -b10;
    A(0, 1) = -b01;
    A(1, 1) = b00;
  }

  /*!
   * \brief Explicit inverse of a nonsymmetric 2x2 matrix. If the matrix is singular an error will
   * be thrown.
   *
   * @param A (in/out) Matrix to be inverted. The inverse will be stored in this variable.
   */
  template <typename T>
  void Inverse2x2(LINALG::Matrix<2, 2, T>& A)
  {
    T det = Determinant2x2(A);
    if (det == 0.0) dserror("Determinant of 2x2 matrix is exactly zero");
    Inverse2x2ReorderMatrixEntries(A);
    A.Scale(1. / det);
  }

  /*!
   * \brief Explicit inverse of a nonsymmetric 2x2 matrix. If the matrix is singular no error will
   * be thrown, and the original matrix will remain untouched.
   *
   * @param A (in/out) Matrix to be inverted. The inverse will be stored in this variable.
   * @param eps (in) If the absolute value of the determinant is smaller than this value, the matrix
   * is considered to be singular.
   * @return True if matrix could be inverted, false if matrix is singular.
   */
  template <typename T>
  bool Inverse2x2DoNotThrowErrorOnZeroDeterminant(LINALG::Matrix<2, 2, T>& A, double eps)
  {
    T det = Determinant2x2(A);
    if (abs(det) < eps) return false;
    Inverse2x2ReorderMatrixEntries(A);
    A.Scale(1. / det);
    return true;
  }

  /**
   * \brief Reorder the matrix entries for the inverse of a nonsymmetric 3x3 matrix.
   *
   * @param A (in/out) Matrix A.
   */
  template <typename T>
  void Inverse3x3ReorderMatrixEntries(LINALG::Matrix<3, 3, T>& A)
  {
    T b00 = A(0, 0);
    T b01 = A(0, 1);
    T b02 = A(0, 2);
    T b10 = A(1, 0);
    T b11 = A(1, 1);
    T b12 = A(1, 2);
    T b20 = A(2, 0);
    T b21 = A(2, 1);
    T b22 = A(2, 2);
    A(0, 0) = b11 * b22 - b21 * b12;
    A(1, 0) = -b10 * b22 + b20 * b12;
    A(2, 0) = b10 * b21 - b20 * b11;
    A(0, 1) = -b01 * b22 + b21 * b02;
    A(1, 1) = b00 * b22 - b20 * b02;
    A(2, 1) = -b00 * b21 + b20 * b01;
    A(0, 2) = b01 * b12 - b11 * b02;
    A(1, 2) = -b00 * b12 + b10 * b02;
    A(2, 2) = b00 * b11 - b10 * b01;
  }

  /*!
   * \brief Explicit inverse of a nonsymmetric 3x3 matrix. If the matrix is singular an error will
   * be thrown.
   *
   * @param A (in/out) Matrix to be inverted. The inverse will be stored in this variable.
   */
  template <typename T>
  void Inverse3x3(LINALG::Matrix<3, 3, T>& A)
  {
    T det = Determinant3x3(A);
    if (det == 0.0) dserror("Determinant of 3x3 matrix is exactly zero");
    Inverse3x3ReorderMatrixEntries(A);
    A.Scale(1. / det);
  }

  /*!
   * \brief Explicit inverse of a nonsymmetric 3x3 matrix. If the matrix is singular no error will
   * be thrown, and the original matrix will remain untouched.
   *
   * @param A (in/out) Matrix to be inverted. The inverse will be stored in this variable.
   * @param eps (in) If the absolute value of the determinant is smaller than this value, the matrix
   * is considered to be singular.
   * @return True if matrix could be inverted, false if matrix is singular.
   */
  template <typename T>
  bool Inverse3x3DoNotThrowErrorOnZeroDeterminant(LINALG::Matrix<3, 3, T>& A, double eps)
  {
    T det = Determinant3x3(A);
    if (abs(det) < eps) return false;
    Inverse3x3ReorderMatrixEntries(A);
    A.Scale(1. / det);
    return true;
  }

  /**
   * \brief Reorder the matrix entries for the inverse of a nonsymmetric 4x4 matrix.
   *
   * @param A (in/out) Matrix A.
   */
  template <typename T>
  void Inverse4x4ReorderMatrixEntries(LINALG::Matrix<4, 4, T>& A)
  {
    T a00 = A(0, 0);
    T a01 = A(0, 1);
    T a02 = A(0, 2);
    T a03 = A(0, 3);

    T a10 = A(1, 0);
    T a11 = A(1, 1);
    T a12 = A(1, 2);
    T a13 = A(1, 3);

    T a20 = A(2, 0);
    T a21 = A(2, 1);
    T a22 = A(2, 2);
    T a23 = A(2, 3);

    T a30 = A(3, 0);
    T a31 = A(3, 1);
    T a32 = A(3, 2);
    T a33 = A(3, 3);

    A(0, 0) = a11 * a33 * a22 - a11 * a32 * a23 - a31 * a13 * a22 + a32 * a21 * a13 +
              a31 * a12 * a23 - a33 * a21 * a12;
    A(0, 1) = -(a31 * a02 * a23 + a33 * a22 * a01 - a33 * a21 * a02 - a31 * a03 * a22 -
                a32 * a23 * a01 + a32 * a21 * a03);
    A(0, 2) = -a11 * a33 * a02 + a11 * a32 * a03 + a33 * a01 * a12 + a31 * a13 * a02 -
              a32 * a01 * a13 - a31 * a12 * a03;
    A(0, 3) = -(-a11 * a02 * a23 + a11 * a03 * a22 - a03 * a21 * a12 - a01 * a13 * a22 +
                a02 * a21 * a13 + a01 * a12 * a23);

    A(1, 0) = -(a33 * a22 * a10 - a32 * a23 * a10 + a30 * a12 * a23 - a33 * a20 * a12 -
                a30 * a13 * a22 + a32 * a20 * a13);
    A(1, 1) = a00 * a33 * a22 - a00 * a32 * a23 + a30 * a02 * a23 - a33 * a20 * a02 +
              a32 * a20 * a03 - a30 * a03 * a22;
    A(1, 2) = -(a00 * a33 * a12 - a00 * a32 * a13 - a33 * a10 * a02 - a30 * a03 * a12 +
                a32 * a10 * a03 + a30 * a02 * a13);
    A(1, 3) = -a00 * a13 * a22 + a00 * a12 * a23 - a10 * a02 * a23 + a13 * a20 * a02 +
              a10 * a03 * a22 - a12 * a20 * a03;

    A(2, 0) = a11 * a30 * a23 - a11 * a33 * a20 - a31 * a10 * a23 + a33 * a21 * a10 +
              a31 * a13 * a20 - a30 * a21 * a13;
    A(2, 1) = -(-a00 * a31 * a23 + a00 * a33 * a21 - a30 * a03 * a21 + a31 * a20 * a03 +
                a30 * a01 * a23 - a33 * a20 * a01);
    A(2, 2) = a33 * a00 * a11 - a33 * a10 * a01 - a30 * a03 * a11 - a31 * a00 * a13 +
              a31 * a10 * a03 + a30 * a01 * a13;
    A(2, 3) = -(a23 * a00 * a11 - a23 * a10 * a01 - a20 * a03 * a11 - a21 * a00 * a13 +
                a21 * a10 * a03 + a20 * a01 * a13);

    A(3, 0) = -(a11 * a30 * a22 - a11 * a32 * a20 - a30 * a21 * a12 - a31 * a10 * a22 +
                a32 * a21 * a10 + a31 * a12 * a20);
    A(3, 1) = -a00 * a31 * a22 + a00 * a32 * a21 + a31 * a20 * a02 + a30 * a01 * a22 -
              a32 * a20 * a01 - a30 * a02 * a21;
    A(3, 2) = -(a32 * a00 * a11 - a32 * a10 * a01 - a30 * a02 * a11 - a31 * a00 * a12 +
                a31 * a10 * a02 + a30 * a01 * a12);
    A(3, 3) = a22 * a00 * a11 - a22 * a10 * a01 - a20 * a02 * a11 - a21 * a00 * a12 +
              a21 * a10 * a02 + a20 * a01 * a12;
  }

  /*!
   * \brief Explicit inverse of a nonsymmetric 4x4 matrix. If the matrix is singular an error will
   * be thrown.
   *
   * @param A (in/out) Matrix to be inverted. The inverse will be stored in this variable.
   */
  template <typename T>
  void Inverse4x4(LINALG::Matrix<4, 4, T>& A)
  {
    T det = Determinant4x4(A);
    if (det == 0.0) dserror("Determinant of 4x4 matrix is exactly zero");
    Inverse4x4ReorderMatrixEntries(A);
    A.Scale(1. / det);
  }

  /*!
   * \brief Explicit inverse of a nonsymmetric 4x4 matrix. If the matrix is singular no error will
   * be thrown, and the original matrix will remain untouched.
   *
   * @param A (in/out) Matrix to be inverted. The inverse will be stored in this variable.
   * @param eps (in) If the absolute value of the determinant is smaller than this value, the matrix
   * is considered to be singular.
   * @return True if matrix could be inverted, false if matrix is singular.
   */
  template <typename T>
  bool Inverse4x4DoNotThrowErrorOnZeroDeterminant(LINALG::Matrix<4, 4, T>& A, double eps)
  {
    T det = Determinant4x4(A);
    if (abs(det) < eps) return false;
    Inverse4x4ReorderMatrixEntries(A);
    A.Scale(1. / det);
    return true;
  }

  /*!
    \brief Solve soe with me*ae^T = de^T (employed for calculating coefficient matrix for dual shape
    functions) farah 07/14

    \param me (in):  me (not transposed!)
    \paramdeA (in):  de (not transposed!)
    \param ae(out):  coef. matrix (not transposed!)
    */
  LINALG::SerialDenseMatrix InvertAndMultiplyByCholesky(
      LINALG::SerialDenseMatrix& me, LINALG::SerialDenseMatrix& de, LINALG::SerialDenseMatrix& ae);

  /*!
   \brief Solve soe with me*ae^T = de^T (employed for calculating coefficient matrix for dual shape
   functions) farah 07/14

   \param me (in):  me (not transposed!)
   \paramdeA (in):  de (not transposed!)
   \param ae(out):  coef. matrix (not transposed!)
   */
  template <const int n>
  void InvertAndMultiplyByCholesky(
      LINALG::Matrix<n, n>& me, LINALG::Matrix<n, n>& de, LINALG::SerialDenseMatrix& ae)
  {
    LINALG::Matrix<n, n> y;
    LINALG::Matrix<n, n> maux = me;

    // calc G with me=G*G^T
    for (int z = 0; z < n; ++z)
    {
      for (int u = 0; u < z + 1; ++u)
      {
        double sum = me(z, u);
        for (int k = 0; k < u; ++k) sum -= me(z, k) * me(u, k);

        if (z > u)
          me(z, u) = sum / me(u, u);
        else if (sum > 0.0)
          me(z, z) = sqrt(sum);
        else
          dserror("matrix is not positive definite!");
      }

      // get y for G*y=De
      const double yfac = 1.0 / me(z, z);
      for (int col = 0; col < n; ++col)
      {
        y(z, col) = yfac * de(col, z);
        for (int u = 0; u < z; ++u) y(z, col) -= yfac * me(z, u) * y(u, col);
      }
    }

    // get y for G^T*x=y
    for (int z = n - 1; z > -1; --z)
    {
      const double xfac = 1.0 / me(z, z);
      for (int col = 0; col < n; ++col)
      {
        ae(col, z) = xfac * y(z, col);
        for (int u = n - 1; u > z; --u) ae(col, z) -= xfac * me(u, z) * ae(col, u);
      }
    }

    me = maux;
    return;
  }

  /*!
   \brief Solve soe with me*ae^T = de^T (employed for calculating coefficient matrix for dual shape
   functions) farah 07/14

   \param me (in):  me (not transposed!)
   \paramdeA (in):  de (not transposed!)
   \param ae(out):  coef. matrix (not transposed!)
   */
  template <const int n>
  void InvertAndMultiplyByCholesky(
      LINALG::Matrix<n, n>& me, LINALG::Matrix<n, n>& de, LINALG::Matrix<n, n>& ae)
  {
    LINALG::Matrix<n, n> y;

    // calc G with me=G*G^T
    for (int z = 0; z < n; ++z)
    {
      for (int u = 0; u < z + 1; ++u)
      {
        double sum = me(z, u);
        for (int k = 0; k < u; ++k) sum -= me(z, k) * me(u, k);

        if (z > u)
          me(z, u) = sum / me(u, u);
        else if (sum > 0.0)
          me(z, z) = sqrt(sum);
        else
          dserror("matrix is not positive definite!");
      }

      // get y for G*y=de^T
      const double yfac = 1.0 / me(z, z);
      for (int col = 0; col < n; ++col)
      {
        y(z, col) = yfac * de(col, z);
        for (int u = 0; u < z; ++u) y(z, col) -= yfac * me(z, u) * y(u, col);
      }
    }

    // get x for G^T*x=y
    for (int z = n - 1; z > -1; --z)
    {
      const double xfac = 1.0 / me(z, z);
      for (int col = 0; col < n; ++col)
      {
        ae(col, z) = xfac * y(z, col);
        for (int u = n - 1; u > z; --u) ae(col, z) -= xfac * me(u, z) * ae(col, u);
      }
    }

    return;
  }

  /*!
   \brief Invert a symmetric dim*dim square matrix (positive definite)
   farah 07/14
   \param A (in/out): Matrix to be inverted
   */
  template <const int n>
  void SymmetricPositiveDefiniteInverse(LINALG::Matrix<n, n>& A)
  {
    LINALG::Matrix<n, n> y(true);
    LINALG::Matrix<n, n> ae(true);

    // calc G with me=G*G^T
    for (int z = 0; z < n; ++z)
    {
      for (int u = 0; u < z + 1; ++u)
      {
        double sum = A(z, u);
        for (int k = 0; k < u; ++k) sum -= A(z, k) * A(u, k);

        if (z > u)
          A(z, u) = sum / A(u, u);
        else if (sum > 0.0)
          A(z, z) = sqrt(sum);
        else
          dserror("matrix is not positive definite!");
      }

      // get y for G*y=de^T
      const double yfac = 1.0 / A(z, z);
      for (int col = 0; col < n; ++col)
      {
        if (col == z) y(z, col) = yfac;
        for (int u = 0; u < z; ++u) y(z, col) -= yfac * A(z, u) * y(u, col);
      }
    }

    // get x for G^T*x=y
    for (int z = n - 1; z > -1; --z)
    {
      const double xfac = 1.0 / A(z, z);
      for (int col = 0; col < n; ++col)
      {
        ae(z, col) = xfac * y(z, col);
        for (int u = n - 1; u > z; --u) ae(z, col) -= xfac * A(u, z) * ae(u, col);
      }
    }

    // get result
    A = ae;

    return;
  }

  /*!
   \brief Invert a symmetric dim*dim square matrix

   \param A (in/out): Matrix to be inverted
   \param dim (in) :  Dimension of matrix
   */
  void SymmetricInverse(Epetra_SerialDenseMatrix& A, const int dim);

}  // namespace LINALG

#endif
