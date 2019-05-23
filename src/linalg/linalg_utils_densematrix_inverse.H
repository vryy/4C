/*!----------------------------------------------------------------------
\file linalg_utils_densematrix_inverse.H

\brief Inverte dense matrices up to 4x4.

<pre>
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*-----------------------------------------------------------------------*/
#ifndef LINALG_UTILS_DENSEMATRIX_INVERSE_H
#define LINALG_UTILS_DENSEMATRIX_INVERSE_H


#include "linalg_utils_densematrix_determinant.H"


namespace LINALG
{
  /**
   * \brief Reorder the matrix entries for the inverse of a nonsymmetric 2x2 matrix.
   *
   * @param A (in/out) Matrix A.
   */
  template <typename T>
  void Inverse2x2ReorderMatrixEntries(LINALG::TMatrix<T, 2, 2>& A)
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
  void Inverse2x2(LINALG::TMatrix<T, 2, 2>& A)
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
  bool Inverse2x2DoNotThrowErrorOnZeroDeterminant(LINALG::TMatrix<T, 2, 2>& A, double eps)
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
  void Inverse3x3ReorderMatrixEntries(LINALG::TMatrix<T, 3, 3>& A)
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
  void Inverse3x3(LINALG::TMatrix<T, 3, 3>& A)
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
  bool Inverse3x3DoNotThrowErrorOnZeroDeterminant(LINALG::TMatrix<T, 3, 3>& A, double eps)
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
  void Inverse4x4ReorderMatrixEntries(LINALG::TMatrix<T, 4, 4>& A)
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
  void Inverse4x4(LINALG::TMatrix<T, 4, 4>& A)
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
  bool Inverse4x4DoNotThrowErrorOnZeroDeterminant(LINALG::TMatrix<T, 4, 4>& A, double eps)
  {
    T det = Determinant4x4(A);
    if (abs(det) < eps) return false;
    Inverse4x4ReorderMatrixEntries(A);
    A.Scale(1. / det);
    return true;
  }

}  // namespace LINALG

#endif
