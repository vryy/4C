/*----------------------------------------------------------------------*/
/*! \file
\file linalg_utils_densematrix_determinant.H

\brief Determinant functions for dense matrices up to 4x4.

<pre>
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*-----------------------------------------------------------------------*/
#ifndef LINALG_UTILS_DENSEMATRIX_DETERMINANT_H
#define LINALG_UTILS_DENSEMATRIX_DETERMINANT_H


#include "linalg_fixedsizematrix.H"


namespace LINALG
{
  /**
   *\brief Explicit determinant of a nonsymmetric 2x2 matrix.
   *
   * @param A (in) Matrix A.
   * @return Determinant of matrix A.
   */
  template <typename T>
  T Determinant2x2(const LINALG::Matrix<2, 2, T>& A)
  {
    T b00 = A(0, 0);
    T b01 = A(0, 1);
    T b10 = A(1, 0);
    T b11 = A(1, 1);
    T det = b00 * b11 - b01 * b10;
    return det;
  }

  /**
   *\brief Explicit determinant of a nonsymmetric 3x3 matrix.
   *
   * @param A (in) Matrix A.
   * @return Determinant of matrix A.
   */
  template <typename T>
  T Determinant3x3(const LINALG::Matrix<3, 3, T>& A)
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
    T a = b11 * b22 - b21 * b12;
    T b = -b10 * b22 + b20 * b12;
    T c = b10 * b21 - b20 * b11;
    T det = b00 * a + b01 * b + b02 * c;
    return det;
  }

  /**
   *\brief Explicit determinant of a nonsymmetric 4x4 matrix.
   *
   * @param A (in) Matrix A.
   * @return Determinant of matrix A.
   */
  template <typename T>
  T Determinant4x4(const LINALG::Matrix<4, 4, T>& A)
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

    T det = a00 * a11 * a33 * a22 - a00 * a11 * a32 * a23 - a00 * a31 * a13 * a22 +
            a00 * a32 * a21 * a13 + a00 * a31 * a12 * a23 - a00 * a33 * a21 * a12 +
            a11 * a30 * a02 * a23 - a11 * a33 * a20 * a02 + a11 * a32 * a20 * a03 -
            a11 * a30 * a03 * a22 - a31 * a10 * a02 * a23 - a33 * a22 * a10 * a01 +
            a33 * a21 * a10 * a02 + a33 * a20 * a01 * a12 + a30 * a03 * a21 * a12 +
            a31 * a13 * a20 * a02 + a31 * a10 * a03 * a22 + a30 * a01 * a13 * a22 +
            a32 * a23 * a10 * a01 - a32 * a21 * a10 * a03 - a32 * a20 * a01 * a13 -
            a30 * a02 * a21 * a13 - a31 * a12 * a20 * a03 - a30 * a01 * a12 * a23;
    return det;
  }
}  // namespace LINALG

#endif
