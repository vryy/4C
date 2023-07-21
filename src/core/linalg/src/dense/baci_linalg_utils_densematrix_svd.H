/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of singular value decomposition (SVD) methods for namespace CORE::LINALG

\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_LINALG_UTILS_DENSEMATRIX_SVD_H
#define BACI_LINALG_UTILS_DENSEMATRIX_SVD_H

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include <Teuchos_LAPACK.hpp>
#include "baci_lib_exporter.H"

namespace CORE::LINALG
{
  /*!
   \brief Compute singular value decomposition (SVD) of a real M-by-N matrix A
   A = U * SIGMA * transpose(V)

   \param A (in/out):    Matrix to be decomposed
   \param U (in/out):    M-by-M orthogonal matrix
   \param SIGMA (in/out):M-by-N matrix which is zero except for its min(m,n) diagonal elements
   \param Vt (in/out):   V is a N-by-N orthogonal matrix, actually returned is V^T
   */
  void SVD(const Epetra_SerialDenseMatrix& A, CORE::LINALG::SerialDenseMatrix& U,
      CORE::LINALG::SerialDenseMatrix& SIGMA, CORE::LINALG::SerialDenseMatrix& Vt);

  /*!
   \brief Singular value decomposition (SVD) of a real M-by-N matrix in fixed
   size format

   A = Q * S * VT

   \tparam row Number of rows
   \tparam col Number of columns

   \param A (in):        M-by-N matrix to be decomposed
   \param Q (out):       M-by-M orthogonal matrix
   \param S (out):       M-by-N matrix which is zero except for its min(m,n) diagonal elements
   \param VT (out):      N-by-N orthogonal matrix (transpose of V)
   */
  template <unsigned int rows, unsigned int cols>
  void SVD(const CORE::LINALG::Matrix<rows, cols>& A, CORE::LINALG::Matrix<rows, rows>& Q,
      CORE::LINALG::Matrix<rows, cols>& S, CORE::LINALG::Matrix<cols, cols>& VT)
  {
    Matrix<rows, cols> tmp(A.A(), false);  // copy, because content of matrix is destroyed
    const char jobu = 'A';                 // compute and return all M columns of U
    const char jobvt = 'A';                // compute and return all N rows of V^T
    std::vector<double> s(std::min(rows, cols));
    int info;
    int lwork = std::max(3 * std::min(rows, cols) + std::max(rows, cols), 5 * std::min(rows, cols));
    std::vector<double> work(lwork);
    double rwork;

    Teuchos::LAPACK<int, double> lapack;
    lapack.GESVD(jobu, jobvt, rows, cols, tmp.A(), tmp.M(), s.data(), Q.A(), Q.M(), VT.A(), VT.M(),
        work.data(), lwork, &rwork, &info);

    if (info) dserror("Lapack's dgesvd returned %d", info);

    for (unsigned int i = 0; i < std::min(rows, cols); ++i)
    {
      for (unsigned int j = 0; j < std::min(rows, cols); ++j)
      {
        S(i, j) = (i == j) * s[i];  // 0 for off-diagonal, otherwise s
      }
    }
    return;
  }

}  // namespace CORE::LINALG

#endif  // LINALG_UTILS_DENSEMATRIX_SVD_H
