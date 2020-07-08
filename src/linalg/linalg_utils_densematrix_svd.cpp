/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of singular value decomposition (SVD) methods for namespace LINALG

\level 0
*/
/*----------------------------------------------------------------------*/

#include "../headers/compiler_definitions.h" /* access to fortran routines */
#include "linalg_utils_densematrix_svd.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 |  singular value decomposition (SVD) of a real M-by-N matrix A.       |
 |  Wrapper for Lapack/Epetra_Lapack                           maf 05/08|
 *----------------------------------------------------------------------*/
void LINALG::SVD(const Epetra_SerialDenseMatrix& A, LINALG::SerialDenseMatrix& Q,
    LINALG::SerialDenseMatrix& S, LINALG::SerialDenseMatrix& VT)
{
  Epetra_SerialDenseMatrix tmp(A);  // copy, because content of A is destroyed
  Epetra_LAPACK lapack;
  const char jobu = 'A';   // compute and return all M columns of U
  const char jobvt = 'A';  // compute and return all N rows of V^T
  const int n = tmp.N();
  const int m = tmp.M();
  std::vector<double> s(std::min(n, m));
  int info;
  int lwork = std::max(3 * std::min(m, n) + std::max(m, n), 5 * std::min(m, n));
  std::vector<double> work(lwork);

  lapack.GESVD(jobu, jobvt, m, n, tmp.A(), tmp.LDA(), &s[0], Q.A(), Q.LDA(), VT.A(), VT.LDA(),
      &work[0], &lwork, &info);

  if (info) dserror("Lapack's dgesvd returned %d", info);

  for (int i = 0; i < std::min(n, m); ++i)
  {
    for (int j = 0; j < std::min(n, m); ++j)
    {
      S(i, j) = (i == j) * s[i];  // 0 for off-diagonal, otherwise s
    }
  }
  return;
}
