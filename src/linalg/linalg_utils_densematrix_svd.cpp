/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of singular value decomposition (SVD) methods for namespace LINALG

\level 0
*/
/*----------------------------------------------------------------------*/

#include <Teuchos_ScalarTraits.hpp>

#include "linalg_utils_densematrix_svd.H"
#include "dserror.H"

void LINALG::SVD(const Epetra_SerialDenseMatrix& A, LINALG::SerialDenseMatrix& Q,
    LINALG::SerialDenseMatrix& S, LINALG::SerialDenseMatrix& VT)
{
  Epetra_SerialDenseMatrix tmp(A);  // copy, because content of A is destroyed
  const char jobu = 'A';            // compute and return all M columns of U
  const char jobvt = 'A';           // compute and return all N rows of V^T
  const int n = tmp.N();
  const int m = tmp.M();
  std::vector<double> s(std::min(n, m));
  int info;
  int lwork = std::max(3 * std::min(m, n) + std::max(m, n), 5 * std::min(m, n));
  std::vector<double> work(lwork);
  double rwork;

  Teuchos::LAPACK<int, double> lapack;
  lapack.GESVD(jobu, jobvt, m, n, tmp.A(), tmp.LDA(), &s[0], Q.A(), Q.LDA(), VT.A(), VT.LDA(),
      &work[0], lwork, &rwork, &info);

  if (info) dserror("Lapack's dgesvd returned %d", info);

  // 0 for off-diagonal, otherwise s
  S.Zero();
  for (int i = 0; i < std::min(n, m); ++i)
  {
    S(i, i) = s[i];
  }
}
