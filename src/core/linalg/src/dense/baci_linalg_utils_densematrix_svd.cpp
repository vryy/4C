/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of singular value decomposition (SVD) methods for namespace CORE::LINALG

\level 0
*/
/*----------------------------------------------------------------------*/

#include "baci_linalg_utils_densematrix_svd.hpp"

#include "baci_utils_exceptions.hpp"

#include <Teuchos_ScalarTraits.hpp>

FOUR_C_NAMESPACE_OPEN

void CORE::LINALG::SVD(const CORE::LINALG::SerialDenseMatrix::Base& A,
    CORE::LINALG::SerialDenseMatrix& Q, CORE::LINALG::SerialDenseMatrix& S,
    CORE::LINALG::SerialDenseMatrix& VT)
{
  CORE::LINALG::SerialDenseMatrix tmp(A);  // copy, because content of A is destroyed
  const char jobu = 'A';                   // compute and return all M columns of U
  const char jobvt = 'A';                  // compute and return all N rows of V^T
  const int n = tmp.numCols();
  const int m = tmp.numRows();
  std::vector<double> s(std::min(n, m));
  int info;
  int lwork = std::max(3 * std::min(m, n) + std::max(m, n), 5 * std::min(m, n));
  std::vector<double> work(lwork);
  double rwork;

  Teuchos::LAPACK<int, double> lapack;
  lapack.GESVD(jobu, jobvt, m, n, tmp.values(), tmp.stride(), s.data(), Q.values(), Q.stride(),
      VT.values(), VT.stride(), work.data(), lwork, &rwork, &info);

  if (info) FOUR_C_THROW("Lapack's dgesvd returned %d", info);

  // 0 for off-diagonal, otherwise s
  S.putScalar(0.0);
  for (int i = 0; i < std::min(n, m); ++i)
  {
    S(i, i) = s[i];
  }
}

FOUR_C_NAMESPACE_CLOSE
