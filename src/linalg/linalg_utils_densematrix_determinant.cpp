/*----------------------------------------------------------------------*/
/*! \file

\brief Determinant functions for dense matrices up to 4x4 and LU determinant

\level 0
\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/

#include "../headers/compiler_definitions.h" /* access to fortran routines */
#include "linalg_utils_densematrix_determinant.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double LINALG::DeterminantLU(const Epetra_SerialDenseMatrix& A)
{
#ifdef DEBUG
  if (A.M() != A.N()) dserror("Matrix is not square");
#endif
  Epetra_SerialDenseMatrix tmp(A);
  Epetra_LAPACK lapack;
  const int n = tmp.N();
  const int m = tmp.M();
  std::vector<int> ipiv(n);
  int info;
  lapack.GETRF(m, n, tmp.A(), tmp.LDA(), &ipiv[0], &info);
  if (info < 0)
    dserror("Lapack's dgetrf returned %d", info);
  else if (info > 0)
    return 0.0;
  double d = tmp(0, 0);
  for (int i = 1; i < n; ++i) d *= tmp(i, i);
  // swapping rows of A changes the sign of the determinant, so we have to
  // undo lapack's permutation w.r.t. the determinant
  // note the fortran indexing convention in ipiv
  for (int i = 0; i < n; ++i)
    if (ipiv[i] != i + 1) d *= -1.0;
  return d;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double LINALG::DeterminantSVD(const Epetra_SerialDenseMatrix& A)
{
  // todo: the sign of the determinant is wrong
  // once this method is fixed activate unit tests!
  dserror("the sign of the determinant is wrong due to a bug in this method!");

#ifdef DEBUG
  if (A.M() != A.N()) dserror("Matrix is not square");
#endif
  Epetra_SerialDenseMatrix tmp(A);
  Epetra_LAPACK lapack;
  const int n = tmp.N();
  const int m = tmp.M();
  std::vector<double> s(std::min(n, m));
  int info;
  int lwork = std::max(3 * std::min(m, n) + std::max(m, n), 5 * std::min(m, n));
  std::vector<double> work(lwork);
  lapack.GESVD('N', 'N', m, n, tmp.A(), tmp.LDA(), &s[0], NULL, tmp.LDA(), NULL, tmp.LDA(),
      &work[0], &lwork, &info);
  if (info) dserror("Lapack's dgesvd returned %d", info);
  double d = s[0];
  for (int i = 1; i < n; ++i) d *= s[i];
  return d;
}
