/*----------------------------------------------------------------------*/
/*! \file

\brief Determinant functions for dense matrices up to 4x4 and LU determinant

\level 0
*/
/*----------------------------------------------------------------------*/

#include <Teuchos_LAPACK.hpp>
#include "linalg_utils_densematrix_determinant.H"
#include "utils_exceptions.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double CORE::LINALG::DeterminantLU(const Epetra_SerialDenseMatrix& A)
{
#ifdef DEBUG
  if (A.M() != A.N()) dserror("Matrix is not square");
#endif
  Epetra_SerialDenseMatrix tmp(A);
  const int n = tmp.N();
  const int m = tmp.M();
  std::vector<int> ipiv(n);
  int info;

  Teuchos::LAPACK<int, double> lapack;
  lapack.GETRF(m, n, tmp.A(), tmp.LDA(), ipiv.data(), &info);

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
