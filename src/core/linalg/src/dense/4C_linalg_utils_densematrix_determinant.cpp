/*----------------------------------------------------------------------*/
/*! \file

\brief Determinant functions for dense matrices up to 4x4 and LU determinant

\level 0
*/
/*----------------------------------------------------------------------*/

#include "4C_linalg_utils_densematrix_determinant.hpp"

#include "4C_utils_exceptions.hpp"

#include <Teuchos_LAPACK.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Core::LinAlg::DeterminantLU(const Core::LinAlg::SerialDenseMatrix& A)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (A.numRows() != A.numCols()) FOUR_C_THROW("Matrix is not square");
#endif
  Core::LinAlg::SerialDenseMatrix tmp(A);
  const int n = tmp.numCols();
  const int m = tmp.numRows();
  std::vector<int> ipiv(n);
  int info;

  Teuchos::LAPACK<int, double> lapack;
  lapack.GETRF(m, n, tmp.values(), tmp.stride(), ipiv.data(), &info);

  if (info < 0)
    FOUR_C_THROW("Lapack's dgetrf returned %d", info);
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

FOUR_C_NAMESPACE_CLOSE
