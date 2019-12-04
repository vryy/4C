/*----------------------------------------------------------------------*/
/*! \file

\brief Inverse dense matrices up to 4x4 and other inverse methods.

\level 0
\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/

#include "../headers/compiler_definitions.h" /* access to fortran routines */
#include "linalg_utils_densematrix_inverse.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*
 |  invert a dense symmetric matrix         )                mwgee 12/06|
 *----------------------------------------------------------------------*/
void LINALG::SymmetricInverse(Epetra_SerialDenseMatrix& A, const int dim)
{
  if (A.M() != A.N()) dserror("Matrix is not square");
  if (A.M() != dim) dserror("Dimension supplied does not match matrix");

  double* a = A.A();
  char uplo[5];
  strncpy(uplo, "L ", 2);
  std::vector<int> ipiv(dim);
  int lwork = 10 * dim;
  std::vector<double> work(lwork);
  int info = 0;
  int n = dim;
  int m = dim;

  dsytrf(uplo, &m, a, &n, &(ipiv[0]), &(work[0]), &lwork, &info);
  if (info) dserror("dsytrf returned info=%d", info);

  dsytri(uplo, &m, a, &n, &(ipiv[0]), &(work[0]), &info);
  if (info) dserror("dsytri returned info=%d", info);

  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < i; ++j) A(j, i) = A(i, j);
  return;
}

/*----------------------------------------------------------------------*
 |  Solve soe with me*ae^T = de^T and return me^-1           farah 07/14|
 *----------------------------------------------------------------------*/
LINALG::SerialDenseMatrix LINALG::InvertAndMultiplyByCholesky(
    LINALG::SerialDenseMatrix& me, LINALG::SerialDenseMatrix& de, LINALG::SerialDenseMatrix& ae)
{
  const int n = me.N();

  LINALG::SerialDenseMatrix y(n, n);
  LINALG::SerialDenseMatrix y_identity(n, n, true);
  LINALG::SerialDenseMatrix meinv(n, n, true);

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
      if (col == z) y_identity(z, col) = yfac;
      for (int u = 0; u < z; ++u)
      {
        y(z, col) -= yfac * me(z, u) * y(u, col);
        y_identity(z, col) -= yfac * me(z, u) * y_identity(u, col);
      }
    }
  }

  // get y for G^T*x=y
  for (int z = n - 1; z > -1; --z)
  {
    const double xfac = 1.0 / me(z, z);
    for (int col = 0; col < n; ++col)
    {
      ae(col, z) = xfac * y(z, col);
      meinv(col, z) = xfac * y_identity(z, col);
      for (int u = n - 1; u > z; --u)
      {
        ae(col, z) -= xfac * me(u, z) * ae(col, u);
        meinv(col, z) -= xfac * me(u, z) * meinv(col, u);
      }
    }
  }

  return meinv;
}
