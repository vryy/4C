/*----------------------------------------------------------------------*/
/*! \file

\brief Inverse dense matrices up to 4x4 and other inverse methods.

\level 0
*/
/*----------------------------------------------------------------------*/

#include "4C_linalg_utils_densematrix_inverse.hpp"

#include "4C_linalg_fortran_definitions.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  invert a dense symmetric matrix         )                mwgee 12/06|
 *----------------------------------------------------------------------*/
void Core::LinAlg::SymmetricInverse(Core::LinAlg::SerialDenseMatrix& A, const int dim)
{
  if (A.numRows() != A.numCols()) FOUR_C_THROW("Matrix is not square");
  if (A.numRows() != dim) FOUR_C_THROW("Dimension supplied does not match matrix");

  double* a = A.values();
  char uplo[5];
  strcpy(uplo, "L ");
  std::vector<int> ipiv(dim);
  int lwork = 10 * dim;
  std::vector<double> work(lwork);
  int info = 0;
  int n = dim;
  int m = dim;

  dsytrf(uplo, &m, a, &n, ipiv.data(), work.data(), &lwork, &info);
  if (info) FOUR_C_THROW("dsytrf returned info=%d", info);

  dsytri(uplo, &m, a, &n, ipiv.data(), work.data(), &info);
  if (info) FOUR_C_THROW("dsytri returned info=%d", info);

  for (int i = 0; i < dim; ++i)
    for (int j = 0; j < i; ++j) A(j, i) = A(i, j);
  return;
}

/*----------------------------------------------------------------------*
 |  Solve soe with me*ae^T = de^T and return me^-1           farah 07/14|
 *----------------------------------------------------------------------*/
Core::LinAlg::SerialDenseMatrix Core::LinAlg::InvertAndMultiplyByCholesky(
    Core::LinAlg::SerialDenseMatrix& me, Core::LinAlg::SerialDenseMatrix& de,
    Core::LinAlg::SerialDenseMatrix& ae)
{
  const int n = me.numCols();

  Core::LinAlg::SerialDenseMatrix y(n, n);
  Core::LinAlg::SerialDenseMatrix y_identity(n, n, true);
  Core::LinAlg::SerialDenseMatrix meinv(n, n, true);

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
        FOUR_C_THROW("matrix is not positive definite!");
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

FOUR_C_NAMESPACE_CLOSE
