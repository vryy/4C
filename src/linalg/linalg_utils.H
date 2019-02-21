/*!----------------------------------------------------------------------
\file linalg_utils.H

\brief A collection of helper methods for namespace LINALG

<pre>
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*-----------------------------------------------------------------------*/
#ifndef LINALG_UTILS_H
#define LINALG_UTILS_H

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_LAPACK.h"
#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_exporter.H"
#include "linalg_blocksparsematrix.H"

// forward declarations
class Epetra_IntVector;

namespace IO
{
  class DiscretizationWriter;
}

namespace LINALG
{
  // forward declaration
  // class Solver;

  /*!
   \brief Create a new Epetra_CrsMatrix and return RefcountPtr to it

   \param rowmap (in): row map of matrix
   \param npr (in): estimated number of entries per row.
   (need not be exact, better should be too big rather then too small)
   */
  Teuchos::RCP<Epetra_CrsMatrix> CreateMatrix(const Epetra_Map& rowmap, const int npr);

  /*!
   \brief Create a new Epetra_Vector and return RefcountPtr to it

   \param rowmap (in): row map of vector
   \param init (in): initializa vector to zero upon construction
   */
  Teuchos::RCP<Epetra_Vector> CreateVector(const Epetra_BlockMap& rowmap, const bool init = true);

  /*!
   \brief Create a new Epetra_MultiVector and return RefcountPtr to it

   \param rowmap (in): row map of vector
   \param rowmap (in): number of vectors
   \param init (in): initializa vector to zero upon construction
   */
  Teuchos::RCP<Epetra_MultiVector> CreateMultiVector(
      const Epetra_BlockMap& rowmap, const int numrows, const bool init = true);

  /*!
   \brief Export a vector to a different map

   Values of source are copied to target where maps don't have to match.
   Prerequisite: Either the map of source OR the map of target has to be unique
   (will be tested)
   \warning When source is overlapping (and therefore target is unique), values
   in the overlapping region are inserted into the target on a first come
   first serve basis, meaning they should be equal in the source to
   be deterministic
   \param source (in) : source vector values are taken from
   \param target (out): target vector values will be inserted in
   */
  void Export(const Epetra_MultiVector& source, Epetra_MultiVector& target);

  /*!
   \brief Export a vector to a different map

   Values of source are copied to target where maps don't have to match.
   Prerequisite: Either the map of source OR the map of target has to be unique
   (will be tested)
   \warning When source is overlapping (and therefore target is unique), values
   in the overlapping region are inserted into the target on a first come
   first serve basis, meaning they should be equal in the source to
   be deterministic
   \param source (in) : source vector values are taken from
   \param target (out): target vector values will be inserted in
   */
  void Export(const Epetra_IntVector& source, Epetra_IntVector& target);

  /*
   \brief Insert a list of elements in a given global row of the matrix with
   out-of-bound and error check

   \param A (in):           Considered matrix
   \param GlobalRow(in):    Row number (in global coordinates) to put elements.
   \param NumEntries(in):   Number of entries.
   \param Values(in):       Values to enter.
   \param Indices(in):      Global column indices corresponding to values.
   */
  void InsertGlobalValues(Teuchos::RCP<Epetra_CrsMatrix> mat, int GlobalRow, int NumEntries,
      double* Values, int* Indices);

  /*!
   \brief Make Epetra_SerialDenseMatrix symmetric by averaging upper and lower traingular indices

   \param A (in/out): Matrix to be symmetrised
   */
  void SymmetriseMatrix(Epetra_SerialDenseMatrix& A);

  /*!
   \brief Compute all eigenvalues of a real symmetric matrix A

   \param A (in):        Matrix to be analysed
   \param L (out):       Vector of eigenvalues in ascending order
   \param postproc (in): flag indicating whether we are using this
   routine for postprocessing only (in that
   case dserror is replaced with a warning)
   */
  void SymmetricEigenValues(
      Epetra_SerialDenseMatrix& A, Epetra_SerialDenseVector& L, const bool postproc = false);

  /*!
   \brief Compute all eigenvalues and eigenvectors of a real symmetric matrix A

   \param A (in/out):    in: Matrix to be analysed, out: eigenvectors
   (i.e. original matrix is destroyed!!!)
   \param L (out):       Vector of eigenvalues in ascending order
   \param postproc (in): flag indicating whether we are using this
   routine for postprocessing only (in that
   case dserror is replaced with a warning)
   */
  void SymmetricEigenProblem(
      Epetra_SerialDenseMatrix& A, Epetra_SerialDenseVector& L, const bool postproc = false);
#if 0
void SymmetricEigenProblem(IO::DiscretizationWriter& output,
    LINALG::Solver& solver,
    Teuchos::RCP<LINALG::SparseMatrix> A,
    Teuchos::RCP<LINALG::SparseMatrix> M,
    Epetra_Vector& invtoggle,
    int nev,
    int step=0,
    double time=0.0);
#endif
  /*!
   \brief Compute all eigenvalues and, optionally, eigenvectors
   of a real symmetric matrix A

   \param A (in/out):    Matrix to be analysed, if eigv=true stores eigenvectors
   \param L (in/out):    Vector of eigenvalues in ascending order
   \param eigv (in):     flag to evaluate also eigenvectors ('N'=no, 'V'=yes)
   \param postproc (in): flag indicating whether we are using this
   routine for postprocessing only (in that
   case dserror is replaced with a warning)
   */
  void SymmetricEigen(Epetra_SerialDenseMatrix& A, Epetra_SerialDenseVector& L, const char eigv,
      const bool postproc = false);

  /*!
   \brief Compute all eigenvalues the generalized Eigenvalue problem
   Ax = lambda Bx via QZ-algorithm (B is singular) and returns the
   maximum Eigenvalue.

   \param A (in):    A Matrix
   \param B (in):    B Matrix

   */
  double GeneralizedEigen(Epetra_SerialDenseMatrix& A, Epetra_SerialDenseMatrix& B);

  /*!
   \brief Compute all eigenvalues and eigenvectors of a real symmetric matrix A

   A = V * S * VT

   \param A (in):        M-by-M matrix to be decomposed
   \param S (out):       M-by-N matrix which is zero except for its diagonal entries holding the
   eigenvalues \param V (out):       M-by-M orthonormal matrix of eigenvectors
   */
  template <unsigned int dim>
  void SYEV(LINALG::Matrix<dim, dim>& A, LINALG::Matrix<dim, dim>& S, LINALG::Matrix<dim, dim>& V)
  {
    const char jobz = 'V';               // Compute eigenvalues and eigenvectors.
    const char uplo = 'U';               // Upper triangle of A is stored;
    const int N = dim;                   // The order of the matrix A.  N >= 0.
    Matrix<dim, dim> tmp(A.A(), false);  // copy, because content of matrix is destroyed
    const int lda = dim;                 // The leading dimension of the array A.  LDA >=max(1,N).
    std::vector<double> w(dim);
    const int lwork = 2 * dim * dim + 6 * dim + 1;
    std::vector<double> work(lwork);
    const int liwork = 3 + 5 * dim;
    std::vector<int> iwork(liwork);
    int info;

    Epetra_LAPACK lapack;
    lapack.SYEVD(
        jobz, uplo, N, tmp.A(), lda, &(w[0]), &(work[0]), lwork, &(iwork[0]), liwork, &info);

    if (info) dserror("Lapack's SYEVD returned %d", info);

    // return eigenvectors
    V.Update(tmp);

    // return eigenvalues
    S.Clear();
    for (unsigned int i = 0; i < dim; ++i) S(i, i) = w[i];

    return;
  }

  /*!
   \brief Compute singular value decomposition (SVD) of a real M-by-N matrix A
   A = U * SIGMA * transpose(V)

   \param A (in/out):    Matrix to be decomposed
   \param U (in/out):    M-by-M orthogonal matrix
   \param SIGMA (in/out):M-by-N matrix which is zero except for its min(m,n) diagonal elements
   \param Vt (in/out):   V is a N-by-N orthogonal matrix, actually returned is V^T
   */
  void SVD(const Epetra_SerialDenseMatrix& A, LINALG::SerialDenseMatrix& U,
      LINALG::SerialDenseMatrix& SIGMA, LINALG::SerialDenseMatrix& Vt);

  /*!
   \brief Singular value decomposition (SVD) of a real M-by-N matrix in fixed
   size format

   A = Q * S * VT

   \param A (in):        M-by-N matrix to be decomposed
   \param Q (out):       M-by-M orthogonal matrix
   \param S (out):       M-by-N matrix which is zero except for its min(m,n) diagonal elements
   \param VT (out):      N-by-N orthogonal matrix (transpose of V)
   */
  template <unsigned int rows, unsigned int cols>
  void SVD(LINALG::Matrix<rows, cols>& A, LINALG::Matrix<rows, rows>& Q,
      LINALG::Matrix<rows, cols>& S, LINALG::Matrix<cols, cols>& VT)
  {
    Matrix<rows, cols> tmp(A.A(), false);  // copy, because content of matrix is destroyed
    Epetra_LAPACK lapack;
    const char jobu = 'A';   // compute and return all M columns of U
    const char jobvt = 'A';  // compute and return all N rows of V^T
    std::vector<double> s(std::min(rows, cols));
    int info;
    int lwork = std::max(3 * std::min(rows, cols) + std::max(rows, cols), 5 * std::min(rows, cols));
    std::vector<double> work(lwork);

    lapack.GESVD(jobu, jobvt, rows, cols, tmp.A(), tmp.M(), &s[0], Q.A(), Q.M(), VT.A(), VT.M(),
        &work[0], &lwork, &info);

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

  /*!
   \brief Explicit inverse and determinant of nonsymmetric 3x3 matrix

   \param A (in/out): Matrix to be inverted
   \param dim (in) :  Dimension of matrix
   */
  double NonsymInverse3x3(Epetra_SerialDenseMatrix& A);

  /*!
   \brief Determinant of a nonsymmetric matrix using SVD

   \return the determinant of the matrix A if nonzero
   */
  double DeterminantSVD(const Epetra_SerialDenseMatrix& A);

  /*!
   \brief Determinant of a nonsymmetric matrix using LU

   \return the determinant of the matrix A
   */
  double DeterminantLU(const Epetra_SerialDenseMatrix& A);

  /*!
   \brief Explicit determinant of a nonsymmetric 3x3 matrix

   \return the determinant of the matrix A
   */
  template <typename T>
  T Determinant3x3(LINALG::TMatrix<T, 3, 3>& A)
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

  template <typename T>
  T Determinant2x2(LINALG::TMatrix<T, 2, 2>& A)
  {
    T b00 = A(0, 0);
    T b01 = A(0, 1);
    T b10 = A(1, 0);
    T b11 = A(1, 1);
    T det = b00 * b11 - b01 * b10;
    return det;
  }

  /*!
   * \brief Explicit inverse of a nonsymmetric 4x4 matrix
   *
   * If the determinant is zero and the flag fail_on_zero_determinant is false, the function will
   * return without an error, but the values in matrix A may not be used in further calculations.
   * The user has to make sure this case is handled correctly.
   *
   * @param A (in/out) Matrix to be inverted. The inverse will be stored in this variable.
   * @param fail_on_zero_det (flag) If true, a dserror will be raised if the determinant of A is
   * zero.
   * @return the determinant of the matrix A
   */
  template <typename T>
  T Inverse4x4(LINALG::TMatrix<T, 4, 4>& A, bool fail_on_zero_determinant = true)
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

    A(0, 0) = a11 * a33 * a22 - a11 * a32 * a23 - a31 * a13 * a22 + a32 * a21 * a13 +
              a31 * a12 * a23 - a33 * a21 * a12;
    A(0, 1) = -(a31 * a02 * a23 + a33 * a22 * a01 - a33 * a21 * a02 - a31 * a03 * a22 -
                a32 * a23 * a01 + a32 * a21 * a03);
    A(0, 2) = -a11 * a33 * a02 + a11 * a32 * a03 + a33 * a01 * a12 + a31 * a13 * a02 -
              a32 * a01 * a13 - a31 * a12 * a03;
    A(0, 3) = -(-a11 * a02 * a23 + a11 * a03 * a22 - a03 * a21 * a12 - a01 * a13 * a22 +
                a02 * a21 * a13 + a01 * a12 * a23);

    A(1, 0) = -(a33 * a22 * a10 - a32 * a23 * a10 + a30 * a12 * a23 - a33 * a20 * a12 -
                a30 * a13 * a22 + a32 * a20 * a13);
    A(1, 1) = a00 * a33 * a22 - a00 * a32 * a23 + a30 * a02 * a23 - a33 * a20 * a02 +
              a32 * a20 * a03 - a30 * a03 * a22;
    A(1, 2) = -(a00 * a33 * a12 - a00 * a32 * a13 - a33 * a10 * a02 - a30 * a03 * a12 +
                a32 * a10 * a03 + a30 * a02 * a13);
    A(1, 3) = -a00 * a13 * a22 + a00 * a12 * a23 - a10 * a02 * a23 + a13 * a20 * a02 +
              a10 * a03 * a22 - a12 * a20 * a03;

    A(2, 0) = a11 * a30 * a23 - a11 * a33 * a20 - a31 * a10 * a23 + a33 * a21 * a10 +
              a31 * a13 * a20 - a30 * a21 * a13;
    A(2, 1) = -(-a00 * a31 * a23 + a00 * a33 * a21 - a30 * a03 * a21 + a31 * a20 * a03 +
                a30 * a01 * a23 - a33 * a20 * a01);
    A(2, 2) = a33 * a00 * a11 - a33 * a10 * a01 - a30 * a03 * a11 - a31 * a00 * a13 +
              a31 * a10 * a03 + a30 * a01 * a13;
    A(2, 3) = -(a23 * a00 * a11 - a23 * a10 * a01 - a20 * a03 * a11 - a21 * a00 * a13 +
                a21 * a10 * a03 + a20 * a01 * a13);

    A(3, 0) = -(a11 * a30 * a22 - a11 * a32 * a20 - a30 * a21 * a12 - a31 * a10 * a22 +
                a32 * a21 * a10 + a31 * a12 * a20);
    A(3, 1) = -a00 * a31 * a22 + a00 * a32 * a21 + a31 * a20 * a02 + a30 * a01 * a22 -
              a32 * a20 * a01 - a30 * a02 * a21;
    A(3, 2) = -(a32 * a00 * a11 - a32 * a10 * a01 - a30 * a02 * a11 - a31 * a00 * a12 +
                a31 * a10 * a02 + a30 * a01 * a12);
    A(3, 3) = a22 * a00 * a11 - a22 * a10 * a01 - a20 * a02 * a11 - a21 * a00 * a12 +
              a21 * a10 * a02 + a20 * a01 * a12;

    T det = a00 * a11 * a33 * a22 - a00 * a11 * a32 * a23 - a00 * a31 * a13 * a22 +
            a00 * a32 * a21 * a13 + a00 * a31 * a12 * a23 - a00 * a33 * a21 * a12 +
            a11 * a30 * a02 * a23 - a11 * a33 * a20 * a02 + a11 * a32 * a20 * a03 -
            a11 * a30 * a03 * a22 - a31 * a10 * a02 * a23 - a33 * a22 * a10 * a01 +
            a33 * a21 * a10 * a02 + a33 * a20 * a01 * a12 + a30 * a03 * a21 * a12 +
            a31 * a13 * a20 * a02 + a31 * a10 * a03 * a22 + a30 * a01 * a13 * a22 +
            a32 * a23 * a10 * a01 - a32 * a21 * a10 * a03 - a32 * a20 * a01 * a13 -
            a30 * a02 * a21 * a13 - a31 * a12 * a20 * a03 - a30 * a01 * a12 * a23;

    if (det == 0.0)
    {
      if (fail_on_zero_determinant) dserror("Determinant of 4x4 matrix is exactly zero");
    }
    else
    {
      // Only scale the matrix if the determinant is not 0.
      A.Scale(1. / det);
    }
    return det;
  }  // namespace LINALG

  /*!
   * \brief Explicit inverse of a nonsymmetric 3x3 matrix
   *
   * If the determinant is zero and the flag fail_on_zero_determinant is false, the function will
   * return without an error, but the values in matrix A may not be used in further calculations.
   * The user has to make sure this case is handled correctly.
   *
   * @param A (in/out) Matrix to be inverted. The inverse will be stored in this variable.
   * @param fail_on_zero_det (flag) If true, a dserror will be raised if the determinant of A is
   * zero.
   * @return the determinant of the matrix A
   */
  template <typename T>
  T Inverse3x3(LINALG::TMatrix<T, 3, 3>& A, bool fail_on_zero_determinant = true)
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
    A(0, 0) = b11 * b22 - b21 * b12;
    A(1, 0) = -b10 * b22 + b20 * b12;
    A(2, 0) = b10 * b21 - b20 * b11;
    A(0, 1) = -b01 * b22 + b21 * b02;
    A(1, 1) = b00 * b22 - b20 * b02;
    A(2, 1) = -b00 * b21 + b20 * b01;
    A(0, 2) = b01 * b12 - b11 * b02;
    A(1, 2) = -b00 * b12 + b10 * b02;
    A(2, 2) = b00 * b11 - b10 * b01;
    T det = b00 * A(0, 0) + b01 * A(1, 0) + b02 * A(2, 0);
    if (det == 0.0)
    {
      if (fail_on_zero_determinant) dserror("Determinant of 3x3 matrix is exactly zero");
    }
    else
    {
      // Only scale the matrix if the determinant is not 0.
      A.Scale(1. / det);
    }
    return det;
  }

  /*!
   * \brief Explicit inverse of a nonsymmetric 2x2 matrix
   *
   * If the determinant is zero and the flag fail_on_zero_determinant is true, the function will
   * return without an error, but the values in matrix A may not be used in further calculations.
   * The user has to make sure this case is handled correctly.
   *
   * @param A (in/out) Matrix to be inverted. The inverse will be stored in this variable.
   * @param fail_on_zero_det (flag) If true, a dserror will be raised if the determinant of A is
   * zero.
   * @return the determinant of the matrix A
   */
  template <typename T>
  T Inverse2x2(LINALG::TMatrix<T, 2, 2>& A, bool fail_on_zero_determinant = true)
  {
    T b00 = A(0, 0);
    T b01 = A(0, 1);
    T b10 = A(1, 0);
    T b11 = A(1, 1);
    A(0, 0) = b11;
    A(1, 0) = -b10;
    A(0, 1) = -b01;
    A(1, 1) = b00;
    T det = b00 * b11 - b01 * b10;
    if (det == 0.0)
    {
      if (fail_on_zero_determinant) dserror("Determinant of 2x2 matrix is exactly zero");
    }
    else
    {
      // Only scale the matrix if the determinant is not 0.
      A.Scale(1. / det);
    }
    return det;
  }

  /*!
   \brief Solve soe with me*ae^T = de^T (employed for calculating coefficient matrix for dual shape
   functions) farah 07/14

   \param me (in):  me (not transposed!)
   \paramdeA (in):  de (not transposed!)
   \param ae(out):  coef. matrix (not transposed!)
   */
  LINALG::SerialDenseMatrix InvertAndMultiplyByCholesky(
      LINALG::SerialDenseMatrix& me, LINALG::SerialDenseMatrix& de, LINALG::SerialDenseMatrix& ae);

  /*!
   \brief Solve soe with me*ae^T = de^T (employed for calculating coefficient matrix for dual shape
   functions) farah 07/14

   \param me (in):  me (not transposed!)
   \paramdeA (in):  de (not transposed!)
   \param ae(out):  coef. matrix (not transposed!)
   */
  template <const int n>
  void InvertAndMultiplyByCholesky(
      LINALG::Matrix<n, n>& me, LINALG::Matrix<n, n>& de, LINALG::SerialDenseMatrix& ae)
  {
    LINALG::Matrix<n, n> y;
    LINALG::Matrix<n, n> maux = me;

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
        for (int u = 0; u < z; ++u) y(z, col) -= yfac * me(z, u) * y(u, col);
      }
    }

    // get y for G^T*x=y
    for (int z = n - 1; z > -1; --z)
    {
      const double xfac = 1.0 / me(z, z);
      for (int col = 0; col < n; ++col)
      {
        ae(col, z) = xfac * y(z, col);
        for (int u = n - 1; u > z; --u) ae(col, z) -= xfac * me(u, z) * ae(col, u);
      }
    }

    me = maux;
    return;
  }

  /*!
   \brief Solve soe with me*ae^T = de^T (employed for calculating coefficient matrix for dual shape
   functions) farah 07/14

   \param me (in):  me (not transposed!)
   \paramdeA (in):  de (not transposed!)
   \param ae(out):  coef. matrix (not transposed!)
   */
  template <const int n>
  void InvertAndMultiplyByCholesky(
      LINALG::Matrix<n, n>& me, LINALG::Matrix<n, n>& de, LINALG::Matrix<n, n>& ae)
  {
    LINALG::Matrix<n, n> y;

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

      // get y for G*y=de^T
      const double yfac = 1.0 / me(z, z);
      for (int col = 0; col < n; ++col)
      {
        y(z, col) = yfac * de(col, z);
        for (int u = 0; u < z; ++u) y(z, col) -= yfac * me(z, u) * y(u, col);
      }
    }

    // get x for G^T*x=y
    for (int z = n - 1; z > -1; --z)
    {
      const double xfac = 1.0 / me(z, z);
      for (int col = 0; col < n; ++col)
      {
        ae(col, z) = xfac * y(z, col);
        for (int u = n - 1; u > z; --u) ae(col, z) -= xfac * me(u, z) * ae(col, u);
      }
    }

    return;
  }

  /*!
   \brief Invert a symmetric dim*dim square matrix (positive definite)
   farah 07/14
   \param A (in/out): Matrix to be inverted
   */
  template <const int n>
  void SymmetricPositiveDefiniteInverse(LINALG::Matrix<n, n>& A)
  {
    LINALG::Matrix<n, n> y(true);
    LINALG::Matrix<n, n> ae(true);

    // calc G with me=G*G^T
    for (int z = 0; z < n; ++z)
    {
      for (int u = 0; u < z + 1; ++u)
      {
        double sum = A(z, u);
        for (int k = 0; k < u; ++k) sum -= A(z, k) * A(u, k);

        if (z > u)
          A(z, u) = sum / A(u, u);
        else if (sum > 0.0)
          A(z, z) = sqrt(sum);
        else
          dserror("matrix is not positive definite!");
      }

      // get y for G*y=de^T
      const double yfac = 1.0 / A(z, z);
      for (int col = 0; col < n; ++col)
      {
        if (col == z) y(z, col) = yfac;
        for (int u = 0; u < z; ++u) y(z, col) -= yfac * A(z, u) * y(u, col);
      }
    }

    // get x for G^T*x=y
    for (int z = n - 1; z > -1; --z)
    {
      const double xfac = 1.0 / A(z, z);
      for (int col = 0; col < n; ++col)
      {
        ae(z, col) = xfac * y(z, col);
        for (int u = n - 1; u > z; --u) ae(z, col) -= xfac * A(u, z) * ae(u, col);
      }
    }

    // get result
    A = ae;

    return;
  }

  /*!
   \brief Invert a symmetric dim*dim square matrix

   \param A (in/out): Matrix to be inverted
   \param dim (in) :  Dimension of matrix
   */
  void SymmetricInverse(Epetra_SerialDenseMatrix& A, const int dim);

  /*!
   \brief Invert a nonsymmetric dim*dim square matrix

   \param A (in/out): Matrix to be inverted
   \param dim (in) :  Dimension of matrix
   */
  void NonSymmetricInverse(Epetra_SerialDenseMatrix& A, const int dim);

  /*!
   \brief Assemble an Epetra_SerialDenseMatrix into an Epetra_CrsMatrix

   This is an individual call.
   Will only assemble locally and will never do any commmunication.
   All values that can not be assembled locally will be ignored.
   Will use the communicator and rowmap from matrix A to determine ownerships.
   Local matrix Aele may be \b square or \b rectangular.

   This version of 'Assemble' does not work for a matrix A that is already
   Filled()! If matrix A is not Filled(), it will be enlarged as required.

   \note The user must provide an \b additional input vector 'lmcol'
   containing the column gids for assembly seperately!

   \param A (out)         : Sparse matrix to be assembled on
   \param Aele (in)       : dense matrix to be assembled
   \param lmrow (in)      : vector with row gids
   \param lmrowowner (in) : vector with owner procs of row gids
   \param lmcol (in)      : vector with column gids
   */
  void Assemble(Epetra_CrsMatrix& A, const Epetra_SerialDenseMatrix& Aele,
      const std::vector<int>& lmrow, const std::vector<int>& lmrowowner,
      const std::vector<int>& lmcol);

  /*!
   \brief Assemble an Epetra_SerialDenseVector into an Epetra_Vector

   This is an individual call.
   Will only assemble locally and will never do any commmunication.
   All values that can not be assembled locally will be ignored.
   Will use the communicator from vector V to determine ownerships.

   \param V (out)   : Vector to be assembled on
   \param Vele (in) : dense vector to be assembled
   \param lm (in) : vector with gids
   \param lmowner (in) : vector with owner procs of gids
   */
  void Assemble(Epetra_Vector& V, const Epetra_SerialDenseVector& Vele, const std::vector<int>& lm,
      const std::vector<int>& lmowner);

  /*!
   \brief Assemble a LINALG::Matrix<3,1> into an Epetra_Vector

   This is a wrapper that converts LINALG::Matrix<3,1> to
   Epetra_SerialDenseVector (view version)

   \param V (out)   : Vector to be assembled on
   \param Vele (in) : dense vector to be assembled
   \param lm (in) : vector with gids
   \param lmowner (in) : vector with owner procs of gids
   */
  void Assemble(Epetra_Vector& V, LINALG::Matrix<3, 1>& Vele, const std::vector<int>& lm,
      const std::vector<int>& lmowner);

  /*!
   \brief Assemble a LINALG::Matrix<3,1> into an Epetra_Vector in case the
   owner is always the same

   This is a wrapper that converts LINALG::Matrix<3,1> to
   Epetra_SerialDenseVector (view version)

   \param V (out)   : Vector to be assembled on
   \param Vele (in) : dense vector to be assembled
   \param lm (in) : vector with gids
   \param lmowner (in) : vector with owner procs of gids
   */
  void Assemble(
      Epetra_Vector& V, LINALG::Matrix<3, 1>& Vele, const std::vector<int>& lm, const int& lmowner);

  /*!
   \brief Assemble a double into an Epetra_Vector

   This is an wrapper.

   \param V (out)   : Vector to be assembled on
   \param Vele (in) : dense vector to be assembled
   \param lm (in) : vector with gids
   \param lmowner (in) : vector with owner procs of gids
   */
  void Assemble(Epetra_Vector& V, double& Vele, const int& lm, const int& lmowner);

  ///////////////////////////
  /*!
   \brief Assemble a Epetra_SerialDenseVector into a Epetra_MultiVector

   This is an individual call.
   Will only assemble locally and will never do any commmunication.
   All values that can not be assembled locally will be ignored.
   Will use the communicator from vector V to determine ownerships.

   \param V (out)   : Vector to be assembled on
   \param n (in)   : column index of MultiVector to be assembled on
   \param Vele (in) : dense vector to be assembled
   \param lm (in) : vector with gids
   \param lmowner (in) : vector with owner procs of gids
   */
  void Assemble(Epetra_MultiVector& V, const int n, const Epetra_SerialDenseVector& Vele,
      const std::vector<int>& lm, const std::vector<int>& lmowner);

  /*! \brief Assemble a source Epetra_Vector into a target Epetra_Vector
   *
   *  The map of the source vector has to be a sub-map of the target vector and
   *  the maps must have the same processor distribution. This method does not
   *  build up any communication between different processors!
   *  The entries of the source vector are added to the target vector:
   *
   *      target_vector[GID] = source_scalar * source[GID] + target_scalar * target[GID]
   *
   *  \note The remaining GIDs, which are no part of the source map stay untouched!
   *
   *  \param scalar_target (in) : scale the target entries by this factor
   *  \param target        (out): target vector (part of the source)
   *  \param scalar_source (in) : scale the source entries by this factor
   *  \param source        (in) : source vector
   *
   *  \author hiermeier \date 03/17 */
  void AssembleMyVector(double scalar_target, Epetra_Vector& target, double scalar_source,
      const Epetra_Vector& source);

  /*! \brief Extract a partial Epetra_Vector from a given source vector
   *         on each proc without communication
   *
   *  This methods uses a given partial map to create the partial target vector.
   *
   *  \param source     (in) : source vector ( read-only )
   *  \param target_map (in) : map of the new target vector ( read-only )
   *
   *  \return the extracted partial Epetra_Vector as RCP
   *
   *  \author hiermeier \date 03/17 */
  Teuchos::RCP<Epetra_Vector> ExtractMyVector(
      const Epetra_Vector& source, const Epetra_Map& target_map);

  /*! \brief Extract a partial Eptra_Vector from a given source vector
   *         on each proc without communication
   *
   *  \param source (in) : source vector ( read-only )
   *  \param target (out): this target vector is going to be filled
   *
   *  \author hiermeier \date 03/17 */
  void ExtractMyVector(const Epetra_Vector& source, Epetra_Vector& target);

  /*! \brief Extract a partial Epetra_Vector from a given source vector
   *         and scale the entries on each proc without communication
   *
   *  \param scalar_source (in) : scale the target entries by this factor
   *  \param source        (in) : source vector ( read-only )
   *  \param scalar_target (in) : scale the source entries by this factor
   *  \param target        (out): this target vector is going to be filled with
   *                              the combined scaled values of the source and
   *                              target vectors
   *
   *  \author hiermeier \date 03/17 */
  void ExtractMyVector(double scalar_source, const Epetra_Vector& source, double scalar_target,
      Epetra_Vector& target);

  /*!
   \brief Call FillComplete on a Epetra_CrsMatrix (for square matrices only!)
   */
  void Complete(Epetra_CrsMatrix& A);

  /*!
   \brief Call FillComplete on a Epetra_CrsMatrix (for square matrices only!)

   This is the Teuchos::RCP wrapper of the above method.
   */
  inline void Complete(Teuchos::RCP<Epetra_CrsMatrix> A)
  {
    LINALG::Complete(*A);
    return;
  }

  /*!
   \brief Call FillComplete on a Epetra_CrsMatrix (for rectangular and square matrices)
   */
  void Complete(Epetra_CrsMatrix& A, const Epetra_Map& domainmap, const Epetra_Map& rangemap);

  /*!
   \brief Call FillComplete on a Epetra_CrsMatrix (for rectangular and square matrices)

   This is the Teuchos::RCP wrapper of the above method.
   */
  inline void Complete(
      Teuchos::RCP<Epetra_CrsMatrix> A, const Epetra_Map& domainmap, const Epetra_Map& rangemap)
  {
    LINALG::Complete(*A, domainmap, rangemap);
    return;
  }

  /*!
   \brief Add a (transposed) Epetra_CrsMatrix to another: B = B*scalarB + A(^T)*scalarA

   Add one matrix to another.

   The matrix B may or may not be completed. If B is completed, no new elements can be
   inserted and the addition only succeeds in case the sparsity pattern of B is a superset of
   the sparsity pattern of A (otherwise: dserror).

   Performance characterization: If B is filled (completed), this function is pretty fast,
   typically on the order of two to four matrix-vector products with B. The case where B is
   un-filled runs much slower (on the order of up to 100 matrix-vector products).

   Sparsity patterns of A and B need not match and A and B can be
   nonsymmetric in value and pattern.

   Row map of A has to be a processor-local subset of the row map of B.

   Note that this is a true parallel add, even in the transposed case!

   \param A          (in)     : Matrix to add to B (must have Filled()==true)
   \param transposeA (in)     : flag indicating whether transposed of A should be used
   \param scalarA    (in)     : scaling factor for A
   \param B          (in/out) : Matrix to be added to (must have Filled()==false)
   \param scalarB    (in)     : scaling factor for B
   */
  void Add(const Epetra_CrsMatrix& A, const bool transposeA, const double scalarA,
      Epetra_CrsMatrix& B, const double scalarB);

  /*!
   \brief Add a (transposed) Epetra_CrsMatrix to another: B = B*scalarB + A(^T)*scalarA

   Add one matrix to another.

   The matrix B may or may not be completed. If B is completed, no new elements can be
   inserted and the addition only succeeds in case the sparsity pattern of B is a superset of
   the sparsity pattern of A (otherwise: dserror).

   Performance characterization: If B is filled (completed), this function is pretty fast,
   typically on the order of two to four matrix-vector products with B. The case where B is
   un-filled runs much slower (on the order of up to 100 matrix-vector products).

   Sparsity patterns of A and B need not match and A and B can be
   nonsymmetric in value and pattern.

   Row map of A has to be a processor-local subset of the row map of B.


   Note that this is a true parallel add, even in the transposed case!
   This is the Teuchos::RCP wrapper of the above method.

   \param A          (in)     : Matrix to add to B (must have Filled()==true)
   \param transposeA (in)     : flag indicating whether transposed of A should be used
   \param scalarA    (in)     : scaling factor for A
   \param B          (in/out) : Matrix to be added to (must have Filled()==false)
   \param scalarB    (in)     : scaling factor for B
   */
  inline void Add(const Teuchos::RCP<Epetra_CrsMatrix> A, const bool transposeA,
      const double scalarA, Teuchos::RCP<Epetra_CrsMatrix> B, const double scalarB)
  {
    LINALG::Add(*A, transposeA, scalarA, *B, scalarB);
    return;
  }

  /*!
   \brief Add a (transposed) Epetra_CrsMatrix to a LINALG::SparseMatrix: B = B*scalarB +
   A(^T)*scalarA

   Add one matrix to another.

   As opposed to the other Add() functions, this method can handle both the case where
   matrix B is fill-completed (for performance reasons) but does not have to.
   If B is completed and new matrix elements are detected, the matrix is un-completed and
   rebuild internally (expensive).

   Sparsity patterns of A and B need not match and A and B can be
   nonsymmetric in value and pattern.

   Row map of A has to be a processor-local subset of the row map of B.

   Note that this is a true parallel add, even in the transposed case!

   \param A          (in)     : Matrix to add to B (must have Filled()==true)
   \param transposeA (in)     : flag indicating whether transposed of A should be used
   \param scalarA    (in)     : scaling factor for A
   \param B          (in/out) : Matrix to be added to (must have Filled()==false)
   \param scalarB    (in)     : scaling factor for B
   */
  void Add(const Epetra_CrsMatrix& A, const bool transposeA, const double scalarA,
      LINALG::SparseMatrixBase& B, const double scalarB);

  /*!
   \brief Compute transposed matrix of an Epetra_CrsMatrix explicitly

   Returns Teuchos::RCP to the transposed matrix of the input matrix A.

   \param A          (in)     : Matrix to transpose (must have Filled()==true)

   */
  Teuchos::RCP<Epetra_CrsMatrix> Transpose(const Epetra_CrsMatrix& A);

  /*!
   \brief Compute transposed matrix of an Epetra_CrsMatrix explicitly

   Returns Teuchos::RCP to the transposed matrix of the input matrix A.
   This is the Teuchos::RCP wrapper of the above method.

   \param A          (in)     : Matrix to transpose (must have Filled()==true)

   */
  inline Teuchos::RCP<Epetra_CrsMatrix> Transpose(const Teuchos::RCP<Epetra_CrsMatrix> A)
  {
    return LINALG::Transpose(*A);
  }

  /*!
   \brief Multiply a (transposed) Epetra_CrsMatrix with another (transposed): C = A(^T)*B(^T)

   Multiply one matrix with another. Both matrices must be completed. Sparsity
   Respective Range, Row and Domain maps of A(^T) and B(^T) have to match.

   Note that this is a true parallel multiplication, even in the transposed case!

   \param A          (in)     : Matrix to multiply with B (must have Filled()==true)
   \param transA     (in)     : flag indicating whether transposed of A should be used
   \param B          (in)     : Matrix to multiply with A (must have Filled()==true)
   \param transB     (in)     : flag indicating whether transposed of B should be used
   \param complete   (in)     : flag indicating whether FillComplete should be called on C upon
   exit, (defaults to true) \return Matrix product A(^T)*B(^T)
   */
  Teuchos::RCP<Epetra_CrsMatrix> Multiply(const Epetra_CrsMatrix& A, bool transA,
      const Epetra_CrsMatrix& B, bool transB, bool complete = true);

  /*!
   \brief Multiply a (transposed) Epetra_CrsMatrix with another (transposed): C = A(^T)*B(^T)

   Multiply one matrix with another. Both matrices must be completed. Sparsity
   Respective Range, Row and Domain maps of A(^T) and B(^T) have to match.

   Note that this is a true parallel multiplication, even in the transposed case!
   This is the Teuchos::RCP wrapper of the above method.

   \param A          (in)     : Matrix to multiply with B (must have Filled()==true)
   \param transA     (in)     : flag indicating whether transposed of A should be used
   \param B          (in)     : Matrix to multiply with A (must have Filled()==true)
   \param transB     (in)     : flag indicating whether transposed of B should be used
   \param complete   (in)     : flag indicating whether FillComplete should be called on C upon
   exit, (defaults to true) \return Matrix product A(^T)*B(^T)
   */
  inline Teuchos::RCP<Epetra_CrsMatrix> Multiply(const Teuchos::RCP<Epetra_CrsMatrix>& A,
      bool transA, const Teuchos::RCP<Epetra_CrsMatrix>& B, bool transB, bool complete = true)
  {
    return Multiply(*A, transA, *B, transB, complete);
  }

  /*!
   \brief Triple matrix product: D = A(^T)*B(^T)*C(^T)

   Multiply one matrix with another. All input matrices must be completed. Sparsity
   Respective Range, Row and Domain maps of A(^T) and B(^T) C(^T) have to match.

   Note that this is a true parallel multiplication, even in the transposed case!

   \param A          (in)     : Matrix to multiply with B (must have Filled()==true)
   \param transA     (in)     : flag indicating whether transposed of A should be used
   \param B          (in)     : Matrix to multiply with C (must have Filled()==true)
   \param transB     (in)     : flag indicating whether transposed of B should be used
   \param C          (in)     : Matrix C (must have Filled()==true)
   \param transC     (in)     : flag indicating whether transposed of C should be used
   \param complete   (in)     : flag indicating whether FillComplete should be called on C upon
   exit, (defaults to true) \return Matrix product A(^T)*B(^T)*C(^T)
   */
  Teuchos::RCP<Epetra_CrsMatrix> Multiply(const Epetra_CrsMatrix& A, bool transA,
      const Epetra_CrsMatrix& B, bool transB, const Epetra_CrsMatrix& C, bool transC,
      bool complete = true);

  /*!
   \brief Triple matrix product: D = A(^T)*B(^T)*C(^T)

   Multiply one matrix with another. All input matrices must be completed. Sparsity
   Respective Range, Row and Domain maps of A(^T) and B(^T) C(^T) have to match.

   Note that this is a true parallel multiplication, even in the transposed case!
   This is the Teuchos::RCP wrapper of the above method.

   \param A          (in)     : Matrix to multiply with B (must have Filled()==true)
   \param transA     (in)     : flag indicating whether transposed of A should be used
   \param B          (in)     : Matrix to multiply with C (must have Filled()==true)
   \param transB     (in)     : flag indicating whether transposed of B should be used
   \param C          (in)     : Matrix C (must have Filled()==true)
   \param transC     (in)     : flag indicating whether transposed of C should be used
   \param complete   (in)     : flag indicating whether FillComplete should be called on C upon
   exit, (defaults to true) \return Matrix product A(^T)*B(^T)*C(^T)
   */
  inline Teuchos::RCP<Epetra_CrsMatrix> Multiply(const Teuchos::RCP<Epetra_CrsMatrix>& A,
      bool transA, const Teuchos::RCP<Epetra_CrsMatrix>& B, bool transB,
      const Teuchos::RCP<Epetra_CrsMatrix>& C, bool transC, bool complete = true)
  {
    return Multiply(*A, transA, *B, transB, *C, transC, complete);
  }

  /*!
   \brief Apply dirichlet boundary condition to a linear system of equations

   Modifies a system of equations such that dirichlet boundary conditions are enforced.
   Prescribed dirichlet BC values are supplied in dbcval and dbctoggle, where
   a prescribed value is dbcval[i] and dbctoggle[i] = 1.0. No BC is enforced in
   all places where dbctoggle[i] != 1.0.<br>
   Let us denote the \f$ A_{2 \times 2} \f$ blocks of \f$A\f$ by
   \f$A_{ff}, A_{fD}, A_{Df}, A_{DD}\f$, where \f$f\f$ stands for 'free' and
   \f$D\f$ stands for 'Dirichlet BC'. Then, after a call to this method<br>

   \f$ A_{ff} = A_{ff}, \f$<br>
   \f$ A_{fD} = A_{fD}, \f$<br>
   \f$ A_{Df} = 0_{Df}, \f$<br>
   \f$ A_{DD} = I_{DD}, \f$<br>
   \f$ x_{D} = dbcval_{D}, \f$<br>
   \f$ b_{D} = dbcval_{D} \f$<br>

   and

   \f$ A_{ff} x_f + A_{fD} x_D = b_f \f$<br>
   \f$ 0 x_f + I_{DD} x_D = x_D \f$.<br>

   \note The matrix is then nonsymmetric. When using iterative methods on this
   linear system of equations that depend on the symmetry of the matrix (such as e.g. CG),
   the initial guess supplied to the solver has to be exact at the
   Dirichlet BCs. This should be easy, as the values at the Dirichlet BCs
   are known.

   \note The mask of matrix \f$A\f$ is not modified. That is the
   entries in \f$A_{Df}\f$ and \f$A_{DD}\f$ are set to zero, not
   removed. This way the matrix can be reused in the next step.

   \param A         (in/out) : Matrix of Ax=b
   \param x         (in/out) : initial guess vector x of Ax=b
   \param b         (in/out) : rhs vector b of Ax=b
   \param dbcval    (in)     : vector holding prescribed dirichlet values
   \param dbctoggle (in)     : vector holding 1.0 where dirichlet should be applied
   and 0.0 everywhere else
   */
  void ApplyDirichlettoSystem(Teuchos::RCP<SparseOperator> A, Teuchos::RCP<Epetra_Vector>& x,
      Teuchos::RCP<Epetra_Vector>& b, const Teuchos::RCP<const Epetra_Vector> dbcval,
      const Teuchos::RCP<const Epetra_Vector> dbctoggle);

  /*!
   \brief Apply dirichlet boundary condition to a linear system of equations

   This is a flexible routine. The vectors x and dbcval might have different
   maps. The map does not need to contain all Dirichlet dofs.

   The purpose is to set Dirichlet values at a subset of all Dirichlet
   boundaries.

   \param A (in/out)         : Matrix of Ax=b
   \param x (in/out)         : vector x of Ax=b
   \param b (in/out)         : vector b of Ax=b
   \param dbcval (in)        : vector holding values that are supposed to be prescribed
   \param dbcmap (in)        : unique map of all dofs that should be constrained

   \pre The map dbcmap must be subset of the maps of the vectors.
   */
  void ApplyDirichlettoSystem(Teuchos::RCP<LINALG::SparseOperator> A,
      Teuchos::RCP<Epetra_Vector>& x, Teuchos::RCP<Epetra_Vector>& b,
      const Teuchos::RCP<const Epetra_Vector>& dbcval, const Epetra_Map& dbcmap);

  /*!
   \brief Apply dirichlet boundary condition to a linear system of equations

   This is a flexible routine. The vectors x and dbcval might have different
   maps. The map does not need to contain all Dirichlet dofs.

   The purpose is to set Dirichlet values at a subset of all Dirichlet
   boundaries.

   Special in this routine is the ability to insert rows of general rotation
   matrices (stored in #trafo) rather than simply put ones and zeros
   at the rows associated Dirichlet DOFs.

   \param A (in/out)         : Matrix of Ax=b
   \param x (in/out)         : vector x of Ax=b
   \param b (in/out)         : vector b of Ax=b
   \param trafo (in)         : global matrix holding rotation matrices to convert
   from global to local co-ordinate systems
   \param dbcval (in)        : vector holding values that are supposed to be prescribed
   \param dbcmap (in)        : unique map of all dofs that should be constrained

   \pre The map dbcmap must be subset of the maps of the vectors.
   */
  void ApplyDirichlettoSystem(Teuchos::RCP<LINALG::SparseOperator> A,
      Teuchos::RCP<Epetra_Vector>& x, Teuchos::RCP<Epetra_Vector>& b,
      Teuchos::RCP<const LINALG::SparseMatrix> trafo,
      const Teuchos::RCP<const Epetra_Vector>& dbcval, const Epetra_Map& dbcmap);

  /*!
   \brief Apply dirichlet boundary condition to a linear system of equations

   This is a flexible routine. The vectors b and dbcval might have different
   maps. The map does not need to contain all Dirichlet dofs.

   The purpose is to set Dirichlet values at a subset of all Dirichlet
   boundaries.

   Special in this routine is the ability to insert rows of general rotation
   matrices (stored in #trafo) rather than simply put ones and zeros
   at the rows associated Dirichlet DOFs.

   \param A (in/out)         : Matrix of Ax=b
   \param b (in/out)         : vector b of Ax=b
   \param trafo (in)         : global matrix holding rotation matrices to convert
   from global to local co-ordinate systems
   \param dbcval (in)        : vector holding values that are supposed to be prescribed
   \param dbcmap (in)        : unique map of all dofs that should be constrained

   \pre The map dbcmap must be subset of the maps of the vectors.
   */
  void ApplyDirichlettoSystem(Teuchos::RCP<LINALG::SparseOperator> A,
      Teuchos::RCP<Epetra_Vector>& b, Teuchos::RCP<const LINALG::SparseMatrix> trafo,
      const Teuchos::RCP<const Epetra_Vector>& dbcval, const Epetra_Map& dbcmap);

  /*!
   \brief Apply dirichlet boundary condition to a linear system of equations


   \param x (in/out)         : vector x of Ax=b
   \param b (in/out)         : vector b of Ax=b
   \param dbcval (in)        : vector holding values that are supposed to be prescribed
   \param dbctoggle (in)     : vector holding 1.0 where dirichlet should be applied
   and 0.0 everywhere else
   */
  void ApplyDirichlettoSystem(Teuchos::RCP<Epetra_Vector>& x, Teuchos::RCP<Epetra_Vector>& b,
      const Teuchos::RCP<const Epetra_Vector> dbcval,
      const Teuchos::RCP<const Epetra_Vector> dbctoggle);

  /*!
   \brief Apply dirichlet boundary condition to a linear system of equations

   This is a flexible routine. The vectors x and dbcval might have different
   maps. The dbcmap does not need to contain all Dirichlet dofs, but the vectors
   all dofs defined in it.

   The purpose is to set Dirichlet values at a subset of all Dirichlet
   boundaries.

   \param x (in/out)         : vector x of Ax=b
   \param b (in/out)         : vector b of Ax=b
   \param dbcval (in)        : vector holding values that are supposed to be prescribed
   \param dbcmap (in)        : unique map of all dofs that should be constrained

   \pre The map dbcmap must be subset of the maps of the vectors.
   */
  void ApplyDirichlettoSystem(Teuchos::RCP<Epetra_Vector>& x, Teuchos::RCP<Epetra_Vector>& b,
      Teuchos::RCP<const Epetra_Vector> dbcval, const Epetra_Map& dbcmap);

  /*!
   \brief Apply dirichlet boundary condition to a linear system of equations

   This is a flexible routine. The vectors x and dbcval might have different
   maps. The map does not need to contain all Dirichlet dofs.

   NOTE: Vector b does not need to contain all Dirichlet dofs defined in dbcmap

   The purpose is to set Dirichlet values at a subset of all Dirichlet
   boundaries.

   \param b (in/out)         : vector b of Ax=b
   \param dbcval (in)        : vector holding values that are supposed to be prescribed
   \param dbcmap (in)        : unique map of all dofs that should be constrained

   \pre The map dbcmap must be subset of the maps of the vectors.
   */
  void ApplyDirichlettoSystem(Teuchos::RCP<Epetra_Vector>& b,
      const Teuchos::RCP<const Epetra_Vector> dbcval, const Epetra_Map& dbcmap);

  /*!
   \brief Convert a Dirichlet toggle vector in a Dirichlet map

   The purpose of the routine is a smooth transition from Dirichlet toggle vectors
   to Dirichlet condition maps. Eventually, this method should be removed.

   A Dirichlet toogle vector is a real vector which holds a 1.0 at DOF subjected
   to Dirichlet boundary conditions and a 0.0 at every remaining/free DOF.

   \param dbctoggle (in)     : the Dirichlet toggle vector
   \return MapExtractor object which stores the Dirichlet condition and remaining (other) DOF map

   \author bborn
   \date 10/08
   */
  Teuchos::RCP<LINALG::MapExtractor> ConvertDirichletToggleVectorToMaps(
      const Teuchos::RCP<const Epetra_Vector>& dbctoggle);

  /*!
   \brief split a matrix into a 2x2 block system where the rowmap of one of the blocks is given

   Splits a given matrix into a 2x2 block system where the rowmap of one of the blocks is given
   on input. Blocks A11 and A22 are assumed to be square.
   All values on entry have to be Teuchos::null except the given rowmap and matrix A.
   Note that either A11rowmap or A22rowmap or both have to be nonzero. In case
   both rowmaps are supplied they have to be an exact and nonoverlapping split of A->RowMap().
   Matrix blocks are FillComplete() on exit.

   \param A         : Matrix A on input
   \param A11rowmap : rowmap of A11 or null
   \param A22rowmap : rowmap of A22 or null
   \param A11       : on exit matrix block A11
   \param A12       : on exit matrix block A12
   \param A21       : on exit matrix block A21
   \param A22       : on exit matrix block A22
   */
  bool SplitMatrix2x2(Teuchos::RCP<Epetra_CrsMatrix> A, Teuchos::RCP<Epetra_Map>& A11rowmap,
      Teuchos::RCP<Epetra_Map>& A22rowmap, Teuchos::RCP<Epetra_CrsMatrix>& A11,
      Teuchos::RCP<Epetra_CrsMatrix>& A12, Teuchos::RCP<Epetra_CrsMatrix>& A21,
      Teuchos::RCP<Epetra_CrsMatrix>& A22);

  /*!
   \brief split a matrix into a 2x2 block system

   Splits a given matrix into a 2x2 block system. All values on entry have to be
   Teuchos::null except the given rowmap(s) / domainmap(s) and matrix A.
   Note that either A11rowmap or A22rowmap or both have to be nonzero!
   Note that either A11domainmap or A22domainmap or both have to be nonzero!
   In case both rowmaps / domainmaps are supplied they have to be an exact and
   nonoverlapping split of A->RowMap() / A->DomainMap().
   Matrix blocks are FillComplete() on exit.

   \param A            : Matrix A on input
   \param A11rowmap    : rowmap of A11 or null
   \param A22rowmap    : rowmap of A22 or null
   \param A11domainmap : domainmap of A11 or null
   \param A22domainmap : domainmap of A22 or null
   \param A11          : on exit matrix block A11
   \param A12          : on exit matrix block A12
   \param A21          : on exit matrix block A21
   \param A22          : on exit matrix block A22
   */
  bool SplitMatrix2x2(Teuchos::RCP<LINALG::SparseMatrix> A, Teuchos::RCP<Epetra_Map>& A11rowmap,
      Teuchos::RCP<Epetra_Map>& A22rowmap, Teuchos::RCP<Epetra_Map>& A11domainmap,
      Teuchos::RCP<Epetra_Map>& A22domainmap, Teuchos::RCP<LINALG::SparseMatrix>& A11,
      Teuchos::RCP<LINALG::SparseMatrix>& A12, Teuchos::RCP<LINALG::SparseMatrix>& A21,
      Teuchos::RCP<LINALG::SparseMatrix>& A22);

  /** \brief Insert a diagonal row vector into a unfilled SparseMatrix
   *         on each proc without communication
   *
   *  \param mat (out) : Unfilled matrix
   *  \param diag (in) : Given diagonal (row-layout)
   *
   *  Return 0, if successful. If the given matrix is already filled, the method
   *  returns -1. In this case you should use ReplaceDiagonalValues(), instead.
   *
   *  \author hiermeier \date 03/17 */
  int InsertMyRowDiagonalIntoUnfilledMatrix(LINALG::SparseMatrix& mat, const Epetra_Vector& diag);

  /*!
   \brief split a rowmap of matrix A

   splits A->RowMap() into 2 maps, where one is given on input and the other map
   is created as complementary map. The complementary map is returned.

   \param Amap      : Map to split on input
   \param Agiven    : on entry submap that is given and part of Amap
   \return the remainder map of Amap that is not overlapping with Agiven
   */
  Teuchos::RCP<Epetra_Map> SplitMap(const Epetra_Map& Amap, const Epetra_Map& Agiven);

  /*!
   \brief merges two given Epetra_Maps

   merges input map1 and input map2, both of which have to be unique,
   but may be overlapping, to a new map and returns Teuchos::RCP to it.

   \param map1         : one map to be merged
   \param map2         : the other map to be merged
   \param allowoverlap : when set to false, an error is thrown if the result
   map is overlapping (default = true, overlap allowed)
   \return the (sorted) merged map of input maps map1 and map2
   */
  Teuchos::RCP<Epetra_Map> MergeMap(
      const Epetra_Map& map1, const Epetra_Map& map2, bool overlap = true);

  /*!
   \brief find the intersection set of two given Epetra_Maps

   Find the insection set of input map1 and input map2.

   \param map1         : first map
   \param map2         : second map
   \return the (sorted) intersection map of input maps map1 and map2
   */
  Teuchos::RCP<Epetra_Map> IntersectMap(const Epetra_Map& map1, const Epetra_Map& map2);


  /*!
   \brief merges two given Epetra_Maps

   merges input map1 and input map2 (given as Teuchos::RCP), both of which
   have to be unique, but may be overlapping, to a new map and returns
   Teuchos::RCP to it. The case that one or both input Teuchos::RCPs are null is
   detected and handled appropriately.

   \param map1         : one map to be merged
   \param map2         : the other map to be merged
   \param allowoverlap : when set to false, an error is thrown if the result
   map is overlapping (default = true, overlap allowed)
   \return the (sorted) merged map of input maps map1 and map2
   */
  Teuchos::RCP<Epetra_Map> MergeMap(const Teuchos::RCP<const Epetra_Map>& map1,
      const Teuchos::RCP<const Epetra_Map>& map2, bool overlap = true);

  /*! \brief compute a set of dof maps from a provided set of node maps
   *
   *  \param(in) dofset_id  considered dof set id. In the default case this value
   *                        can be set to zero.
   *  \param(in) node_maps  set of provided node maps
   *  \param(in) discret    underlying discretization
   *  \param(out) dof_maps  set of corresponding dof maps
   *
   *  \author hiermeier \date 03/18 */
  void ComputeDofMapsFromNodeMaps(const int dofset_id,
      const std::vector<Teuchos::RCP<Epetra_Map>>& node_maps,
      const DRT::DiscretizationInterface& discret, std::vector<Teuchos::RCP<Epetra_Map>>& dof_maps);

  /*! \brief Compute a dof map corresponding to the given node map
   *
   *  \param(in) dofset_id  considered dof set id. In the default case this value
   *                        can be set to zero.
   *  \param(in) node_map   provided node map
   *  \param(in) discret    underlying discretization
   *  \return corresponding dof map
   *
   *  \author hiermeier \date 03/18 */
  Teuchos::RCP<Epetra_Map> ComputeDofMapFromNodeMap(
      const int dofset_id, const Epetra_Map& node_map, const DRT::DiscretizationInterface& discret);

  /*! \brief Extract a locally overlapping sub-map from the source map
   *
   *  Find all gids in the source map which are also part of the reference map.
   *  This function works only locally on each proc. No communication is
   *  performed.
   *
   *  \param(in) src_map  source map which shall be considered
   *  \param(in) ref_map  reference map
   *  \return overlapping sub-map
   *
   *  \author hiermeier \date 03/18 */
  Teuchos::RCP<Epetra_Map> ExtractMyOverlappingSubMap(
      const Epetra_BlockMap& src_map, const Epetra_BlockMap& ref_map);

  /// create an Epetra_Map from a set of gids
  /*!
   This is one of the basic operations that is needed every so often.

   \param gids The local gids of this map
   \param comm The map's communicator
   */
  Teuchos::RCP<Epetra_Map> CreateMap(const std::set<int>& gids, const Epetra_Comm& comm);

  /// create an Epetra_Map from a vector of gids
  /*!
   This is one of the basic operations that is needed every so often.

   \param gids The local gids of this map
   \param comm The map's communicator
   */
  Teuchos::RCP<Epetra_Map> CreateMap(const std::vector<int>& gids, const Epetra_Comm& comm);

  /*!
   \brief split a vector into 2 non-overlapping pieces (Teuchos::RCP version)

   \param xmap    : map of vector to be split
   \param x       : vector to be split
   \param x1map   : map of first vector to be extracted
   \param x1      : first vector to be extracted
   \param x2map   : map of second vector to be extracted
   \param x2      : second vector to be extracted

   */
  bool SplitVector(const Epetra_Map& xmap, const Epetra_Vector& x, Teuchos::RCP<Epetra_Map>& x1map,
      Teuchos::RCP<Epetra_Vector>& x1, Teuchos::RCP<Epetra_Map>& x2map,
      Teuchos::RCP<Epetra_Vector>& x2);

  /*!
   \brief split a vector into 2 non-overlapping pieces (Teuchos::RCP version)

   \param xmap    : map of vector to be split
   \param x       : vector to be split
   \param x1map   : map of first vector to be extracted
   \param x1      : first vector to be extracted
   \param x2map   : map of second vector to be extracted
   \param x2      : second vector to be extracted

   */
  bool SplitVector(const Epetra_Map& xmap, const Epetra_Vector& x,
      Teuchos::RCP<const Epetra_Map>& x1map, Teuchos::RCP<Epetra_Vector>& x1,
      Teuchos::RCP<const Epetra_Map>& x2map, Teuchos::RCP<Epetra_Vector>& x2);

  /*!
   \brief Print sparsity pattern of a matrix to postscript file

   creates a file Epetra::Matrix.ps in current directory where the exact
   name of file depends on the exact type of class.

   \note works in parallel and serial!
   */
  void PrintSparsityToPostscript(const Epetra_RowMatrix& A);

  /*!
   \brief Print content of a sparse matrix in Matlab format

   The sparse matrix is printed in Matlab sparse format to
   a file with specified file name.

   \note - Works in parallel and serial!
   \note - When reading this in Matlab, entries are located according to
   their global row and column ID, even if corresponding DOFs are condensed
   from global system.
   */
  void PrintMatrixInMatlabFormat(
      std::string fname, const Epetra_CrsMatrix& A, const bool newfile = true);

  /*!
   \brief Print content of SerialDenseMatrix in Matlab format

   */
  void PrintSerialDenseMatrixInMatlabFormat(
      std::string fname, const Epetra_SerialDenseMatrix& A, const bool newfile = true);

  /*!
   \brief Print content of a block sparse matrix in Matlab format

   The block sparse matrix is printed in Matlab sparse format to
   a file with specified file name.

   \note - Works in parallel and serial!
   \note - When reading this in Matlab, entries are located according to
   their global row and column ID, even if corresponding DOFs are condensed
   from global system.
   */
  void PrintBlockMatrixInMatlabFormat(std::string fname, const BlockSparseMatrixBase& A);

  /*!
   \brief Print content of a vector in Matlab format

   The Vector V is printed in Matlab sparse format to
   a file with specified file name.

   */
  void PrintVectorInMatlabFormat(
      std::string fname, const Epetra_Vector& V, const bool newfile = true);

  /*!
   \brief Print content of a map in Matlab format

   The map is printed in Matlab numbering (+1) to
   a file with specified file name.

   */
  void PrintMapInMatlabFormat(std::string fname, const Epetra_Map& map, const bool newfile = true);

  /*!
   \brief Gather information of type vector<T> on a subset of processors

   This template gathers information provided in sdata on a subset of
   processors tprocs, where the length of the array tprocs is ntargetprocs.
   The redistributed data is returned in rdata which has appropiate size
   on output (size of rdata on input is arbitrary). ntargetprocs can be
   one to reduce data to one proc, it also can be equal to the total number
   of processors to make sdata redundant on all procs.

   \note Functionality of this method is equal to that of Epetra_Comm::GatherAll
   except for that the Epetra version demands the data to be of constant
   size over all procs which this method does not require!

   \param sdata (in) : Information to be gathered on tprocs.
   Length of sdata can be different on every proc.
   \param rdata (out): Information from sdata gathered on a subset of procs.
   size of rdata on input is arbitrary, it is exact on output.
   \param ntargetprocs (in): length of tprocs
   \param tprocs (in): vector of procs ids the information in sdata shall be
   gathered on.
   \param comm (in):   communicator to be used.


   */
  template <typename T>
  void Gather(std::vector<T>& sdata, std::vector<T>& rdata, const int ntargetprocs,
      const int* tprocs, const Epetra_Comm& comm)
  {
    const int myrank = comm.MyPID();
    const int numproc = comm.NumProc();
    if (numproc == 1)
    {
      rdata = sdata;
      return;  // nothing to do in serial
    }
    // build a map of data
    std::map<int, std::vector<T>> datamap;
    datamap[myrank] = sdata;
    // build a source map
    Epetra_Map source(numproc, 1, &myrank, 0, comm);
    // build a target map which is redundant on all target procs and zero everywhere else
    bool iamtarget = false;
    for (int i = 0; i < ntargetprocs; ++i)
      if (tprocs[i] == myrank)
      {
        iamtarget = true;
        break;
      }
    std::vector<int> targetvec(0);
    if (iamtarget)
    {
      targetvec.resize(numproc);
      for (int i = 0; i < numproc; ++i) targetvec[i] = i;
    }
    const int tnummyelements = (int)targetvec.size();
    Epetra_Map target(-1, tnummyelements, &targetvec[0], 0, comm);
    // build an exporter and export data
    DRT::Exporter exporter(source, target, comm);
    exporter.Export(datamap);
    // put data from map in rdata
    rdata.clear();
    int count = 0;
    typename std::map<int, std::vector<T>>::const_iterator curr;
    for (curr = datamap.begin(); curr != datamap.end(); ++curr)
    {
      const std::vector<T>& current = curr->second;
      const int size = (int)current.size();
      rdata.resize((int)rdata.size() + size);
      for (int i = 0; i < size; ++i) rdata[count + i] = current[i];
      count += size;
    }
    return;
  }

  /*!
   \brief Gather information of type set<T> on a subset of processors

   This template gathers information provided in sdata on a subset of
   processors tprocs, where the length of the array tprocs is ntargetprocs.
   The redistributed data is returned in rdata which has appropiate size
   on output (size of rdata on input is arbitrary). ntargetprocs can be
   one to reduce data to one proc, it also can be equal to the total number
   of processors to make sdata redundant on all procs.

   \note Functionality of this method is equal to that of Epetra_Comm::GatherAll
   except for that the Epetra version demands the data to be of constant
   size over all procs which this method does not require!

   \param sdata (in) : Information to be gathered on tprocs.
   Length of sdata can be different on every proc.
   \param rdata (out): Information from sdata gathered on a subset of procs.
   size of rdata on input is arbitrary, it is exact on output.
   \param ntargetprocs (in): length of tprocs
   \param tprocs (in): vector of procs ids the information in sdata shall be
   gathered on.
   \param comm (in):   communicator to be used.


   */
  template <typename T>
  void Gather(std::set<T>& sdata, std::set<T>& rdata, const int ntargetprocs, const int* tprocs,
      const Epetra_Comm& comm)
  {
    const int myrank = comm.MyPID();
    const int numproc = comm.NumProc();
    if (numproc == 1)
    {
      rdata = sdata;
      return;  // nothing to do in serial
    }
    // build a map of data
    std::map<int, std::set<T>> datamap;
    datamap[myrank] = sdata;
    // build a source map
    Epetra_Map source(numproc, 1, &myrank, 0, comm);
    // build a target map which is redundant on all target procs and zero everywhere else
    bool iamtarget = false;
    for (int i = 0; i < ntargetprocs; ++i)
      if (tprocs[i] == myrank)
      {
        iamtarget = true;
        break;
      }
    std::vector<int> targetvec(0);
    if (iamtarget)
    {
      targetvec.resize(numproc);
      for (int i = 0; i < numproc; ++i) targetvec[i] = i;
    }
    const int tnummyelements = (int)targetvec.size();
    Epetra_Map target(-1, tnummyelements, &targetvec[0], 0, comm);
    // build an exporter and export data
    DRT::Exporter exporter(source, target, comm);
    exporter.Export(datamap);
    // put data from map in rdata
    rdata.clear();
    typename std::map<int, std::set<T>>::const_iterator curr;
    for (curr = datamap.begin(); curr != datamap.end(); ++curr)
    {
      const std::set<T>& current = curr->second;
      typename std::set<T>::const_iterator setiter;
      for (setiter = current.begin(); setiter != current.end(); ++setiter)
      {
        rdata.insert(*setiter);
      }
    }
    return;
  }

  /*!
   \brief Gather information of type std::map<int, std::set<T> > on a subset of processors

   This template gathers information provided in sdata on a subset of
   processors tprocs, where the length of the array tprocs is ntargetprocs.
   The redistributed data is returned in rdata which has appropiate size
   on output (size of rdata on input is arbitrary). ntargetprocs can be
   one to reduce data to one proc, it also can be equal to the total number
   of processors to make sdata redundant on all procs.

   \note Functionality of this method is equal to that of Epetra_Comm::GatherAll
   except for that the Epetra version demands the data to be of constant
   size over all procs which this method does not require!

   \param sdata (in) : Information to be gathered on tprocs.
   Length of sdata can be different on every proc.
   \param rdata (out): Information from sdata gathered on a subset of procs.
   size of rdata on input is arbitrary, it is exact on output.
   \param ntargetprocs (in): length of tprocs
   \param tprocs (in): vector of procs ids the information in sdata shall be
   gathered on.
   \param comm (in):   communicator to be used.


   */
  template <typename T>
  void Gather(std::map<int, std::set<T>>& sdata, std::map<int, std::set<T>>& rdata,
      const int ntargetprocs, const int* tprocs, const Epetra_Comm& comm)
  {
    const int myrank = comm.MyPID();
    const int numproc = comm.NumProc();
    if (numproc == 1)
    {
      rdata = sdata;
      return;  // nothing to do in serial
    }
    // build a map of data
    std::map<int, std::map<int, std::set<T>>> datamap;
    datamap[myrank] = sdata;
    // build a source map
    Epetra_Map source(numproc, 1, &myrank, 0, comm);
    // build a target map which is redundant on all target procs and zero everywhere else
    bool iamtarget = false;
    for (int i = 0; i < ntargetprocs; ++i)
      if (tprocs[i] == myrank)
      {
        iamtarget = true;
        break;
      }
    std::vector<int> targetvec(0);
    if (iamtarget)
    {
      targetvec.resize(numproc);
      for (int i = 0; i < numproc; ++i) targetvec[i] = i;
    }
    const int tnummyelements = (int)targetvec.size();
    Epetra_Map target(-1, tnummyelements, &targetvec[0], 0, comm);
    // build an exporter and export data
    DRT::Exporter exporter(source, target, comm);
    exporter.Export(datamap);
    // put data from map in rdata
    rdata.clear();
    typename std::map<int, std::map<int, std::set<T>>>::const_iterator curr;
    for (curr = datamap.begin(); curr != datamap.end(); ++curr)
    {
      const std::map<int, std::set<T>>& innercurr = curr->second;
      typename std::map<int, std::set<T>>::const_iterator inneriter;
      for (inneriter = innercurr.begin(); inneriter != innercurr.end(); ++inneriter)
      {
        rdata[inneriter->first].insert(inneriter->second.begin(), inneriter->second.end());
      }
    }
    return;
  }

  /*!
   \brief Gather information of type std::map<int, std::vector<T> > on a subset of processors

   This template gathers information provided in sdata on a subset of
   processors tprocs, where the length of the array tprocs is ntargetprocs.
   The redistributed data is returned in rdata which has appropiate size
   on output (size of rdata on input is arbitrary). ntargetprocs can be
   one to reduce data to one proc, it also can be equal to the total number
   of processors to make sdata redundant on all procs.

   \note Functionality of this method is equal to that of Epetra_Comm::GatherAll
   except for that the Epetra version demands the data to be of constant
   size over all procs which this method does not require!

   \param sdata (in) : Information to be gathered on tprocs.
   Length of sdata can be different on every proc.
   \param rdata (out): Information from sdata gathered on a subset of procs.
   size of rdata on input is arbitrary, it is exact on output.
   \param ntargetprocs (in): length of tprocs
   \param tprocs (in): vector of procs ids the information in sdata shall be
   gathered on.
   \param comm (in):   communicator to be used.


   */
  template <typename T>
  void Gather(std::map<int, std::vector<T>>& sdata, std::map<int, std::vector<T>>& rdata,
      const int ntargetprocs, const int* tprocs, const Epetra_Comm& comm)
  {
    const int myrank = comm.MyPID();
    const int numproc = comm.NumProc();
    if (numproc == 1)
    {
      rdata = sdata;
      return;  // nothing to do in serial
    }
    // build a map of data
    std::map<int, std::map<int, std::vector<T>>> datamap;
    datamap[myrank] = sdata;
    // build a source map
    Epetra_Map source(numproc, 1, &myrank, 0, comm);
    // build a target map which is redundant on all target procs and zero everywhere else
    bool iamtarget = false;
    for (int i = 0; i < ntargetprocs; ++i)
      if (tprocs[i] == myrank)
      {
        iamtarget = true;
        break;
      }
    std::vector<int> targetvec(0);
    if (iamtarget)
    {
      targetvec.resize(numproc);
      for (int i = 0; i < numproc; ++i) targetvec[i] = i;
    }
    const int tnummyelements = (int)targetvec.size();
    Epetra_Map target(-1, tnummyelements, &targetvec[0], 0, comm);
    // build an exporter and export data
    DRT::Exporter exporter(source, target, comm);
    exporter.Export(datamap);
    // put data from map in rdata
    rdata.clear();
    typename std::map<int, std::map<int, std::vector<T>>>::const_iterator curr1;
    typename std::map<int, std::vector<T>>::const_iterator curr2;
    typename std::vector<T>::const_iterator curr3;
    for (curr1 = datamap.begin(); curr1 != datamap.end(); ++curr1)
    {
      const std::map<int, std::vector<T>>& data = curr1->second;
      for (curr2 = data.begin(); curr2 != data.end(); ++curr2)
      {
        const std::vector<T>& vectordata = curr2->second;
        for (curr3 = vectordata.begin(); curr3 != vectordata.end(); ++curr3)
        {
          rdata[curr2->first].push_back(*curr3);
        }
      }
    }
    return;
  }

  /*!
   \brief Gather information of type map<T,U> on a subset of processors

   This template gathers information provided in sdata on a subset of
   processors tprocs, where the length of the array tprocs is ntargetprocs.
   The redistributed data is returned in rdata which has appropiate size
   on output (size of rdata on input is arbitrary). ntargetprocs can be
   one to reduce data to one proc, it also can be equal to the total number
   of processors to make sdata redundant on all procs.

   \param sdata (in) : Information to be gathered on tprocs.
   Length of sdata can be different on every proc.
   \param rdata (out): Information from sdata gathered on a subset of procs.
   size of rdata on input is arbitrary, it is exact on output.
   \param ntargetprocs (in): length of tprocs
   \param tprocs (in): vector of procs ids the information in sdata shall be
   gathered on.
   \param comm (in):   communicator to be used.

   */

  template <typename T, typename U>
  void Gather(std::map<T, U>& sdata, std::map<T, U>& rdata, const int ntargetprocs,
      const int* tprocs, const Epetra_Comm& comm)
  {
    const int myrank = comm.MyPID();
    const int numproc = comm.NumProc();
    if (numproc == 1)
    {
      rdata = sdata;
      return;  // nothing to do in serial
    }

    // build a map of data
    std::map<int, std::map<T, U>> datamap;
    datamap[myrank] = sdata;

    // build a source map
    Epetra_Map source(numproc, 1, &myrank, 0, comm);
    // build a target map which is redundant on all target procs and zero everywhere else
    bool iamtarget = false;
    for (int i = 0; i < ntargetprocs; ++i)
      if (tprocs[i] == myrank)
      {
        iamtarget = true;
        break;
      }
    std::vector<int> targetvec(0);
    if (iamtarget)
    {
      targetvec.resize(numproc);
      for (int i = 0; i < numproc; ++i) targetvec[i] = i;
    }
    const int tnummyelements = (int)targetvec.size();
    Epetra_Map target(-1, tnummyelements, &targetvec[0], 0, comm);
    // build an exporter and export data
    DRT::Exporter exporter(source, target, comm);
    exporter.Export(datamap);
    // put data from map in rdata
    rdata.clear();
    typename std::map<int, std::map<T, U>>::const_iterator curr;
    for (curr = datamap.begin(); curr != datamap.end(); ++curr)
    {
      const std::map<T, U>& current = curr->second;
      typename std::map<T, U>::const_iterator mapiter;
      for (mapiter = current.begin(); mapiter != current.end(); ++mapiter)
      {
        rdata.insert(std::make_pair(mapiter->first, mapiter->second));
      }
    }
    return;
  }

  /*!
   \brief Gather information of type set<T> from all processors

   This template gathers information provided in data on all processors.
   The redistributed data is returned in data which has appropriate size
   on output.

   \note Functionality of this method is equal to that of Epetra_Comm::GatherAll
   except for that the Epetra version demands the data to be of constant
   size over all procs which this method does not require!

   \param data (in/out) : Information to be gathered.
   Length of data can be different on every proc.
   \param comm (in):   communicator to be used.

   */
  template <typename T>
  void GatherAll(std::set<T>& data, const Epetra_Comm& comm)
  {
    // ntargetprocs is equal to the total number of processors to make data redundant on all procs
    const int numprocs = comm.NumProc();
    int allproc[numprocs];
    for (int i = 0; i < numprocs; ++i) allproc[i] = i;

    Gather<T>(data, data, numprocs, allproc, comm);
    return;
  }

  // nagler 07/2012
  /*!
   \brief Gather information of type map<T,U> from all processors

   This template gathers information provided in data on all processors.
   The redistributed data is returned in data which has appropriate size
   on output.

   \param data (in/out) : Information to be gathered.
   Length of data can be different on every proc.
   \param comm (in):   communicator to be used.

   */

  template <typename T, typename U>
  void GatherAll(std::map<T, U>& data, const Epetra_Comm& comm)
  {
    const int numprocs = comm.NumProc();
    int allproc[numprocs];
    for (int i = 0; i < numprocs; ++i) allproc[i] = i;

    Gather<T, U>(data, data, numprocs, allproc, comm);

    return;
  }

  //                                                              sudhakar 02/2014
  /*!
   \brief Gather information of type map<int,vector<T> > from all processors

   This template gathers information provided in data on all processors.
   The redistributed data is returned in data which has appropriate size
   on output.

   \param data (in/out) : Information to be gathered.
   Length of data can be different on every proc.
   \param comm (in):   communicator to be used.

   */

  template <typename T>
  void GatherAll(std::map<int, std::vector<T>>& data, const Epetra_Comm& comm)
  {
    const int numprocs = comm.NumProc();
    int allproc[numprocs];
    for (int i = 0; i < numprocs; ++i) allproc[i] = i;

    Gather<T>(data, data, numprocs, allproc, comm);

    return;
  }

  /*!
   \brief Gather information of type vector<T> from all processors

   This template gathers information provided in data on all processors.
   The redistributed data is returned in data which has appropriate size
   on output.

   \note Functionality of this method is equal to that of Epetra_Comm::GatherAll
   except for that the Epetra version demands the data to be of constant
   size over all procs which this method does not require!

   \param data (in/out) : Information to be gathered.
   Length of data can be different on every proc.
   \param comm (in):   communicator to be used.

   */
  template <typename T>
  void GatherAll(std::vector<T>& data, const Epetra_Comm& comm)
  {
    // ntargetprocs is equal to the total number of processors to make data redundant on all procs
    const int numprocs = comm.NumProc();
    int allproc[numprocs];
    for (int i = 0; i < numprocs; ++i) allproc[i] = i;

    Gather<T>(data, data, numprocs, allproc, comm);
    return;
  }

  /// Create an allreduced vector of gids from the given Epetra_Map
  /*!
   We have nodes and elements arbitrary global ids. On rare occasions, however,
   we need to allreduce a particular map to one or more processors.
   This is a building block for such occasions. We allreduce the gids
   of the given Epetra_Map into a vector ordered by processor number.

   \note You are not supposed to use redundant vectors in normal
   situations. If you happen to need this method you are probably
   about to do something illegal.

   \param rredundant (o) redundant vector of global ids
   \param emap (i) unique distributed Epetra_Map

   \author u.kue
   \date 05/07
   */
  void AllreduceEMap(std::vector<int>& rredundant, const Epetra_Map& emap);

  /// Create an allreduced gid to index map from the given Epetra_Map
  /*!
   We have nodes and elements with unique but otherwise arbitrary
   global ids. But unfortunately we need an allreduced vector of dof
   numbers during the dof assignment phase. In order to use such a
   vector we need to map from global ids to vector indexes. Here we
   provide that map.

   \note You are not supposed to use redundant vectors in normal
   situations. If you happen to need this method you are probably
   about to do something illegal.

   \param idxmap (o) map from global ids to (redundant) vector indexes
   \param emap (i) unique distributed Epetra_Map

   \author u.kue
   \date 05/07
   */
  void AllreduceEMap(std::map<int, int>& idxmap, const Epetra_Map& emap);

  /// Create an allreduced gid to index map from the given Epetra_Map
  /// on a distinct processor, all other procs create empty maps instead.
  /*!
   This method is currently used within the parallel post_drt_ensight
   filter in order to import all values stored in a distributed Epetra_Vector
   to processor 0 for writing them into file.

   \note see also documentation for the usual AllreduceEMap methods

   \param emap (i) any distributed Epetra_Map
   \param pid (i)  processor id where you want to have the allreduced map
   exclusively
   \author gjb
   \date 11/07
   */
  Teuchos::RCP<Epetra_Map> AllreduceEMap(const Epetra_Map& emap, const int pid);

  /// Create an allreduced Epetra_Map from the given Epetra_Map
  /// and give it to all processors.
  /*!
   This method is currently used within the constraint management, since
   current values of constraint values and langrange multipliers are distributed
   uniquely for computation. At some places we need the full information of these
   values on every processor, so this method has to be used.

   \note You are not supposed to use redundant vectors in normal
   situations. If you happen to need this method you are probably
   about to do something illegal.

   \param emap (i) any distributed Epetra_Map

   \author tk
   \date 04/08
   */
  Teuchos::RCP<Epetra_Map> AllreduceEMap(const Epetra_Map& emap);

  /// Create an allreduced Epetra_Map from the given Epetra_Map
  /// and give it to all processors.
  /*!
   Here, we have a overlapping source map and still want to have a fully
   redundant map on all processors without duplicated entries.

   \author u.kue
   \date 08/09
   */
  Teuchos::RCP<Epetra_Map> AllreduceOverlappingEMap(const Epetra_Map& emap);

  /// Create an allreduced Epetra_Map from the given Epetra_Map
  /// on a distinct processor, all other procs create empty maps instead.
  /*!
   \param emap (i) any distributed overlapping Epetra_Map
   \param pid (i)  processor id where you want to have the allreduced and
    sorted map exclusively
   \author ghamm
   \date 10/14
   */
  Teuchos::RCP<Epetra_Map> AllreduceOverlappingEMap(const Epetra_Map& emap, const int pid);

  /// find position of my map elements in a consecutive vector
  /*!
   The idea is to put the entries of a given map into a redundant
   vector, ordered by processor number. The map is assumed to be
   nonoverlapping. Here we figure out the index of our first entry in
   that vector.

   \note You are not supposed to use redundant vectors in normal
   situations. If you happen to need this method you are probably
   about to do something illegal.

   \param nummyelements (i) number of elements on this proc
   \param comm (i) communicator

   \return vector position of first entry on each processor

   \author u.kue
   \date 05/07
   */
  int FindMyPos(int nummyelements, const Epetra_Comm& comm);

  /// create an allreduced sorted copy of the source vectors
  void AllreduceVector(
      const std::vector<int>& src, std::vector<int>& dest, const Epetra_Comm& comm);

  /// Communication between all pairs of processes, with distinct data for each.
  /*!
   Sends a different vector<int> to each processes. The size of each vector may
   be different, zero-length vectors are allowed.
   Communication is implemented with the MPI function MPI_Alltoallv.

   \param comm (i) communicator
   \param send (i) vector of length comm.NumProc(), j-th element to be send to
   j-th processor.
   \param recv (o) vector of length comm.NumProc(), j-th element received from
   j-th processor.

   \author h.kue
   \date 09/07
   */
  void AllToAllCommunication(const Epetra_Comm& comm, const std::vector<std::vector<int>>& send,
      std::vector<std::vector<int>>& recv);

  /// return the first slot of a pair
  /*!
   To be used with stl algorithms.

   This should be part of stl but is not. So we define our own version.
   */
  template <typename pair_type>
  struct select1st
      : public std::unary_function<const pair_type&, const typename pair_type::first_type&>
  {
    const typename pair_type::first_type& operator()(const pair_type& v) const { return v.first; }
  };

  /// return the second slot of a pair
  /*!
   To be used with stl algorithms.

   This should be part of stl but is not. So we define our own version.
   */
  template <typename pair_type>
  struct select2nd
      : public std::unary_function<const pair_type&, const typename pair_type::second_type&>
  {
    const typename pair_type::second_type& operator()(const pair_type& v) const { return v.second; }
  };

  /// mem_fun_t version from gcc adapted to Teuchos::RCP
  /*!
   \note This is an internal class you never use directly. Use the rcp_fun
   method instead.
   */
  template <class _Ret, class _Tp>
  class rcp_fun_t : public std::unary_function<const _Tp*, _Ret>
  {
   public:
    explicit rcp_fun_t(_Ret (_Tp::*__pf)()) : _M_f(__pf) {}

    _Ret operator()(const Teuchos::RCP<_Tp> p) const { return ((&*p)->*_M_f)(); }

   private:
    _Ret (_Tp::*_M_f)();
  };

  /// mem_fun_t version from gcc adapted to Teuchos::RCP
  /*!
   \note This is an internal class you never use directly. Use the rcp_fun
   method instead.
   */
  template <class _Ret, class _Tp>
  class const_rcp_fun_t : public std::unary_function<const _Tp*, _Ret>
  {
   public:
    explicit const_rcp_fun_t(_Ret (_Tp::*__pf)() const) : _M_f(__pf) {}

    _Ret operator()(const Teuchos::RCP<_Tp> p) const { return ((&*p)->*_M_f)(); }

   private:
    _Ret (_Tp::*_M_f)() const;
  };

  /// mem_fun version from gcc adapted to Teuchos::RCP
  /*!
   This is needed if you want to use std containers (vector, set, list, ...)
   that contain Teuchos::RCPed objects with std algorithms. Use it in place of the
   usual std::mem_fun method.
   */
  template <class _Ret, class _Tp>
  inline rcp_fun_t<_Ret, _Tp> rcp_fun(_Ret (_Tp::*__f)())
  {
    return rcp_fun_t<_Ret, _Tp>(__f);
  }

  /// mem_fun version from gcc adapted to Teuchos::RCP
  /*!
   This is needed if you want to use std containers (vector, set, list, ...)
   that contain Teuchos::RCPed objects with std algorithms. Use it in place of the
   usual std::mem_fun method.
   */
  template <class _Ret, class _Tp>
  inline const_rcp_fun_t<_Ret, _Tp> rcp_fun(_Ret (_Tp::*__f)() const)
  {
    return const_rcp_fun_t<_Ret, _Tp>(__f);
  }

}  // namespace LINALG

#endif  // #ifndef LINALG_UTILS_H
