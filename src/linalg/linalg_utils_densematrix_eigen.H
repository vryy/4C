/*----------------------------------------------------------------------*/
/*! \file
 * \file linalg_utils_densematrix_eigen.H

\brief A collection of eigenvalue methods for namespace LINALG

\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*-----------------------------------------------------------------------*/

#ifndef SRC_LINALG_LINALG_UTILS_DENSEMATRIX_EIGEN_H_
#define SRC_LINALG_LINALG_UTILS_DENSEMATRIX_EIGEN_H_

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_LAPACK.h"
#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_exporter.H"
#include "linalg_blocksparsematrix.H"

namespace LINALG
{
  // forward declaration
  // class Solver;

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

}  // namespace LINALG

#endif /* SRC_LINALG_LINALG_UTILS_DENSEMATRIX_EIGEN_H_ */
