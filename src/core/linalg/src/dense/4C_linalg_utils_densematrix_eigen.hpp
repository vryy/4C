/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of eigenvalue methods for namespace Core::LinAlg

\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_LINALG_UTILS_DENSEMATRIX_EIGEN_HPP
#define FOUR_C_LINALG_UTILS_DENSEMATRIX_EIGEN_HPP

#include "4C_config.hpp"

#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"

#include <Teuchos_LAPACK.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  /*!
   \brief Compute all eigenvalues of a real symmetric matrix A

   \param A (in):        Matrix to be analysed
   \param L (out):       Vector of eigenvalues in ascending order
   \param postproc (in): flag indicating whether we are using this
   routine for postprocessing only (in that
   case FOUR_C_THROW is replaced with a warning)
   */
  void SymmetricEigenValues(Core::LinAlg::SerialDenseMatrix& A, Core::LinAlg::SerialDenseVector& L,
      const bool postproc = false);

  /*!
   \brief Compute all eigenvalues and eigenvectors of a real symmetric matrix A

   \param A (in/out):    in: Matrix to be analysed, out: eigenvectors
   (i.e. original matrix is destroyed!!!)
   \param L (out):       Vector of eigenvalues in ascending order
   \param postproc (in): flag indicating whether we are using this
   routine for postprocessing only (in that
   case FOUR_C_THROW is replaced with a warning)
   */
  void SymmetricEigenProblem(Core::LinAlg::SerialDenseMatrix& A, Core::LinAlg::SerialDenseVector& L,
      const bool postproc = false);

  /*!
   \brief Compute all eigenvalues and, optionally, eigenvectors
   of a real symmetric matrix A

   \param A (in/out):    Matrix to be analysed, if eigv=true stores eigenvectors
   \param L (in/out):    Vector of eigenvalues in ascending order
   \param eval_eigenvectors (in):     flag to evaluate also eigenvectors
   \param postproc (in): flag indicating whether we are using this
   routine for postprocessing only (in that
   case FOUR_C_THROW is replaced with a warning)
   */
  void SymmetricEigen(Core::LinAlg::SerialDenseMatrix& A, Core::LinAlg::SerialDenseVector& L,
      bool eval_eigenvectors, bool postproc = false);

  /*!
   \brief Compute all eigenvalues the generalized Eigenvalue problem
   Ax = lambda Bx via QZ-algorithm (B is singular) and returns the
   maximum Eigenvalue.

   \param A (in):    A Matrix
   \param B (in):    B Matrix

   */
  double GeneralizedEigen(
      Core::LinAlg::SerialDenseMatrix::Base& A, Core::LinAlg::SerialDenseMatrix::Base& B);

  /*!
   \brief Compute all eigenvalues and eigenvectors of a real symmetric matrix A

   A = V * S * VT

   \param A (in):        M-by-M matrix to be decomposed
   \param S (out):       M-by-N matrix which is zero except for its diagonal entries holding the
   eigenvalues \param V (out):       M-by-M orthonormal matrix of eigenvectors
   */
  template <unsigned int dim>
  void SYEV(Core::LinAlg::Matrix<dim, dim>& A, Core::LinAlg::Matrix<dim, dim>& S,
      Core::LinAlg::Matrix<dim, dim>& V)
  {
    const char jobz = 'V';               // Compute eigenvalues and eigenvectors.
    const char uplo = 'U';               // Upper triangle of A is stored;
    const int N = dim;                   // The order of the matrix A.  N >= 0.
    Matrix<dim, dim> tmp(A.A(), false);  // copy, because content of matrix is destroyed
    const int lda = dim;                 // The leading dimension of the array A.  LDA >=max(1,N).
    std::vector<double> w(dim);
    const int lwork = 2 * dim * dim + 6 * dim + 1;
    std::vector<double> work(lwork);
    int info;

    Teuchos::LAPACK<int, double> lapack;
    lapack.SYEV(jobz, uplo, N, tmp.A(), lda, w.data(), work.data(), lwork, &info);

    if (info) FOUR_C_THROW("Lapack's SYEV returned %d", info);

    // return eigenvectors
    V.Update(tmp);

    // return eigenvalues
    S.Clear();
    for (unsigned int i = 0; i < dim; ++i) S(i, i) = w[i];

    return;
  }

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
