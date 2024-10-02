/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for the sparse algebra manipulation utils

\level 0
*/
/*----------------------------------------------------------------------*/

#include <gtest/gtest.h>

#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <Epetra_MpiComm.h>
#include <EpetraExt_CrsMatrixIn.h>

FOUR_C_NAMESPACE_OPEN

namespace
{
  class SparseAlgebraManipulationTest : public testing::Test
  {
   public:
    //! Testing parameters
    Teuchos::RCP<Epetra_Comm> comm_;

   protected:
    SparseAlgebraManipulationTest() { comm_ = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD)); }
  };

  /** The test setup is based on a simple 1d poisson problem with the given matrix "poisson1d.mm"
   * constructed with MATLAB, having the well known [1 -2 1] tri-diagonal entries.
   *
   * Due to the threshold tol=1.1, the off-diagonal values are dropped resulting in a diagonal
   * matrix with -2 values on the main diagonal.
   */
  TEST_F(SparseAlgebraManipulationTest, ThresholdMatrix1)
  {
    Epetra_CrsMatrix* A;
    const char* filename = "test_matrices/poisson1d.mm";

    int err = EpetraExt::MatrixMarketFileToCrsMatrix(filename, *comm_, A);
    if (err != 0) FOUR_C_THROW("Matrix read failed.");
    Teuchos::RCP<Epetra_CrsMatrix> A_crs = Teuchos::rcpFromRef(*A);
    Core::LinAlg::SparseMatrix A_sparse(A_crs, Core::LinAlg::Copy);

    const double tol = 1.1;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> A_thresh =
        Core::LinAlg::threshold_matrix(A_sparse, tol);

    // Check for global entries
    const int A_thresh_nnz = A_thresh->epetra_matrix()->NumGlobalNonzeros();
    EXPECT_EQ(A_thresh_nnz, 20);

    // Check for overall norm of matrix
    EXPECT_NEAR(A_thresh->norm_frobenius(), 8.944271909999159, 1e-12);
  }

  /** The test setup is based on the simple matrix  "filter.mm" constructed with MATLAB.
   *
   * Due to the threshold all values smaller than 1e-5 are dropped.
   */
  TEST_F(SparseAlgebraManipulationTest, ThresholdMatrix2)
  {
    Epetra_CrsMatrix* A;
    const char* filename = "test_matrices/filter.mm";

    int err = EpetraExt::MatrixMarketFileToCrsMatrix(filename, *comm_, A);
    if (err != 0) FOUR_C_THROW("Matrix read failed.");
    Teuchos::RCP<Epetra_CrsMatrix> A_crs = Teuchos::rcpFromRef(*A);
    Core::LinAlg::SparseMatrix A_sparse(A_crs, Core::LinAlg::Copy);

    const double tol = 1e-5;
    Teuchos::RCP<Core::LinAlg::SparseMatrix> A_thresh =
        Core::LinAlg::threshold_matrix(A_sparse, tol);

    // Check for global entries
    const int A_thresh_nnz = A_thresh->epetra_matrix()->NumGlobalNonzeros();
    EXPECT_EQ(A_thresh_nnz, 13);

    // Check for overall norm of matrix
    EXPECT_NEAR(A_thresh->norm_frobenius(), 7.549834435270750, 1e-12);
  }

  /** The test setup is based on the simple matrix  "filter.mm" constructed with MATLAB.
   *
   * Due to the threshold all indices based on values smaller than 1e-5 are dropped from the graph.
   */
  TEST_F(SparseAlgebraManipulationTest, ThresholdMatrixGraph)
  {
    Epetra_CrsMatrix* A;
    const char* filename = "test_matrices/filter.mm";

    int err = EpetraExt::MatrixMarketFileToCrsMatrix(filename, *comm_, A);
    if (err != 0) FOUR_C_THROW("Matrix read failed.");
    Teuchos::RCP<Epetra_CrsMatrix> A_crs = Teuchos::rcpFromRef(*A);
    Core::LinAlg::SparseMatrix A_sparse(A_crs, Core::LinAlg::Copy);

    const double tol = 1e-5;
    Teuchos::RCP<Epetra_CrsGraph> G = Core::LinAlg::threshold_matrix_graph(A_sparse, tol);

    // Check for global entries
    const int A_thresh_nnz = G->NumGlobalNonzeros();
    EXPECT_EQ(A_thresh_nnz, 13);
  }
}  // namespace

FOUR_C_NAMESPACE_CLOSE
