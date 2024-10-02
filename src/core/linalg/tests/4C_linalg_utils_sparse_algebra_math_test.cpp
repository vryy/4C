/*----------------------------------------------------------------------*/
/*! \file

\brief Unit tests for the Vector wrapper

\level 0
*/


#include <gtest/gtest.h>

#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include "4C_linalg_utils_sparse_algebra_print.hpp"

// Epetra related headers
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <Teuchos_RCPDecl.hpp>

FOUR_C_NAMESPACE_OPEN

namespace
{
  class SparseAlgebraMathTest : public testing::Test
  {
   public:
    //! Testing parameters
    Teuchos::RCP<Epetra_Comm> comm_;
    static constexpr double TOL = 1.0e-14;

   protected:
    SparseAlgebraMathTest() { comm_ = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD)); }
  };

  /** The test setup is based on a simple 1d poisson problem with the given matrix "poisson1d.mm"
   * constructed with MATLAB, having the well known [1 -2 1] tri-diagonal entries.
   *
   * All test results can be reconstructed loading the matrix into MATLAB and using the command
   * A_inverse = matrix_sparse_inverse(A, A), by using the given MATLAB script.
   */
  TEST_F(SparseAlgebraMathTest, MatrixSparseInverse1)
  {
    Epetra_CrsMatrix* A;
    const char* filename = "test_matrices/poisson1d.mm";

    int err = EpetraExt::MatrixMarketFileToCrsMatrix(filename, *comm_, A);
    if (err != 0) FOUR_C_THROW("Matrix read failed.");
    Teuchos::RCP<Epetra_CrsMatrix> A_crs = Teuchos::rcpFromRef(*A);
    Core::LinAlg::SparseMatrix A_sparse(A_crs, Core::LinAlg::Copy);

    Teuchos::RCP<Core::LinAlg::SparseMatrix> A_inverse = Core::LinAlg::matrix_sparse_inverse(
        A_sparse, Teuchos::rcp(new Epetra_CrsGraph(A->Graph())));

    // Check for global entries
    const int A_sparse_nnz = A_sparse.epetra_matrix()->NumGlobalNonzeros();
    const int A_inverse_nnz = A_inverse->epetra_matrix()->NumGlobalNonzeros();
    EXPECT_EQ(A_sparse_nnz, A_inverse_nnz);

    // Check for overall norm of matrix inverse
    EXPECT_NEAR(A_inverse->norm_frobenius(), 3.037251711528645, SparseAlgebraMathTest::TOL);
  }

  /** The test setup is based on the given nonsymmetric matrix "nonsym.mm" constructed with MATLAB.
   *  All test results can be reconstructed loading the matrix into MATLAB and using the command
   *  A_inverse = inv(A).
   *  A is given as:
   *    10     0     5     0     1
   *     0    20     0    10     0
   *     0     0    30     0    20
   *     0     0     0    40     0
   *     0     0     0     0    59
   *  With it's sparse inverse A_inverse:
   *  0.1000         0   -0.0167         0    0.0040
   *      0    0.0500         0   -0.0125         0
   *      0         0    0.0333         0   -0.0113
   *      0         0         0    0.0250         0
   *      0         0         0         0    0.0169
   */
  TEST_F(SparseAlgebraMathTest, MatrixSparseInverse2)
  {
    Epetra_CrsMatrix* A;
    const char* filename = "test_matrices/nonsym.mm";

    int err = EpetraExt::MatrixMarketFileToCrsMatrix(filename, *comm_, A);
    if (err != 0) FOUR_C_THROW("Matrix read failed.");
    Teuchos::RCP<Epetra_CrsMatrix> A_crs = Teuchos::rcpFromRef(*A);
    Core::LinAlg::SparseMatrix A_sparse(A_crs, Core::LinAlg::Copy);

    Teuchos::RCP<Core::LinAlg::SparseMatrix> A_inverse = Core::LinAlg::matrix_sparse_inverse(
        A_sparse, Teuchos::rcp(new Epetra_CrsGraph(A->Graph())));

    // Check for global entries
    const int A_sparse_nnz = A_sparse.epetra_matrix()->NumGlobalNonzeros();
    const int A_inverse_nnz = A_inverse->epetra_matrix()->NumGlobalNonzeros();
    EXPECT_EQ(A_sparse_nnz, A_inverse_nnz);

    // Check for overall norm of matrix inverse
    EXPECT_NEAR(A_inverse->norm_frobenius(), 0.1235706050986417, SparseAlgebraMathTest::TOL);

    // Check fist matrix row of inverse
    if (comm_->MyPID() == 0)
    {
      double* values;
      int length;
      A_inverse->epetra_matrix()->ExtractMyRowView(0, length, values);

      EXPECT_NEAR(values[0], 0.1, SparseAlgebraMathTest::TOL);
      EXPECT_NEAR(values[1], -0.016666666666666673, SparseAlgebraMathTest::TOL);
      EXPECT_NEAR(values[2], 0.0046666666666666688, SparseAlgebraMathTest::TOL);
    }
  }

}  // namespace

FOUR_C_NAMESPACE_CLOSE
