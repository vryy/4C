// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_linalg_utils_sparse_algebra_assemble.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_linalg_utils_sparse_algebra_io.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_unittest_utils_support_files_test.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  class SparseAlgebraAssembleTest : public testing::Test
  {
   public:
    //! Testing parameters
    MPI_Comm comm_;

   protected:
    SparseAlgebraAssembleTest() { comm_ = MPI_COMM_WORLD; }
  };

  /** The test setup is based on a simple 1d poisson problem with the given matrix "poisson1d.mm"
   * constructed with MATLAB, having the well known [1 -2 1] tri-diagonal entries.
   *
   * The matrix is missing Dirichlet boundary conditions.
   */
  TEST_F(SparseAlgebraAssembleTest, HasDirichletBoundaryCondition)
  {
    bool has_dirichlet;
    int rank = Core::Communication::my_mpi_rank(comm_);

    Core::LinAlg::SparseMatrix A = Core::LinAlg::read_matrix_market_file_as_sparse_matrix(
        TESTING::get_support_file_path("test_matrices/poisson1d.mm").c_str(), comm_);

    has_dirichlet = Core::LinAlg::has_dirichlet_boundary_condition(A);
    EXPECT_EQ(has_dirichlet, false);

    // Explicitly modify a row to be Dirichlet
    if (rank == 0)
    {
      int num_entries;
      int* indices;
      double* values;

      A.extract_my_row_view(0, num_entries, values, indices);

      for (int i = 0; i < num_entries; i++)
      {
        if (i == 0)
          values[i] = 1.0;
        else
          values[i] = 0.0;
      }

      A.replace_global_values(0, num_entries, values, indices);
    }

    has_dirichlet = Core::LinAlg::has_dirichlet_boundary_condition(A);
    EXPECT_EQ(has_dirichlet, true);
  }
}  // namespace

FOUR_C_NAMESPACE_CLOSE
