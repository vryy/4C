// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_linalg_sparsematrix.hpp"

#include "4C_linalg_map.hpp"

#include <mpi.h>

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace
{
  constexpr int num_rows = 5;

  // Tridiagonal matrix with an optional additional entry in the upper right corner
  std::unique_ptr<Core::LinAlg::SparseMatrix> make_matrix(
      const Core::LinAlg::Map& map, double diagonal, bool corner, bool optimize_data_storage)
  {
    auto matrix = std::make_unique<Core::LinAlg::SparseMatrix>(map, 3);
    for (int gid = 0; gid < num_rows; ++gid)
    {
      matrix->set_value(diagonal, gid, gid);
      if (gid > 0) matrix->set_value(-1.0, gid, gid - 1);
      if (gid < num_rows - 1) matrix->set_value(-1.0, gid, gid + 1);
    }
    if (corner) matrix->set_value(0.5, 0, num_rows - 1);

    matrix->complete({.optimize_data_storage = optimize_data_storage});
    return matrix;
  }

  class LocalSparsityPatternTest : public testing::TestWithParam<bool>
  {
   protected:
    const Core::LinAlg::Map map_{num_rows, 0, MPI_COMM_SELF};
  };

  TEST_P(LocalSparsityPatternTest, PatternIsInCrsForm)
  {
    const auto matrix = make_matrix(map_, 2.0, false, GetParam());

    const Core::LinAlg::LocalSparsityPattern pattern = matrix->local_sparsity_pattern();

    EXPECT_EQ(pattern.row_offsets, (std::vector<std::size_t>{0, 2, 5, 8, 11, 13}));
    EXPECT_EQ(pattern.column_indices, (std::vector<int>{0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4}));
  }

  TEST_P(LocalSparsityPatternTest, ChangedValuesKeepPattern)
  {
    const auto matrix = make_matrix(map_, 2.0, false, GetParam());
    const auto other_values = make_matrix(map_, 7.0, false, GetParam());

    EXPECT_TRUE(other_values->has_local_sparsity_pattern(matrix->local_sparsity_pattern()));
    EXPECT_TRUE(matrix->has_local_sparsity_pattern(matrix->local_sparsity_pattern()));
  }

  TEST_P(LocalSparsityPatternTest, AdditionalEntryChangesPattern)
  {
    const auto matrix = make_matrix(map_, 2.0, false, GetParam());
    const auto with_corner = make_matrix(map_, 2.0, true, GetParam());

    EXPECT_FALSE(with_corner->has_local_sparsity_pattern(matrix->local_sparsity_pattern()));
    EXPECT_FALSE(matrix->has_local_sparsity_pattern(with_corner->local_sparsity_pattern()));
  }

  TEST_P(LocalSparsityPatternTest, DifferentNumberOfRowsChangesPattern)
  {
    const auto matrix = make_matrix(map_, 2.0, false, GetParam());

    Core::LinAlg::LocalSparsityPattern pattern = matrix->local_sparsity_pattern();
    pattern.row_offsets.pop_back();

    EXPECT_FALSE(matrix->has_local_sparsity_pattern(pattern));
  }

  INSTANTIATE_TEST_SUITE_P(OptimizedAndNotOptimizedStorage, LocalSparsityPatternTest,
      testing::Values(true, false), [](const testing::TestParamInfo<bool>& info)
      { return info.param ? "OptimizedStorage" : "NotOptimizedStorage"; });
}  // namespace

FOUR_C_NAMESPACE_CLOSE
