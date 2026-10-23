// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_linalg_sparsematrix.hpp"

#include "4C_linalg_graph.hpp"
#include "4C_linalg_map.hpp"

#include <mpi.h>

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace
{
  constexpr int num_rows = 5;

  // Tridiagonal matrix with an optional additional entry in the upper right corner. The global ids
  // start at first_gid.
  std::unique_ptr<Core::LinAlg::SparseMatrix> make_matrix(
      double diagonal, bool corner, bool optimize_data_storage, int first_gid = 0)
  {
    const Core::LinAlg::Map map(num_rows, first_gid, MPI_COMM_SELF);

    auto matrix = std::make_unique<Core::LinAlg::SparseMatrix>(map, 3);
    for (int i = 0; i < num_rows; ++i)
    {
      matrix->set_value(diagonal, first_gid + i, first_gid + i);
      if (i > 0) matrix->set_value(-1.0, first_gid + i, first_gid + i - 1);
      if (i < num_rows - 1) matrix->set_value(-1.0, first_gid + i, first_gid + i + 1);
    }
    if (corner) matrix->set_value(0.5, first_gid, first_gid + num_rows - 1);

    matrix->complete({.optimize_data_storage = optimize_data_storage});
    return matrix;
  }

  class SparseMatrixGraphTest : public testing::TestWithParam<bool>
  {
  };

  TEST_P(SparseMatrixGraphTest, ChangedValuesKeepGraph)
  {
    const auto matrix = make_matrix(2.0, false, GetParam());
    const auto other_values = make_matrix(7.0, false, GetParam());

    const Core::LinAlg::Graph graph = matrix->get_current_graph();

    EXPECT_TRUE(graph.is_identical_to(other_values->get_current_graph()));
    EXPECT_TRUE(other_values->get_current_graph().is_identical_to(graph));
    EXPECT_TRUE(graph.is_identical_to(matrix->get_current_graph()));
  }

  TEST_P(SparseMatrixGraphTest, AdditionalEntryChangesGraph)
  {
    const auto matrix = make_matrix(2.0, false, GetParam());
    const auto with_corner = make_matrix(2.0, true, GetParam());

    EXPECT_FALSE(with_corner->get_current_graph().is_identical_to(matrix->get_current_graph()));
    EXPECT_FALSE(matrix->get_current_graph().is_identical_to(with_corner->get_current_graph()));
  }

  TEST_P(SparseMatrixGraphTest, DifferentMapsChangeGraph)
  {
    const auto matrix = make_matrix(2.0, false, GetParam());
    const auto shifted = make_matrix(2.0, false, GetParam(), 10);

    EXPECT_FALSE(shifted->get_current_graph().is_identical_to(matrix->get_current_graph()));
  }

  TEST_P(SparseMatrixGraphTest, GraphOutlivesMatrix)
  {
    auto matrix = make_matrix(2.0, false, GetParam());
    const Core::LinAlg::Graph graph = matrix->get_current_graph();
    const auto same_graph = make_matrix(3.0, false, GetParam());

    matrix.reset();

    EXPECT_TRUE(graph.is_identical_to(same_graph->get_current_graph()));
  }

  INSTANTIATE_TEST_SUITE_P(OptimizedAndNotOptimizedStorage, SparseMatrixGraphTest,
      testing::Values(true, false), [](const testing::TestParamInfo<bool>& info)
      { return info.param ? "OptimizedStorage" : "NotOptimizedStorage"; });
}  // namespace

FOUR_C_NAMESPACE_CLOSE
