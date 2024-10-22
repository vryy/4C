// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_config.hpp"

#include <mpi.h>

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  // The following flag makes DeathTests thread-safe at the cost of a severely increased runtime
  // (only of the DeathTests). See also
  // https://github.com/google/googletest/blob/main/docs/advanced.md#death-test-styles
  (void)(::testing::GTEST_FLAG(death_test_style) = "threadsafe");

  MPI_Init(&argc, &argv);

  const int result = RUN_ALL_TESTS();

  MPI_Finalize();

  return result;
}
