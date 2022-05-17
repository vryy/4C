/*----------------------------------------------------------------------*/
/*! \file
\brief This file provide the entry point for any GoogleTest test
       executable running with MPI Initialized.
\level 0
*/



#include <gtest/gtest.h>

#include <mpi.h>

/*!
 * Release and delete Google Test's test listeners for all MPI ranks > 0.
 */
void DeleteDuplicateTestListeners()
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank > 0)
  {
    testing::TestEventListeners& listeners = testing::UnitTest::GetInstance()->listeners();
    delete listeners.Release(listeners.default_result_printer());
  }
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  // The following flag makes DeathTests thread-safe at the cost of a severely increased runtime
  // (only of the DeathTests). See also
  // https://github.com/google/googletest/blob/main/docs/advanced.md#death-test-styles
  (void)(::testing::GTEST_FLAG(death_test_style) = "threadsafe");

  MPI_Init(&argc, &argv);

  DeleteDuplicateTestListeners();
  const int result = RUN_ALL_TESTS();

  MPI_Finalize();

  return result;
}
