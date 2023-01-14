/*---------------------------------------------------------------------------*/
/*! \file
\brief Unittests for demangle utility
\level 1
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>
#include "lib_utils_gid_vector.H"
#include <Epetra_MpiComm.h>


namespace
{
  TEST(BroadcastMapVector, Vector)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // at least two procs required
    ASSERT_GT(comm.NumProc(), 1);

    const int myPID = comm.MyPID();
    std::vector<double> vec_in(myPID + 1, static_cast<double>(myPID));
    std::vector<double> vec_out;
    DRT::UTILS::BroadcastVector(vec_in, vec_out, comm);

    std::vector<double> vec_expected;
    for (int pid = 0; pid < comm.NumProc(); ++pid)
    {
      for (int i = 0; i < pid + 1; ++i) vec_expected.emplace_back(pid);
    }

    EXPECT_EQ(vec_out, vec_expected);
  }

  TEST(BroadcastMapVector, PairVector)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // at least two procs required
    ASSERT_GT(comm.NumProc(), 1);

    const int myPID = comm.MyPID();
    std::vector<std::pair<int, double>> pair_vec_in(
        myPID + 1, std::make_pair(myPID, static_cast<double>(myPID)));
    std::vector<std::pair<int, double>> pair_vec_out;
    DRT::UTILS::BroadcastPairVector(pair_vec_in, pair_vec_out, comm);

    std::vector<std::pair<int, double>> pair_vec_expected;
    for (int pid = 0; pid < comm.NumProc(); ++pid)
    {
      for (int i = 0; i < pid + 1; ++i)
      {
        pair_vec_expected.emplace_back(std::make_pair(pid, static_cast<double>(pid)));
      }
    }

    EXPECT_EQ(pair_vec_out, pair_vec_expected);
  }

  TEST(BroadcastMapVector, Map)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // at least two procs required
    ASSERT_GT(comm.NumProc(), 1);

    const int myPID = comm.MyPID();

    std::map<int, double> map_in;
    map_in.insert(std::make_pair(myPID, static_cast<double>(myPID)));

    std::map<int, double> map_out;
    DRT::UTILS::BroadcastMap(map_in, map_out, comm);

    std::map<int, double> map_expected;
    for (int pid = 0; pid < comm.NumProc(); ++pid)
    {
      map_expected.insert(std::make_pair(pid, static_cast<double>(pid)));
    }

    EXPECT_EQ(map_out, map_expected);
  }
}  // namespace
