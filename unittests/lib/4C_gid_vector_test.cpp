/*---------------------------------------------------------------------------*/
/*! \file
\brief Unittests for demangle utility
\level 1
*/
/*----------------------------------------------------------------------*/
#include <gtest/gtest.h>

#include "4C_lib_utils_gid_vector.hpp"

#include <Epetra_MpiComm.h>


namespace
{
  using namespace FourC;

  TEST(BroadcastMapVector, Vector)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // at least two procs required
    ASSERT_GT(comm.NumProc(), 1);

    const int myPID = comm.MyPID();
    const std::vector<double> vec_in(myPID + 1, static_cast<double>(myPID));
    const auto vec_out = DRT::UTILS::BroadcastVector(vec_in, comm);

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
    const std::vector<std::pair<int, double>> pair_vec_in(
        myPID + 1, std::make_pair(myPID, static_cast<double>(myPID)));
    const auto pair_vec_out = DRT::UTILS::BroadcastPairVector(pair_vec_in, comm);

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
    const std::map<int, double> map_in = {std::make_pair(myPID, static_cast<double>(myPID))};
    const auto map_out = DRT::UTILS::BroadcastMap(map_in, comm);

    std::map<int, double> map_expected;
    for (int pid = 0; pid < comm.NumProc(); ++pid)
    {
      map_expected.insert(std::make_pair(pid, static_cast<double>(pid)));
    }

    EXPECT_EQ(map_out, map_expected);
  }

  TEST(BroadcastMapVector, Set)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // at least two procs required
    ASSERT_GT(comm.NumProc(), 1);

    const int myPID = comm.MyPID();
    const std::set<int> set_in = {myPID};
    const std::set<int> set_out = DRT::UTILS::BroadcastSet(set_in, comm);

    std::set<int> set_expected;
    for (int pid = 0; pid < comm.NumProc(); ++pid)
    {
      set_expected.insert(pid);
    }

    EXPECT_EQ(set_out, set_expected);
  }
}  // namespace
