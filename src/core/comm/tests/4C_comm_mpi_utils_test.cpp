// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_comm_mpi_utils.hpp"

#include <Epetra_MpiComm.h>


namespace
{
  using namespace FourC;

  TEST(AllGather, Vector)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // at least two procs required
    ASSERT_GT(comm.NumProc(), 1);

    const int myPID = comm.MyPID();
    const std::vector<double> vec_in(myPID + 1, static_cast<double>(myPID));
    const auto vec_out = Core::Communication::all_gather(vec_in, comm);

    std::vector<double> vec_expected;
    for (int pid = 0; pid < comm.NumProc(); ++pid)
    {
      for (int i = 0; i < pid + 1; ++i) vec_expected.emplace_back(pid);
    }

    EXPECT_EQ(vec_out, vec_expected);
  }

  TEST(AllGather, PairVector)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // at least two procs required
    ASSERT_GT(comm.NumProc(), 1);

    const int myPID = comm.MyPID();
    const std::vector<std::pair<int, double>> pair_vec_in(
        myPID + 1, std::make_pair(myPID, static_cast<double>(myPID)));
    const auto pair_vec_out = Core::Communication::all_gather(pair_vec_in, comm);

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

  TEST(AllGather, Map)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // at least two procs required
    ASSERT_GT(comm.NumProc(), 1);

    const int myPID = comm.MyPID();
    const std::map<int, double> map_in = {std::make_pair(myPID, static_cast<double>(myPID))};
    const auto map_out = Core::Communication::all_gather(map_in, comm);

    std::map<int, double> map_expected;
    for (int pid = 0; pid < comm.NumProc(); ++pid)
    {
      map_expected.insert(std::make_pair(pid, static_cast<double>(pid)));
    }

    EXPECT_EQ(map_out, map_expected);
  }

  TEST(AllGather, UnorderedMap)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // at least two procs required
    ASSERT_GT(comm.NumProc(), 1);

    const int myPID = comm.MyPID();
    const std::unordered_map<int, double> map_in = {
        std::make_pair(myPID, static_cast<double>(myPID))};
    const auto map_out = Core::Communication::all_gather(map_in, comm);

    std::unordered_map<int, double> map_expected;
    for (int pid = 0; pid < comm.NumProc(); ++pid)
    {
      map_expected.insert(std::make_pair(pid, static_cast<double>(pid)));
    }

    EXPECT_EQ(map_out, map_expected);
  }

  TEST(AllGather, UnorderedMultiMap)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // at least two procs required
    ASSERT_GT(comm.NumProc(), 1);

    const int myPID = comm.MyPID();
    const std::unordered_multimap<int, int> map_in = {std::make_pair(1, myPID)};
    const auto map_out = Core::Communication::all_gather(map_in, comm);

    std::unordered_multimap<int, int> map_expected;
    for (int pid = 0; pid < comm.NumProc(); ++pid)
    {
      map_expected.insert(std::make_pair(1, pid));
    }

    EXPECT_EQ(map_out, map_expected);
  }

  TEST(AllGather, Set)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // at least two procs required
    ASSERT_GT(comm.NumProc(), 1);

    const int myPID = comm.MyPID();
    const std::set<int> set_in = {myPID};
    const std::set<int> set_out = Core::Communication::all_gather(set_in, comm);

    std::set<int> set_expected;
    for (int pid = 0; pid < comm.NumProc(); ++pid)
    {
      set_expected.insert(pid);
    }

    EXPECT_EQ(set_out, set_expected);
  }


  TEST(AllGather, VectorNonTrivial)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // at least two procs required
    ASSERT_GT(comm.NumProc(), 1);

    const int myPID = comm.MyPID();
    const std::vector<std::string> vec_in(myPID + 1, std::to_string(myPID));
    const auto vec_out = Core::Communication::all_gather(vec_in, comm);

    std::vector<std::string> vec_expected;
    for (int pid = 0; pid < comm.NumProc(); ++pid)
    {
      for (int i = 0; i < pid + 1; ++i) vec_expected.emplace_back(std::to_string(pid));
    }

    EXPECT_EQ(vec_out, vec_expected);
  }

  TEST(AllGather, ComplicatedMap)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // at least two procs required
    ASSERT_GT(comm.NumProc(), 1);

    const int myPID = comm.MyPID();
    const std::map<std::string, std::pair<int, double>> in = {
        {std::to_string(myPID), {myPID, 2.0}}};
    const auto out = Core::Communication::all_gather(in, comm);

    std::map<std::string, std::pair<int, double>> expected;
    for (int pid = 0; pid < comm.NumProc(); ++pid)
    {
      expected.insert({std::to_string(pid), {pid, 2.0}});
    }

    EXPECT_EQ(out, expected);
  }
}  // namespace
