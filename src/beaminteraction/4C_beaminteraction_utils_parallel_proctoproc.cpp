// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_utils_parallel_proctoproc.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_comm_mpi_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_MpiComm.h>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::Utils::i_send_receive_any(Core::FE::Discretization& discret,
    std::map<int, std::vector<std::pair<int, std::vector<int>>>> const& toranktosenddata,
    std::vector<std::pair<int, std::vector<int>>>& recvdata)
{
  // build exporter
  Core::Communication::Exporter exporter(discret.get_comm());
  int const numproc = Core::Communication::num_mpi_ranks(discret.get_comm());
  int const myrank = Core::Communication::my_mpi_rank(discret.get_comm());

  // -----------------------------------------------------------------------
  // send
  // -----------------------------------------------------------------------
  // ---- pack data for sending -----
  std::map<int, std::vector<char>> sdata;
  std::vector<int> targetprocs(numproc, 0);
  std::map<int, std::vector<std::pair<int, std::vector<int>>>>::const_iterator p;
  for (p = toranktosenddata.begin(); p != toranktosenddata.end(); ++p)
  {
    std::vector<std::pair<int, std::vector<int>>>::const_iterator iter;
    for (iter = p->second.begin(); iter != p->second.end(); ++iter)
    {
      Core::Communication::PackBuffer data;
      add_to_pack(data, *iter);
      sdata[p->first].insert(sdata[p->first].end(), data().begin(), data().end());
    }
    targetprocs[p->first] = 1;
  }

  // ---- send ----
  const int length = sdata.size();
  std::vector<MPI_Request> request(length);
  int tag = 0;
  for (std::map<int, std::vector<char>>::const_iterator p = sdata.begin(); p != sdata.end(); ++p)
  {
    exporter.i_send(
        myrank, p->first, (p->second).data(), (int)(p->second).size(), 1234, request[tag]);
    ++tag;
  }
  if (tag != length) FOUR_C_THROW("Number of messages is mixed up");

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  // ---- prepare receiving procs -----
  std::vector<int> summedtargets(numproc, 0);
  discret.get_comm().SumAll(targetprocs.data(), summedtargets.data(), numproc);

  // ---- receive ----
  for (int rec = 0; rec < summedtargets[myrank]; ++rec)
  {
    std::vector<char> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.receive_any(from, tag, rdata, length);
    if (tag != 1234)
      FOUR_C_THROW("Received on proc %i data with wrong tag from proc %i", myrank, from);

    // ---- unpack ----
    {
      // Put received nodes into discretization
      Core::Communication::UnpackBuffer buffer(rdata);
      while (!buffer.at_end())
      {
        std::pair<int, std::vector<int>> pair;
        extract_from_pack(buffer, pair);
        recvdata.push_back(pair);
      }
    }
  }

  // wait for all communications to finish
  for (int i = 0; i < length; ++i) exporter.wait(request[i]);

  // safety, should be a no time operation if everything works fine before
  discret.get_comm().Barrier();
}

FOUR_C_NAMESPACE_CLOSE
