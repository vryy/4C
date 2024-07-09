/*-----------------------------------------------------------*/
/*! \file

\brief utils for proc to proc communication


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_beaminteraction_utils_parallel_proctoproc.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_MpiComm.h>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Discret::UTILS::ISendReceiveAny(Teuchos::RCP<Core::FE::Discretization> const& discret,
    std::map<int, std::vector<std::pair<int, std::vector<int>>>> const& toranktosenddata,
    std::vector<std::pair<int, std::vector<int>>>& recvdata)
{
  // build exporter
  Core::Communication::Exporter exporter(discret->get_comm());
  int const numproc = discret->get_comm().NumProc();
  int const myrank = discret->get_comm().MyPID();

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
      Core::Communication::ParObject::add_to_pack(data, *iter);
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
  discret->get_comm().SumAll(targetprocs.data(), summedtargets.data(), numproc);

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
      std::vector<char>::size_type index = 0;
      while (index < rdata.size())
      {
        std::pair<int, std::vector<int>> pair;
        Core::Communication::ParObject::extract_from_pack(index, rdata, pair);
        recvdata.push_back(pair);
      }
      if (index != rdata.size())
        FOUR_C_THROW("Mismatch in size of data %d <-> %d", static_cast<int>(rdata.size()), index);
    }
  }

  // wait for all communications to finish
  for (int i = 0; i < length; ++i) exporter.wait(request[i]);

  // safety, should be a no time operation if everything works fine before
  discret->get_comm().Barrier();
}

FOUR_C_NAMESPACE_CLOSE
