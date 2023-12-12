/*-----------------------------------------------------------*/
/*! \file

\brief utils for proc to proc communication


\level 3

*/
/*-----------------------------------------------------------*/

#include "baci_beaminteraction_utils_parallel_proctoproc.H"

#include "baci_comm_exporter.H"
#include "baci_lib_discret.H"
#include "baci_utils_exceptions.H"

#include <Epetra_MpiComm.h>

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::UTILS::ISendReceiveAny(Teuchos::RCP<BACI::DRT::Discretization> const& discret,
    std::map<int, std::vector<std::pair<int, std::vector<int>>>> const& toranktosenddata,
    std::vector<std::pair<int, std::vector<int>>>& recvdata)
{
  // build exporter
  CORE::COMM::Exporter exporter(discret->Comm());
  int const numproc = discret->Comm().NumProc();
  int const myrank = discret->Comm().MyPID();

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
      CORE::COMM::PackBuffer data;
      CORE::COMM::ParObject::AddtoPack(data, *iter);
      data.StartPacking();
      CORE::COMM::ParObject::AddtoPack(data, *iter);
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
    exporter.ISend(
        myrank, p->first, (p->second).data(), (int)(p->second).size(), 1234, request[tag]);
    ++tag;
  }
  if (tag != length) dserror("Number of messages is mixed up");

  // -----------------------------------------------------------------------
  // receive
  // -----------------------------------------------------------------------
  // ---- prepare receiving procs -----
  std::vector<int> summedtargets(numproc, 0);
  discret->Comm().SumAll(targetprocs.data(), summedtargets.data(), numproc);

  // ---- receive ----
  for (int rec = 0; rec < summedtargets[myrank]; ++rec)
  {
    std::vector<char> rdata;
    int length = 0;
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny(from, tag, rdata, length);
    if (tag != 1234) dserror("Received on proc %i data with wrong tag from proc %i", myrank, from);

    // ---- unpack ----
    {
      // Put received nodes into discretization
      std::vector<char>::size_type index = 0;
      while (index < rdata.size())
      {
        std::pair<int, std::vector<int>> pair;
        CORE::COMM::ParObject::ExtractfromPack(index, rdata, pair);
        recvdata.push_back(pair);
      }
      if (index != rdata.size())
        dserror("Mismatch in size of data %d <-> %d", static_cast<int>(rdata.size()), index);
    }
  }

  // wait for all communications to finish
  for (int i = 0; i < length; ++i) exporter.Wait(request[i]);

  // safety, should be a no time operation if everything works fine before
  discret->Comm().Barrier();
}

BACI_NAMESPACE_CLOSE
