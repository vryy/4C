/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of communication methods for namespace CORE::LINALG

\level 0
*/
/*----------------------------------------------------------------------*/

#include "baci_linalg_utils_densematrix_communication.H"

#include "baci_utils_exceptions.H"

#include <numeric>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int CORE::LINALG::FindMyPos(int nummyelements, const Epetra_Comm& comm)
{
  const int myrank = comm.MyPID();
  const int numproc = comm.NumProc();

  std::vector<int> snum(numproc, 0);
  std::vector<int> rnum(numproc);
  snum[myrank] = nummyelements;

  comm.SumAll(snum.data(), rnum.data(), numproc);

  return std::accumulate(rnum.data(), rnum.data() + myrank, 0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::LINALG::AllreduceVector(
    const std::vector<int>& src, std::vector<int>& dest, const Epetra_Comm& comm)
{
  // communicate size
  int localsize = static_cast<int>(src.size());
  int globalsize;
  comm.SumAll(&localsize, &globalsize, 1);

  // communicate values
  int pos = FindMyPos(localsize, comm);
  std::vector<int> sendglobal(globalsize, 0);
  dest.resize(globalsize);
  std::copy(src.begin(), src.end(), sendglobal.data() + pos);
  comm.SumAll(sendglobal.data(), dest.data(), globalsize);

  // sort & unique
  std::sort(dest.begin(), dest.end());
  std::vector<int>::iterator i = std::unique(dest.begin(), dest.end());
  dest.erase(i, dest.end());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::LINALG::AllreduceEMap(std::vector<int>& rredundant, const Epetra_Map& emap)
{
  const int mynodepos = FindMyPos(emap.NumMyElements(), emap.Comm());

  std::vector<int> sredundant(emap.NumGlobalElements(), 0);

  int* gids = emap.MyGlobalElements();
  std::copy(gids, gids + emap.NumMyElements(), sredundant.data() + mynodepos);

  rredundant.resize(emap.NumGlobalElements());
  emap.Comm().SumAll(sredundant.data(), rredundant.data(), emap.NumGlobalElements());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void CORE::LINALG::AllreduceEMap(std::map<int, int>& idxmap, const Epetra_Map& emap)
{
#ifdef DEBUG
  if (not emap.UniqueGIDs()) dserror("works only for unique Epetra_Maps");
#endif

  idxmap.clear();

  std::vector<int> rredundant;
  AllreduceEMap(rredundant, emap);

  for (std::size_t i = 0; i < rredundant.size(); ++i)
  {
    idxmap[rredundant[i]] = i;
  }
}

/*----------------------------------------------------------------------*
 |  create an allreduced map on a distinct processor (public)  gjb 12/07|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> CORE::LINALG::AllreduceEMap(const Epetra_Map& emap, const int pid)
{
#ifdef DEBUG
  if (not emap.UniqueGIDs()) dserror("works only for unique Epetra_Maps");
#endif
  std::vector<int> rv;
  AllreduceEMap(rv, emap);
  Teuchos::RCP<Epetra_Map> rmap;

  if (emap.Comm().MyPID() == pid)
  {
    rmap = Teuchos::rcp(new Epetra_Map(-1, rv.size(), rv.data(), 0, emap.Comm()));
    // check the map
    dsassert(rmap->NumMyElements() == rmap->NumGlobalElements(),
        "Processor with pid does not get all map elements");
  }
  else
  {
    rv.clear();
    rmap = Teuchos::rcp(new Epetra_Map(-1, 0, NULL, 0, emap.Comm()));
    // check the map
    dsassert(rmap->NumMyElements() == 0, "At least one proc will keep a map element");
  }
  return rmap;
}

/*----------------------------------------------------------------------*
 |  create an allreduced map on EVERY processor (public)        tk 12/07|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> CORE::LINALG::AllreduceEMap(const Epetra_Map& emap)
{
#ifdef DEBUG
  if (not emap.UniqueGIDs()) dserror("works only for unique Epetra_Maps");
#endif
  std::vector<int> rv;
  AllreduceEMap(rv, emap);
  Teuchos::RCP<Epetra_Map> rmap;

  rmap = Teuchos::rcp(new Epetra_Map(-1, rv.size(), rv.data(), 0, emap.Comm()));

  return rmap;
}

/*----------------------------------------------------------------------*
|  create an allreduced map on EVERY processor (public)                 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> CORE::LINALG::AllreduceOverlappingEMap(const Epetra_Map& emap)
{
  std::vector<int> rv;
  AllreduceEMap(rv, emap);

  // remove duplicates
  std::set<int> rs(rv.begin(), rv.end());
  rv.assign(rs.begin(), rs.end());

  return Teuchos::rcp(new Epetra_Map(-1, rv.size(), rv.data(), 0, emap.Comm()));
}

/*----------------------------------------------------------------------*
| create an allreduced map on a distinct processor (public)  ghamm 10/14|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> CORE::LINALG::AllreduceOverlappingEMap(
    const Epetra_Map& emap, const int pid)
{
  std::vector<int> rv;
  AllreduceEMap(rv, emap);
  Teuchos::RCP<Epetra_Map> rmap;

  if (emap.Comm().MyPID() == pid)
  {
    // remove duplicates only on proc pid
    std::set<int> rs(rv.begin(), rv.end());
    rv.assign(rs.begin(), rs.end());

    rmap = Teuchos::rcp(new Epetra_Map(-1, rv.size(), rv.data(), 0, emap.Comm()));
    // check the map
    dsassert(rmap->NumMyElements() == rmap->NumGlobalElements(),
        "Processor with pid does not get all map elements");
  }
  else
  {
    rv.clear();
    rmap = Teuchos::rcp(new Epetra_Map(-1, 0, NULL, 0, emap.Comm()));
    // check the map
    dsassert(rmap->NumMyElements() == 0, "At least one proc will keep a map element");
  }
  return rmap;
}

/*----------------------------------------------------------------------*
 |  Send and receive lists of ints.  (heiner 09/07)                     |
 *----------------------------------------------------------------------*/
void CORE::LINALG::AllToAllCommunication(const Epetra_Comm& comm,
    const std::vector<std::vector<int>>& send, std::vector<std::vector<int>>& recv)
{
  if (comm.NumProc() == 1)
  {
    dsassert(send.size() == 1, "there has to be just one entry for sending");

    // make a copy
    recv.clear();
    recv.push_back(send[0]);
  }
  else
  {
    const Epetra_MpiComm& mpicomm = dynamic_cast<const Epetra_MpiComm&>(comm);

    std::vector<int> sendbuf;
    std::vector<int> sendcounts;
    sendcounts.reserve(comm.NumProc());
    std::vector<int> sdispls;
    sdispls.reserve(comm.NumProc());

    int displacement = 0;
    sdispls.push_back(0);
    for (std::vector<std::vector<int>>::const_iterator iter = send.begin(); iter != send.end();
         ++iter)
    {
      sendbuf.insert(sendbuf.end(), iter->begin(), iter->end());
      sendcounts.push_back(iter->size());
      displacement += iter->size();
      sdispls.push_back(displacement);
    }

    std::vector<int> recvcounts(comm.NumProc());

    // initial communication: Request. Send and receive the number of
    // ints we communicate with each process.

    int status = MPI_Alltoall(
        sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, mpicomm.GetMpiComm());

    if (status != MPI_SUCCESS) dserror("MPI_Alltoall returned status=%d", status);

    std::vector<int> rdispls;
    rdispls.reserve(comm.NumProc());

    displacement = 0;
    rdispls.push_back(0);
    for (std::vector<int>::const_iterator iter = recvcounts.begin(); iter != recvcounts.end();
         ++iter)
    {
      displacement += *iter;
      rdispls.push_back(displacement);
    }

    std::vector<int> recvbuf(rdispls.back());

    // transmit communication: Send and get the data.

    status = MPI_Alltoallv(sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT,
        recvbuf.data(), recvcounts.data(), rdispls.data(), MPI_INT, mpicomm.GetMpiComm());
    if (status != MPI_SUCCESS) dserror("MPI_Alltoallv returned status=%d", status);

    recv.clear();
    for (int proc = 0; proc < comm.NumProc(); ++proc)
    {
      recv.push_back(
          std::vector<int>(recvbuf.data() + rdispls[proc], recvbuf.data() + rdispls[proc + 1]));
    }
  }
}

/*----------------------------------------------------------------------*
 |  Send and receive lists of ints.                                     |
 *----------------------------------------------------------------------*/
void CORE::LINALG::AllToAllCommunication(
    const Epetra_Comm& comm, const std::vector<std::vector<int>>& send, std::vector<int>& recv)
{
  if (comm.NumProc() == 1)
  {
    dsassert(send.size() == 1, "there has to be just one entry for sending");

    // make a copy
    recv.clear();
    recv = send[0];
  }
  else
  {
    const Epetra_MpiComm& mpicomm = dynamic_cast<const Epetra_MpiComm&>(comm);

    std::vector<int> sendbuf;
    std::vector<int> sendcounts;
    sendcounts.reserve(comm.NumProc());
    std::vector<int> sdispls;
    sdispls.reserve(comm.NumProc());

    int displacement = 0;
    sdispls.push_back(0);
    for (std::vector<std::vector<int>>::const_iterator iter = send.begin(); iter != send.end();
         ++iter)
    {
      sendbuf.insert(sendbuf.end(), iter->begin(), iter->end());
      sendcounts.push_back(iter->size());
      displacement += iter->size();
      sdispls.push_back(displacement);
    }

    std::vector<int> recvcounts(comm.NumProc());

    // initial communication: Request. Send and receive the number of
    // ints we communicate with each process.

    int status = MPI_Alltoall(
        sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, mpicomm.GetMpiComm());

    if (status != MPI_SUCCESS) dserror("MPI_Alltoall returned status=%d", status);

    std::vector<int> rdispls;
    rdispls.reserve(comm.NumProc());

    displacement = 0;
    rdispls.push_back(0);
    for (std::vector<int>::const_iterator iter = recvcounts.begin(); iter != recvcounts.end();
         ++iter)
    {
      displacement += *iter;
      rdispls.push_back(displacement);
    }

    std::vector<int> recvbuf(rdispls.back());

    // transmit communication: Send and get the data.

    recv.clear();
    recv.resize(rdispls.back());

    status = MPI_Alltoallv(sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_INT, recv.data(),
        recvcounts.data(), rdispls.data(), MPI_INT, mpicomm.GetMpiComm());
    if (status != MPI_SUCCESS) dserror("MPI_Alltoallv returned status=%d", status);
  }
}
