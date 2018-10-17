/*---------------------------------------------------------------------------*/
/*!
\file particle_communication_utils.cpp

\brief communication utils for particle engine

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_communication_utils.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | communicate data via non-buffered send from proc to proc   sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEENGINE::COMMUNICATION::ImmediateRecvBlockingSend(const Epetra_Comm& comm,
    std::map<int, std::vector<char>>& sdata, std::map<int, std::vector<char>>& rdata)
{
  // number of processors
  int const numproc = comm.NumProc();

  // processor id
  int const myrank = comm.MyPID();

  // number of processors receiving data from this processor
  int const numsendtoprocs = sdata.size();
  int counter = 0;

  // mpi communicator
  const Epetra_MpiComm* mpicomm = dynamic_cast<const Epetra_MpiComm*>(&comm);
  if (!mpicomm) dserror("dynamic cast to Epetra_MpiComm failed!");

  // ---- communicate target processors to all processors ----
  std::vector<int> targetprocs(numproc, 0);
  std::vector<int> summedtargets(numproc, 0);

  for (auto& p : sdata) targetprocs[p.first] = 1;

  comm.SumAll(targetprocs.data(), summedtargets.data(), numproc);

  // number of processors this processor receives data from
  int const numrecvfromprocs = summedtargets[myrank];

  // ---- send size of messages to receiving processors ----
  std::vector<MPI_Request> sizerequest(numsendtoprocs);
  counter = 0;
  for (auto& p : sdata)
  {
    int torank = p.first;
    if (myrank == torank) dserror("processor should not send messages to itself!");
    if (torank == -1) dserror("processor can not send messages to processor -1!");

    int msgsizetosend = static_cast<int>((p.second).size());

    // perform non-blocking send operation
    MPI_Isend(&msgsizetosend, 1, MPI_INT, torank, 1234, mpicomm->Comm(), &sizerequest[counter]);

    ++counter;
  }
  if (counter != numsendtoprocs) dserror("number of messages is mixed up!");

  // ---- receive size of messages ----
  std::vector<MPI_Request> recvrequest(numrecvfromprocs);
  std::set<int> recvrunning;
  for (int rec = 0; rec < numrecvfromprocs; ++rec)
  {
    int msgsize = -1;
    int msgtag = -1;
    int msgsource = -1;
    int msgsizetorecv = -1;

    // probe for any message to come
    MPI_Status status;
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, mpicomm->Comm(), &status);

    // get sender, tag and size
    msgsource = status.MPI_SOURCE;
    msgtag = status.MPI_TAG;
    MPI_Get_count(&status, MPI_INT, &msgsize);

    // check message tag
    if (msgtag != 1234)
      dserror("Received on proc %i data with wrong tag from proc %i", myrank, msgsource);

    // check size of message
    if (msgsize != 1) dserror("message size not correct (one int expected)!");

    // perform blocking receive operation
    MPI_Recv(&msgsizetorecv, msgsize, MPI_INT, msgsource, msgtag, mpicomm->Comm(), &status);

    // check received size of message
    if (msgsizetorecv < 0) dserror("received message size is negative!");

    // resize receiving buffer to received size
    (rdata[msgsource]).resize(msgsizetorecv);

    // perform non-blocking receive operation
    MPI_Irecv((void*)(&(rdata[msgsource])[0]), msgsizetorecv, MPI_CHAR, msgsource, msgtag,
        mpicomm->Comm(), &recvrequest[rec]);
    recvrunning.insert(rec);
  }

  // wait until every processor is ready to receive all data
  comm.Barrier();

  // ---- send data to already waiting processors ----
  counter = 0;
  for (auto& p : sdata)
  {
    // perform blocking send operation
    MPI_Send((void*)(&((p.second)[0])), static_cast<int>((p.second).size()), MPI_CHAR, p.first,
        1234, mpicomm->Comm());

    ++counter;
  }
  if (counter != numsendtoprocs) dserror("number of messages is mixed up!");

  // clear send buffer after all successful communication
  sdata.clear();

  // ---- wait for completion of running receive operations ----
  for (int rec : recvrunning)
  {
    // wait for non-blocking receive operation
    MPI_Status status;
    MPI_Wait(&recvrequest[rec], &status);
  }
}
