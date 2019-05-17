/*---------------------------------------------------------------------*/
/*!

\brief Implementation of exporter class

\level 0

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/


/*!----------------------------------------------------------------------
 *----------------------------------------------------------------------*/

#include "vector"
#include "drt_exporter.H"
#include "drt_utils.H"
#include "drt_dserror.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Exporter::Exporter(const Epetra_Comm& comm)
    : dummymap_(0, 0, comm),
      frommap_(dummymap_),
      tomap_(dummymap_),
      comm_(comm),
      myrank_(comm.MyPID()),
      numproc_(comm.NumProc())
{
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Exporter::Exporter(const Epetra_Map& frommap, const Epetra_Map& tomap, const Epetra_Comm& comm)
    : dummymap_(0, 0, comm),
      frommap_(frommap),
      tomap_(tomap),
      comm_(comm),
      myrank_(comm.MyPID()),
      numproc_(comm.NumProc())
{
  ConstructExporter();
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Exporter::Exporter(const DRT::Exporter& old)
    : dummymap_(old.dummymap_),
      frommap_(old.frommap_),
      tomap_(old.tomap_),
      comm_(old.comm_),
      myrank_(old.myrank_),
      numproc_(old.numproc_),
      sendplan_(old.sendplan_)
#if 0
recvplan_(old.recvplan_)
#endif
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Exporter::~Exporter() { return; }


/*----------------------------------------------------------------------*
 |  do a send of data (public)                               mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Exporter::ISend(const int frompid, const int topid, const char* data, const int dsize,
    const int tag, MPI_Request& request)
{
  if (MyPID() != frompid) return;
  const Epetra_MpiComm* comm = dynamic_cast<const Epetra_MpiComm*>(&(Comm()));
  if (!comm) dserror("Comm() is not a Epetra_MpiComm\n");
  MPI_Isend((void*)data, dsize, MPI_CHAR, topid, tag, comm->Comm(), &request);
  return;
}

/*----------------------------------------------------------------------*
 |  do a send of data (public)                               mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Exporter::ISend(const int frompid, const int topid, const int* data, const int dsize,
    const int tag, MPI_Request& request)
{
  if (MyPID() != frompid) return;
  const Epetra_MpiComm* comm = dynamic_cast<const Epetra_MpiComm*>(&(Comm()));
  if (!comm) dserror("Comm() is not a Epetra_MpiComm\n");
  MPI_Isend((void*)data, dsize, MPI_INT, topid, tag, comm->Comm(), &request);
  return;
}

/*----------------------------------------------------------------------*
 |  do a send of data (public)                               mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Exporter::ISend(const int frompid, const int topid, const double* data, const int dsize,
    const int tag, MPI_Request& request)
{
  if (MyPID() != frompid) return;
  const Epetra_MpiComm* comm = dynamic_cast<const Epetra_MpiComm*>(&(Comm()));
  if (!comm) dserror("Comm() is not a Epetra_MpiComm\n");
  MPI_Isend((void*)data, dsize, MPI_DOUBLE, topid, tag, comm->Comm(), &request);
  return;
}

/*----------------------------------------------------------------------*
 |  receive anything (public)                                mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Exporter::ReceiveAny(int& source, int& tag, std::vector<char>& recvbuff, int& length)
{
  const Epetra_MpiComm* comm = dynamic_cast<const Epetra_MpiComm*>(&(Comm()));
  if (!comm) dserror("Comm() is not a Epetra_MpiComm\n");
  MPI_Status status;
  // probe for any message to come
  MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm->Comm(), &status);
  // get sender, tag and length
  source = status.MPI_SOURCE;
  tag = status.MPI_TAG;
  MPI_Get_count(&status, MPI_CHAR, &length);
  if (length > (int)recvbuff.size()) recvbuff.resize(length);
  // receive the message
  MPI_Recv(&recvbuff[0], length, MPI_CHAR, source, tag, comm->Comm(), &status);
  return;
}

/*----------------------------------------------------------------------*
 |  receive specific (public)                                mwgee 03/07|
 *----------------------------------------------------------------------*/
void DRT::Exporter::Receive(
    const int source, const int tag, std::vector<char>& recvbuff, int& length)
{
  const Epetra_MpiComm* comm = dynamic_cast<const Epetra_MpiComm*>(&(Comm()));
  if (!comm) dserror("Comm() is not a Epetra_MpiComm\n");
  MPI_Status status;
  // probe for any message to come
  MPI_Probe(source, tag, comm->Comm(), &status);
  MPI_Get_count(&status, MPI_INT, &length);
  if (length > (int)recvbuff.size()) recvbuff.resize(length);
  // receive the message
  MPI_Recv(&recvbuff[0], length, MPI_CHAR, source, tag, comm->Comm(), &status);
  return;
}

/*----------------------------------------------------------------------*
 |  receive anything (public)                                mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Exporter::ReceiveAny(int& source, int& tag, std::vector<int>& recvbuff, int& length)
{
  const Epetra_MpiComm* comm = dynamic_cast<const Epetra_MpiComm*>(&(Comm()));
  if (!comm) dserror("Comm() is not a Epetra_MpiComm\n");
  MPI_Status status;
  // probe for any message to come
  MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm->Comm(), &status);
  // get sender, tag and length
  source = status.MPI_SOURCE;
  tag = status.MPI_TAG;
  MPI_Get_count(&status, MPI_INT, &length);
  if (length > (int)recvbuff.size()) recvbuff.resize(length);
  // receive the message
  MPI_Recv(&recvbuff[0], length, MPI_INT, source, tag, comm->Comm(), &status);
  return;
}

/*----------------------------------------------------------------------*
 |  receive specific (public)                                mwgee 03/07|
 *----------------------------------------------------------------------*/
void DRT::Exporter::Receive(
    const int source, const int tag, std::vector<int>& recvbuff, int& length)
{
  const Epetra_MpiComm* comm = dynamic_cast<const Epetra_MpiComm*>(&(Comm()));
  if (!comm) dserror("Comm() is not a Epetra_MpiComm\n");
  MPI_Status status;
  // probe for any message to come
  MPI_Probe(source, tag, comm->Comm(), &status);
  MPI_Get_count(&status, MPI_INT, &length);
  if (length > (int)recvbuff.size()) recvbuff.resize(length);
  // receive the message
  MPI_Recv(&recvbuff[0], length, MPI_INT, source, tag, comm->Comm(), &status);
  return;
}

/*----------------------------------------------------------------------*
 |  receive anything (public)                                mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Exporter::ReceiveAny(int& source, int& tag, std::vector<double>& recvbuff, int& length)
{
  const Epetra_MpiComm* comm = dynamic_cast<const Epetra_MpiComm*>(&(Comm()));
  if (!comm) dserror("Comm() is not a Epetra_MpiComm\n");
  MPI_Status status;
  // probe for any message to come
  MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm->Comm(), &status);
  // get sender, tag and length
  source = status.MPI_SOURCE;
  tag = status.MPI_TAG;
  MPI_Get_count(&status, MPI_DOUBLE, &length);
  if (length > (int)recvbuff.size()) recvbuff.resize(length);
  // receive the message
  MPI_Recv(&recvbuff[0], length, MPI_DOUBLE, source, tag, comm->Comm(), &status);
  return;
}

/*----------------------------------------------------------------------*
 |  receive specific (public)                                mwgee 03/07|
 *----------------------------------------------------------------------*/
void DRT::Exporter::Receive(
    const int source, const int tag, std::vector<double>& recvbuff, int& length)
{
  const Epetra_MpiComm* comm = dynamic_cast<const Epetra_MpiComm*>(&(Comm()));
  if (!comm) dserror("Comm() is not a Epetra_MpiComm\n");
  MPI_Status status;
  // probe for any message to come
  MPI_Probe(source, tag, comm->Comm(), &status);
  MPI_Get_count(&status, MPI_DOUBLE, &length);
  if (length > (int)recvbuff.size()) recvbuff.resize(length);
  // receive the message
  MPI_Recv(&recvbuff[0], length, MPI_DOUBLE, source, tag, comm->Comm(), &status);
  return;
}

/*----------------------------------------------------------------------*
 |  reduce all (public)                                       umay 10/07|
 *----------------------------------------------------------------------*/
void DRT::Exporter::Allreduce(std::vector<int>& sendbuff, std::vector<int>& recvbuff, MPI_Op mpi_op)
{
  const Epetra_MpiComm* comm = dynamic_cast<const Epetra_MpiComm*>(&(Comm()));
  if (!comm) dserror("Comm() is not a Epetra_MpiComm\n");

  int length = (int)sendbuff.size();
  if (length > (int)recvbuff.size()) recvbuff.resize(length);

  MPI_Allreduce(&(sendbuff[0]), &(recvbuff[0]), length, MPI_INT, mpi_op, comm->Comm());
  return;
}

/*----------------------------------------------------------------------*
 |  the actual exporter constructor (private)                mwgee 05/07|
 *----------------------------------------------------------------------*/
void DRT::Exporter::ConstructExporter()
{
  if (SourceMap().SameAs(TargetMap())) return;

  // allocate a sendplan array and init to zero
  // SendPlan():
  // SendPlan()(lid,proc)    = 1 for data with local id lid needs sending to proc
  // SendPlan()(lid,proc)    = 0 otherwise
  // SendPlan()(lid,MyPID()) = 0 always! (I never send to myself)
  SendPlan().resize(NumProc());

#if 0
  // allocate a receive plan
  // RecvPlan():
  // RecvPlan()(lid,proc) = 1 data with local id lid will be received from proc
  // RecvPlan()(lid,proc) = 0 otherwise
  // RecvPlan()(lid,MyPID()) = 0 always! (I never receive from myself)
  RecvPlan().resize(NumProc());
#endif

  // To build these plans, everybody has to communicate what he has and wants:
  // bundle this info to save on communication:
  int sizes[2];
  sizes[0] = SourceMap().NumMyElements();
  sizes[1] = TargetMap().NumMyElements();
  const int sendsize = sizes[0] + sizes[1];
  std::vector<int> sendbuff;
  sendbuff.reserve(sendsize);
  std::copy(SourceMap().MyGlobalElements(),
      SourceMap().MyGlobalElements() + SourceMap().NumMyElements(), std::back_inserter(sendbuff));
  std::copy(TargetMap().MyGlobalElements(),
      TargetMap().MyGlobalElements() + TargetMap().NumMyElements(), std::back_inserter(sendbuff));

  for (int proc = 0; proc < NumProc(); ++proc)
  {
    int recvsizes[2];
    recvsizes[0] = sizes[0];
    recvsizes[1] = sizes[1];
    Comm().Broadcast(recvsizes, 2, proc);
    const int recvsize = recvsizes[0] + recvsizes[1];
    std::vector<int> recvbuff(recvsize);
    if (proc == MyPID()) std::copy(sendbuff.begin(), sendbuff.end(), &recvbuff[0]);
    Comm().Broadcast(&recvbuff[0], recvsize, proc);
    // const int* have = &recvbuff[0];            // this is what proc has
    const int* want = &recvbuff[recvsizes[0]];  // this is what proc needs

#if 0
    // Loop what proc has and what I need (RecvPlan)
    // (I do not receive from myself)
    if (proc != MyPID())
      for (int i=0; i<recvsizes[0]; ++i)
      {
        const int gid = have[i];
        if (TargetMap().MyGID(gid))
        {
          const int lid = TargetMap().LID(gid);
          RecvPlan()[proc].insert(lid);
        }
      }
#endif

    // Loop what proc wants and what I have (SendPlan)
    if (proc != MyPID())
      for (int i = 0; i < recvsizes[1]; ++i)
      {
        const int gid = want[i];
        if (SourceMap().MyGID(gid))
        {
          const int lid = SourceMap().LID(gid);
          SendPlan()[proc].insert(lid);
        }
      }
    Comm().Barrier();
  }  // for (int proc=0; proc<NumProc(); ++proc)



#if 0
  // make test print of RecvPlan
  for (int proc=0; proc<NumProc(); ++proc)
  {
    if (MyPID()==proc)
    {
      for (int i=0; i<RecvPlan().M(); ++i)
        for (int j=0; j<RecvPlan().N(); ++j)
          if (RecvPlan()(i,j))
            cout << "Proc " << MyPID() << " wants gid "
                 << TargetMap().MyGlobalElements()[i]
                 << " from proc " << j << endl;
    }
    fflush(stdout);
    Comm().Barrier();
  }
  cout << endl;
  // make test print of SendPlan
  for (int proc=0; proc<NumProc(); ++proc)
  {
    if (MyPID()==proc)
    {
      for (int i=0; i<SendPlan().M(); ++i)
        for (int j=0; j<SendPlan().N(); ++j)
          if (SendPlan()(i,j))
          cout << "Proc " << MyPID() << " sends gid "
               << SourceMap().MyGlobalElements()[i]
               << " to proc " << j << endl;
    }
    fflush(stdout);
    Comm().Barrier();
  }
#endif

  return;
}


/*----------------------------------------------------------------------*
 |  communicate objects (public)                             mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Exporter::GenericExport(ExporterHelper& helper)
{
  if (SendPlan().size() == 0) return;
  // if (SourceMap().SameAs(TargetMap())) return;

  helper.PreExportTest(this);

  //------------------------------------------------ do the send/recv loop
  for (int i = 0; i < NumProc() - 1; ++i)
  {
    int tproc = MyPID() + 1 + i;
    int sproc = MyPID() - 1 - i;
    if (tproc < 0) tproc += NumProc();
    if (sproc < 0) sproc += NumProc();
    if (tproc > NumProc() - 1) tproc -= NumProc();
    if (sproc > NumProc() - 1) sproc -= NumProc();
    // cout << "Proc " << MyPID() << " tproc " << tproc << " sproc " << sproc << endl;
    // fflush(stdout);

    //------------------------------------------------ do sending to tproc
    // gather all objects to be send
    DRT::PackBuffer sendblock;
    std::vector<int> sendgid;
    sendgid.reserve(SendPlan()[tproc].size());

    // count

    for (std::set<int>::iterator i = SendPlan()[tproc].begin(); i != SendPlan()[tproc].end(); ++i)
    {
      const int lid = *i;
      const int gid = SourceMap().GID(lid);
      helper.PackObject(gid, sendblock);
    }

    // pack

    sendblock.StartPacking();

    for (std::set<int>::iterator i = SendPlan()[tproc].begin(); i != SendPlan()[tproc].end(); ++i)
    {
      const int lid = *i;
      const int gid = SourceMap().GID(lid);
      if (helper.PackObject(gid, sendblock)) sendgid.push_back(gid);
    }

    // send tproc no. of chars tproc must receive
    std::vector<int> snmessages(2);
    snmessages[0] = sendblock().size();
    snmessages[1] = sendgid.size();

    MPI_Request sizerequest;
    ISend(MyPID(), tproc, &snmessages[0], 2, 1, sizerequest);

    // do the sending of the objects
    MPI_Request sendrequest;
    ISend(MyPID(), tproc, &sendblock()[0], sendblock().size(), 2, sendrequest);

    MPI_Request sendgidrequest;
    ISend(MyPID(), tproc, &sendgid[0], sendgid.size(), 3, sendgidrequest);

    //---------------------------------------- do the receiving from sproc
    // receive how many messages I will receive from sproc
    std::vector<int> rnmessages(2);
    int source = sproc;
    int length = 0;
    int tag = 1;
    // do a blocking specific receive
    Receive(source, tag, rnmessages, length);
    if (length != 2 or tag != 1) dserror("Messages got mixed up");

    // receive the objects
    std::vector<char> recvblock(rnmessages[0]);
    tag = 2;
    ReceiveAny(source, tag, recvblock, length);
    if (tag != 2) dserror("Messages got mixed up");

    // receive the gids
    std::vector<int> recvgid(rnmessages[1]);
    tag = 3;
    ReceiveAny(source, tag, recvgid, length);
    if (tag != 3) dserror("Messages got mixed up");

    std::vector<char>::size_type index = 0;
    int j = 0;
    while (index < recvblock.size())
    {
      int gid = recvgid[j];
      helper.UnpackObject(gid, index, recvblock);
      j += 1;
    }

    //----------------------------------- do waiting for messages to tproc to leave
    Wait(sizerequest);
    Wait(sendrequest);
    Wait(sendgidrequest);

    // make sure we do not get mixed up messages as we use wild card receives here
    Comm().Barrier();
  }  // for (int i=0; i<NumProc()-1; ++i)

  helper.PostExportCleanup(this);

  return;
}


/*----------------------------------------------------------------------*
 |  communicate objects (public)                             mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Exporter::Export(std::map<int, int>& data)
{
  PODExporterHelper<int> helper(data);
  GenericExport(helper);
}


/*----------------------------------------------------------------------*
 |  communicate objects (public)                             mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Exporter::Export(std::map<int, double>& data)
{
  PODExporterHelper<double> helper(data);
  GenericExport(helper);
}


/*----------------------------------------------------------------------*
 |  communicate objects (public)                             u.kue 07/09|
 *----------------------------------------------------------------------*/
void DRT::Exporter::Export(std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>& data)
{
  AnyObjectExporterHelper<Epetra_SerialDenseMatrix> helper(data);
  GenericExport(helper);
}
