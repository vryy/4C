/*!----------------------------------------------------------------------
\file exporter.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "vector"
#include "drt_exporter.H"
#include "drt_utils.H"
#include "drt_dserror.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Exporter::Exporter(const Epetra_Map& frommap,const Epetra_Map& tomap, 
                        const Epetra_Comm& comm) :
frommap_(frommap),
tomap_(tomap),
comm_(comm),
myrank_(comm.MyPID()),
numproc_(comm.NumProc())
{
  // allocate a sendplan array and init to zero
  // SendPlan():
  // SendPlan()(lid,proc)    = 1 for data with local id lid needs sending to proc
  // SendPlan()(lid,proc)    = 0 otherwise
  // SendPlan()(lid,MyPID()) = 0 always! (I never send to myself)
  SendPlan().Shape(SourceMap().NumMyElements(),NumProc());
  
  // allocate a receive plan
  // RecvPlan():
  // RecvPlan()(lid,proc) = 1 data with local id lid will be received from proc
  // RecvPlan()(lid,proc) = 0 otherwise
  // RecvPlan()(lid,MyPID()) = 0 always! (I never receive from myself)
  RecvPlan().Shape(TargetMap().NumMyElements(),NumProc());
  
  // allocate a send buffer for ParObject packs
  SendBuff().resize(SourceMap().NumMyElements());
  SendSize().resize(SourceMap().NumMyElements());
  for (int i=0; i<(int)SendSize().size(); ++i) SendSize()[i] = 0;
  
  // To build these plans, everybody has to communicate what he has and wants:
  // bundle this info to save on communication:
  int sizes[2];
  sizes[0] = SourceMap().NumMyElements(); 
  sizes[1] = TargetMap().NumMyElements(); 
  const int sendsize = sizes[0]+sizes[1];
  vector<int> sendbuff(sendsize);
  for (int i=0; i<sizes[0]; ++i)
    sendbuff[i] = SourceMap().MyGlobalElements()[i]; 
  for (int i=0; i<sizes[1]; ++i)
    sendbuff[i+sizes[0]] = TargetMap().MyGlobalElements()[i]; 
  
  for (int proc=0; proc<NumProc(); ++proc)
  {
    int recvsizes[2];
    recvsizes[0] = sizes[0];
    recvsizes[1] = sizes[1];
    Comm().Broadcast(recvsizes,2,proc);
    const int recvsize = recvsizes[0]+recvsizes[1];
    vector<int> recvbuff(recvsize);
    if (proc==MyPID()) 
      for (int i=0; i<recvsize; ++i) recvbuff[i] = sendbuff[i];
    Comm().Broadcast(&recvbuff[0],recvsize,proc);
    const int* have = &recvbuff[0];            // this is what proc has
    const int* want = &recvbuff[recvsizes[0]]; // this is what proc needs
    
    // Loop what proc has and what I need (RecvPlan)
    // (I do not receive from myself)
    if (proc != MyPID())
      for (int i=0; i<recvsizes[0]; ++i)
      {
        const int gid = have[i];
        if (TargetMap().MyGID(gid))
        {
          const int lid = TargetMap().LID(gid);
          RecvPlan()(lid,proc) = 1;
        }
      }
    
    // Loop what proc wants and what I have (SendPlan)
    if (proc != MyPID())
      for (int i=0; i<recvsizes[1]; ++i)
      {
        const int gid = want[i];
        if (SourceMap().MyGID(gid))
        {
          const int lid = SourceMap().LID(gid);
          SendPlan()(lid,proc) = 1;
        }
      }
    Comm().Barrier();
  } // for (int proc=0; proc<NumProc(); ++proc)



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
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Exporter::Exporter(const DRT::Exporter& old) :
frommap_(old.frommap_),
tomap_(old.tomap_),
comm_(old.comm_),
myrank_(old.myrank_),
numproc_(old.numproc_),
sendplan_(old.sendplan_),
recvplan_(old.recvplan_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Exporter::~Exporter()
{
  return;
}


#ifdef PARALLEL
/*----------------------------------------------------------------------*
 |  do a send of data (public)                               mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Exporter::ISend(const int frompid, 
                          const int topid, 
                          const char* data, 
                          const int dsize,
                          const int tag, 
                          MPI_Request& request)
{
  if (MyPID()!=frompid) return;
  const Epetra_MpiComm* comm = dynamic_cast<const Epetra_MpiComm*>(&(Comm()));
  if (!comm) dserror("Comm() is not a Epetra_MpiComm\n");
  MPI_Isend((void*)data,dsize,MPI_CHAR,topid,tag,comm->Comm(),&request);
  return;
}
#endif  

#ifdef PARALLEL
/*----------------------------------------------------------------------*
 |  receive anything (public)                                mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Exporter::ReceiveAny(int& source, int&tag, 
                               vector<char>& recvbuff,int& length)
{
  const Epetra_MpiComm* comm = dynamic_cast<const Epetra_MpiComm*>(&(Comm()));
  if (!comm) dserror("Comm() is not a Epetra_MpiComm\n");
  MPI_Status status;
  // probe for any message to come
  MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,comm->Comm(),&status);
  // get sender, tag and length
  source = status.MPI_SOURCE;
  tag    = status.MPI_TAG;
  MPI_Get_count(&status,MPI_CHAR,&length);
  if (length>(int)recvbuff.size()) recvbuff.resize(length);
  // receive the message
  MPI_Recv(&recvbuff[0],length,MPI_CHAR,source,tag,comm->Comm(),&status);
  return;
}
#endif  


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
