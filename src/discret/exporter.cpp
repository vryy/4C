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

#include "exporter.H"
#include "dserror.H"
#include "vector"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Exporter::Exporter(const Epetra_Map& frommap, 
                                      const Epetra_Map& tomap, 
                                      const Epetra_Comm& comm) :
frommap_(frommap),
tomap_(tomap),
comm_(comm),
myrank_(comm.MyPID()),
numproc_(comm.NumProc())
{
  // allocate the send plan vectors
  SendPID().resize(SourceMap().NumMyElements());
  SendLID().resize(SourceMap().NumMyElements());
  
  // allocate the recv plan vectors
  RecvPID().resize(TargetMap().NumMyElements());
  RecvLID().resize(TargetMap().NumMyElements());
  
  int ierr=0;
  // fill the send information
  ierr = TargetMap().RemoteIDList(SourceMap().NumMyElements(),
                                  SourceMap().MyGlobalElements(),
                                  &SendPID()[0],
                                  &SendLID()[0]);
  if (ierr) dserror("TargetMap().RemoteIDList returned err=%d",ierr);

                           
  // fill the receive information
  ierr = SourceMap().RemoteIDList(TargetMap().NumMyElements(),
                                  TargetMap().MyGlobalElements(),
                                  &RecvPID()[0],
                                  &RecvLID()[0]);                           
  if (ierr) dserror("SourceMap().RemoteIDList returned err=%d",ierr);
  
  cout << "Proc " << MyPID() << " SourceMap().NumMyElements() " << SourceMap().NumMyElements() << endl;
  cout << "Proc " << MyPID() << " TargetMap().NumMyElements() " << TargetMap().NumMyElements() << endl;
  
  for (int proc=0; proc<NumProc(); ++proc)
  {
    if (proc==MyPID())
    {
      cout << "Proc " << proc << endl
           << "SendPID()\n";
      for (int i=0; i<SendPID().size(); ++i)
        cout << SendPID()[i] << endl;
      cout << "RecvPID()\n";
      for (int i=0; i<RecvPID().size(); ++i)
        cout << RecvPID()[i] << endl;
    }
    Comm().Barrier();
  }

#if 0
  // allocate a sendplan array and init to zero
  // SendPlan():
  // SendPlan()(lid,proc) = 1 for data with local id lid needs sending to proc
  // SendPlan()(lid,proc) = 0 otherwise
  sendplan_.Shape(SourceMap().NumMyElements(),NumProc());
  
  // allocate a receive plan and init to MyPID()
  // RecvPlan():
  // RecvPlan()[lid] = proc data with local id lid will be received from proc proc
  // RecvPlan()[lid] = myrank data with local id lid is already owned by myself
  recvplan_.resize(TargetMap().NumMyElements());
  for (int i=0; i<TargetMap().NumMyElements(); ++i) RecvPlan()[i] = MyPID();
  
  // To build these plans, everybody has to communicate what she has and wants:
  // we bundle this info to save on communication:
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
  }
#endif


  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Exporter::Exporter(const CCADISCRETIZATION::Exporter& old) :
frommap_(old.frommap_),
tomap_(old.tomap_),
comm_(old.comm_),
myrank_(old.myrank_),
numproc_(old.numproc_),
sendpid_(old.sendpid_),
sendlid_(old.sendlid_),
recvpid_(old.recvpid_),
recvlid_(old.recvlid_)                                      
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::Exporter::~Exporter()
{
  return;
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
