/*!----------------------------------------------------------------------
\file drt_exporter.H

\class DRT::Exporter

\brief A class to manage explicit mpi communications

The discretization management module uses this class to do most of its
communication. It is used to redistribute grids and to do point-to-point
communications. It is therefore the only place on DRT where explicit calls to MPI
methods are done.<br>
It has strong capabilities in gathering and scattering information in a collective AND
an individual way. Whenever you need explicit communication, check this class first before
implementing your own mpi stuff.

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef EXPORTER_H
#define EXPORTER_H

#include "vector"
#include "map"
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_SerialDenseMatrix.h"
#include "Teuchos_RefCountPtr.hpp"

#include "drt_node.H"
#include "drt_utils.H"
#include "drt_parobject.H"

using namespace std;
using namespace Teuchos;

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{

/*!
\class Exporter

\brief A class to manage explicit mpi communications

The discretization management module uses this class to do most of its
communication. It is used to redistribute grids and to do point-to-point
communications. It is therefore the only place on DRT where explicit calls to MPI
methods are done.<br>
It has strong capabilities in gathering and scattering information in a collective AND
an individual way. Whenever you need explicit communication, check this class first before
implementing your own mpi stuff.

\author gee (gee@lnm.mw.tum.de)
*/
class Exporter
{

public:



  /*!
  \brief Standard Constructor

  this ctor constructs an exporter with no maps. It can than be used to do
  point-to-point communication only, map based exportes are not possible!

  \param comm    (in): Communicator that shall be used in exports
  */
  Exporter(const Epetra_Comm& comm);

  /*!
  \brief Standard Constructor

  \param frommap (in): The source map data shall be exported from
  \param tomap   (in): The target map data shall be exported to
  \param comm    (in): Communicator that shall be used in exports
  */
  Exporter(const Epetra_Map& frommap, const Epetra_Map& tomap, const Epetra_Comm& comm);

  /*!
  \brief Copy Constructor

  */
  Exporter(const DRT::Exporter& old);

  /*!
  \brief Destructor

  */
  virtual ~Exporter();


  //! @name Acess methods

  /*!
  \brief Get communicator
  */
  inline const Epetra_Comm& Comm() const { return comm_; }

  /*!
  \brief Get source map
  */
  inline const Epetra_Map& SourceMap() const { return frommap_; }

  /*!
  \brief Get target map
  */
  inline const Epetra_Map& TargetMap() const { return tomap_; }

  //@}

  //! @name Communication methods

  /*!
  \brief Communicate a map of objects that implement ParObject

  This method takes a map of objects and redistributes them according to
  the send and receive plans. It is implicitly assumed, that the key in
  the map of objects pointwise matches SourceMap(). It is also assumed
  (and tested), that type T implements the ParObject class.

  \param parobjects (in/out): A map of classes T that implement the
                              class ParObject. On input, the map
                              has a distribution matching SourceMap().
                              On output, the map has a distribution of
                              TargetMap().
  */
  template<typename T> void Export(map<int,RefCountPtr<T> >& parobjects);

  /*!
  \brief Communicate a map of vectors of some basic data type T

  This method takes a map of vectors and redistributes them according to
  the send and receive plans. It is implicitly assumed, that the key in
  the map of vectors pointwise matches SourceMap().

  \note T can be int, double or char. The method will not compile will other
        then these basic data types (and will give a rater kryptic error message)

  \param data (in/out): A map of vectors<T>. On input, the map
                        has a distribution matching SourceMap().
                        On output, the map has a distribution of
                        TargetMap().
  */
  template<typename T> void Export(map<int,vector<T> >& data);



#ifdef PARALLEL
  /*!
  \brief Send data from one processor to another (nonblocking)

  The method will send an array of chars in a nonblocking way meaning that this method
  will return immediately on the calling processor - even if communcation has not
  finished yet. The char array must not be altered or destroyed as long as the
  communication might still be in progress. This can be tested for using DRT::Exporter::Wait
  and the request handle returned.
  The receiving processor should call DRT::Exporter::ReceiveAny to receive the message.
  Note that messages from one explicit proc to another explicit proc are non-overtaking meaning
  they will arrive in the order they have been sent.

  \note This is an individual call.

  \param frompid (in)  : sending processors' pid
  \param topid (in)    : target processors' pid
  \param data (in)     : ptr to data to be send
  \param dsize (in)    : size of data (no. of chars)
  \param tag (in)      : tag to be used with message
  \param request (out) : mpi request handle to be used for testing completion of the
                         communication. data may not be altered or destroyed before
                         communcation finalized! One can use DRT::Exporter::Wait for this.

  \note This is an individual call
  */
  void ISend(const int frompid, const int topid, const char* data,
            const int dsize, const int tag, MPI_Request& request);

  /*!
  \brief Send data from one processor to another (nonblocking)

  The method will send an array of ints in a nonblocking way meaning that this method
  will return immediately on the calling processor - even if communcation has not
  finished yet. The int array must not be altered or destroyed as long as the
  communication might still be in progress. This can be tested for using DRT::Exporter::Wait
  and the request handle returned.
  The receiving processor should call DRT::Exporter::ReceiveAny to receive the message.
  Note that messages from one explicit proc to another explicit proc are non-overtaking meaning
  they will arrive in the order they have been sent.

  \note This is an individual call.

  \param frompid (in)  : sending processors' pid
  \param topid (in)    : target processors' pid
  \param data (in)     : ptr to data to be send
  \param dsize (in)    : size of data (no. of integers)
  \param tag (in)      : tag to be used with message
  \param request (out) : mpi request handle to be used for testing completion of the
                         communication. data may not be altered or destroyed before
                         communcation finalized! One can use DRT::Exporter::Wait for this.

  \note This is an individual call
  */
  void ISend(const int frompid, const int topid, const int* data,
            const int dsize, const int tag, MPI_Request& request);

  /*!
  \brief Send data from one processor to another (nonblocking)

  The method will send an array of doubles in a nonblocking way meaning that this method
  will return immediately on the calling processor - even if communcation has not
  finished yet. The double array must not be altered or destroyed as long as the
  communication might still be in progress. This can be tested for using DRT::Exporter::Wait
  and the request handle returned.
  The receiving processor should call DRT::Exporter::ReceiveAny to receive the message.
  Note that messages from one explicit proc to another explicit proc are non-overtaking meaning
  they will arrive in the order they have been sent.

  \note This is an individual call.

  \param frompid (in)  : sending processors' pid
  \param topid (in)    : target processors' pid
  \param data (in)     : ptr to data to be send
  \param dsize (in)    : size of data (no. of doubles)
  \param tag (in)      : tag to be used with message
  \param request (out) : mpi request handle to be used for testing completion of the
                         communication. data may not be altered or destroyed before
                         communcation finalized! One can use DRT::Exporter::Wait for this.

  \note This is an individual call
  */
  void ISend(const int frompid, const int topid, const double* data,
            const int dsize, const int tag, MPI_Request& request);

  /*!
  \brief Receive anything joker (blocking)

  This method receives an MPI_CHAR string message from any source proc with
  any message tag of any length. It simply takes the first message that's
  coming in no matter from which sender of with which tag.
  recvbuff is resized to fit received message.
  the method is blocking for the calling (receiving) proc but for none of the
  other processors.
  It is used together with ISend and Wait to do nonblocking chaotic
  point to point communication.

  \note This is an individual call.

  \warning There is absolutely no guarantee about the order messages are
           received with this method except for one: Messages from the SAME sender
           to the SAME receiver will not overtake each other (which is not a really strong
           statement).

  \param source (output): source the message came from
  \param tag (output): message tag of message received
  \param recvbuff (output): buffer containing received data
  \param length (output): length of message upon receive
  */
  void ReceiveAny(int& source, int& tag, vector<char>& recvbuff, int& length);
  void Receive(const int source,const int tag, vector<char>& recvbuff,int& length);

  /*!
  \brief Receive anything joker (blocking)

  This method receives an MPI_INT message from any source proc with
  any message tag of any length. It simply takes the first message that's
  coming in no matter from which sender of with which tag.
  recvbuff is resized to fit received message.
  the method is blocking for the calling (receiving) proc but for none of the
  other processors.
  It is used together with ISend and Wait to do nonblocking chaotic
  point to point communication.

  \note This is an individual call.

  \warning There is absolutely no guarantee about the order messages are
           received with this method except for one: Messages from the SAME sender
           to the SAME receiver will not overtake each other (which is not a really strong
           statement).

  \param source (output): source the message came from
  \param tag (output): message tag of message received
  \param recvbuff (output): buffer containing received data
  \param length (output): length of message upon receive
  */
  void ReceiveAny(int& source, int& tag, vector<int>& recvbuff, int& length);
  void Receive(const int source,const int tag, vector<int>& recvbuff, int& length);

  /*!
  \brief Receive anything joker (blocking)

  This method receives an MPI_DOUBLE message from any source proc with
  any message tag of any length. It simply takes the first message that's
  coming in no matter from which sender of with which tag.
  recvbuff is resized to fit received message.
  the method is blocking for the calling (receiving) proc but for none of the
  other processors.
  It is used together with ISend and Wait to do nonblocking chaotic
  point to point communication.

  \note This is an individual call.

  \warning There is absolutely no guarantee about the order messages are
           received with this method except for one: Messages from the SAME sender
           to the SAME receiver will not overtake each other (which is not a really strong
           statement).

  \param source (output): source the message came from
  \param tag (output): message tag of message received
  \param recvbuff (output): buffer containing received data
  \param length (output): length of message upon receive
  */
  void ReceiveAny(int& source, int& tag, vector<double>& recvbuff, int& length);

  /*!
  \brief wait for nonblocking send to finish

  The method is used together with Isend and ReceiveAny to guarantee finalization
  of a communication. It is an individual call done by the sending processor to guarantee
  that message was taken from the sendbuffer before destroying the sendbuffer.
  This method is blocking and will return one communication associated with request has
  left the sender.

  \param request (in): mpi request handle

  */
  
  
  /*!
  \brief performs an allreduce operation on all processors 
  		 and sends the result to all processors

  \param sendbuff (input): buffer containing data that has to be sent  
  \param recvbuff (output): buffer containing received data
  \param mpi_op   (input): MPI operation
  */
  void Allreduce(vector<int>& sendbuff, vector<int>& recvbuff, MPI_Op mpi_op);
  
  
  void Wait(MPI_Request& request) { MPI_Status status; MPI_Wait(&request,&status); return;}


#endif

  //@}

private:

  /*!
  \brief Do initialization of the exporter
  */
  void ConstructExporter();

  /*!
  \brief Get PID
  */
  inline int MyPID() const { return myrank_; }
  /*!
  \brief Get no. of processors
  */
  inline int NumProc() const { return numproc_; }

  /*!
  \brief Get sendplan_
  */
  inline Epetra_SerialDenseMatrix& SendPlan() { return sendplan_; }

  /*!
  \brief Get recvplan_
  */
  inline Epetra_SerialDenseMatrix& RecvPlan() { return recvplan_; }

  /*!
  \brief Get sendbuff_
  */
  inline vector<vector<char> >& SendBuff() { return sendbuff_; }

  /*!
  \brief Get sendsize_
  */
  inline vector<int>& SendSize() { return sendsize_; }


private:

  //! dummy map in case of empty exporter
  Epetra_Map               dummymap_;
  //! source layout
  const Epetra_Map&        frommap_;
  //! target map
  const Epetra_Map&        tomap_;
  //! communicator
  const Epetra_Comm&       comm_;
  //! PID
  int                      myrank_;
  //! no. of processors
  int                      numproc_;
  //! sending information
  Epetra_SerialDenseMatrix sendplan_;
  //! receiving information
  Epetra_SerialDenseMatrix recvplan_;
  //! a sendbuffer
  vector<vector<char> >    sendbuff_;
  //! sendsize_[i] is length of vector sendbuff_[i]
  vector<int>              sendsize_;


}; // class Exporter
} // namespace DRT



/*----------------------------------------------------------------------*
 |  communicate objects (public)                             mwgee 11/06|
 *----------------------------------------------------------------------*/
template<typename T> void DRT::Exporter::Export(
                                     map<int,RefCountPtr<T> >& parobjects)
{
  if (SendPlan().N()==0) return;
  //if (SourceMap().SameAs(TargetMap())) return;

  // test whether type T implements ParObject
  {
    typename map<int,RefCountPtr<T> >::iterator curr = parobjects.begin();
    if (curr != parobjects.end())
    {
      T* ptr = curr->second.get();
      ParObject* tester = dynamic_cast<ParObject*>(ptr);
      if (!tester) dserror("typename T in template does not implement class ParObject (dynamic_cast failed)");
    }
  }


#ifdef PARALLEL
  // allocate requests for unknown number of sends
  vector<MPI_Request> request(100);
  MPI_Request         sizerequest;

  //------------------------------------------------ do the send/recv loop
  for (int i=0; i<NumProc()-1; ++i)
  {
    int countsend=0; // count how many sends to tproc
    int tproc = MyPID()+1+i;
    int sproc = MyPID()-1-i;
    if (tproc<0) tproc += NumProc();
    if (sproc<0) sproc += NumProc();
    if (tproc>NumProc()-1) tproc -= NumProc();
    if (sproc>NumProc()-1) sproc -= NumProc();
    //cout << "Proc " << MyPID() << " tproc " << tproc << " sproc " << sproc << endl;
    //fflush(stdout);

    //------------------------------------------------ do sending to tproc
    // count how many objects will actually be send to tproc
    int nsend=0;
    for (int i=0; i<SourceMap().NumMyElements(); ++i)
    {
      const int lid = i;
      if (SendPlan()(lid,tproc)!=1) continue;
      const int gid = SourceMap().MyGlobalElements()[lid];
      typename map<int,RefCountPtr<T> >::iterator curr = parobjects.find(gid);
      if (curr==parobjects.end()) continue;
      nsend++;
    }

    // send tproc no. of messages tproc must receive
    vector<int> snmessages(1);
    snmessages[0] = nsend;
    ISend(MyPID(),tproc,&snmessages[0],1,1,sizerequest);

    // do the sending of the objects
    for (int i=0; i<SourceMap().NumMyElements(); ++i)
    {
      const int lid = i;
      if (SendPlan()(lid,tproc)!=1) continue;
      const int gid = SourceMap().MyGlobalElements()[lid];
      typename map<int,RefCountPtr<T> >::iterator curr = parobjects.find(gid);
      if (curr==parobjects.end()) continue;
      RefCountPtr<T> actobject = curr->second;
      // pack the stuff
      actobject->Pack(SendBuff()[lid]);
      SendSize()[lid] = SendBuff()[lid].size();
      if (countsend>=(int)request.size()) request.resize(request.size()+500);
      ISend(MyPID(),tproc,&SendBuff()[lid][0],SendSize()[lid],gid,request[countsend]);
      ++countsend;
    }
    if (countsend!=nsend) dserror("No. of send messages wrong");

    //---------------------------------------- do the receiving from sproc
    // receive how many messages I will receive from sproc
    vector<int> rnmessages(1);
    int source = sproc;
    int tag = 1;
    int length = 0;
    // do a blocking specific receive
    Receive(source,tag,rnmessages,length);
    if (length!=1) dserror("Messages got mixed up");
    int nrecv = rnmessages[0];

    // receive the objects
    vector<char> recvbuff(500);
    for (int i=0; i<nrecv; ++i)
    {
      int source=-1;
      int gid=-1;
      int length=0;
      ReceiveAny(source,gid,recvbuff,length);
      if (source!=sproc) dserror("Messages got mixed up");
      recvbuff.resize(length);
      if (!TargetMap().MyGID(gid)) dserror("Received object with gid that I did not want");
      // check whether I already have this object (might happen)
      // in this case do nothing
      typename map<int,RefCountPtr<T> >::iterator curr = parobjects.find(gid);
      if (curr != parobjects.end()) continue;
      // create object and unpack
      DRT::ParObject* object = DRT::Utils::Factory(recvbuff);
      T* ptr = dynamic_cast<T*>(object);
      if (!ptr) dserror("typename T in template does not implement ParObject (dynamic_cast failed)");
      RefCountPtr<T> refptr = rcp(ptr);
      // add object to my map
      parobjects[gid] = refptr;
    }

    //----------------------------------- do waiting for messages to tproc to leave
    Wait(sizerequest);
    for (int i=0; i<countsend; ++i)
    {
      Wait(request[i]);
    }
    for (int i=0; i<SourceMap().NumMyElements(); ++i)
    {
      SendBuff()[i].clear();
      SendSize()[i] = 0;
    }

    // make sure we do not get mixed up messages as we use wild card receives here
    Comm().Barrier();
  } // for (int i=0; i<NumProc()-1; ++i)

  // loop map and kick out everything that's not in TargetMap()
  map<int,RefCountPtr<T> > newmap;
  typename map<int,RefCountPtr<T> >::iterator fool;
  for (fool=parobjects.begin(); fool != parobjects.end(); ++fool)
    if (TargetMap().MyGID(fool->first))
      newmap[fool->first] = fool->second;
  parobjects.clear();
  for (fool=newmap.begin(); fool != newmap.end(); ++fool)
    parobjects[fool->first] = fool->second;

#endif
  return;
}

#if 0 // old version!
/*----------------------------------------------------------------------*
 |  communicate objects (public)                             mwgee 11/06|
 *----------------------------------------------------------------------*/
template<typename T> void DRT::Exporter::Export(
                                     map<int,RefCountPtr<T> >& parobjects)
{
  if (SendPlan().N()==0) return;
  //if (SourceMap().SameAs(TargetMap())) return;

  // test whether type T implements ParObject
  {
    typename map<int,RefCountPtr<T> >::iterator curr = parobjects.begin();
    if (curr != parobjects.end())
    {
      T* ptr = curr->second.get();
      ParObject* tester = dynamic_cast<ParObject*>(ptr);
      if (!tester) dserror("typename T in template does not implement class ParObject (dynamic_cast failed)");
    }
  }


#ifdef PARALLEL
  // allocate requests for unknown number of sends
  vector<MPI_Request> request(200);
  int nsend=0;

  //-------------------------------------------------------- do the sending
  for (int i=0; i<SourceMap().NumMyElements(); ++i)
  {
    // get the global id
    const int gid = SourceMap().MyGlobalElements()[i];
    const int lid = i;
    // check whether there will be any send to do at all
    bool issend = false;
    for (int j=0; j<NumProc(); ++j)
      if (MyPID() != j)
        if (SendPlan()(lid,j)==1)
        {
          issend = true;
          break;
        }
    if (!issend) continue;


    // get the object to send
    typename map<int,RefCountPtr<T> >::iterator curr = parobjects.find(gid);
    if (curr==parobjects.end()) dserror("Cannot find object with gid %d",gid);
    RefCountPtr<T> actobject = curr->second;
    // pack the stuff
    actobject->Pack(SendBuff()[lid]);
    SendSize()[lid] = SendBuff()[lid].size();
    // do sending
    for (int j=0; j<NumProc(); ++j)
      if (j != MyPID() && SendPlan()(lid,j)==1)
      {
        if (nsend>=(int)request.size())
          request.resize(request.size()+200);
        ISend(MyPID(),j,&SendBuff()[lid][0],SendSize()[lid],gid,request[nsend]);
        ++nsend;
      }
  }

  //-------------------------------------------------------- do receiving
  // count how many receives I want to do
  int nrecv=0;
  for (int i=0; i<RecvPlan().M(); ++i)
    for (int j=0; j<RecvPlan().N(); ++j)
      if (RecvPlan()(i,j)) ++nrecv;
  vector<char> recvbuff(500);
  while (nrecv)
  {
    int source = -1;
    int gid = -1;
    int length = 0;
    ReceiveAny(source,gid,recvbuff,length);
    recvbuff.resize(length);
    --nrecv;
    // check whether message was for me
    if (!TargetMap().MyGID(gid))
      dserror("Received object with gid that I did not want");
    // check whether I already have this object
    // This can happen if I've received it before from someone else
    // In this case, do nothing (keep what I have)
    typename map<int,RefCountPtr<T> >::iterator curr = parobjects.find(gid);
    if (curr != parobjects.end())
      continue;
    // Create an empty object and unpack
    DRT::ParObject* object = DRT::Utils::Factory(recvbuff);
    T* ptr = dynamic_cast<T*>(object);
    if (!ptr) dserror("typename T in template does not implement ParObject (dynamic_cast failed)");
    RefCountPtr<T> refptr = rcp(ptr);
    // add object to the map
    parobjects[gid] = refptr;
  }

  //--------- loop map and kick out everything that's not in TargetMap()
  // (this is done strangely here because of a mysterious bug when
  //  deleting from the parobjects map)
  map<int,RefCountPtr<T> > newmap;
  typename map<int,RefCountPtr<T> >::iterator fool;
  for (fool=parobjects.begin(); fool != parobjects.end(); ++fool)
    if (TargetMap().MyGID(fool->first))
      newmap[fool->first] = fool->second;
  parobjects.clear();
  for (fool=newmap.begin(); fool != newmap.end(); ++fool)
    parobjects[fool->first] = fool->second;

  //---------------------------------------------------------- do waiting
  for (int i=0; i<nsend; ++i) Wait(request[i]);

  //------------------------------------------------------ free sendbuffer
  for (int i=0; i<(int)SendBuff().size(); ++i)
  {
    SendBuff()[i].clear();
    SendSize()[i] = 0;
  }
#endif
  return;
}
#endif // old version!


//static void Huhu() { cout << "Huhu\n"; fflush(stdout); return; }
/*----------------------------------------------------------------------*
 |  communicate objects (public)                             mwgee 11/06|
 *----------------------------------------------------------------------*/
template<typename T> void DRT::Exporter::Export(map<int,vector<T> >& data)
{
  if (SendPlan().N()==0) return;
  //if (SourceMap().SameAs(TargetMap())) return;
#ifdef PARALLEL
  // allocate requests for unknown number of sends
  vector<MPI_Request> request(100);
  MPI_Request         sizerequest;

  // allocate a sendbuffer
  vector<T*> sendbuff(SourceMap().NumMyElements());
  //Huhu();

  //------------------------------------------------ do the send/recv loop
  for (int i=0; i<NumProc()-1; ++i)
  {
    int countsend=0; // count how many sends to tproc
    int tproc = MyPID()+1+i;
    int sproc = MyPID()-1-i;
    if (tproc<0) tproc += NumProc();
    if (sproc<0) sproc += NumProc();
    if (tproc>NumProc()-1) tproc -= NumProc();
    if (sproc>NumProc()-1) sproc -= NumProc();

    //------------------------------------------------ do sending to tproc
    // count how many objects will actually be send to tproc
    int nsend=0;
    for (int i=0; i<SourceMap().NumMyElements(); ++i)
    {
      const int lid = i;
      if (SendPlan()(lid,tproc)!=1) continue;
      const int gid = SourceMap().MyGlobalElements()[lid];
      typename map<int,vector<T> >::iterator curr = data.find(gid);
      if (curr==data.end()) continue;
      nsend++;
    }

    // send tproc no. of messages tproc must receive
    vector<int> snmessages(1);
    snmessages[0] = nsend;
    ISend(MyPID(),tproc,&snmessages[0],1,1,sizerequest);

    // do the sending of the objects
    for (int i=0; i<SourceMap().NumMyElements(); ++i)
    {
      const int lid = i;
      if (SendPlan()(lid,tproc)!=1) continue;
      const int gid = SourceMap().MyGlobalElements()[lid];
      typename map<int,vector<T> >::iterator curr = data.find(gid);
      if (curr==data.end()) continue;
      vector<T>& actobject = curr->second;
      sendbuff[lid]   = &actobject[0];
      SendSize()[lid] = (int)actobject.size();
      if (countsend>=(int)request.size()) request.resize(request.size()+500);
      ISend(MyPID(),tproc,sendbuff[lid],SendSize()[lid],gid,request[countsend]);
      ++countsend;
    }
    if (countsend!=nsend) dserror("No. of send messages wrong");

    //---------------------------------------- do the receiving from sproc
    // receive how many messages I will receive from sproc
    vector<int> rnmessages(1);
    int source = sproc;
    int tag = 1;
    int length = 0;
    // do a blocking specific receive
    Receive(source,tag,rnmessages,length);
    if (length!=1) dserror("Messages got mixed up");
    int nrecv = rnmessages[0];

    // receive the objects
    vector<T> recvbuff;
    for (int i=0; i<nrecv; ++i)
    {
      int source=-1;
      int gid=-1;
      int length=0;
      ReceiveAny(source,gid,recvbuff,length);
      if (source!=sproc) dserror("Messages got mixed up");
      recvbuff.resize(length);
      if (!TargetMap().MyGID(gid)) dserror("Received object with gid that I did not want");
      // check whether I already have this object (might happen)
      // in this case do nothing
      typename map<int,vector<T> >::iterator curr = data.find(gid);
      if (curr != data.end()) continue;
      data[gid] = recvbuff;
    }

    //----------------------------------- do waiting for messages to tproc to leave
    Wait(sizerequest);
    for (int i=0; i<countsend; ++i)
      Wait(request[i]);

    // make sure we do not get mixed up messages as we use wild card receives here
    Comm().Barrier();
  } // for (int i=0; i<NumProc()-1; ++i)

  //--------- loop map and kick out everything that's not in TargetMap()
  // (this is done strangely here because of a mysterious bug when
  //  deleting from the parobjects map)
  map<int,vector<T> > newmap;
  typename map<int,vector<T> >::iterator fool;
  for (fool=data.begin(); fool != data.end(); ++fool)
    if (TargetMap().MyGID(fool->first))
      newmap[fool->first] = fool->second;
  data.clear();
  for (fool=newmap.begin(); fool != newmap.end(); ++fool)
    data[fool->first] = fool->second;

#endif
  return;
}

#if 0 // old version!
/*----------------------------------------------------------------------*
 |  communicate objects (public)                             mwgee 11/06|
 *----------------------------------------------------------------------*/
template<typename T> void DRT::Exporter::Export(map<int,vector<T> >& data)
{
#ifdef PARALLEL
  // allocate requests for unknown number of sends
  vector<MPI_Request> request(200);
  int nsend=0;

  if (SendPlan().N()==0) return;
  //if (SourceMap().SameAs(TargetMap())) return;

  // allocate a sendbuffer
  vector<T*> sendbuff(SourceMap().NumMyElements());

  //-------------------------------------------------------- do the sending
  for (int i=0; i<SourceMap().NumMyElements(); ++i)
  {
    // get the global id
    const int gid = SourceMap().MyGlobalElements()[i];
    const int lid = i;
    // check whether there will be any send to do at all
    bool issend = false;
    for (int j=0; j<NumProc(); ++j)
      if (MyPID() != j)
        if (SendPlan()(lid,j)==1)
        {
          issend = true;
          break;
        }
    if (!issend) continue;


    // get the object to send
    typename map<int,vector<T> >::iterator curr = data.find(gid);
    if (curr==data.end()) dserror("Cannot find object with gid %d",gid);
    vector<T>& actobject = curr->second;
    // pack the stuff
    sendbuff[lid]   = &actobject[0];
    SendSize()[lid] = (int)actobject.size();
    // do sending
    for (int j=0; j<NumProc(); ++j)
      if (j != MyPID() && SendPlan()(lid,j)==1)
      {
        if (nsend>=(int)request.size())
          request.resize(request.size()+200);
        ISend(MyPID(),j,sendbuff[lid],SendSize()[lid],gid,request[nsend]);
        ++nsend;
      }
  }


  //-------------------------------------------------------- do receiving
  // count how many receives I want to do
  int nrecv=0;
  for (int i=0; i<RecvPlan().M(); ++i)
    for (int j=0; j<RecvPlan().N(); ++j)
      if (RecvPlan()(i,j)) ++nrecv;
  vector<T> recvbuff;
  while (nrecv)
  {
    int source = -1;
    int gid = -1;
    int length = 0;
    ReceiveAny(source,gid,recvbuff,length);
    --nrecv;
    // check whether message was for me
    if (!TargetMap().MyGID(gid))
      dserror("Received object with gid that I did not want");
    // check whether I already have this object
    // This can happen if I've received it before from someone else
    // In this case, do nothing (keep what I have)
    typename map<int,vector<T> >::iterator curr = data.find(gid);
    if (curr != data.end())
      continue;
    data[gid].resize(length);
    for (int i=0; i<length; ++i) data[gid][i] = recvbuff[i];
  }

  //--------- loop map and kick out everything that's not in TargetMap()
  // (this is done strangely here because of a mysterious bug when
  //  deleting from the parobjects map)
  map<int,vector<T> > newmap;
  typename map<int,vector<T> >::iterator fool;
  for (fool=data.begin(); fool != data.end(); ++fool)
    if (TargetMap().MyGID(fool->first))
      newmap[fool->first] = fool->second;
  data.clear();
  for (fool=newmap.begin(); fool != newmap.end(); ++fool)
    data[fool->first] = fool->second;

  //---------------------------------------------------------- do waiting
  for (int i=0; i<nsend; ++i) Wait(request[i]);

#endif
  return;
}
#endif


#endif  // #ifndef EXPORTER_H
#endif  // #ifdef CCADISCRET
