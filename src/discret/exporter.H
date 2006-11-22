/*!----------------------------------------------------------------------
\file exporter.H

\class CCADISCRETIZATION::Exporter

\brief A class to manage communcation patterns

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
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

#include "node.H"
#include "utils.H"

using namespace std;
using namespace Teuchos;

/*!
\brief CCADISCRET: namespace of the ccarat discretization module

*/
namespace CCADISCRETIZATION
{

/*!
\class Exporter

\brief A class to manage communcation patterns

\author gee (gee@lnm.mw.tum.de)
*/
class Exporter
{
public:


  // @{ \name Constructors and destructors

  /*!
  \brief Standard Constructor
  
  */
  Exporter(const Epetra_Map& frommap, const Epetra_Map& tomap, const Epetra_Comm& comm);

  /*!
  \brief Copy Constructor

  */
  Exporter(const CCADISCRETIZATION::Exporter& old);

  /*!
  \brief Destructor

  */
  virtual ~Exporter();

  //@}

  // @{ \name Acess methods

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

  // @{ \name Communication methods

  /*!
  \brief Communicate objects
  
  this method takes a map of objects and redistributes them according to
  the send and receive plans. It is implicitly assumed, that the map of nodes
  pointwise matches SourceMap(). It is also assumed (and tested), that type T 
  implements the ParObject class.
  
  */
  template<typename T> void Export(map<int,RefCountPtr<T> >& parobjects);

  //@}


private:

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
  inline vector<const char*>& SendBuff() { return sendbuff_; }

  /*!
  \brief Get sendsize_
  */
  inline vector<int>& SendSize() { return sendsize_; }

#ifdef PARALLEL
  /*!
  \brief Send data from one processor to another (nonblocking)

  \note This is an individual call
  */
  void ISend(const int frompid, const int topid, const char* data, 
            const int dsize, const int tag, MPI_Request& request);
  /*!
  \brief Receive anything joker (blocking)
  
  This method receives an MPI char string message from any source proc with
  any message tag of any length. It simply takes the first message that's 
  coming in. recvbuff is resized to fit received message.
  the method is blocking for the calling proc.
  It is used together with ISend and Wait to do nonblocking chaotic 
  point to point communication.
  
  \param source (output): source the message came from
  \param tag (output): message tag
  \param recvbuff (output): buffer containing received data
  \param length (output): length of message upon receive
  */
  void ReceiveAny(int& source, int&tag, vector<char>& recvbuff, int& length);
  
  /*!
  \brief wait for nonblocking send to finish
  */
  void Wait(MPI_Request& request) { MPI_Status status; MPI_Wait(&request,&status); return;}

            
#endif

private:

  const Epetra_Map&        frommap_;   // source layout
  const Epetra_Map&        tomap_;     // target map 
  const Epetra_Comm&       comm_;      // communicator
  int                      myrank_;    // PID
  int                      numproc_;   // no. of processors
  
  Epetra_SerialDenseMatrix sendplan_;  // sending information
  Epetra_SerialDenseMatrix recvplan_;  // receiving information

  vector<const char*>      sendbuff_;  // a sendbuffer
  vector<int>              sendsize_;  // sendsize_[i] is length of sendbuff_[i]
  
  
}; // class Exporter
} // namespace CCADISCRETIZATION




/*----------------------------------------------------------------------*
 |  communicate objects (public)                             mwgee 11/06|
 *----------------------------------------------------------------------*/
template<typename T> void CCADISCRETIZATION::Exporter::Export(
                                     map<int,RefCountPtr<T> >& parobjects)
{
  // test whether type T implements ParObject
  {
    typename map<int,RefCountPtr<T> >::iterator curr = parobjects.begin();
    if (curr != parobjects.end())
    {
      T* ptr = curr->second.get();
      ParObject* tester = dynamic_cast<ParObject*>(ptr);
      if (!tester) dserror("typename T in template does not implement ParObject (dynamic_cast failed)");
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
    SendBuff()[lid] = actobject->Pack(SendSize()[lid]);
    // do sending
    for (int j=0; j<NumProc(); ++j)
      if (j != MyPID() && SendPlan()(lid,j)==1)
      {
        if (nsend>=(int)request.size()) 
          request.resize(request.size()+200);
        ISend(MyPID(),j,SendBuff()[lid],SendSize()[lid],gid,request[nsend]);
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
    CCADISCRETIZATION::ParObject* object = 
                                 CCADISCRETIZATION::Factory(&(recvbuff[0]));
    T* ptr = dynamic_cast<T*>(object);
    if (!ptr) dserror("typename T in template does not implement ParObject (dynamic_cast failed)");
    RefCountPtr<T> refptr = rcp(ptr);
    // add object to the map
    parobjects[gid] = refptr;
  }
  
  //--------- loop map and kick out everything that's not in TargetMap()
  typename map<int,RefCountPtr<T> >::iterator curr;
  for (curr=parobjects.begin(); curr != parobjects.end(); ++curr)
  {
    if (TargetMap().MyGID(curr->first)) continue;
    parobjects.erase(curr->first);
  }
  
  //---------------------------------------------------------- do waiting
  for (int i=0; i<nsend; ++i) Wait(request[i]);

  //------------------------------------------------------ free sendbuffer
  for (int i=0; i<(int)SendBuff().size(); ++i)
    if (SendBuff()[i]) 
    {
      delete [] SendBuff()[i];
      SendBuff()[i] = NULL;
      SendSize()[i] = 0;
    }
    
#endif  
  return;
}







#endif  // #ifndef EXPORTER_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
