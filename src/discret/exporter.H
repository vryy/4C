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
  \brief Communicate nodes
  
  this method takes a map of nodes and redistirbutes them according to
  the send and receive plans. It is implicitly assumed, that the map of nodes
  pointwise matches SourceMap().
  
  */
  void Export(map<int,RefCountPtr<CCADISCRETIZATION::Node> >& nodes);

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











#endif  // #ifndef EXPORTER_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
