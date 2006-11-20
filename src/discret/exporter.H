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
#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialDenseMatrix.h"

using namespace std;

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
  \brief Get sendpid_
  */
  inline vector<int>& SendPID() { return sendpid_; }
  /*!
  \brief Get sendlid_
  */
  inline vector<int>& SendLID() { return sendlid_; }
  /*!
  \brief Get recvpid_
  */
  inline vector<int>& RecvPID() { return recvpid_; }
  /*!
  \brief Get recvlid_
  */
  inline vector<int>& RecvLID() { return recvlid_; }

private:

  const Epetra_Map&        frommap_;   // source layout
  const Epetra_Map&        tomap_;     // target map 
  const Epetra_Comm&       comm_;      // communicator
  int                      myrank_;    // PID
  int                      numproc_;   // no. of processors
  
  /*! 
    length SourceMap().NumMyElements()
    For each element in SourceMap, the target processor is stored.
    Target processor might be herself as well.
  */
  vector <int>             sendpid_;   

  /*! 
    length SourceMap().NumMyElements()
    For each element in SourceMap, the local index lid on the target processor 
    is stored. Target processor might be herself as well.
  */
  vector <int>             sendlid_;   // length SourceMap().NumMyElements()
  
  /*! 
    length TargetMap().NumMyElements()
    For each element in TargetMap, the source processor is stored.
    Source processor might be herself as well.
  */
  vector <int>             recvpid_;   // length TargetMap().NumMyElements()

  /*! 
    length TargetMap().NumMyElements()
    For each element in TargetMap, the local index lid on the source 
    processor is stored. Source processor might be herself as well.
  */
  vector <int>             recvlid_;   // length TargetMap().NumMyElements()

}; // class Exporter
} // namespace CCADISCRETIZATION











#endif  // #ifndef EXPORTER_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
