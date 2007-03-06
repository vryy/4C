/*!----------------------------------------------------------------------
\file drt_utils.H
\brief A collection of helper methods for namespace DRT

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef UTILS_H
#define UTILS_H


#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Teuchos_RefCountPtr.hpp"

#include "drt_parobject.H"

using namespace std;
using namespace Teuchos;

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{
namespace Utils
{
  
  /*!
  \brief Create an instance of a ParObject depending on the type stored in data
  
  An instance of ParObject is allocated and returned. The type of instance
  depends on the first entry in data and is an int (defined at the top of
  ParObject.H)
  
  \param data (in): A char vector used for communication or io
  
  \warning All instances of ParObject have to store the parobject id 
           at the vey beginning of the char vector
           
  \return Allocates and returns the correct instance of ParObject where data is
          already unpacked in the instance. The calling method is responsible for
          freeing this instance!
  */
  DRT::ParObject* Factory(const vector<char>& data);

  /*!
  \brief Partition a graph with metis
  
  Take an input graph with any row map and a vector of point weights with the same row
  map and redistribute the graph using the partitioning library metis.
  The output graph has a row map according to the new distribution.
  
  \param graph : (in) any graph
  \param weights: (int) a vector of weight factors for each row of graph. 
                  Map of weights must be pointwise equal to graph.RowMap().
                  
  \return graph with new distribution                  
  */
  RefCountPtr<Epetra_CrsGraph> PartGraphUsingMetis(const Epetra_CrsGraph& graph,
                                                   const Epetra_Vector& weights);

  /*!
  \brief Locally extract a subset of values from an Epetra_Vector
  
  Extracts lm.size() values from a distributed epetra vector and stores them into local.
  this is NOT a parallel method, meaning that all values to be extracted on a processor
  must be present in global on that specific processor. This usually means that global
  has to be in column map style.
  
  \param global (in): global distributed vector with values to be extracted
  \param local (out): vector holding values etxracted from global
  \param lm     (in): vector containing global ids to be extracted. Size of lm
                      determines number of values to be extracted.
  */
  void ExtractMyValues(const Epetra_Vector& global, vector<double>& local, const vector<int> lm);

} // namespace Utils
} // namespace DRT


#endif  // #ifndef UTILS_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
