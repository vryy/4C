/*---------------------------------------------------------------------*/
/*!

\brief A collection of helper methods for namespace DRT

\level 0

\maintainer Martin Kronbichler

*/
/*---------------------------------------------------------------------*/

#ifndef UTILS_METIS_H
#define UTILS_METIS_H

#include "Teuchos_RCP.hpp"
#include "Epetra_CrsGraph.h"
#include "Epetra_Vector.h"

namespace DRT
{
  namespace UTILS
  {
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
    Teuchos::RCP<Epetra_CrsGraph> PartGraphUsingMetis(
        const Epetra_CrsGraph& graph, const Epetra_Vector& weights);


  }  // namespace UTILS
}  // namespace DRT


#endif  // #ifndef UTILS_PARMETIS_H
