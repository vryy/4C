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
#ifndef UTILS_H
#define UTILS_H


#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Teuchos_RefCountPtr.hpp"

#include "drt_parobject.H"
#include "drt_element.H"
#include "drt_utils_fem_shapefunctions.H"
#include "drt_utils_dgfem_basisfunctions.H"
#include "drt_utils_integration.H"
#include "drt_utils_local_connectivity_matrices.H"

using namespace std;
using namespace Teuchos;

namespace LINALG
{
  class MapExtractor;
}


/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{
namespace UTILS
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
  \brief Create an instance of a finie element depending on the type of element

  A DRT::Element derived class is allocated and returned. The type of element
  allocated depends on the input parameter eletype.

  \param eletype (in): A string containing the type of element
  \param id (in)     : id of the new element to be created
  \param owner (in)  : owner of the new element

  */
  RefCountPtr<DRT::Element> Factory(const string eletype, const int id, const int owner);

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
  \param local (out): vector holding values extracted from global
  \param lm     (in): vector containing global ids to be extracted. Size of lm
                      determines number of values to be extracted.
  */
  void ExtractMyValues(const Epetra_Vector& global, vector<double>& local, const vector<int> lm);

  /// find all local nodes from discretization marked with condition
  /*!

    Loop all conditions of the given discretization, find the ones with the
    specified name and collect the locally owned node ids in the suppied
    set. The nodes vector is unique and ordered on output.

    \param dis : (in) Discretization
    \param condname : (in) name of condition in question
    \param nodes : (out) empty set on input, filled with nodal gids of local nodes

    \author u.kue
    \date 06/07
   */
  void FindConditionedNodes(const DRT::Discretization& dis,
                            std::string condname,
                            std::vector<int>& nodes);

  /// Split dof row map of discretization into velocity and pressure maps
  /*!

    We assume that each nodes possesses ndim velocity dofs and (optionally)
    one pressure dof. The dof row map is thus split in two.

    \param dis : (in) Discretization
    \param ndim : (in) dimensions of problem (2 or 3)
    \param extractor : (out) ready made map splitter

    \author u.kue
    \date 02/08
   */
  void SetupFluidSplit(const DRT::Discretization& dis,
                       int ndim,
                       LINALG::MapExtractor& extractor);

} // namespace UTILS
} // namespace DRT


#endif  // #ifndef UTILS_H
#endif  // #ifdef CCADISCRET
