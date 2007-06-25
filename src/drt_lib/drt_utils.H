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
#include "drt_element.H"

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

  /// Create an allreduced vector of gids from the given Epetra_Map
  /*!
    We have nodes and elements with unique but otherwise arbitrary
    global ids. On rare occations, however, we need to allreduce a
    particular map to one or more processors. This is a building block
    for such occations. We allreduce the gids of the given Epetra_Map
    into a vector ordered by processor number.

    It is assumed that the given map does not overlap!

    \note You are not supposed to use redundant vectors in normal
    situations. If you happen to need this method you are probably
    about to do something illegal.

    \param rredundant (o) redundant vector of global ids
    \param emap (i) distributed Epetra_Map

    \author u.kue
    \date 05/07
   */
  void AllreduceEMap(vector<int>& rredundant, const Epetra_Map& emap);

  /// Create an allreduced gid to index map from the given Epetra_Map
  /*!
    We have nodes and elements with unique but otherwise arbitrary
    global ids. But unfortunately we need an allreduced vector of dof
    numbers during the dof assignment phase. In order to use such a
    vector we need to map from global ids to vector indexes. Here we
    provide that map.

    \note You are not supposed to use redundant vectors in normal
    situations. If you happen to need this method you are probably
    about to do something illegal.

    \param idxmap (o) map from global ids to (redundant) vector indexes
    \param emap (i) distributed Epetra_Map

    \author u.kue
    \date 05/07
   */
  void AllreduceEMap(map<int,int>& idxmap, const Epetra_Map& emap);

  /// find position of my map elements in a consecutive vector
  /*!
    The idea is to put the entries of a given map into a redundant
    vector, ordered by processor number. The map is assumed to be
    nonoverlapping. Here we figure out the index of our first entry in
    that vector.

    \note You are not supposed to use redundant vectors in normal
    situations. If you happen to need this method you are probably
    about to do something illegal.

    \param myrank (i) rank
    \param numproc (i) number of processors
    \param emap (i) distributed Epetra_Map

    \return vector position of first entry on each processor

    \author u.kue
    \date 05/07
   */
  int FindMyPos(int myrank, int numproc, const Epetra_Map& emap);

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


  /*!
  \brief Most common nodal mappings describing the relation between volume, surface and line
         node numbering

  They should be used as the only reference for such relationships.
  
  The corresponding graphics and a detailed description can be found
  in the Baci guide in the Convention chapter.
  
  The numbering of lower order elements is included in the higher order element, such that
  e.g. the hex8 volume element uses only the first 8 nodes of the hex27 mapping
  
  
  */
  
  // 6 Surfaces of a Hex27 element with 9 nodes per surface
  const int     eleNodeNumbering_hex27_surfaces[6][9] = {{0,  3,  2,  1, 11, 10,  9,  8, 20},
                                                         {0,  1,  5,  4,  8, 13, 16, 12, 21},
                                                         {1,  2,  6,  5,  9, 14, 17, 13, 22},
                                                         {2,  3,  7,  6, 10, 15, 18, 14, 23},
                                                         {0,  4,  7,  3, 12, 19, 15, 11, 24},
                                                         {4,  5,  6,  7, 16, 17, 18, 19, 25}};
  
  // 12 Lines of a Hex27 element with 3 nodes per line
  const int     eleNodeNumbering_hex27_lines[12][3] = {{0,  1,  8},
                                                       {1,  2,  9},
                                                       {2,  3, 10},
                                                       {0,  3, 11},
                                                       {0,  4, 12},
                                                       {1,  5, 13},
                                                       {2,  6, 14},
                                                       {3,  7, 15},
                                                       {4,  5, 16},
                                                       {5,  6, 17},
                                                       {6,  7, 18},
                                                       {4,  7, 19}};

  // 4 Surfaces of a Tet10 element with 6 nodes per surface
  const int     eleNodeNumbering_tet10_surfaces[4][6] = {{0,  1,  3,  4,  8,  7},
                                                         {1,  2,  3,  5,  9,  8},
                                                         {0,  3,  2,  7,  9,  6},
                                                         {0,  2,  1,  6,  5,  4}};

  // 6 Lines of a Tet10 element with 3 nodes per line
  const int     eleNodeNumbering_tet10_lines[6][3] = {{0,  1,  4},
                                                      {1,  2,  5},
                                                      {0,  2,  6},
                                                      {0,  3,  7},
                                                      {1,  3,  8},
                                                      {2,  3,  9}};

  // 4 Lines of a Quad9 element with 3 nodes per line
  const int     eleNodeNumbering_quad9_lines[4][3] = {{0,  1,  4},
                                                      {1,  2,  5},
                                                      {2,  3,  6},
                                                      {0,  3,  7}};

  // 3 Lines of a Tri6 element with 3 nodes per line
  const int     eleNodeNumbering_tri6_lines[3][3] = {{0,  1,  3},
                                                     {1,  2,  4},
                                                     {2,  0,  5}};

} // namespace Utils
} // namespace DRT


#endif  // #ifndef UTILS_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
