/*-----------------------------------------------------------*/
/*! \file

\brief Utility functions for the geometric search.

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_DISCRETIZATION_GEOMETRIC_SEARCH_UTILS_HPP
#define BACI_DISCRETIZATION_GEOMETRIC_SEARCH_UTILS_HPP

#include "baci_config.hpp"

#include "baci_discretization_geometric_search_bounding_volume.hpp"

#include <Epetra_Comm.h>

BACI_NAMESPACE_OPEN

namespace CORE::GEOMETRICSEARCH
{
  /*! \brief Storing information on the geometric search
   */
  struct GeometricSearchInfo
  {
    int primitive_size;
    int predicate_size;
    int coupling_pair_size;
  };

  /*! \brief Prints details on the geometric search algorithm
   */
  void PrintGeometricSearchDetails(const Epetra_Comm &comm, const GeometricSearchInfo info);

  /*! \brief Returns interaction pair indices based on the search output of ArborX
   *
   * @tparam T Data container for stored indice and offset array
   * @param indices Stores the indices of the objects that satisfy predicates
   * @param offset Stores the locations in indices that start a predicate
   * @return Interaction pair indices , i.e., {i_predicate, i_primitive} where i_predicate is
   * satisfied by i_primitive. The indices correspond to the primitive and predicate vectors, NOT
   * element GIDs.
   */
  template <typename T>
  std::vector<std::pair<int, int>> GetPairs(const T &indices, const T &offset)
  {
    std::vector<std::pair<int, int>> pairs;
    for (size_t i_offset = 0; i_offset < offset.size() - 1; i_offset++)
      for (int j = offset[i_offset]; j < offset[i_offset + 1]; j++)
        pairs.emplace_back(std::pair{i_offset, indices[j]});

    return pairs;
  }

  /*! \brief Get the polyhedron representation of a k-DOP
   *
   * @param boundingVolume Bounding volume enclosing the respective element (as k-DOP)
   * @return Points of the polyhedron and connecting polygons
   */
  std::pair<std::vector<LINALG::Matrix<3, 1>>, std::vector<std::vector<int>>>
  GetKDopPolyhedronRepresentation(const BoundingVolume boundingVolume);

}  // namespace CORE::GEOMETRICSEARCH

BACI_NAMESPACE_CLOSE

#endif
