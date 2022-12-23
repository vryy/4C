/*-----------------------------------------------------------*/
/*! \file

\brief Utility functions for the geometric search.

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef GEOMETRIC_SEARCH_UTILS_H
#define GEOMETRIC_SEARCH_UTILS_H

namespace GEOMETRICSEARCH::UTILS
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
  void PrintGeometricSearchDetails(const Epetra_Comm& comm, const UTILS::GeometricSearchInfo info);

}  // namespace GEOMETRICSEARCH::UTILS

#endif
