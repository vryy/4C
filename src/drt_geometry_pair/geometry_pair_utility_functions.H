/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for the geometry pairs.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_UTILITY_FUNCTIONS_H
#define GEOMETRY_PAIR_UTILITY_FUNCTIONS_H


#include <string>


// Forward declarations.
namespace GEOMETRYPAIR
{
  enum class DiscretizationTypeGeometry;
}


namespace GEOMETRYPAIR
{
  /**
   * \brief Convert the enum DiscretizationTypeGeometry to a human readable string.
   * @param discretization_type (in)
   * @return Human readable string representation of the enum.
   */
  std::string DiscretizationTypeGeometryToString(
      const DiscretizationTypeGeometry discretization_type);
}  // namespace GEOMETRYPAIR


#endif
