/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for the geometry pairs.

\level 1
*/
// End doxygen header.



#include "baci_geometry_pair_utility_functions.hpp"

#include "baci_geometry_pair_element_classes.hpp"

BACI_NAMESPACE_OPEN

/**
 *
 */
std::string GEOMETRYPAIR::DiscretizationTypeGeometryToString(
    const DiscretizationTypeGeometry discretization_type)
{
  switch (discretization_type)
  {
    case DiscretizationTypeGeometry::none:
      return "undefined";
    case DiscretizationTypeGeometry::triangle:
      return "triangle";
    case DiscretizationTypeGeometry::quad:
      return "quadrilateral";
    case DiscretizationTypeGeometry::hexahedron:
      return "hexahedron";
    case DiscretizationTypeGeometry::tetraeder:
      return "tetraeder";
    default:
      dserror(
          "GEOMETRYPAIR::DiscretizationTypeGeometryToString: Got unexpected discretization "
          "type.");
      return "";
      break;
  }
}

BACI_NAMESPACE_CLOSE
