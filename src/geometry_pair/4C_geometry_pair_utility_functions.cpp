#include "4C_geometry_pair_utility_functions.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
std::string GEOMETRYPAIR::discretization_type_geometry_to_string(
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
      FOUR_C_THROW(
          "GEOMETRYPAIR::DiscretizationTypeGeometryToString: Got unexpected discretization "
          "type.");
      return "";
      break;
  }
}

FOUR_C_NAMESPACE_CLOSE
