/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of shape function types
\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef BACI_DISCRETIZATION_FEM_GENERAL_SHAPE_FUNCTION_TYPE_H
#define BACI_DISCRETIZATION_FEM_GENERAL_SHAPE_FUNCTION_TYPE_H

#include <map>
#include <string>

BACI_NAMESPACE_OPEN

namespace CORE::FE
{
  /// Type of shape functions used in spatial discretization
  enum class ShapeFunctionType
  {
    undefined,   ///< Undefined
    polynomial,  ///< Polynomial shape functions
    nurbs,       ///< NURBS shape functions
    hdg          ///< Hybridizable Discontinuous Galerkin
  };

  /// Return shape function type enum for a given shape function name
  CORE::FE::ShapeFunctionType StringToShapeFunctionType(std::string name);

  /// Return shape function name for a given shape function type
  std::string ShapeFunctionTypeToString(CORE::FE::ShapeFunctionType shapefunctiontype);

  const std::map<std::string, ShapeFunctionType>& StringToShapeFunctionTypeMap();

}  // namespace CORE::FE

BACI_NAMESPACE_CLOSE

#endif