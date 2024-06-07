/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of shape function types
\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FEM_GENERAL_SHAPE_FUNCTION_TYPE_HPP
#define FOUR_C_FEM_GENERAL_SHAPE_FUNCTION_TYPE_HPP

#include "4C_config.hpp"

#include <map>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
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
  Core::FE::ShapeFunctionType StringToShapeFunctionType(std::string name);

  /// Return shape function name for a given shape function type
  std::string ShapeFunctionTypeToString(Core::FE::ShapeFunctionType shapefunctiontype);

  const std::map<std::string, ShapeFunctionType>& StringToShapeFunctionTypeMap();

}  // namespace Core::FE

FOUR_C_NAMESPACE_CLOSE

#endif