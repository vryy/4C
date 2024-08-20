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
  Core::FE::ShapeFunctionType string_to_shape_function_type(std::string name);

  /// Return shape function name for a given shape function type
  std::string shape_function_type_to_string(Core::FE::ShapeFunctionType shapefunctiontype);

  const std::map<std::string, ShapeFunctionType>& string_to_shape_function_type_map();

}  // namespace Core::FE

FOUR_C_NAMESPACE_CLOSE

#endif