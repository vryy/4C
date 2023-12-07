/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of shape function types
\level 0
*/
/*----------------------------------------------------------------------*/

#include "baci_discretization_fem_general_shape_function_type.H"

#include "baci_utils_exceptions.H"

#include <algorithm>
#include <map>

namespace CORE::FE
{
  namespace
  {
    const std::map<std::string, ShapeFunctionType> string2shapefuntype{
        {"Polynomial", ShapeFunctionType::polynomial}, {"Nurbs", ShapeFunctionType::nurbs},
        {"HDG", ShapeFunctionType::hdg}};
  }  // namespace
  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  ShapeFunctionType StringToShapeFunctionType(std::string name)
  {
    const auto it = string2shapefuntype.find(name);
    if (it != string2shapefuntype.end()) return it->second;
    dserror(
        "'%s' does not name a shape function type. Check for typos or consider adding the shape "
        "function type to the map.",
        name.c_str());
  }

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  std::string ShapeFunctionTypeToString(ShapeFunctionType shapefunctiontype)
  {
    const auto it = std::find_if(string2shapefuntype.begin(), string2shapefuntype.end(),
        [&](const auto& kv) { return kv.second == shapefunctiontype; });


    if (it != string2shapefuntype.end()) return it->first;
    dserror(
        "Could not find the name of the given shape function type or the shapefunction is "
        "undefined.");
  }

  const std::map<std::string, ShapeFunctionType>& StringToShapeFunctionTypeMap()
  {
    return string2shapefuntype;
  }
}  // namespace CORE::FE