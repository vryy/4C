/*----------------------------------------------------------------------*/
/*! \file
\brief Center coordinates of an element
\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_FEM_GENERAL_ELEMENT_CENTER_HPP
#define FOUR_C_DISCRETIZATION_FEM_GENERAL_ELEMENT_CENTER_HPP

#include "4C_config.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Elements
{
  class Element;
}

namespace Core::FE
{
  // return center coordinates of @p ele. Note that the average of all nodal coordinates is
  // computed.
  std::vector<double> element_center_refe_coords(const Core::Elements::Element& ele);
}  // namespace Core::FE

#endif

FOUR_C_NAMESPACE_CLOSE
