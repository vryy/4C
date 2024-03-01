/*! \file

\brief Declaration of a solid-scatra coupling library functions

\level 1
*/

#ifndef BACI_SOLID_SCATRA_3D_ELE_LIB_HPP
#define BACI_SOLID_SCATRA_3D_ELE_LIB_HPP

#include "baci_config.hpp"

#include "baci_inpar_scatra.hpp"
#include "baci_io_linedefinition.hpp"


BACI_NAMESPACE_OPEN

namespace DRT::ELEMENTS
{
  /*!
   * @brief Read the scatra implementation type from the input line definition of the element
   *
   * @param line_definition
   * @return INPAR::SCATRA::ImplType
   */
  INPAR::SCATRA::ImplType ReadScatraImplType(const INPUT::LineDefinition& line_definition);
}  // namespace DRT::ELEMENTS

BACI_NAMESPACE_CLOSE

#endif  // BACI_SOLID_SCATRA_3D_ELE_LIB_HPP
