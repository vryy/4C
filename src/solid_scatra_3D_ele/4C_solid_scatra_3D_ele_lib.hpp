/*! \file

\brief Declaration of a solid-scatra coupling library functions

\level 1
*/

#ifndef FOUR_C_SOLID_SCATRA_3D_ELE_LIB_HPP
#define FOUR_C_SOLID_SCATRA_3D_ELE_LIB_HPP

#include "4C_config.hpp"

#include "4C_inpar_scatra.hpp"
#include "4C_io_linedefinition.hpp"


FOUR_C_NAMESPACE_OPEN

namespace Discret::ELEMENTS
{
  /*!
   * @brief Read the scatra implementation type from the input line definition of the element
   *
   * @param line_definition
   * @return Inpar::ScaTra::ImplType
   */
  Inpar::ScaTra::ImplType ReadScatraImplType(const Input::LineDefinition& line_definition);
}  // namespace Discret::ELEMENTS

FOUR_C_NAMESPACE_CLOSE

#endif
