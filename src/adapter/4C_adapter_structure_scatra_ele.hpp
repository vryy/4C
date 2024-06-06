/*-----------------------------------------------------------*/
/*! \file


\brief factory for structure adapters

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_STRUCTURE_SCATRA_ELE_HPP
#define FOUR_C_ADAPTER_STRUCTURE_SCATRA_ELE_HPP

#include "4C_config.hpp"

#include "4C_inpar_scatra.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Elements
{
  class Element;
}

namespace Adapter
{
  Inpar::ScaTra::ImplType GetScaTraImplType(Core::Elements::Element* ele);
}  // namespace Adapter


FOUR_C_NAMESPACE_CLOSE

#endif
