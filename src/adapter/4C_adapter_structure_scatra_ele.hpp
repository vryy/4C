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

namespace CORE::Elements
{
  class Element;
}

namespace ADAPTER
{
  INPAR::SCATRA::ImplType GetScaTraImplType(CORE::Elements::Element* ele);
}  // namespace ADAPTER


FOUR_C_NAMESPACE_CLOSE

#endif
