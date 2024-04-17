/*-----------------------------------------------------------*/
/*! \file


\brief factory for structure adapters

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_STRUCTURE_SCATRA_ELE_HPP
#define FOUR_C_ADAPTER_STRUCTURE_SCATRA_ELE_HPP

#include "baci_config.hpp"

#include "baci_inpar_scatra.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Element;
}

namespace ADAPTER
{
  INPAR::SCATRA::ImplType GetScaTraImplType(DRT::Element* ele);
}  // namespace ADAPTER


FOUR_C_NAMESPACE_CLOSE

#endif
