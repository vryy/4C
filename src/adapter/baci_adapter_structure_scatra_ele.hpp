/*-----------------------------------------------------------*/
/*! \file


\brief factory for structure adapters

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_ADAPTER_STRUCTURE_SCATRA_ELE_HPP
#define BACI_ADAPTER_STRUCTURE_SCATRA_ELE_HPP

#include "baci_config.hpp"

#include "baci_inpar_scatra.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Element;
}

namespace ADAPTER
{
  INPAR::SCATRA::ImplType GetScaTraImplType(DRT::Element* ele);
}  // namespace ADAPTER


BACI_NAMESPACE_CLOSE

#endif
