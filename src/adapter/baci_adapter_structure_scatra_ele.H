/*-----------------------------------------------------------*/
/*! \file


\brief factory for structure adapters

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_ADAPTER_STRUCTURE_SCATRA_ELE_H
#define BACI_ADAPTER_STRUCTURE_SCATRA_ELE_H

#include "baci_inpar_scatra.H"

namespace DRT
{
  class Element;
}

namespace ADAPTER
{
  INPAR::SCATRA::ImplType GetScaTraImplType(DRT::Element* ele);
}  // namespace ADAPTER


#endif
