/*-----------------------------------------------------------*/
/*! \file


\brief factory for structure adapters

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef ADAPTER_STRUCTURE_SCATRA_ELE_H
#define ADAPTER_STRUCTURE_SCATRA_ELE_H

#include "inpar_scatra.H"

namespace DRT
{
  class Element;
}

namespace ADAPTER
{
  INPAR::SCATRA::ImplType GetScaTraImplType(DRT::Element* ele);
}  // namespace ADAPTER


#endif
