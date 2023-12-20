/*---------------------------------------------------------------------*/
/*! \file

\brief Service for microstructural problems


\level 2

*/
/*---------------------------------------------------------------------*/

#ifndef BACI_STRU_MULTI_MICROSTATIC_NPSUPPORT_H
#define BACI_STRU_MULTI_MICROSTATIC_NPSUPPORT_H

#include "baci_config.H"

BACI_NAMESPACE_OPEN

namespace STRUMULTI
{
  // std::endless loop for supporting procs in multi scale problems
  void np_support_drt();


}  // namespace STRUMULTI

BACI_NAMESPACE_CLOSE

#endif
