/*---------------------------------------------------------------------*/
/*! \file

\brief Service for microstructural problems


\level 2

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_STRU_MULTI_MICROSTATIC_NPSUPPORT_HPP
#define FOUR_C_STRU_MULTI_MICROSTATIC_NPSUPPORT_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace MultiScale
{
  // std::endless loop for supporting procs in multi scale problems
  void np_support_drt();


}  // namespace MultiScale

FOUR_C_NAMESPACE_CLOSE

#endif
