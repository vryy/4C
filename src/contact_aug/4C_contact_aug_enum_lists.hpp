/*----------------------------------------------------------------------------*/
/*! \file
\brief Enumerator lists for the internal use (no input parameters)

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_ENUM_LISTS_HPP
#define FOUR_C_CONTACT_AUG_ENUM_LISTS_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CONTACT
{
  namespace Aug
  {
    enum class MapType : char
    {
      all_slave_nodes,
      active_slave_nodes
    };

    enum class WGapGradientType : char
    {
      vague,
      force_balance,
      constraint_enforcement
    };

    enum class SideType : char
    {
      slave_master,
      slave,
      master
    };
  }  // namespace Aug
}  // namespace CONTACT

FOUR_C_NAMESPACE_CLOSE

#endif
