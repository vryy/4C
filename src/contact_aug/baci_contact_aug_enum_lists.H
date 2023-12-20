/*----------------------------------------------------------------------------*/
/*! \file
\brief Enumerator lists for the internal use (no input parameters)

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef BACI_CONTACT_AUG_ENUM_LISTS_H
#define BACI_CONTACT_AUG_ENUM_LISTS_H

#include "baci_config.H"

BACI_NAMESPACE_OPEN

namespace CONTACT
{
  namespace AUG
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
  }  // namespace AUG
}  // namespace CONTACT

BACI_NAMESPACE_CLOSE

#endif  // CONTACT_AUG_ENUM_LISTS_H
