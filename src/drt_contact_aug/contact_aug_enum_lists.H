/*----------------------------------------------------------------------------*/
/*! \file
\brief Enumerator lists for the internal use (no input parameters)

\level 3

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_ENUM_LISTS_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_ENUM_LISTS_H_

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


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_ENUM_LISTS_H_ */
