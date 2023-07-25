/*---------------------------------------------------------------------*/
/*! \file

\brief registration functionality for parallel objects

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef BACI_MODULE_REGISTRY_PAROBJECTREGISTER_H
#define BACI_MODULE_REGISTRY_PAROBJECTREGISTER_H

namespace DRT
{
  /**
   * This call ensures that every ParObjectType listed within the corresponding source file
   * registers itself with ParObjectRegistry. This call is necessary at the beginning of any
   * executable that uses the parallel communication of data as documented in ParObject.
   */
  void ForceRegistrationOfParObjectTypes();
}  // namespace DRT



#endif
