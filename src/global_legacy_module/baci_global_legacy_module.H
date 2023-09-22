/*---------------------------------------------------------------------*/
/*! \file
\brief This file provides a legacy module which exposes all functionality
       which isn't yet modularized to the main executables.
\level 0
*/
/*---------------------------------------------------------------------*/

#ifndef BACI_GLOBAL_LEGACY_MODULE_H
#define BACI_GLOBAL_LEGACY_MODULE_H

#include "baci_module_registry_callbacks.H"

namespace BACI
{
  /**
   * Expose callbacks for the global legacy module.
   *
   * @note This function exists for historic reasons. Its internals need to be split up to remove
   * forced dependencies on the main apps.
   */
  ModuleCallbacks GlobalLegacyModuleCallbacks();
}  // namespace BACI

#endif
