/*---------------------------------------------------------------------*/
/*! \file
\brief This file provides a legacy module which exposes all functionality
       which isn't yet modularized to the main executables.
\level 0
*/
/*---------------------------------------------------------------------*/

#ifndef BACI_GLOBAL_LEGACY_MODULE_H
#define BACI_GLOBAL_LEGACY_MODULE_H

#include "baci_config.H"

#include "baci_module_registry_callbacks.H"

BACI_NAMESPACE_OPEN

/**
 * Expose callbacks for the global legacy module.
 *
 * @note This function exists for historic reasons. Its internals need to be split up to remove
 * forced dependencies on the main apps.
 */
[[nodiscard]] ModuleCallbacks GlobalLegacyModuleCallbacks();

BACI_NAMESPACE_CLOSE

#endif
