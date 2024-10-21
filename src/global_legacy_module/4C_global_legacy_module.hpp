#ifndef FOUR_C_GLOBAL_LEGACY_MODULE_HPP
#define FOUR_C_GLOBAL_LEGACY_MODULE_HPP

#include "4C_config.hpp"

#include "4C_module_registry_callbacks.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 * Expose callbacks for the global legacy module.
 *
 * @note This function exists for historic reasons. Its internals need to be split up to remove
 * forced dependencies on the main apps.
 */
[[nodiscard]] ModuleCallbacks global_legacy_module_callbacks();

FOUR_C_NAMESPACE_CLOSE

#endif
