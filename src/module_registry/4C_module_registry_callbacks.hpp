/*---------------------------------------------------------------------*/
/*! \file
\brief Central callback structure to register modules in executables
\level 0
*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_MODULE_REGISTRY_CALLBACKS_HPP
#define FOUR_C_MODULE_REGISTRY_CALLBACKS_HPP

#include "4C_config.hpp"

#include "4C_utils_function_manager.hpp"

#include <functional>

FOUR_C_NAMESPACE_OPEN

/**
 * This struct holds callback functions that are populated by modules. The callbacks are invoked
 * by executables if the specific module is configured. All callbacks are optional and do not have
 * to be set inside the modules.
 */
struct ModuleCallbacks
{
  /**
   * Register all types defined in a module that are derived from CORE::COMM::ParObject and can be
   * communicated over MPI.
   *
   * This call ensures that every ParObjectType listed within the corresponding source file
   * registers itself with ParObjectRegistry. This call is necessary at the beginning of any
   * executable that uses the parallel communication of data as documented in ParObject.
   *
   */
  std::function<void()> RegisterParObjectTypes;

  /**
   * Allow the module to attach any custom Function to the @p function_manager object. Inside this
   * callback, a module should call `FunctionManager::AddFunctionDefinition`.
   */
  std::function<void(CORE::UTILS::FunctionManager& function_manager)> AttachFunctionDefinitions;

  /**
   * A callback to attach valid result description lines to.
   * Modules should use the Add() function on the Lines object to add their own valid lines.
   */
  
  std::function<void(INPUT::Lines&)> AttachResultLines;
};

FOUR_C_NAMESPACE_CLOSE

#endif
