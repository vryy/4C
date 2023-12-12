/*---------------------------------------------------------------------*/
/*! \file
\brief Central callback structure to register modules in executables
\level 0
*/
/*---------------------------------------------------------------------*/

#ifndef BACI_MODULE_REGISTRY_CALLBACKS_H
#define BACI_MODULE_REGISTRY_CALLBACKS_H

#include "baci_utils_function_manager.H"

#include <functional>

BACI_NAMESPACE_OPEN

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
};

BACI_NAMESPACE_CLOSE

#endif
