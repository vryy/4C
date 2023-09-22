/*---------------------------------------------------------------------*/
/*! \file
\brief Central callback structure to register modules in executables
\level 0
*/
/*---------------------------------------------------------------------*/

#ifndef BACI_MODULE_REGISTRY_CALLBACKS_H
#define BACI_MODULE_REGISTRY_CALLBACKS_H

#include <functional>

namespace BACI
{
  /**
   * This struct holds callback functions that are populated by modules. The callbacks are invoked
   * by executables if the specific module is configured. All callbacks are optional and do not have
   * to be set inside the modules.
   */
  struct ModuleCallbacks
  {
    /**
     * Register all types defined in a module that are derived from DRT::ParObject and can be
     * communicated over MPI.
     *
     * This call ensures that every ParObjectType listed within the corresponding source file
     * registers itself with ParObjectRegistry. This call is necessary at the beginning of any
     * executable that uses the parallel communication of data as documented in ParObject.
     *
     */
    std::function<void()> RegisterParObjectTypes;
  };

}  // namespace BACI



#endif
