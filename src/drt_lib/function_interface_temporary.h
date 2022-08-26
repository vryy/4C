/*----------------------------------------------------------------------*/
/*! \file
\brief Temporary interface while the function framework is restructured
\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef SRC_DRT_LIB_FUNCTION_INTERFACE_TEMPORARY_H
#define SRC_DRT_LIB_FUNCTION_INTERFACE_TEMPORARY_H

namespace DRT::UTILS
{
  /**
   * This interface only exists to facilitate a restructuring of the function framework.
   *
   * The to be introduced, new function interfaces can temporarily derive from this interface. The
   * FunctionManager stores functions via this interface but gives access only via the newly
   * introduced interfaces.
   */
  class TemporaryFunctionInterface
  {
   public:
    virtual ~TemporaryFunctionInterface() = default;
  };
}  // namespace DRT::UTILS

#endif
