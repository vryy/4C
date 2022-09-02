/*----------------------------------------------------------------------*/
/*! \file
\brief Utility to demangle typeid names
\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef DRT_LIB_DEMANGLE_H
#define DRT_LIB_DEMANGLE_H

#include <string>
#include <memory>
#include <cxxabi.h>

namespace DRT::UTILS
{
  /**
   * Utility function which tries to demangle the given @p mangledString. As a fallback, it returns
   * the input string. You should typically pass the result of typeid().name to this function.
   */
  inline std::string TryDemangle(const char* mangledString)
  {
    int status;
    std::unique_ptr<char[], void (*)(void*)> result(
        abi::__cxa_demangle(mangledString, nullptr, nullptr, &status), std::free);

    return (status == 0 && result) ? std::string(result.get()) : std::string(mangledString);
  }
}  // namespace DRT::UTILS

#endif