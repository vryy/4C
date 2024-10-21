// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_DEMANGLE_HPP
#define FOUR_C_UTILS_DEMANGLE_HPP

#include "4C_config.hpp"

#include <cxxabi.h>

#include <memory>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace Core::Utils
{
  /**
   * Utility function which tries to demangle the given @p mangledString. As a fallback, it returns
   * the input string. You should typically pass the result of typeid().name to this function.
   */
  inline std::string try_demangle(const char* mangledString)
  {
    int status;
    std::unique_ptr<char[], void (*)(void*)> result(
        abi::__cxa_demangle(mangledString, nullptr, nullptr, &status), std::free);

    return (status == 0 && result) ? std::string(result.get()) : std::string(mangledString);
  }
}  // namespace Core::Utils

FOUR_C_NAMESPACE_CLOSE

#endif