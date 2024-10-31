// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_utils_demangle.hpp"

#include <cxxabi.h>

#include <memory>

FOUR_C_NAMESPACE_OPEN


std::string Core::Utils::try_demangle(const char* mangled_name)
{
  int status;
  std::unique_ptr<char[], void (*)(void*)> result(
      abi::__cxa_demangle(mangled_name, nullptr, nullptr, &status), std::free);

  return (status == 0 && result) ? std::string(result.get()) : std::string(mangled_name);
}

FOUR_C_NAMESPACE_CLOSE