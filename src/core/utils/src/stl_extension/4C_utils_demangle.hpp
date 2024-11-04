// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_DEMANGLE_HPP
#define FOUR_C_UTILS_DEMANGLE_HPP

#include "4C_config.hpp"

#include <string>
#include <typeinfo>

FOUR_C_NAMESPACE_OPEN

namespace Core::Utils
{
  /**
   * Utility function which tries to demangle the given @p mangled_name. As a fallback, it returns
   * the input string. You should typically pass the result of a call to typeid().name to this
   * function.
   */
  std::string try_demangle(const char* mangled_name);


  /**
   * This function is a wrapper around typeid(). It returns the dynamic type of the given @p expr.
   * This function avoids a common misunderstanding of the typeid operator. Unlike the decltype()
   * operator, the typeid() operator actually needs to evaluate the expression to determine
   * the dynamic type. This may have side-effects. Wrapping the operation in a function makes it
   * clear to programmers and the compiler that the expression is evaluated.
   */
  template <typename T>
  const std::type_info& get_dynamic_type(const T& expr)
  {
    return typeid(expr);
  }

  /**
   * Utility function which tries to give the demangled type name of the given @p expression. If the
   * @p expression has a polymorphic type, the demangled dynamic type name is returned.
   * As a fallback, if demangling fails, it returns the mangled type name.
   */
  template <typename T>
  std::string get_dynamic_type_name(const T& expr)
  {
    return try_demangle(get_dynamic_type(expr).name());
  }

  /**
   * Utility function which tries to give the demangled type name of the given type @p T. As a
   * fallback, if demangling fails, it returns the mangled type name.
   */
  template <typename T>
  std::string get_type_name()
  {
    return try_demangle(typeid(T).name());
  }
}  // namespace Core::Utils

FOUR_C_NAMESPACE_CLOSE

#endif