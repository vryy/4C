// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_MATHOPERATIONS_SACADO_HPP
#define FOUR_C_UTILS_MATHOPERATIONS_SACADO_HPP

#include "4C_config.hpp"

#include "4C_utils_mathoperations.hpp"

#include <Sacado.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core
{
  template <typename T>
    requires(Sacado::IsFad<T>::value || Sacado::IsExpr<T>::value)
  struct MathOperations<T>
  {
    template <typename T2>
    static constexpr std::decay_t<T2> abs(T2&& t)
    {
      // Note: This only works since Sacado exports their math functions to namespace std and not
      // one of their own.
      return std::abs(std::forward<T2>(t));
    }
    template <typename T2>
    static constexpr std::decay_t<T2> sqrt(T2&& t)
    {
      // Note: This only works since Sacado exports their math functions to namespace std and not
      // one of their own.
      return std::sqrt(std::forward<T2>(t));
    }
    template <typename T2>
    static constexpr std::decay_t<T2> pow(T2&& base, int exponent)
    {
      // Note: This only works since Sacado exports their math functions to namespace std and not
      // one of their own.
      return std::pow(std::forward<T2>(base), exponent);
    }
  };
}  // namespace Core
FOUR_C_NAMESPACE_CLOSE

#endif
