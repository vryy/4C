// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_SHARED_PTR_FROM_REF_HPP
#define FOUR_C_UTILS_SHARED_PTR_FROM_REF_HPP

#include "4C_config.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::Utils
{
  /*!
   * \brief Create a shared pointer from a reference.
   *
   * In an ideal world, we would not need this function. However, there are cases where we need to
   * pass an argument to a function that expects a shared pointer, but we only have a reference to
   * the object. This function creates a shared pointer from a reference and uses a custom deleter
   * that does nothing. This way, the shared pointer does not take ownership of the object and does
   * not delete it when it goes out of scope.
   *
   * @note Make sure that the object the reference points to is still alive when the shared pointer
   * is used.
   */
  template <typename T>
  std::shared_ptr<T> shared_ptr_from_ref(T& ref)
  {
    // Use the aliasing constructor as suggested here: https://stackoverflow.com/a/36691828
    return std::shared_ptr<T>(std::shared_ptr<T>{}, &ref);
  }
}  // namespace Core::Utils

FOUR_C_NAMESPACE_CLOSE

#endif
