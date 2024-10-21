// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_LAZY_PTR_HPP
#define FOUR_C_UTILS_LAZY_PTR_HPP

#include "4C_config.hpp"

#include <functional>
#include <memory>
#include <utility>

FOUR_C_NAMESPACE_OPEN

namespace Core::Utils
{
  /**
   * A class which lazily allocates an object of type T. This class behaves mostly like a
   * std::shared_ptr, except that the object is only allocated when it is accessed for the first
   * time.
   *
   * @tparam T The type of the object to be allocated.
   */
  template <typename T>
  class LazyPtr
  {
   public:
    /**
     * Default constructor. The constructed LazyPtr will not allocate any object. Accessing the
     * content before copy or move assigning another LazyPtr will fail at runtime.
     */
    LazyPtr() = default;

    /**
     * Construct from a std::shared_ptr<T>.
     */
    LazyPtr(std::shared_ptr<T> ptr) : ptr_(std::move(ptr)) {}

    /**
     * Construct from a std::unique_ptr<T>.
     */
    LazyPtr(std::unique_ptr<T>&& ptr) : ptr_(std::move(ptr)) {}

    /**
     * Construct a LazyPtr with a @p constructor function which is called when the object is
     * accessed for the first time. The constructor function must return a type compatible with
     * std::shared_ptr<T>.
     */
    template <typename Fn>
    LazyPtr(Fn constructor) : constructor_(std::move(constructor))
    {
    }

    T& operator*() const { return *construct_if_necessary(); }

    T* operator->() const { return construct_if_necessary(); }

    T* get() const { return construct_if_necessary(); }

   private:
    decltype(auto) construct_if_necessary() const
    {
      if (!ptr_) ptr_ = constructor_();
      return ptr_.get();
    }

    //! The actual pointer storage is mutable because is is allocated lazily. The information
    //! is encoded in the constructor_.
    mutable std::shared_ptr<T> ptr_;
    std::function<std::shared_ptr<T>()> constructor_;
  };

}  // namespace Core::Utils

FOUR_C_NAMESPACE_CLOSE

#endif
