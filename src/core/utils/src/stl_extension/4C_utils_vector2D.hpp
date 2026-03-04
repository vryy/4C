// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_VECTOR2D_HPP
#define FOUR_C_UTILS_VECTOR2D_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

#include <algorithm>
#include <memory>
#include <ranges>
#include <span>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::Utils
{
  /**
   * @brief A 2D vector container that stores a variable length vector in a vector with contiguous
   * memory
   *
   * This class provides a convenient interface for managing 2D arrays of arbitrary type T.
   * Internally, data is stored as a single contiguous vector for memory efficiency and cache
   * locality.
   *
   * @tparam T The element type stored in the 2D vector.
   *
   * @details
   * - Access: Use operator()(item, component) for element access
   * - Iteration: Use items() to iterate over items as spans
   * - Capacity: Use reserve() and push_back() for dynamic item insertion
   */
  template <typename T>
  class Vector2D
  {
   private:
    std::unique_ptr<T[]> data_ = nullptr;
    std::size_t size_ = 0;
    std::size_t num_components_ = 0;
    std::size_t capacity_ = 0;

   public:
    using value_type = T;

    Vector2D() = default;

    /**
     * @brief Constructs a 2D vector with a specified number of components per item.
     *
     * @param num_components The number of components for the 2D vector.
     */
    Vector2D(std::size_t num_components) : num_components_(num_components) {}

    /**
     * @brief Constructs a 2D vector with specified dimensions, initialized with a default value.
     *
     * @param items The number of items in the 2D vector.
     * @param num_components The number of components of the entries in the vector.
     * @param default_val The value to initialize all elements with. Defaults to
     * T{}(default-constructed value).
     */
    Vector2D(std::size_t items, std::size_t num_components, const T& default_val = T{})
        : data_(std::make_unique<T[]>(items * num_components)),
          size_(items),
          num_components_(num_components),
          capacity_(items * num_components)
    {
      std::fill_n(data_.get(), size_ * num_components_, default_val);
    }

    // Rule of five
    Vector2D(const Vector2D& other)
        : data_(std::make_unique<T[]>(other.size_ * other.num_components_)),
          size_(other.size_),
          num_components_(other.num_components_),
          capacity_(other.size_ * other.num_components_)
    {
      std::copy_n(other.data_.get(), size_ * num_components_, data_.get());
    }

    Vector2D& operator=(const Vector2D& other)
    {
      if (this != &other)
      {
        auto new_data = std::make_unique<T[]>(other.size_ * other.num_components_);
        std::copy_n(other.data_.get(), other.size_ * other.num_components_, new_data.get());
        data_ = std::move(new_data);
        size_ = other.size_;
        num_components_ = other.num_components_;
        capacity_ = other.size_ * other.num_components_;
      }

      return *this;
    }

    Vector2D(Vector2D&&) noexcept = default;
    Vector2D& operator=(Vector2D&&) noexcept = default;
    ~Vector2D() = default;

    /**
     * @brief Reserves storage space for the 2D vector.
     *
     * Increases the capacity of the underlying data storage to accommodate
     * at least the specified number of items, while keeping the current number
     * of components unchanged. If the requested size is less than or equal to the
     * current size, this function has no effect.
     *
     * @param new_size The number of items to reserve storage for.
     *
     * @pre num_components_ must not be zero. An assertion will fail if this precondition
     *      is violated.
     *
     * @note This operation reallocates the underlying storage and invalidates
     *       any pointers or references to existing elements.
     *
     * @complexity O(n) where n is new_size * num_components_
     */
    void reserve(std::size_t new_size)
    {
      FOUR_C_ASSERT_ALWAYS(
          num_components_ > 0, "Cannot reserve storage for a 2D vector with zero components.");
      if (new_size * num_components_ <= capacity_) return;

      auto new_data = std::make_unique<T[]>(new_size * num_components_);
      if (data_)
      {
        std::move(data_.get(), data_.get() + (size_ * num_components_), new_data.get());
      }
      data_ = std::move(new_data);
      capacity_ = new_size * num_components_;
    }

    /**
     * @brief Adds a new item of data to the 2D vector.
     *
     * Appends an item of elements to the end of the 2D vector. If this is the first item being
     * added, the number of components is automatically set to the size of the provided item data.
     * Subsequent items must have the same size as the first item.
     *
     * @param item_data The item of data to append. Its size must match the number of components
     *                 (num_components_) if the vector is not empty.
     *
     * @throw Raises an error if item_data.size() does not equal num_components_ when
     * num_components_ is already set (i.e., when adding an item to a non-empty 2D vector).
     */
    void push_back(const std::vector<T>& data)
    {
      if (num_components_ == 0)
      {
        FOUR_C_ASSERT_ALWAYS(data.size() > 0,
            "Cannot add an empty item to the 2D vector. Please provide at least one component.");
        num_components_ = data.size();
      }
      FOUR_C_ASSERT_ALWAYS(data.size() == num_components_,
          "Size of item data ({}) does not match number of components ({}).", data.size(),
          num_components_);

      // Exponential growth strategy (similar to std::vector)
      if ((size_ + 1) * num_components_ > capacity_)
      {
        reserve(size_ == 0 ? 1 : (size_ * 2));
      }

      std::copy(data.begin(), data.end(), data_.get() + (size_ * num_components_));
      ++size_;
    }


    /**
     * @brief Access a specific component of an entry
     *
     * @param item Item index (0-based)
     * @param component Component index (0-based)
     * @return Reference to the element at position [item, component]
     * + component)
     */
    T& operator()(std::size_t item, std::size_t component)
    {
      return data_[item * num_components_ + component];
    }

    const T& operator()(std::size_t item, std::size_t component) const
    {
      return data_[item * num_components_ + component];
    }

    /// @brief Get raw access to the underlying data storage as a pointer.
    T* data() { return data_.get(); }

    /// @brief Get raw access to the underlying data storage as a pointer.
    const T* data() const { return data_.get(); }

    /**
     * @brief Access an item of the 2D vector as a span (with a bound-check).
     *
     * @param item The item index to access.
     * @return A span view of the requested item containing @c num_components_ elements.
     * @throws Assertion failure if @c item is greater than or equal to the number of items.
     */
    std::span<T> at(std::size_t item)
    {
      FOUR_C_ASSERT_ALWAYS(
          item < size_, "Item index {} is out of bounds for Vector2D with {} items.", item, size_);
      return std::span<T>(data() + item * num_components_, num_components_);
    }

    /**
     * @brief Access an item of the 2D vector as a span (without a bound-check).
     *
     * @param item The item index to access.
     * @return A span view of the requested item containing @c num_components_ elements.
     * @throws Assertion failure if @c item is greater than or equal to the number of items.
     */
    std::span<T> operator[](std::size_t item)
    {
      return std::span<T>(data() + item * num_components_, num_components_);
    }

    /**
     * @brief Access an item of the 2D vector as a span (with a bound-check).
     *
     * @param item The item index to access.
     * @return A span view of the requested item containing @c num_components_ elements.
     * @throws Assertion failure if @c item is greater than or equal to the number of items.
     */
    std::span<const T> at(std::size_t item) const
    {
      FOUR_C_ASSERT_ALWAYS(
          item < size_, "Item index {} is out of bounds for Vector2D with {} items.", item, size_);
      return std::span<const T>(data() + item * num_components_, num_components_);
    }

    /**
     * @brief Access an item of the 2D vector as a span (without a bound-check).
     *
     * @param item The item index to access.
     * @return A span view of the requested item containing @c num_components_ elements.
     * @throws Assertion failure if @c item is greater than or equal to the number of items.
     */
    std::span<const T> operator[](std::size_t item) const
    {
      return std::span<const T>(data() + item * num_components_, num_components_);
    }

    /**
     * @brief Get the number of items in the 2D vector.
     * @return The number of items.
     */
    [[nodiscard]] std::size_t size() const { return size_; }

    /**
     * @brief Get the number of components per item in the 2D vector.
     * @return The number of components.
     */
    [[nodiscard]] std::size_t num_components() const { return num_components_; }


    /**
     * @brief Returns a range view of all items in the 2D vector as spans.
     *
     * @details Provides a lazy-evaluated range that iterates over each item of the 2D vector,
     * transforming item indices into non-owning span objects that reference the underlying data.
     * This allows for convenient iteration and access to individual items without copying.
     *
     * @return A range view where each element is a `std::span<T>` representing an item,
     * containing `num_components_` consecutive elements starting from the appropriate offset.
     *
     * @note The returned view is lazy-evaluated; iteration occurs on-demand.
     * @note The underlying data must remain valid while the view is in use.
     */
    auto items()
    {
      return std::ranges::iota_view<std::size_t, std::size_t>{0, size_} |
             std::views::transform([this](std::size_t item)
                 { return std::span<T>(data() + item * num_components_, num_components_); });
    }

    /**
     * @brief Returns a range view of all items in the 2D vector as spans (const version).
     *
     * @details Provides a lazy-evaluated range that iterates over each item of the 2D vector,
     * transforming item indices into non-owning span objects that reference the underlying data.
     * This allows for convenient iteration and access to individual items without copying.
     *
     * @return A range view where each element is a `std::span<const T>` representing an item,
     * containing `num_components_` consecutive elements starting from the appropriate offset.
     *
     * @note The returned view is lazy-evaluated; iteration occurs on-demand.
     * @note The underlying data must remain valid while the view is in use.
     */
    auto items() const
    {
      return std::ranges::iota_view<std::size_t, std::size_t>{0, size_} |
             std::views::transform([this](std::size_t item)
                 { return std::span<const T>(data() + item * num_components_, num_components_); });
    }
  };
}  // namespace Core::Utils

FOUR_C_NAMESPACE_CLOSE

#endif
