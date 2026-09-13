// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_TENSOR_META_UTILS_HPP
#define FOUR_C_LINALG_TENSOR_META_UTILS_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <sstream>
#include <tuple>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg::Internal
{
  enum class TensorBoundCheck : std::uint8_t
  {
    no_check,
    check_with_assertions,
    check
  };

  enum class OrderType : bool
  {
    row_major,
    column_major
  };

  template <typename... Indices>
  constexpr std::array<std::size_t, sizeof...(Indices)> to_array(Indices... i)
  {
    return {static_cast<std::size_t>(i)...};
  }

  template <OrderType order_type, std::size_t... n>
  consteval std::array<std::size_t, sizeof...(n)> get_index_offset()
  {
    constexpr std::array<std::size_t, sizeof...(n)> shape = {n...};
    std::array<std::size_t, sizeof...(n)> index_offset = {1};

    if constexpr (order_type == OrderType::column_major)
    {
      index_offset[0] = 1;
      std::partial_sum(
          shape.begin(), shape.end() - 1, index_offset.begin() + 1, std::multiplies<>());
    }
    else
    {
      index_offset[sizeof...(n) - 1] = 1;

      std::partial_sum(
          shape.rbegin(), shape.rend() - 1, index_offset.rbegin() + 1, std::multiplies<>());
    }
    return index_offset;
  }

  template <std::size_t... n, typename... Indices>
  void check_bounds(Indices... i)
  {
    static_assert(
        sizeof...(Indices) == sizeof...(n), "number of indices must match number of dimensions");

    constexpr std::array<std::size_t, sizeof...(n)> shape = {n...};

    std::array<std::size_t, sizeof...(n)> index = to_array(i...);

    auto join = [](const auto& array) -> std::string
    {
      std::ostringstream joined_string;
      auto begin = array.begin();
      if (array.begin() != array.end())
      {
        joined_string << *begin++;
        for (; begin != array.end(); ++begin) joined_string << ',' << *begin;
      }
      return joined_string.str();
    };

    FOUR_C_ASSERT_ALWAYS(std::inner_product(index.begin(), index.end(), shape.begin(), 0,
                             std::plus<>(), std::greater_equal<>()) == 0,
        "Given index is out of bounds. Given index is ({}) but Tensor has shape ({}).", join(index),
        join(shape));
  }

  template <OrderType order_type, TensorBoundCheck bound_check, std::size_t... n,
      typename... Indices>
  constexpr std::size_t get_flat_index(Indices... i)
  {
    static_assert(
        sizeof...(Indices) == sizeof...(n), "number of indices must match number of dimensions");

    if constexpr (bound_check == TensorBoundCheck::check) check_bounds<n...>(i...);
    if constexpr (bound_check == TensorBoundCheck::check_with_assertions)
    {
#ifdef FOUR_C_ENABLE_ASSERTIONS
      check_bounds<n...>(i...);
#endif
    }

    constexpr std::array<std::size_t, sizeof...(n)> index_offset =
        get_index_offset<order_type, n...>();

    std::array<std::size_t, sizeof...(n)> index = to_array(i...);

    return std::inner_product(index.begin(), index.end(), index_offset.begin(), std::size_t(0));
  }

  template <OrderType order_type, TensorBoundCheck bound_check, std::size_t... n>
  constexpr std::tuple<decltype(n)...> get_md_index(std::size_t flat_index)
  {
    if constexpr (bound_check == TensorBoundCheck::check)
    {
      FOUR_C_ASSERT_ALWAYS(
          flat_index < (n * ...), "flat index must be smaller than the tensor size!");
    }
    if constexpr (bound_check == TensorBoundCheck::check_with_assertions)
    {
      FOUR_C_ASSERT(flat_index < (n * ...), "flat index must be smaller than the tensor size!");
    }
    static_assert(
        order_type == OrderType::column_major, "Currently only column major order is supported");

    std::array<std::size_t, sizeof...(n)> multi_index{};
    constexpr std::array<std::size_t, sizeof...(n)> index_offset =
        get_index_offset<order_type, n...>();

    for (std::size_t i = sizeof...(n); i > 0; i--)
    {
      std::size_t index = flat_index / index_offset[i - 1];
      flat_index = flat_index % index_offset[i - 1];
      multi_index[i - 1] = index;
    }

    return std::tuple_cat(multi_index);
  }

  template <typename Number, std::size_t... n>
  struct TensorInitializerList;

  template <typename Number, std::size_t n1>
  struct TensorInitializerList<Number, n1>
  {
    // Note, we use C-style arrays here so that we get size checks when initializing with brace
    // initializers
    using type = Number[n1];
  };

  template <typename Number, std::size_t n1, std::size_t... n>
  struct TensorInitializerList<Number, n1, n...>
  {
    // Note, we use C-style arrays here so that we get size checks when initializing with brace
    // initializers
    using type = typename TensorInitializerList<Number, n...>::type[n1];
  };

  template <typename Number, std::size_t... n>
  const Number* get_view_to_first_element(
      const typename TensorInitializerList<Number, n...>::type& lst)
  {
    // This is to get a view to the first element in the nested array
    return static_cast<const Number*>(static_cast<const void*>(&lst));
  }

  // A helper class to loop over each dimension of a multi-rank tensor in row major order
  template <std::size_t, std::size_t... n>
  struct MultiDimensionalForRowMajor;

  template <std::size_t rank, std::size_t n1>
  struct MultiDimensionalForRowMajor<rank, n1>
  {
    static constexpr void multi_for(auto function, const auto& index)
    {
      for (std::size_t i = 0; i < n1; ++i)
      {
        function(std::tuple_cat(index, std::make_tuple(i)));
      }
    }
  };

  template <std::size_t rank, std::size_t n1, std::size_t... n>
  struct MultiDimensionalForRowMajor<rank, n1, n...>
  {
    static constexpr void multi_for(auto function, const auto& index)
    {
      for (std::size_t i = 0; i < n1; ++i)
      {
        MultiDimensionalForRowMajor<rank, n...>::multi_for(
            function, std::tuple_cat(index, std::make_tuple(i)));
      }
    }
  };

  template <std::size_t... n>
  consteval std::array<std::size_t, (n * ...)> order_type_mapping()
  {
    std::array<std::size_t, (n * ...)> order_type_mapping{};

    std::size_t i = 0;
    MultiDimensionalForRowMajor<sizeof...(n), n...>::multi_for(
        [&](const auto& index)
        {
          order_type_mapping[i] = std::apply(
              [](auto... indices) -> std::size_t
              {
                return get_flat_index<OrderType::column_major, TensorBoundCheck::no_check, n...>(
                    indices...);
              },
              index);

          ++i;
        },
        std::make_tuple());

    return order_type_mapping;
  }

  template <std::size_t... n>
  consteval std::array<std::tuple<decltype(n)...>, (n * ...)> get_array_of_indices()
  {
    std::array<std::tuple<decltype(n)...>, (n * ...)> index_list{};

    std::size_t i = 0;
    MultiDimensionalForRowMajor<sizeof...(n), n...>::multi_for(
        [&](const auto& index)
        {
          index_list[i] = index;
          ++i;
        },
        std::make_tuple());

    return index_list;
  }
}  // namespace Core::LinAlg::Internal

FOUR_C_NAMESPACE_CLOSE

#endif