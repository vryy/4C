// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_STRING_CONVERTER_HPP
#define FOUR_C_IO_STRING_CONVERTER_HPP

#include "4C_config.hpp"

#include "4C_utils_exceptions.hpp"
#include "4C_utils_string.hpp"

#include <algorithm>
#include <array>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

FOUR_C_NAMESPACE_OPEN


namespace Core::IO
{
  /*!
   * @brief Convert a string matching a predefined pattern into an object of type T.
   *
   * StringConverter<T> can read strings following a predefined pattern. The power of this
   * class comes form the fact that this pattern is determined at compile-time based on the template
   * type T. Many data types T can be expressed by such a pattern.
   *
   * For example, this class allows you to read a string as follows:
   *
   * @code
   *   std::string str = "1: 1.0, 2.0, 3.0; 2: 4.0, 5.0, 6.0; 3: 7.0, 8.0, 9.0"
   *   auto data = Core::IO::StringConverter<std::map<int, std::vector<double>>>::Parse(str);
   * @endcode
   *
   * We distinguish map-compatible and list-compatible data types and the associated map-like
   * {':', '=', '@', '#'} and list-like separators {',', ';', '|', '%'}. Each non-arithmetic type T
   * is either map- or list-compatible, and thus can be described using those separators.
   * - map-compatible types are std::map, std::unordered_map, and std::multimap
   * - list-compatible types are std::vector, std::set, std::list, std::pair, std::array,
   * std::tuple
   *
   * A maximum of four nested map- and list-compatible types is currently supported.
   *
   *
   * ## Implementation details
   *
   * The StringConverter<T> uses the Internal::Rank<T> class to decide how many different
   * separators are required to convert the string to the given (nested) type T.
   * Arithmetic types have rank 0, for list-compatible types the list-rank increments by one for
   * each level, for map-compatible types, both, the map- and list-rank increment by one for each
   * nested level. To determine the string pattern compatible with type T, start from the innermost
   * type level. The assigned rank (and thus the required separator) increases with each level up.
   *
   * As an example, assume the type T is a map with keys of type integer and values of type vector
   * of doubles. Starting at the innermost level, the arithmetic types' (int, double) ranks are 0.
   * Moving one level up, the vector is list-compatible, its list-rank is incremented by 1,
   * resulting in rank 1 (e.g. "VecEle1, VecEle2, VecEle3"). Hence the appropriate separator is ","
   * (the list-compatible separator for rank 1). Again, moving one level up, the map is
   * map-compatible, both, its map- and list-rank is incremented by 1, resulting in list-rank 2 and
   * map-rank 1. The list-like collection of key-value pairs thus needs to be separated by ";" (rank
   * 2), while the key and the value are separated by ":" (rank 1). We thus arrive at the pattern
   * "Key1: VecEle1, VecEle2, VecEle3; Key2: VecEle1, VecEle2, VecEle3; ..."
   *
   *
   * @note More supported types may be added via template specialisation.
   */
  template <typename T, class Enable = void>
  struct StringConverter
  {
    static T Parse(const std::string& str) = delete;
  };

  namespace Internal
  {
    template <typename TypeVec>
    void check_dimension(const std::vector<TypeVec>& vec, size_t expectedSize)
    {
      if (vec.size() != expectedSize)
      {
        FOUR_C_THROW("Parsed %d values but expected %d", vec.size(), expectedSize);
      }
    }

    template <std::size_t rank, typename TypeArr, std::size_t dim_arr>
    constexpr char get_separator_at_rank(const std::array<TypeArr, dim_arr>& sep_list)
    {
      static_assert(rank <= dim_arr,
          "Requesting a separator of too high rank. Your data structure is too deeply nested.");
      return sep_list[rank - 1];
    }

    //! @name Default separators for list and map compatible types for the ranks 1 to 4.
    //! The separator for a given object is selected according to the rank of the list/map object.
    //! The separator with the lowest index is used for the innermost structure. No separator is
    //! required for rank 0.
    //@{
    constexpr std::array<const char, 4> default_list_separator{',', ';', '|', '%'};
    constexpr std::array<const char, 4> default_map_separator{':', '=', '@', '#'};
    //@}

    /**
     * A helper struct to figure out whether a type behaves like a list or a map and determine
     * the associated rank.
     *
     * @tparam T The relevant type.
     * @tparam Enable A parameter to leverage SFINAE and selectively enable specializations.
     */
    template <typename T, typename Enable = void>
    struct StringPatternTraits;

    /**
     * @brief Recursively determine the max of the list rank of the given types
     */
    template <class... Types>
    constexpr int max_list_rank()
    {
      return std::max({StringPatternTraits<Types>::list_rank...});
    }

    /**
     * @brief Recursively determine the max of the map rank of the given types
     */
    template <class... Types>
    constexpr int max_map_rank()
    {
      return std::max({StringPatternTraits<Types>::map_rank...});
    }

    /**
     * If a type has the ListTrait, it behaves like a list and increases the overall list rank by 1
     * compared to its contents.
     */
    template <typename... ContainedTypes>
    struct ListTrait
    {
      static constexpr int list_rank = 1 + max_list_rank<ContainedTypes...>();
      static constexpr int map_rank = max_map_rank<ContainedTypes...>();
      static constexpr bool is_list_compatible = true;
      static constexpr bool is_map_compatible = false;
    };

    /**
     * If a type has the MapTrait, it behaves like a map and increases the overall list _and_ map
     * rank by 1 compared to its contents.
     */
    template <typename... ContainedTypes>
    struct MapTrait
    {
      static constexpr int list_rank = 1 + max_list_rank<ContainedTypes...>();
      static constexpr int map_rank = 1 + max_map_rank<ContainedTypes...>();
      static constexpr bool is_list_compatible = false;
      static constexpr bool is_map_compatible = true;
    };

    /**
     * Arithmetic types are neither lists nor maps and have no associated rank.
     */
    template <typename T>
    struct StringPatternTraits<T, std::enable_if_t<std::is_arithmetic_v<T>>>
    {
      static constexpr int list_rank = 0;
      static constexpr int map_rank = 0;
      static constexpr bool is_list_compatible = false;
      static constexpr bool is_map_compatible = false;
    };

    template <typename T, typename Alloc>
    struct StringPatternTraits<std::vector<T, Alloc>> : ListTrait<T>
    {
    };

    template <typename T, typename Alloc>
    struct StringPatternTraits<std::list<T, Alloc>> : ListTrait<T>
    {
    };

    template <typename T, std::size_t n>
    struct StringPatternTraits<std::array<T, n>> : ListTrait<T>
    {
    };

    template <typename... Ts>
    struct StringPatternTraits<std::tuple<Ts...>> : ListTrait<Ts...>
    {
    };

    template <typename First, typename Second>
    struct StringPatternTraits<std::pair<First, Second>> : ListTrait<First, Second>
    {
    };

    template <typename Key, typename Value, typename Compare, typename Alloc>
    struct StringPatternTraits<std::map<Key, Value, Compare, Alloc>> : MapTrait<Key, Value>
    {
    };

    template <typename Key, typename Value, typename Compare, typename Alloc>
    struct StringPatternTraits<std::unordered_map<Key, Value, Compare, Alloc>>
        : MapTrait<Key, Value>
    {
    };

    template <typename Key, typename Value, typename Compare, typename Alloc>
    struct StringPatternTraits<std::multimap<Key, Value, Compare, Alloc>> : MapTrait<Key, Value>
    {
    };

    /**
     * @brief Parse the split string into an std::vector or std::list.
     */
    template <class T,
        std::enable_if_t<Internal::StringPatternTraits<T>::is_list_compatible, int> = 0>
    void parse_split_string(T& t, const std::vector<std::string>& split_str)
    {
      for (const auto& str : split_str)
        t.insert(t.end(), StringConverter<typename T::value_type>::parse(str));
    };

    /**
     * @brief Parse the split string into an std::array
     *
     * A check is performed to ensure the number of elements found in the split string is N.
     */
    template <typename ValueType, std::size_t n>
    void parse_split_string(std::array<ValueType, n>& t, const std::vector<std::string>& split_str)
    {
      check_dimension(split_str, n);
      for (unsigned int i = 0; i < n; ++i)
        t[i] = Core::IO::StringConverter<ValueType>::parse(split_str[i]);
    }

    /**
     * @brief Parse the split string into an std::pair
     *
     * A check is performed to ensure the number of elements found in the split string is 2.
     */
    template <typename Key, typename Value>
    void parse_split_string(std::pair<Key, Value>& t, const std::vector<std::string>& split_str)
    {
      check_dimension(split_str, 2);
      t = std::make_pair(Core::IO::StringConverter<Key>::parse(split_str[0]),
          StringConverter<Value>::parse(split_str[1]));
    }

    template <std::size_t index, typename Tuple>
    void parse_split_string_helper(Tuple& t, const std::vector<std::string>& split_str)
    {
      if constexpr (index < std::tuple_size<Tuple>::value)
      {
        std::get<index>(t) =
            StringConverter<std::decay_t<decltype(std::get<index>(t))>>::parse(split_str[index]);
        parse_split_string_helper<index + 1>(t, split_str);
      }
    };

    /**
     * @brief Parse the split string into an std::tuple
     *
     * A check is performed to ensure the number of elements found in the split string is equal to
     * the tuple size.
     */
    template <typename... Args>
    void parse_split_string(std::tuple<Args...>& t, const std::vector<std::string>& split_str)
    {
      check_dimension(split_str, sizeof...(Args));
      parse_split_string_helper<0>(t, split_str);
    };
  }  // namespace Internal

  /**
   * @brief Convert a string into an integer
   */
  template <>
  struct StringConverter<int>
  {
    static int parse(const std::string& str)
    {
      std::istringstream is(str);
      int i;
      is >> i;
      return i;
    }
  };

  /**
   * @brief Convert a string into a double.
   */
  template <>
  struct StringConverter<double>
  {
    static double parse(const std::string& str)
    {
      std::istringstream is(str);
      double d;
      is >> d;
      return d;
    }
  };

  /**
   * @brief Convert a string into a char.
   */
  template <>
  struct StringConverter<char>
  {
    static char parse(const std::string& str)
    {
      std::istringstream is(str);
      char c;
      is >> c;
      return c;
    }
  };

  /**
   * @brief Convert a string into a bool.
   */
  template <>
  struct StringConverter<bool>
  {
    static bool parse(const std::string& str)
    {
      std::istringstream is(str);
      bool d;

      if (!(is >> std::boolalpha >> d))
      {
        FOUR_C_THROW("String %s cannot be converted to a boolean", str.c_str());
      }

      return d;
    }
  };

  /**
   * @brief Convert a string into a list-compatible type T.
   *
   * For list-compatible types (e.g., std::vector, std::set, std::list, std::array, std::pair,
   * std::tuple), the elements need to be separated by a list-separator (example: 1,2,3,4,5). The
   * size of the created list-compatible object is either determined dynamically depending on the
   * string definition (std::vector, std::set, std::list) or given by the size of a fixed-sized type
   * T (std::array, std::pair, std::tuple).
   */
  template <class T>
  struct StringConverter<T, std::enable_if_t<Internal::StringPatternTraits<T>::is_list_compatible>>
  {
    static T parse(const std::string& str)
    {
      const char sep = Internal::get_separator_at_rank<Internal::StringPatternTraits<T>::list_rank>(
          Internal::default_list_separator);
      auto split_str = Core::Utils::split_string_list(str, sep);

      T t;
      Internal::parse_split_string(t, split_str);
      return t;
    }
  };

  /**
   * @brief Convert a string into a map-compatible type T.
   *
   * For map-compatible types, the key-value pairs need to be separated by a list-separator
   * (example: 1:"abc",2:"def",3:"hij"), while the key and the value are separated by a
   * map-separator (example: 1:"abc"). A checks are performed to ensure the number of string
   * elements after splitting for the map-separator is 2, since we create one key and one value.
   * The number of key-value pairs in the map is determined dynamically depending on the string
   * definition.
   */
  template <class T>
  struct StringConverter<T, std::enable_if_t<Internal::StringPatternTraits<T>::is_map_compatible>>
  {
    static T parse(const std::string& str)
    {
      T t;

      const char sep_map =
          Internal::get_separator_at_rank<Internal::StringPatternTraits<T>::map_rank>(
              Internal::default_map_separator);
      const char sep_list =
          Internal::get_separator_at_rank<Internal::StringPatternTraits<T>::list_rank>(
              Internal::default_list_separator);

      auto split_str = Core::Utils::split_string_list(str, sep_list);

      for (const auto& split_str_i : split_str)
      {
        auto key_val = Core::Utils::split_string_list(split_str_i, sep_map);
        Internal::check_dimension(key_val, 2);

        t.insert(std::make_pair(StringConverter<typename T::key_type>::parse(key_val[0]),
            StringConverter<typename T::mapped_type>::parse(key_val[1])));
      }

      return t;
    }
  };
};  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE

#endif
