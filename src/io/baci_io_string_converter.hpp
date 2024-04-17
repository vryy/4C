
/*---------------------------------------------------------------------*/
/*! \file

\brief Functionalities to convert strings to arbitrary types

\level 0

*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_IO_STRING_CONVERTER_HPP
#define FOUR_C_IO_STRING_CONVERTER_HPP

#include "baci_config.hpp"

#include "baci_io_utils_reader.hpp"
#include "baci_utils_exceptions.hpp"

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


namespace IO
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
   *   auto data = IO::StringConverter<std::map<int, std::vector<double>>>::Parse(str);
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
   * The StringConverter<T> uses the INTERNAL::Rank<T> class to decide how many different
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
    static T Parse(const std::string &str) = delete;
  };

  namespace INTERNAL
  {
    template <typename TypeVec>
    void CheckDimension(const std::vector<TypeVec> &vec, size_t expectedSize)
    {
      if (vec.size() != expectedSize)
      {
        dserror("Parsed %d values but expected %d", vec.size(), expectedSize);
      }
    }

    template <std::size_t rank, typename TypeArr, std::size_t DimArr>
    constexpr char GetSeparatorAtRank(const std::array<TypeArr, DimArr> &sep_list)
    {
      static_assert(rank <= DimArr,
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
    constexpr int MaxListRank()
    {
      return std::max({StringPatternTraits<Types>::list_rank...});
    }

    /**
     * @brief Recursively determine the max of the map rank of the given types
     */
    template <class... Types>
    constexpr int MaxMapRank()
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
      static constexpr int list_rank = 1 + MaxListRank<ContainedTypes...>();
      static constexpr int map_rank = MaxMapRank<ContainedTypes...>();
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
      static constexpr int list_rank = 1 + MaxListRank<ContainedTypes...>();
      static constexpr int map_rank = 1 + MaxMapRank<ContainedTypes...>();
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

    template <typename T, std::size_t N>
    struct StringPatternTraits<std::array<T, N>> : ListTrait<T>
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
        std::enable_if_t<INTERNAL::StringPatternTraits<T>::is_list_compatible, int> = 0>
    void ParseSplitString(T &t, const std::vector<std::string> &split_str)
    {
      for (const auto &str : split_str)
        t.insert(t.end(), StringConverter<typename T::value_type>::Parse(str));
    };

    /**
     * @brief Parse the split string into an std::array
     *
     * A check is performed to ensure the number of elements found in the split string is N.
     */
    template <typename ValueType, std::size_t N>
    void ParseSplitString(std::array<ValueType, N> &t, const std::vector<std::string> &split_str)
    {
      CheckDimension(split_str, N);
      for (unsigned int i = 0; i < N; ++i)
        t[i] = IO::StringConverter<ValueType>::Parse(split_str[i]);
    }

    /**
     * @brief Parse the split string into an std::pair
     *
     * A check is performed to ensure the number of elements found in the split string is 2.
     */
    template <typename Key, typename Value>
    void ParseSplitString(std::pair<Key, Value> &t, const std::vector<std::string> &split_str)
    {
      CheckDimension(split_str, 2);
      t = std::make_pair(IO::StringConverter<Key>::Parse(split_str[0]),
          StringConverter<Value>::Parse(split_str[1]));
    }

    template <std::size_t index, typename Tuple>
    void ParseSplitStringHelper(Tuple &t, const std::vector<std::string> &split_str)
    {
      if constexpr (index < std::tuple_size<Tuple>::value)
      {
        std::get<index>(t) =
            StringConverter<std::decay_t<decltype(std::get<index>(t))>>::Parse(split_str[index]);
        ParseSplitStringHelper<index + 1>(t, split_str);
      }
    };

    /**
     * @brief Parse the split string into an std::tuple
     *
     * A check is performed to ensure the number of elements found in the split string is equal to
     * the tuple size.
     */
    template <typename... Args>
    void ParseSplitString(std::tuple<Args...> &t, const std::vector<std::string> &split_str)
    {
      CheckDimension(split_str, sizeof...(Args));
      ParseSplitStringHelper<0>(t, split_str);
    };
  }  // namespace INTERNAL

  /**
   * @brief Convert a string into an integer
   */
  template <>
  struct StringConverter<int>
  {
    static int Parse(const std::string &str)
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
    static double Parse(const std::string &str)
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
    static char Parse(const std::string &str)
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
    static bool Parse(const std::string &str)
    {
      std::istringstream is(str);
      bool d;

      if (!(is >> std::boolalpha >> d))
      {
        dserror("String %s cannot be converted to a boolean", str.c_str());
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
  struct StringConverter<T, std::enable_if_t<INTERNAL::StringPatternTraits<T>::is_list_compatible>>
  {
    static T Parse(const std::string &str)
    {
      const char sep = INTERNAL::GetSeparatorAtRank<INTERNAL::StringPatternTraits<T>::list_rank>(
          INTERNAL::default_list_separator);
      auto split_str = DRT::UTILS::SplitStringList(str, sep);

      T t;
      INTERNAL::ParseSplitString(t, split_str);
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
  struct StringConverter<T, std::enable_if_t<INTERNAL::StringPatternTraits<T>::is_map_compatible>>
  {
    static T Parse(const std::string &str)
    {
      T t;

      const char sep_map = INTERNAL::GetSeparatorAtRank<INTERNAL::StringPatternTraits<T>::map_rank>(
          INTERNAL::default_map_separator);
      const char sep_list =
          INTERNAL::GetSeparatorAtRank<INTERNAL::StringPatternTraits<T>::list_rank>(
              INTERNAL::default_list_separator);

      auto split_str = DRT::UTILS::SplitStringList(str, sep_list);

      for (const auto &split_str_i : split_str)
      {
        auto key_val = DRT::UTILS::SplitStringList(split_str_i, sep_map);
        INTERNAL::CheckDimension(key_val, 2);

        t.insert(std::make_pair(StringConverter<typename T::key_type>::Parse(key_val[0]),
            StringConverter<typename T::mapped_type>::Parse(key_val[1])));
      }

      return t;
    }
  };
};  // namespace IO

FOUR_C_NAMESPACE_CLOSE

#endif
