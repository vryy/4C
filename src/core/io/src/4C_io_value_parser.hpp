// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_VALUE_PARSER_HPP
#define FOUR_C_IO_VALUE_PARSER_HPP

#include "4C_config.hpp"

#include "4C_utils_demangle.hpp"
#include "4C_utils_exceptions.hpp"

#include <array>
#include <string>
#include <string_view>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  namespace Internal
  {
    void parse(std::string_view string, int& value);
    void parse(std::string_view string, double& value);
    void parse(std::string_view string, std::string& value);
  }  // namespace Internal

  /**
   * A helper to parse values as defined in the .dat file format into C++ data. This
   * is a low-level class intended for use inside more user-friendly input mechanisms. Based on a
   * input string_view, it allows to read values of different types in sequence and validate that
   * the format matches. Note that the whole class works on string_views, so the original data must
   * outlive the parser.
   */
  class ValueParser
  {
   public:
    /**
     * Set up the ValueParser and give an optional additional scope message. This information is
     * prepended to all error messages. Example:
     *
     * @code
     *   ValueParser parser(line, "While reading section MY PARAMETERS: ");
     * @endcode
     *
     * The input @p line is a string_view onto the data to parse. The parser will not copy the data,
     * so the original data must outlive the parser.
     */
    ValueParser(std::string_view line, std::string user_scope_message);

    //! Read the next string and ensure it matches the expectation.
    void consume(const std::string& expected);

    //! Read a single value of given type.
    template <typename T>
    T read()
    {
      std::string_view read_string = advance_token();
      T read_object;
      Internal::parse(read_string, read_object);
      return read_object;
    }

    //! Read an array of a given type.
    template <typename T, size_t n>
    std::array<T, n> read_array()
    {
      std::array<T, n> array;
      for (size_t i = 0; i < n; ++i) array[i] = read<T>();

      return array;
    }

    /**
     * Return the next token without consuming it. The return token may be empty if the parser has
     * reached the end of the input string.
     */
    [[nodiscard]] std::string_view peek() const;

    //! Check if this parser reached the end of the input string.
    [[nodiscard]] bool at_end() const;

    /**
     * Get anything that hasn't been parsed by previous calls to consume() or read(). Note that
     * the returned string_view is a view into the original input string, so it is only valid as
     * long as the original string is valid. The returned string_view is empty if the parser has
     * reached the end of the input string.
     *
     * @note The returned string_view might contain leading whitespace.
     */
    [[nodiscard]] std::string_view get_unparsed_remainder() const;

   private:
    std::string_view advance_token();

    //! The data to parse from.
    std::string_view line_;

    //! The current position in the line_.
    std::size_t current_index_{0};

    //! Prepend a user message for better error messages.
    std::string user_scope_{};
  };
}  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE

#endif
