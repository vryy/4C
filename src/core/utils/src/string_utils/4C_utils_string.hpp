// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_UTILS_STRING_HPP
#define FOUR_C_UTILS_STRING_HPP

#include "4C_config.hpp"

#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN


namespace Core::Utils
{
  /*!
   * @brief Remove all leading and trailing whitespaces from a string.
   *
   * Note: consecutive whitespaces inside the std::string will be reduced to a single space.
   */
  FOUR_C_API(FOUR_C_CORE) std::string trim(const std::string& line);

  /*!
   * @brief Split the @p input string into multiple substrings between the @p delimiter.
   */
  FOUR_C_API(FOUR_C_CORE)
  std::vector<std::string> split(const std::string& input, const std::string& delimiter);

  /*!
   * @brief Remove comments, trailing and leading whitespaces
   *
   * @param[in] line arbitrary string
   * @param[in] comment_marker The marker that indicates a comment string
   * @result same string stripped
   */
  FOUR_C_API(FOUR_C_CORE)
  std::string strip_comment(const std::string& line, const std::string& comment_marker = "//");

  /*!
   * @brief Convert all characters in a string into lower case (wrapper to the corresponding routine
   * in boost::algorithm)
   *
   *  @param[in] line arbitrary string
   *  @result same string in all lower case
   */
  FOUR_C_API(FOUR_C_CORE) std::string to_lower(const std::string& line);

  /*!
   * @brief Split a string into a vector of strings by a given separator
   *
   *  @param[in] str input string
   *  @param[in] separator separator string the string is split by
   *  @result vector of strings
   */
  FOUR_C_API(FOUR_C_CORE)
  std::vector<std::string> split_string_list(const std::string& str, const std::string& separator);

  /*!
   * @brief Split a string into a vector of strings by a given separator
   *
   *  @param[in] str input string
   *  @param[in] separator separator character the string is split by
   *  @result vector of strings
   */
  FOUR_C_API(FOUR_C_CORE)
  std::vector<std::string> split_string_list(const std::string& str, const char separator);
}  // namespace Core::Utils

FOUR_C_NAMESPACE_CLOSE

#endif
