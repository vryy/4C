/*---------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for readers

\level 0


*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_UTILS_STRING_HPP
#define FOUR_C_UTILS_STRING_HPP

#if (BOOST_MAJOR_VERSION == 1) && (BOOST_MINOR_VERSION >= 47)
#include "4C_config.hpp"

#include <boost/algorithm/string/case_conv.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim_all.hpp>

#include <cctype>
#include <cstring>
#endif

FOUR_C_NAMESPACE_OPEN


namespace CORE::UTILS
{
  /*!
   * @brief Remove all leading and trailing whitespaces from a string.
   *
   * Note: consecutive whitespaces inside the std::string will be reduced to a single space.
   */
  std::string Trim(const std::string& line);

  /*!
   * @brief Split the @p input string into multiple substrings between the @p delimiter.
   */
  std::vector<std::string> Split(const std::string& input, const std::string& delimiter);

  /*!
   * @brief Remove comments, trailing and leading whitespaces
   *
   * @param[in] line arbitrary string
   * @result same string stripped
   */
  std::string StripComment(const std::string& line);

  /*!
   * @brief Convert all characters in a string into lower case (wrapper to the corresponding routine
   * in boost::algorithm)
   *
   *  @param[in] line arbitrary string
   *  @result same string in all lower case
   */
  std::string ToLower(const std::string& line);

  /*!
   * @brief Split a string into a vector of strings by a given separator
   *
   *  @param[in] str input string
   *  @param[in] separator separator string the string is split by
   *  @result vector of strings
   */
  std::vector<std::string> SplitStringList(const std::string& str, const std::string& separator);

  /*!
   * @brief Split a string into a vector of strings by a given separator
   *
   *  @param[in] str input string
   *  @param[in] separator separator character the string is split by
   *  @result vector of strings
   */
  std::vector<std::string> SplitStringList(const std::string& str, const char separator);
}  // namespace CORE::UTILS

FOUR_C_NAMESPACE_CLOSE

#endif
