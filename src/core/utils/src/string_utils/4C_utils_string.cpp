// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_utils_string.hpp"

#include <algorithm>
#include <iterator>

FOUR_C_NAMESPACE_OPEN

namespace Core::Utils
{
  std::string trim(const std::string& line)
  {
    // Replace all whitespace character sequences with a single space. Remove leading and trailing
    // spaces.
    std::string result;
    for (auto c : line)
    {
      if (std::isspace(c))
      {
        if (!result.empty() && result.back() != ' ') result.push_back(' ');
      }
      else
      {
        result.push_back(c);
      }
    }

    // Remove trailing space
    if (!result.empty() && result.back() == ' ') result.pop_back();

    return result;
  }

  std::vector<std::string> split(const std::string& input, const std::string& delimiter)
  {
    // Split input string into substrings using the specified delimiter
    std::vector<std::string> result;
    std::string::size_type start = 0;
    std::string::size_type end = input.find(delimiter, start);
    while (end != std::string::npos)
    {
      result.push_back(input.substr(start, end - start));
      start = end + delimiter.length();
      end = input.find(delimiter, start);
    }
    result.push_back(input.substr(start, end));
    return result;
  }

  std::string strip_comment(const std::string& line, const std::string& comment_marker)
  {
    // remove comments
    std::string::size_type loc = line.find(comment_marker);
    std::string newline = line.substr(0, loc);

    // remove trailing and leading whitespaces
    // compact internal whitespaces
    newline = trim(newline);

    return newline;
  }

  std::string to_lower(const std::string& line)
  {
    std::string lower_line;
    lower_line.reserve(line.size());
    std::ranges::transform(
        line, std::back_inserter(lower_line), [](unsigned char c) { return std::tolower(c); });
    return lower_line;
  }

  std::vector<std::string> split_string_list(const std::string& str, const std::string& separator)
  {
    // Keep the currently remaining part of the input string in 'tmp' and
    // keep chopping elements of the list off the front
    std::string tmp = str;

    // Remove whitespace from the end of the string
    while (tmp.size() != 0 && tmp.back() == ' ') tmp.erase(tmp.size() - 1, 1);

    // Split the input list until it is empty. In every iteration, 'tmp' represents the remaining
    // portion of the string after the next separator. Since trailing spaces have already been
    // removed, 'tmp' will eventually be empty if 'str' ended with a separator, even if there was
    // space after the last separator.
    std::vector<std::string> split_list;
    while (tmp.size() != 0)
    {
      std::string name;
      name = tmp;

      if (name.find(separator) != std::string::npos)
      {
        name.erase(name.find(separator), std::string::npos);
        tmp.erase(0, tmp.find(separator) + separator.size());
      }
      else
        tmp = "";

      // Strip spaces from this element's front and end
      while ((name.size() != 0) && (name[0] == ' ')) name.erase(0, 1);
      while (name.size() != 0 && name.back() == ' ') name.erase(name.size() - 1, 1);

      split_list.push_back(name);
    }

    return split_list;
  }

  std::vector<std::string> split_string_list(const std::string& str, const char separator)
  {
    return split_string_list(str, std::string(1, separator));
  }
}  // namespace Core::Utils


FOUR_C_NAMESPACE_CLOSE
