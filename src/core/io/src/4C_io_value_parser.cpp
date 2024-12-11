// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_io_value_parser.hpp"

FOUR_C_NAMESPACE_OPEN

namespace
{
  //! Helper to find the next token and update the given @p index into the line.
  std::string_view advance_token_impl(std::string_view line, std::size_t& index)
  {
    // Skip whitespace
    while (index < line.size() && std::isspace(line[index])) ++index;

    // Find the end of the token
    std::size_t start_of_token = index;
    while (index < line.size() && !std::isspace(line[index])) ++index;

    return line.substr(start_of_token, index - start_of_token);
  }
}  // namespace


void Core::IO::Internal::parse(std::string_view string, int& value)
{
  std::string string_copy(string);
  std::size_t end;
  value = std::stoi(string_copy.data(), &end);
  if (end != string_copy.size())
    FOUR_C_THROW("Could not parse '%s' as an integer value.", string_copy.c_str());
}

void Core::IO::Internal::parse(std::string_view string, double& value)
{
  std::string string_copy(string);
  std::size_t end;
  value = std::stod(string_copy.data(), &end);
  if (end != string_copy.size())
    FOUR_C_THROW("Could not parse '%s' as a double value.", string_copy.c_str());
}

void Core::IO::Internal::parse(std::string_view string, std::string& value)
{
  value = std::string(string);
}

Core::IO::ValueParser::ValueParser(std::string_view line, std::string user_scope_message)
    : line_(line), user_scope_(std::move(user_scope_message))
{
}


void Core::IO::ValueParser::consume(const std::string& expected)
{
  std::string_view read_string = advance_token();
  if (read_string != std::string_view(expected))
    FOUR_C_THROW("%sCould not read expected string '%s'.", user_scope_.c_str(), expected.c_str());
}


std::string_view Core::IO::ValueParser::peek() const
{
  // Copy the current index to avoid modifying the parser state
  std::size_t temp_index = current_index_;
  return advance_token_impl(line_, temp_index);
}


bool Core::IO::ValueParser::at_end() const { return current_index_ == line_.size(); }


std::string_view Core::IO::ValueParser::get_unparsed_remainder() const
{
  return line_.substr(current_index_);
}


std::string_view Core::IO::ValueParser::advance_token()
{
  auto token = advance_token_impl(line_, current_index_);

  FOUR_C_ASSERT_ALWAYS(!token.empty(), "%sExpected more tokens, but reached the end of the line.",
      user_scope_.c_str());

  return token;
}


FOUR_C_NAMESPACE_CLOSE
