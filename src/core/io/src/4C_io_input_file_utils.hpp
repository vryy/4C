// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_INPUT_FILE_UTILS_HPP
#define FOUR_C_IO_INPUT_FILE_UTILS_HPP

#include "4C_config.hpp"

#include "4C_io_linedefinition.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <ostream>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO::InputFileUtils
{

  /**
   * Print a section header padded with dashes to 67 characters.
   */
  void print_section_header(std::ostream& out, const std::string& header);

  /**
   * Print the values of parameter entries from a dat file based on the provided parameter list.
   *
   * This function prints key-value pairs in the dat file format, including parameter names and
   * their values. It processes sublists and parameters, formatting them appropriately.
   *
   * @param stream      The output stream to which the information will be printed.
   * @param list        The parameter list containing the parameters and sublists to be printed.
   * @param comment     A flag indicating whether to print comments (default is true).
   */
  void print_dat(std::ostream& stream, const Teuchos::ParameterList& list, bool comment = true);

  /**
   * Return true if the @p list contains any parameter that has whitespace in the key name.
   *
   * @note This is needed for the NOX parameters whose keywords and value have white spaces and
   * thus '=' are inserted to distinguish them.
   */
  bool need_to_print_equal_sign(const Teuchos::ParameterList& list);

  /**
   * Print all @p possible_lines into a dat file section with given @p header.
   */
  void print_section(std::ostream& out, const std::string& header,
      const std::vector<Input::LineDefinition>& possible_lines);

  /**
   * Read all lines in a @p section of @p input that match the @p possible_lines.
   * Every line in the @p section must be readable as one of the @p possible_lines. Otherwise, an
   * exception is thrown.
   *
   * @see read_matching_lines_in_section()
   */
  std::vector<Input::LineDefinition> read_all_lines_in_section(Core::IO::InputFile& input,
      const std::string& section, const std::vector<Input::LineDefinition>& possible_lines);


  /**
   * Read only lines in a @p section of @p input that match the @p possible_lines. This implies
   * that, potentially, no lines are read at all, resulting in an empty returned vector. In
   * addition to the vector of parsed lines, the second returned value contains all unparsed input
   * lines.
   *
   * @see read_all_lines_in_section()
   */
  std::pair<std::vector<Input::LineDefinition>, std::vector<std::string>>
  read_matching_lines_in_section(Core::IO::InputFile& input, const std::string& section,
      const std::vector<Input::LineDefinition>& possible_lines);

}  // namespace Core::IO::InputFileUtils


FOUR_C_NAMESPACE_CLOSE

#endif
