/*----------------------------------------------------------------------*/
/*! \file
 * \brief Helpers to read and write dat files
\level 0
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_IO_DAT_FILE_UTILS_HPP
#define FOUR_C_IO_DAT_FILE_UTILS_HPP

#include "4C_config.hpp"

#include "4C_io_linedefinition.hpp"

#include <ostream>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO::DatFileUtils
{

  /**
   * Print a section header padded with dashes to 67 characters.
   */
  void print_section_header(std::ostream& out, const std::string& header);

  /**
   * Print all @p possible_lines into a dat file section with given @p header.
   */
  void print_section(std::ostream& out, const std::string& header,
      const std::vector<Input::LineDefinition>& possible_lines);

  /**
   * Read all lines in a @p section of @p reader that match the @p possible_lines.
   * Every line in the @p section must be readable as one of the @p possible_lines. Otherwise, an
   * exception is thrown.
   *
   * @see read_matching_lines_in_section()
   */
  std::vector<Input::LineDefinition> read_all_lines_in_section(Core::IO::DatFileReader& reader,
      const std::string& section, const std::vector<Input::LineDefinition>& possible_lines);


  /**
   * Read only lines in a @p section of @p reader that match the @p possible_lines. This implies
   * that, potentially, no lines are read at all, resulting in an empty returned vector. In
   * addition to the vector of parsed lines, the second returned value contains all unparsed input
   * lines.
   *
   * @see read_all_lines_in_section()
   */
  std::pair<std::vector<Input::LineDefinition>, std::vector<std::string>>
  read_matching_lines_in_section(Core::IO::DatFileReader& reader, const std::string& section,
      const std::vector<Input::LineDefinition>& possible_lines);

}  // namespace Core::IO::DatFileUtils


FOUR_C_NAMESPACE_CLOSE

#endif
