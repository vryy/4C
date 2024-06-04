/*----------------------------------------------------------------------*/
/*! \file
 * \brief Helpers to read and write dat files
\level 0
*/
/*----------------------------------------------------------------------*/

#include "4C_io_dat_file_utils.hpp"

FOUR_C_NAMESPACE_OPEN

void IO::DatFileUtils::print_section_header(std::ostream& out, const std::string& header)
{
  constexpr std::size_t max_line_width = 65ul;
  FOUR_C_THROW_UNLESS(header.length() <= max_line_width, "Header '%s' too long", header.c_str());

  out << "--";
  out << std::string(std::max(max_line_width - header.length(), 0ul), '-');
  out << header << '\n';
}



void IO::DatFileUtils::print_section(std::ostream& out, const std::string& header,
    const std::vector<INPUT::LineDefinition>& possible_lines)
{
  print_section_header(out, header);

  for (const auto& line : possible_lines)
  {
    out << "// ";
    line.Print(out);
    out << '\n';
  }
}


std::vector<INPUT::LineDefinition> IO::DatFileUtils::read_all_lines_in_section(
    INPUT::DatFileReader& reader, const std::string& section,
    const std::vector<INPUT::LineDefinition>& possible_lines)
{
  auto [parsed_lines, unparsed_lines] =
      read_matching_lines_in_section(reader, section, possible_lines);

  // In this function, encountering unparsed lines is an error, so construct a nice message.
  if (unparsed_lines.size() > 0)
  {
    std::stringstream out;
    out << "Read failed in section " << std::quoted(section) << '\n';
    for (const auto& unparsed : unparsed_lines)
    {
      out << "  " << std::quoted(unparsed) << '\n';
    }
    out << "Valid lines are:\n";
    std::for_each(possible_lines.begin(), possible_lines.end(),
        [&](const INPUT::LineDefinition& def)
        {
          def.Print(out);
          out << '\n';
        });
    FOUR_C_THROW(out.str().c_str());
  }

  return parsed_lines;
}


std::pair<std::vector<INPUT::LineDefinition>, std::vector<std::string>>
IO::DatFileUtils::read_matching_lines_in_section(INPUT::DatFileReader& reader,
    const std::string& section, const std::vector<INPUT::LineDefinition>& possible_lines)
{
  const std::vector<const char*> lines_in_section = reader.Section("--" + section);

  std::vector<std::string> unparsed_lines;
  std::vector<INPUT::LineDefinition> parsed_lines;

  const auto process_line = [&](const std::string& input_line)
  {
    for (const auto& definition : possible_lines)
    {
      std::stringstream l{input_line};

      // Make a copy that potentially gets filled by the Read.
      auto parsed_definition = definition;
      if (parsed_definition.Read(l))
      {
        parsed_lines.emplace_back(std::move(parsed_definition));
        return;
      }
    }
    unparsed_lines.emplace_back(input_line);
  };

  for (const auto& input_line : lines_in_section)
  {
    process_line(input_line);
  }

  FOUR_C_ASSERT(
      unparsed_lines.size() + parsed_lines.size() == lines_in_section.size(), "Internal error.");

  return {parsed_lines, unparsed_lines};
}
FOUR_C_NAMESPACE_CLOSE
