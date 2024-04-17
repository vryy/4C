/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief Read csv input
\level 1
*/
/*-----------------------------------------------------------------------------------------------*/

#include "baci_io_csv_reader.hpp"

#include "baci_utils_exceptions.hpp"

#include <algorithm>
#include <fstream>

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*/
std::vector<std::vector<double>> IO::ReadCsvAsColumns(
    const int number_of_columns, std::istream& csv_stream)
{
  // prepare variables
  std::vector<std::vector<double>> values(number_of_columns);
  std::string line, cell;

  // read lines of csv file
  while (std::getline(csv_stream, line))
  {
    // check if we only have #number_of_columns commas, otherwise throw an error
    if (std::count(line.begin(), line.end(), ',') != number_of_columns - 1)
    {
      dserror(
          "\nInvalid csv file!\n"
          "The csv file has to consist of only %d columns separated by commas but without a "
          "trailing comma! Might also be that your columns do not have the same length!",
          number_of_columns);
    }

    // do not read in line if it is a header
    if (line[0] == '#') continue;

    try
    {
      std::stringstream line_stream(line);
      for (int i = 0; i < number_of_columns; ++i)
      {
        std::getline(line_stream, cell, ',');
        values[i].emplace_back(std::stod(cell));
      }
    }
    catch (...)
    {
      dserror(
          "\nInvalid csv file!\n"
          "Besides a recommended header line starting with '#' the csv file must only consist of "
          "numbers in %d columns separated by commas.",
          number_of_columns);
    }
  }

  return values;
}

/*-----------------------------------------------------------------------------------------------*/
std::vector<std::vector<double>> IO::ReadCsvAsColumns(
    const int number_of_columns, const std::string& csv_file_path)
{
  std::ifstream csv_file_stream(csv_file_path);
  if (csv_file_stream.fail()) dserror("Invalid csv file!");

  return ReadCsvAsColumns(number_of_columns, csv_file_stream);
}

FOUR_C_NAMESPACE_CLOSE
