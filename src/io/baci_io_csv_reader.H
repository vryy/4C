/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief Read csv input
\level 1
*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef BACI_IO_CSV_READER_H
#define BACI_IO_CSV_READER_H

/*-----------------------------------------------------------------------------------------------*/
#include "baci_config.H"

#include "baci_io_string_converter.H"
#include "baci_utils_demangle.H"

#include <fstream>
#include <sstream>
#include <vector>

BACI_NAMESPACE_OPEN

namespace IO
{
  /*!
   * @brief Reads and processes csv file such that a vector of column vectors is returned
   *
   * @param[in] number_of_columns  number of columns in the csv file
   * @param[in] csv_file_path      absolute path to csv file
   * @return vector of column vectors read from csv file
   */
  std::vector<std::vector<double>> ReadCsvAsColumns(
      int number_of_columns, const std::string& csv_file_path);

  /*!
   * @brief Processes csv stream such that a vector of column vectors is returned
   *
   * @param[in] number_of_columns  number of columns in the csv stream
   * @param[in] csv_stream         csv input stream
   * @return vector of column vectors read from csv stream
   */
  std::vector<std::vector<double>> ReadCsvAsColumns(
      int number_of_columns, std::istream& csv_stream);

  /*!
   * @brief Read a csv stream line by line. Parse each line into an object of type T and append the
   * result to the data vector.
   *
   * @param[in] csv_stream csv input stream
   * @param[out] data vector the csv's lines are read into
   * @tparam type T of object one line is read into
   */
  template <typename T>
  void ReadCsvAsLines(std::istream& csv_stream, std::vector<T>& data);

  template <typename T>
  void ReadCsvAsLines(std::istream& csv_stream, std::vector<T>& data)
  {
    // prepare variables
    T line_data;
    std::string line_str;

    // read lines of csv file
    while (std::getline(csv_stream, line_str))
    {
      // do not read in line if it is a header
      if (line_str[0] == '#') continue;

      try
      {
        line_data = IO::StringConverter<T>::Parse(line_str);
        data.emplace_back(line_data);
      }
      catch (...)
      {
        dserror(
            "Could not read line '%s' from csv file. Likely the string's pattern is not "
            "convertible to an object of type %s",
            line_str.c_str(), CORE::UTILS::TryDemangle(typeid(T).name()).c_str());
      }
    }
  }
}  // namespace IO

BACI_NAMESPACE_CLOSE

#endif
