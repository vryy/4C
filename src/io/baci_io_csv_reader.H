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
   * @brief Read a @p csv_stream line by line and parse each line into an object of type @p T using
   * `IO::StringConverter<T>::Parse(line_string)`. Return a vector containing all those objects.
   *
   * @param[in] csv_stream csv input stream
   * @tparam T type of object one line is read into
   */
  template <typename T>
  std::vector<T> ReadCsvAsLines(std::istream& csv_stream);

  /*!
   * @brief Read a @p csv_stream line by line and parse each line into an object of type @p T using
   * `IO::StringConverter<T>::Parse(line_string)`. The parsed objects are then reduced into another
   * object of @p ReturnType. This process is also known as a `fold` over the data.
   * The user can specify which @p operation should be performed by supplying a callable that takes
   * the already accumulated data of type @p ReturnType and the result of parsing a single CSV line
   * into a type @p T.
   *
   * Assume you have a csv stream, where each line follows the pattern `"key:val_1,val_2,val_3"`.
   * Those lines can be parsed into objects of type `T = std::map<int, std::array<int, 3>>`.
   * You want to create an `std::map<int,int>` containing the sum of values for each key.
   * Hence, you need an operation (e.g., a lambda function) that creates a `std::map<int, int>`
   * (ReturnType) from objects of type `std::map<int, std::array<int, 3>>` (T) by summing up the
   * array entries:
   *
   * @code {.cpp}
   * auto operation = [](ReducedType &&acc, T &&next)
   * {
   *   for (const auto &[key, value] : next)
   *   {
   *     acc[key] = value[0] + value[1] + value[2];
   *   }
   *   return acc;
   * };
   * @endcode
   *
   * The desired map could then be read from the csv_stream:
   *
   * @code {.cpp}
   * using ReducedType = std::map<int, int>;
   * using T = std::map<int, std::array<int, 3>>;
   * ReducedType read_data = IO::ReadCsvAsLines<T, ReducedType>(csv_stream, operator);
   * @endcode
   *
   * @param[in] csv_stream csv input stream
   * @param[in] operation Binary operation function object that is apply to create the operated data
   *                      from the parsed data. Its signature must be:
   *                      @code {.cpp}
   *                      ReturnType operation(ReturnType&& a, T&& b)
   *                      @endcode
   * @tparam T type of object one line is read into
   * @tparam ReturnType type of the result created through the binary operation
   */
  template <typename T, typename ReturnType, typename BinaryOperation>
  ReturnType ReadCsvAsLines(std::istream& csv_stream, BinaryOperation operation);

  template <typename T>
  std::vector<T> ReadCsvAsLines(std::istream& csv_stream)
  {
    return ReadCsvAsLines<T, std::vector<T>>(csv_stream,
        [](std::vector<T>&& accumulator, T&& next)
        {
          accumulator.emplace_back(std::move(next));
          return accumulator;
        });
  }

  template <typename T, typename ReturnType, typename BinaryOperation>
  ReturnType ReadCsvAsLines(std::istream& csv_stream, BinaryOperation operation)
  {
    std::string line_str;
    ReturnType operated_data;

    // read lines of csv file
    while (std::getline(csv_stream, line_str))
    {
      // do not read in line if it is a header
      if (line_str[0] == '#') continue;

      try
      {
        // parse line string and apply the specified operation on the parsed data
        T parsed_data = IO::StringConverter<T>::Parse(line_str);
        operated_data = operation(std::move(operated_data), std::move(parsed_data));
      }
      catch (...)
      {
        dserror(
            "Could not read line '%s' from csv file. Likely the string's pattern is not "
            "convertible to an object of type %s",
            line_str.c_str(), CORE::UTILS::TryDemangle(typeid(T).name()).c_str());
      }
    }
    return operated_data;
  }
}  // namespace IO

BACI_NAMESPACE_CLOSE

#endif