/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief Read files
\level 1
*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef FOUR_C_IO_FILE_READER_HPP
#define FOUR_C_IO_FILE_READER_HPP

/*-----------------------------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_io_string_converter.hpp"
#include "4C_utils_demangle.hpp"

#include <fstream>
#include <sstream>
#include <vector>

FOUR_C_NAMESPACE_OPEN

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
   * @brief Read a @p input_stream line by line and parse each line into an object of type @p T
   * using `IO::StringConverter<T>::Parse(line_string)`. Return a vector containing all those
   * objects.
   *
   * @param[in] input_stream input stream
   * @tparam T type of object one line is read into
   */
  template <typename T>
  std::vector<T> convert_lines(std::istream& input_stream);

  /*!
   * @brief Read an @p input_stream line by line and parse each line into an object of type @p T
   * using `IO::StringConverter<T>::Parse(line_string)`. The parsed objects are then reduced into
   * another object of @p ReturnType. This process is also known as a `fold` over the data. You
   * can specify which @p operation should be performed by supplying a callable that takes the
   * already accumulated data of type @p ReturnType and the result of parsing a single line into a
   * type @p T.
   *
   * Assume you have an input stream, where each line follows the pattern `"key:val_1,val_2,val_3"`.
   * Those lines can be parsed into objects of type `T = std::map<int, std::array<int, 3>>`.
   * You want to create an `std::map<int,int>` containing the sum of values for each key.
   * Hence, you need an operation (e.g., a lambda function) that creates a `std::map<int, int>`
   * (ReturnType) from objects of type `std::map<int, std::array<int, 3>>` (T) by summing up the
   * array entries:
   *
   * @code {.cpp}
   * auto operation = [](ReducedType acc, T &&next)
   * {
   *   for (const auto &[key, value] : next)
   *   {
   *     acc[key] = value[0] + value[1] + value[2];
   *   }
   *   return acc;
   * };
   * @endcode
   *
   * The desired map could then be read from the input_stream:
   *
   * @code {.cpp}
   * using ReducedType = std::map<int, int>;
   * using T = std::map<int, std::array<int, 3>>;
   * ReducedType converted_data = IO::convert_lines<T, ReducedType>(input_stream, operator);
   * @endcode
   *
   * @param[in] input_stream input stream
   * @param[in] operation Binary operation function object that is apply to create the operated data
   *                      from the parsed data. Its signature must be:
   *                      @code {.cpp}
   *                      ReturnType operation(ReturnType a, T&& b)
   *                      @endcode
   * @tparam T type of object one line is read into
   * @tparam ReturnType type of the result created through the binary operation
   */
  template <typename T, typename ReturnType, typename BinaryOperation>
  ReturnType convert_lines(std::istream& input_stream, BinaryOperation operation);

  template <typename T>
  std::vector<T> convert_lines(std::istream& input_stream)
  {
    return convert_lines<T, std::vector<T>>(input_stream,
        [](std::vector<T> accumulator, T&& next)
        {
          accumulator.emplace_back(std::move(next));
          return accumulator;
        });
  }

  template <typename T, typename ReturnType, typename BinaryOperation>
  ReturnType convert_lines(std::istream& input_stream, BinaryOperation operation)
  {
    std::string line_str;
    ReturnType operated_data;

    // read the input stream line by line
    while (std::getline(input_stream, line_str))
    {
      // do not read in line if it is a header
      if (line_str[0] == '#') continue;

      try
      {
        // parse line string and apply the specified operation on the parsed data
        T parsed_data = IO::StringConverter<T>::Parse(line_str);
        operated_data = operation(std::forward<ReturnType>(operated_data), std::move(parsed_data));
      }
      catch (...)
      {
        FOUR_C_THROW(
            "Could not read line '%s' from input stream. Likely the string's pattern is not "
            "convertible to an object of type %s",
            line_str.c_str(), CORE::UTILS::TryDemangle(typeid(T).name()).c_str());
      }
    }
    return operated_data;
  }
}  // namespace IO

FOUR_C_NAMESPACE_CLOSE

#endif