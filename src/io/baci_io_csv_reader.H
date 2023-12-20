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
  std::vector<std::vector<double>> ReadCsv(int number_of_columns, const std::string& csv_file_path);

  /*!
   * @brief Processes csv stream such that a vector of column vectors is returned
   *
   * @param[in] number_of_columns  number of columns in the csv stream
   * @param[in] csv_stream         csv input stream
   * @return vector of column vectors read from csv stream
   */
  std::vector<std::vector<double>> ReadCsv(int number_of_columns, std::istream& csv_stream);

}  // namespace IO

BACI_NAMESPACE_CLOSE

#endif
