// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_FILE_READER_HPP
#define FOUR_C_IO_FILE_READER_HPP

/*-----------------------------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_utils_demangle.hpp"

#include <fstream>
#include <sstream>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  /*!
   * @brief Reads and processes csv file such that a vector of column vectors is returned
   *
   * @param[in] number_of_columns  number of columns in the csv file
   * @param[in] csv_file_path      absolute path to csv file
   * @return vector of column vectors read from csv file
   */
  FOUR_C_API(FOUR_C_CORE)
  std::vector<std::vector<double>> read_csv_as_columns(
      int number_of_columns, const std::string& csv_file_path);

  /*!
   * @brief Processes csv stream such that a vector of column vectors is returned
   *
   * @param[in] number_of_columns  number of columns in the csv stream
   * @param[in] csv_stream         csv input stream
   * @return vector of column vectors read from csv stream
   */
  FOUR_C_API(FOUR_C_CORE)
  std::vector<std::vector<double>> read_csv_as_columns(
      int number_of_columns, std::istream& csv_stream);

}  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE

#endif