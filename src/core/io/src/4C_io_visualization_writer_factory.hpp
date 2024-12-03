// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_VISUALIZATION_WRITER_FACTORY_HPP
#define FOUR_C_IO_VISUALIZATION_WRITER_FACTORY_HPP


#include "4C_config.hpp"

#include "4C_io_visualization_parameters.hpp"
#include "4C_io_visualization_writer_base.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  /**
   * @brief Creates the visualization writer that is specified in the parameters object
   */
  [[nodiscard]] std::unique_ptr<VisualizationWriterBase> visualization_writer_factory(
      const VisualizationParameters& parameters, MPI_Comm comm,
      const std::string& visualization_data_name);
}  // namespace Core::IO

FOUR_C_NAMESPACE_CLOSE

#endif