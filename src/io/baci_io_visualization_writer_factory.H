/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Factory that creates the visualization writer that is specified in the input file

\level 0

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BACI_IO_VISUALIZATION_WRITER_FACTORY_H
#define BACI_IO_VISUALIZATION_WRITER_FACTORY_H


#include "baci_config.H"

#include "baci_io_visualization_parameters.H"
#include "baci_io_visualization_writer_base.H"

#include <Epetra_Comm.h>

#include <memory>

BACI_NAMESPACE_OPEN

namespace IO
{
  /**
   * @brief Creates the visualization writer that is specified in the parameters object
   */
  [[nodiscard]] std::unique_ptr<VisualizationWriterBase> VisualizationWriterFactory(
      VisualizationParameters parameters, const Epetra_Comm& comm,
      const std::string& visualization_data_name);
}  // namespace IO

BACI_NAMESPACE_CLOSE

#endif