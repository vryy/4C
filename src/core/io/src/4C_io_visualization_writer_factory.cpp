/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Define a container that holds the visualization data and the visualization writer

\level 0

*/
/*-----------------------------------------------------------------------------------------------*/

#include "4C_io_visualization_writer_factory.hpp"

#include "4C_io_visualization_writer_vtu_per_rank.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
std::unique_ptr<CORE::IO::VisualizationWriterBase> CORE::IO::VisualizationWriterFactory(
    const VisualizationParameters& parameters, const Epetra_Comm& comm,
    const std::string& visualization_data_name)
{
  if (parameters.writer_ == OutputWriter::vtu_per_rank)
  {
    return std::make_unique<VisualizationWriterVtuPerRank>(
        parameters, comm, visualization_data_name);
  }
  else
  {
    FOUR_C_THROW("You have to select a valid visualization writer in the input file");
  }
}
FOUR_C_NAMESPACE_CLOSE
