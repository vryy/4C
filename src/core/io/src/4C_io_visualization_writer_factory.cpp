#include "4C_io_visualization_writer_factory.hpp"

#include "4C_io_visualization_writer_vtu_per_rank.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
std::unique_ptr<Core::IO::VisualizationWriterBase> Core::IO::visualization_writer_factory(
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
