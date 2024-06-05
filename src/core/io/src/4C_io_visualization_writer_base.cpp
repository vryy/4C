/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Define a container that holds the visualization data and the visualization writer

\level 0

*/
/*-----------------------------------------------------------------------------------------------*/

#include "4C_io_visualization_writer_base.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
CORE::IO::VisualizationWriterBase::VisualizationWriterBase(
    CORE::IO::VisualizationParameters parameters, const Epetra_Comm& comm,
    std::string visualization_data_name)
    : parameters_(std::move(parameters)),
      comm_(comm),
      visualization_data_name_(std::move(visualization_data_name))
{
}

/**
 *
 */
void CORE::IO::VisualizationWriterBase::write_visualization_data_to_disk(
    const VisualizationData& visualization_data, const double visualziation_time,
    const int visualization_step)
{
  visualization_data.ConsistencyCheck();
  InitializeTimeStep(visualziation_time, visualization_step);
  write_field_data_to_disk(visualization_data.GetFieldDataMap());
  WriteGeometryToDisk(visualization_data.GetPointCoordinates(),
      visualization_data.GetCellConnectivity(), visualization_data.GetCellOffsets(),
      visualization_data.GetCellTypes(), visualization_data.GetFaceConnectivity(),
      visualization_data.GetFaceOffsets());
  for (const auto& data_name : visualization_data.GetPointDataNames())
  {
    write_point_data_vector_to_disk(visualization_data.GetPointDataVariant(data_name),
        visualization_data.get_point_data_dimension(data_name), data_name);
  }
  for (const auto& data_name : visualization_data.GetCellDataNames())
  {
    write_cell_data_vector_to_disk(visualization_data.GetCellDataVariant(data_name),
        visualization_data.get_cell_data_dimension(data_name), data_name);
  }
  FinalizeTimeStep();
}
FOUR_C_NAMESPACE_CLOSE
