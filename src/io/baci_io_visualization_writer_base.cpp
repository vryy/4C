/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Define a container that holds the visualization data and the visualization writer

\level 0

*/
/*-----------------------------------------------------------------------------------------------*/

#include "baci_io_visualization_writer_base.H"

BACI_NAMESPACE_OPEN


/**
 *
 */
IO::VisualizationWriterBase::VisualizationWriterBase(IO::VisualizationParameters parameters,
    const Epetra_Comm& comm, std::string visualization_data_name)
    : parameters_(parameters),
      comm_(comm),
      visualization_data_name_(std::move(visualization_data_name))
{
}

/**
 *
 */
void IO::VisualizationWriterBase::WriteVisualizationDataToDisk(
    const VisualizationData& visualization_data, const double visualziation_time,
    const int visualization_step)
{
  visualization_data.ConsistencyCheck();
  InitializeTimeStep(visualziation_time, visualization_step);
  WriteFieldDataToDisk(visualization_data.GetFieldDataMap());
  WriteGeometryToDisk(visualization_data.GetPointCoordinates(),
      visualization_data.GetCellConnectivity(), visualization_data.GetCellOffsets(),
      visualization_data.GetCellTypes(), visualization_data.GetFaceConnectivity(),
      visualization_data.GetFaceOffsets());
  for (const auto& data_name : visualization_data.GetPointDataNames())
  {
    WritePointDataVectorToDisk(visualization_data.GetPointDataVariant(data_name),
        visualization_data.GetPointDataDimension(data_name), data_name);
  }
  for (const auto& data_name : visualization_data.GetCellDataNames())
  {
    WriteCellDataVectorToDisk(visualization_data.GetCellDataVariant(data_name),
        visualization_data.GetCellDataDimension(data_name), data_name);
  }
  FinalizeTimeStep();
}
BACI_NAMESPACE_CLOSE
