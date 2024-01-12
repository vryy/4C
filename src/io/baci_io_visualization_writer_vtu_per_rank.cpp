/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Define a container that holds the visualization data and the visualization writer

\level 0

*/
/*-----------------------------------------------------------------------------------------------*/

#include "baci_io_visualization_writer_vtu_per_rank.H"

#include "baci_global_data.H"
#include "baci_inpar_IO_runtime_output.H"
#include "baci_io_control.H"
#include "baci_io_visualization_data.H"

BACI_NAMESPACE_OPEN


/**
 *
 */
IO::VisualizationWriterVtuPerRank::VisualizationWriterVtuPerRank(
    IO::VisualizationParameters parameters, const Epetra_Comm& comm,
    std::string visualization_data_name)
    : VisualizationWriterBase(parameters, comm, std::move(visualization_data_name))
{
  const auto& output_control = GLOBAL::Problem::Instance()->OutputControlFile();
  vtu_writer_.Initialize(comm_.MyPID(), comm_.NumProc(),
      std::pow(10, IO::GetTotalDigitsToReserveInTimeStep(parameters)),
      output_control->DirectoryName(), (output_control->FileNameOnlyPrefix() + "-vtk-files"),
      visualization_data_name_, output_control->RestartName(),
      GLOBAL::Problem::Instance()->Restart(),
      parameters.data_format_ == INPAR::IO_RUNTIME_OUTPUT::OutputDataFormat::binary);
}


/**
 *
 */
void IO::VisualizationWriterVtuPerRank::InitializeTimeStep(
    const double visualziation_time, const int visualization_step)
{
  vtu_writer_.ResetTimeAndTimeStep(visualziation_time, visualization_step);
  vtu_writer_.ResetGeometryName(visualization_data_name_);

  vtu_writer_.InitializeVtkFileStreamsForNewGeometryAndOrTimeStep();
  vtu_writer_.WriteVtkHeaders();
}

/**
 *
 */
void IO::VisualizationWriterVtuPerRank::WriteFieldDataToDisk(
    const std::map<std::string, visualization_vector_type_variant>& field_data_map)
{
  vtu_writer_.WriteVtkFieldDataAndOrTimeAndOrCycle(field_data_map);
}

/**
 *
 */
void IO::VisualizationWriterVtuPerRank::WriteGeometryToDisk(
    const std::vector<double>& point_coordinates,
    const std::vector<IO::index_type>& point_cell_connectivity,
    const std::vector<IO::index_type>& cell_offset, const std::vector<uint8_t>& cell_types,
    const std::vector<IO::index_type>& face_connectivity,
    const std::vector<IO::index_type>& face_offset)
{
  vtu_writer_.WriteGeometryUnstructuredGrid(point_coordinates, point_cell_connectivity, cell_offset,
      cell_types, face_connectivity, face_offset);
}

/**
 *
 */
void IO::VisualizationWriterVtuPerRank::WritePointDataVectorToDisk(
    const visualization_vector_type_variant& data, unsigned int num_components_per_point,
    const std::string& name)
{
  vtu_writer_.WritePointDataVector(data, num_components_per_point, name);
}

/**
 *
 */
void IO::VisualizationWriterVtuPerRank::WriteCellDataVectorToDisk(
    const visualization_vector_type_variant& data, unsigned int num_components_per_point,
    const std::string& name)
{
  vtu_writer_.WriteCellDataVector(data, num_components_per_point, name);
}

/**
 *
 */
void IO::VisualizationWriterVtuPerRank::FinalizeTimeStep()
{
  vtu_writer_.WriteVtkFooters();

  // Write a collection file summarizing all previously written files
  vtu_writer_.WriteVtkCollectionFileForAllWrittenMasterFiles(
      GLOBAL::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix() + "-" +
      visualization_data_name_);
}
BACI_NAMESPACE_CLOSE
