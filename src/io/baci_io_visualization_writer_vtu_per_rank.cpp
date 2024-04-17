/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Define a container that holds the visualization data and the visualization writer

\level 0

*/
/*-----------------------------------------------------------------------------------------------*/

#include "baci_io_visualization_writer_vtu_per_rank.hpp"

#include "baci_inpar_IO_runtime_output.hpp"
#include "baci_io_control.hpp"
#include "baci_io_visualization_data.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
IO::VisualizationWriterVtuPerRank::VisualizationWriterVtuPerRank(
    const IO::VisualizationParameters& parameters, const Epetra_Comm& comm,
    std::string visualization_data_name)
    : VisualizationWriterBase(parameters, comm, std::move(visualization_data_name)),
      vtu_writer_(comm.MyPID(), comm.NumProc(),
          std::pow(10, IO::GetTotalDigitsToReserveInTimeStep(parameters)),
          parameters.directory_name_, (parameters.file_name_prefix_ + "-vtk-files"),
          visualization_data_name_, parameters.restart_from_name_, parameters.restart_time_,
          parameters.data_format_ == INPAR::IO_RUNTIME_OUTPUT::OutputDataFormat::binary)
{
}

/**
 *
 */
void IO::VisualizationWriterVtuPerRank::InitializeTimeStep(
    const double visualziation_time, const int visualization_step)
{
  vtu_writer_.ResetTimeAndTimeStep(visualziation_time, visualization_step);

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
      parameters_.file_name_prefix_ + "-" + visualization_data_name_);
}
FOUR_C_NAMESPACE_CLOSE
