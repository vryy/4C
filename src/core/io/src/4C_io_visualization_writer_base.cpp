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
Core::IO::VisualizationWriterBase::VisualizationWriterBase(
    Core::IO::VisualizationParameters parameters, const Epetra_Comm& comm,
    std::string visualization_data_name)
    : parameters_(std::move(parameters)),
      comm_(comm),
      visualization_data_name_(std::move(visualization_data_name))
{
}

/**
 *
 */
void Core::IO::VisualizationWriterBase::write_visualization_data_to_disk(
    const VisualizationData& visualization_data, const double visualziation_time,
    const int visualization_step)
{
  visualization_data.consistency_check();
  initialize_time_step(visualziation_time, visualization_step);
  write_field_data_to_disk(visualization_data.get_field_data_map());
  write_geometry_to_disk(visualization_data.get_point_coordinates(),
      visualization_data.get_cell_connectivity(), visualization_data.get_cell_offsets(),
      visualization_data.get_cell_types(), visualization_data.get_face_connectivity(),
      visualization_data.get_face_offsets());
  for (const auto& data_name : visualization_data.get_point_data_names())
  {
    write_point_data_vector_to_disk(visualization_data.get_point_data_variant(data_name),
        visualization_data.get_point_data_dimension(data_name), data_name);
  }
  for (const auto& data_name : visualization_data.get_cell_data_names())
  {
    write_cell_data_vector_to_disk(visualization_data.get_cell_data_variant(data_name),
        visualization_data.get_cell_data_dimension(data_name), data_name);
  }
  finalize_time_step();
}
FOUR_C_NAMESPACE_CLOSE
