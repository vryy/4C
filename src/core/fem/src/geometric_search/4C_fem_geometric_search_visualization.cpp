/*-----------------------------------------------------------*/
/*! \file

\brief Visualization writer for geometric search

\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_fem_geometric_search_visualization.hpp"

#include "4C_fem_geometric_search_utils.hpp"
#include "4C_io_visualization_utils.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN

namespace Core::GeometricSearch
{
  GeometricSearchVisualization::GeometricSearchVisualization(
      Core::IO::VisualizationParameters parameters, const Epetra_Comm& comm,
      std::string base_output_name)
      : Core::IO::VisualizationManager(std::move(parameters), comm, std::move(base_output_name))
  {
    my_rank_ = comm.MyPID();
    this->get_visualization_data().RegisterCellData<int>("element_id", 1);
    this->get_visualization_data().RegisterCellData<int>("element_created_on_rank", 1);
    this->get_visualization_data().RegisterCellData<int>("primitive_predicate_flag", 1);
  }

  void GeometricSearchVisualization::write_primitives_and_predicates_to_disk(
      const double visualziation_time, const int visualization_step,
      const std::vector<std::pair<int, Core::GeometricSearch::BoundingVolume>>& primitives,
      const std::vector<std::pair<int, Core::GeometricSearch::BoundingVolume>>& predicates)
  {
#ifndef FOUR_C_WITH_ARBORX
    FOUR_C_THROW(
        "Core::GeometricSearch::write_primitives_and_predicates_to_disk can only be used with "
        "ArborX."
        "To use it, enable ArborX during the configure process.");
#else
    clear_data();

    auto& visualization_data = get_visualization_data();

    // Function that adds bounding volumes to the visualization data
    const auto add_bounding_volumes =
        [&](const auto& bounding_volumes, const int primitive_predicate_flag)
    {
      for (const auto& bounding_volume : bounding_volumes)
      {
        const auto& [points, polygons] = GetKDopPolyhedronRepresentation(bounding_volume.second);
        IO::AppendPolyhedronToVisualizationData(visualization_data, points, polygons);
        visualization_data.GetCellData<int>("element_id").push_back(bounding_volume.first);
        visualization_data.GetCellData<int>("element_created_on_rank").push_back(my_rank_);
        visualization_data.GetCellData<int>("primitive_predicate_flag")
            .push_back(primitive_predicate_flag);
      }
    };

    // All bounding volumes to the output data
    add_bounding_volumes(primitives, 0);
    add_bounding_volumes(predicates, 1);

    // Write the data to disk
    WriteToDisk(visualziation_time, visualization_step);
#endif
  }

}  // namespace Core::GeometricSearch

FOUR_C_NAMESPACE_CLOSE
