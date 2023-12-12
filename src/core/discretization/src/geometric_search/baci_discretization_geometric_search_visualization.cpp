/*-----------------------------------------------------------*/
/*! \file

\brief Visualization writer for geometric search

\level 3

*/
/*-----------------------------------------------------------*/

#include "baci_discretization_geometric_search_visualization.H"

#include "baci_discretization_geometric_search_utils.H"
#include "baci_io_visualization_utils.H"

BACI_NAMESPACE_OPEN

namespace CORE::GEOMETRICSEARCH
{
  GeometricSearchVisualization::GeometricSearchVisualization(
      IO::VisualizationParameters parameters, const Epetra_Comm& comm, std::string base_output_name)
      : IO::VisualizationManager(parameters, comm, base_output_name)
  {
    my_rank_ = comm.MyPID();
    this->GetVisualizationData().RegisterCellData<int>("element_id", 1);
    this->GetVisualizationData().RegisterCellData<int>("element_created_on_rank", 1);
    this->GetVisualizationData().RegisterCellData<int>("primitive_predicate_flag", 1);
  }

  void GeometricSearchVisualization::WritePrimitivesAndPredicatesToDisk(
      const double visualziation_time, const int visualization_step,
      const std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>>& primitives,
      const std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>>& predicates)
  {
#ifndef BACI_WITH_ARBORX
    dserror(
        "CORE::GEOMETRICSEARCH::WritePrimitivesAndPredicatesToDisk can only be used with ArborX."
        "To use it, enable ArborX during the configure process.");
#else
    ClearData();

    auto& visualization_data = GetVisualizationData();

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

}  // namespace CORE::GEOMETRICSEARCH

BACI_NAMESPACE_CLOSE
