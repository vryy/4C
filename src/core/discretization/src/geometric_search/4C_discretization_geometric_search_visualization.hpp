/*-----------------------------------------------------------*/
/*! \file

\brief Visualization writer for geometric search

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_GEOMETRIC_SEARCH_VISUALIZATION_HPP
#define FOUR_C_DISCRETIZATION_GEOMETRIC_SEARCH_VISUALIZATION_HPP

#include "4C_config.hpp"

#include "4C_discretization_geometric_search_bounding_volume.hpp"
#include "4C_io_visualization_manager.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace CORE::GEOMETRICSEARCH
{
  /**
   * @brief Class to write visualization output for geometric search bounding volumes
   */
  class GeometricSearchVisualization : private IO::VisualizationManager
  {
   public:
    /**
     * @brief Constructor (derived from base class), the output data fields are defined here.
     */
    GeometricSearchVisualization(IO::VisualizationParameters parameters, const Epetra_Comm& comm,
        std::string base_output_name);

    /**
     * @brief Write the primitives and predicates of the geometric search to disk.
     *
     * @param visualziation_time (in) Current time that is being visualized
     * @param visualization_step (in) Current time step index
     * @param primitives (in) Vector containing all primitive bounding volumes and their element GID
     * @param predicates (out) Vector containing all predicate bounding volumes and their element
     * GID
     */
    void write_primitives_and_predicates_to_disk(const double visualziation_time,
        const int visualization_step,
        const std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>>& primitives,
        const std::vector<std::pair<int, CORE::GEOMETRICSEARCH::BoundingVolume>>& predicates);

   private:
    //! GID of the current rank
    int my_rank_;
  };
}  // namespace CORE::GEOMETRICSEARCH

FOUR_C_NAMESPACE_CLOSE

#endif
