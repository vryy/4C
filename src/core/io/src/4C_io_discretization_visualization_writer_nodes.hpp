// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_IO_DISCRETIZATION_VISUALIZATION_WRITER_NODES_HPP
#define FOUR_C_IO_DISCRETIZATION_VISUALIZATION_WRITER_NODES_HPP


#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"

#include <memory>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE
namespace Core::IO
{
  class VisualizationManager;

  struct VisualizationParameters;
}  // namespace Core::IO

namespace Core::IO
{
  /**
   * @brief This object allows to write nodal visualization output for a discretization
   */
  class DiscretizationVisualizationWriterNodes
  {
   public:
    /**
     * @brief Constructor
     *
     * @param discretization (in) Pointer to the discretization
     * @param parameters (in)     Visualization parameters
     */
    DiscretizationVisualizationWriterNodes(
        const std::shared_ptr<const Core::FE::Discretization> &discretization,
        VisualizationParameters parameters);

    /**
     * @brief Append result vector with num_dof values per node to output data
     *
     * The results are expected to be stored in the following form (for result_num_dofs_per_node =
     * x, y, z):
     * [node1_x, node1_y, node1_z, node2_x, node2_y, node2_z, node3_x, node3_y, node3_z, ...]
     *
     * @param result_data_dofbased (in) Column based result data vector
     * @param result_num_dofs_per_node  (in) Number of scalar values per node
     * @param read_result_data_from_dofindex (in) Starting DOF index for the nodal DOFs. This is
     * used if not all nodal DOFs should be output, e.g., velocity or pressure in fluid.
     * @param resultname (in) Name of the field to be written to the visualization file
     */
    void append_dof_based_result_data_vector(Core::LinAlg::Vector<double> &result_data_dofbased,
        unsigned int result_num_dofs_per_node, const std::string &resultname);

    /**
     * @brief Append result vector with num_components values per node to output data
     *
     * The results are expected to be stored in the following form (for result_num_dofs_per_node =
     * x, y, z):
     * [
     *   [node1_x, node1_y, node1_z]
     *   [node2_x, node2_y, node2_z],
     *   [node3_x, node3_y, node3_z],
     *   ...
     * ]
     *
     * @param result_data_nodebased (in) Column node based result data vector
     * @param result_num_components_per_node (in) Number of scalar values per node
     * @param resultname (in) Name of the field to be written to the visualization file
     */
    void append_node_based_result_data_vector(
        Core::LinAlg::MultiVector<double> &result_data_nodebased,
        unsigned int result_num_components_per_node, const std::string &resultname);

    /**
     * \brief Write the visualization files to disk
     */
    void write_to_disk(const double visualization_time, const int visualization_step);

    /**
     * @brief Set geometry data from nodes based on reference configuration
     */
    void set_geometry_from_discretization();

   private:
    //! discretization containing nodes of which geometry and result data shall be visualized
    std::shared_ptr<const Core::FE::Discretization> discretization_;

    //! The actual visualization writer object that additionally stores the geometry and result data
    std::shared_ptr<VisualizationManager> visualization_manager_;
  };
}  // namespace Core::IO
FOUR_C_NAMESPACE_CLOSE

#endif
