/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Write nodal coordinates and values at the nodes to a visualization file

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef FOUR_C_IO_DISCRETIZATION_VISUALIZATION_WRITER_NODES_HPP
#define FOUR_C_IO_DISCRETIZATION_VISUALIZATION_WRITER_NODES_HPP


#include "4C_config.hpp"

#include <Epetra_MultiVector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}
namespace CORE::IO
{
  class VisualizationManager;
  struct VisualizationParameters;
}  // namespace CORE::IO

namespace CORE::IO
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
        const Teuchos::RCP<const DRT::Discretization>& discretization,
        VisualizationParameters parameters);

    /**
     * @brief Destructor
     */
    virtual ~DiscretizationVisualizationWriterNodes() = default;

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
    void append_dof_based_result_data_vector(
        const Teuchos::RCP<Epetra_Vector>& result_data_dofbased,
        unsigned int result_num_dofs_per_node, const std::string& resultname);

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
        const Teuchos::RCP<Epetra_MultiVector>& result_data_nodebased,
        unsigned int result_num_components_per_node, const std::string& resultname);

    /**
     * \brief Write the visualization files to disk
     */
    void WriteToDisk(const double visualization_time, const int visualization_step);

    /**
     * @brief Set geometry data from nodes based on reference configuration
     */
    void set_geometry_from_discretization();

   private:
    //! discretization containing nodes of which geometry and result data shall be visualized
    Teuchos::RCP<const DRT::Discretization> discretization_;

    //! The actual visualization writer object that additionally stores the geometry and result data
    Teuchos::RCP<VisualizationManager> visualization_manager_;
  };
}  // namespace CORE::IO
FOUR_C_NAMESPACE_CLOSE

#endif
