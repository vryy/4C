/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Write visualization output for a discretization, i.e., write the mesh and results on the mesh
to disk

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef FOUR_C_IO_DISCRETIZATION_VISUALIZATION_WRITER_MESH_HPP
#define FOUR_C_IO_DISCRETIZATION_VISUALIZATION_WRITER_MESH_HPP

#include "4C_config.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE
namespace Core::Elements
{
  class Element;
}

namespace Core::IO
{
  class VisualizationManager;
  struct VisualizationParameters;
}  // namespace Core::IO

namespace Core::IO
{
  /*!
   * \brief This object allows to write visualization output for a discretization, i.e., write the
   * mesh and results on the mesh to disk
   */
  class DiscretizationVisualizationWriterMesh
  {
   public:
    /**
     * @brief Constructor
     *
     * @param discretization (in)  Pointer to the discretization
     * @param parameters (in)      Visualization parameters
     * @param element_filter (in)  An optional function that returns true for all elements that
     *   should be included in the visualization. By default, all elements in the discretization are
     *   included.
     */
    DiscretizationVisualizationWriterMesh(
        const Teuchos::RCP<const Core::FE::Discretization>& discretization,
        VisualizationParameters parameters,
        std::function<bool(const Core::Elements::Element* element)> element_filter =
            [](const Core::Elements::Element*) { return true; });

    /**
     * @brief Destructor
     */
    virtual ~DiscretizationVisualizationWriterMesh() = default;

    /**
     * @brief Reset state depending if the maps changed or not
     */
    void Reset();

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
        unsigned int result_num_dofs_per_node, unsigned int read_result_data_from_dofindex,
        const std::string& resultname);

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
     * @brief Append result vector with num_components values per element to output data
     *
     * The results are expected to be stored in the following form (for
     * result_num_components_per_element = x, y, z):
     * [ele1_x, ele1_y, ele1_z, ele2_x, ele2_y, ele2_z, ele3_x, ele3_y, ele3_z, ...]
     *
     * @param result_data_elementbased (in) Column based result data vector
     * @param result_num_components_per_element (in) Number of scalar values per element
     * @param resultname (in) Name of the field to be written to the visualization file
     */
    void append_element_based_result_data_vector(
        const Teuchos::RCP<Epetra_MultiVector>& result_data_elementbased,
        unsigned int result_num_components_per_element, const std::string& resultname);

    /**
     * \brief Write owner id for each element
     *
     * In order to process all solid elements (beam elements are handled by a separate output
     * writer), it is sufficient that each processor writes all elements in his row map. This can
     * easily be achieved by appending the PID of the writing processor.
     *
     * @param resultname (in) Name of the owner field in the visualization file
     */
    void AppendElementOwner(const std::string resultname);

    /**
     * \brief Write the 4C internal element GIDs for each element
     *
     * @param resultname (in) Name of the field in the visualization file
     */
    void AppendElementGID(const std::string& resultname);

    /**
     * \brief Write ghosting information to the elements
     *
     * For more details look at the documentation to Core::IO::append_element_ghosting_information
     */
    void append_element_ghosting_information();

    /**
     * \brief Write the 4C internal node GIDs for each node
     *
     * @param resultname (in) Name of the field in the visualization file
     */
    void AppendNodeGID(const std::string& resultname);

    /**
     * \brief Write the visualization files to disk
     */
    void WriteToDisk(const double visualization_time, const int visualization_step);

   private:
    /** \brief Determine and set geometry data from elements based on reference configuration
     *
     * To simplify parallel ouput, we loop over each row element and let the element write its
     * topology to the the output file. This results in a "discontinuous" visualization, i.e., the
     * nodes of adjacent elements are not connected in the output file. In ParaView the global
     * nodal connectivity can be restored with the CleanToGrid filter.
     */
    void set_geometry_from_discretization();


   private:
    //! discretization containing elements of which geometry and result data shall be visualized
    Teuchos::RCP<const Core::FE::Discretization> discretization_;

    //! The actual visualization writer object that additionally stores the geometry and result data
    Teuchos::RCP<VisualizationManager> visualization_manager_;

    //! A filter function that returns true for all elements that should be visualized.
    std::function<bool(const Core::Elements::Element* element)> element_filter_;

    //! Node row and col maps the geometry of visualization writer is based on
    Teuchos::RCP<Epetra_Map> noderowmap_last_geometry_set_;
    Teuchos::RCP<Epetra_Map> nodecolmap_last_geometry_set_;
  };

  /**
   * \brief Add a vector with the length of num_proc to each element, which contains a 1 for the
   * ranks that ghost the element.
   * @param discretization (in) discretization
   * @param visualization_manager (in/out) Visualization writer object
   * @param element_predicate (in) A predicate function which returns whether a given element
   * should be included in the output.
   */
  void append_element_ghosting_information(const Core::FE::Discretization& discretization,
      VisualizationManager& visualization_manager,
      const std::function<bool(const Core::Elements::Element* ele)>& element_predicate);

}  // namespace Core::IO
FOUR_C_NAMESPACE_CLOSE

#endif
