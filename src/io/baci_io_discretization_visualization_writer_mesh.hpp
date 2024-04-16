/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Write visualization output for a discretization, i.e., write the mesh and results on the mesh
to disk

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef FOUR_C_IO_DISCRETIZATION_VISUALIZATION_WRITER_MESH_HPP
#define FOUR_C_IO_DISCRETIZATION_VISUALIZATION_WRITER_MESH_HPP

#include "baci_config.hpp"

#include "baci_io_control.hpp"

#include <Teuchos_RCP.hpp>

// forward declarations
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_Map;

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}
namespace IO
{
  class VisualizationManager;
  struct VisualizationParameters;
}  // namespace IO

namespace IO
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
     */
    DiscretizationVisualizationWriterMesh(
        const Teuchos::RCP<const DRT::Discretization>& discretization,
        VisualizationParameters parameters);

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
    void AppendDofBasedResultDataVector(const Teuchos::RCP<Epetra_Vector>& result_data_dofbased,
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
    void AppendNodeBasedResultDataVector(
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
    void AppendElementBasedResultDataVector(
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
     * \brief Write the BACI internal element GIDs for each element
     *
     * @param resultname (in) Name of the field in the visualization file
     */
    void AppendElementGID(const std::string& resultname);

    /**
     * \brief Write ghosting information to the elements
     *
     * For more details look at the documentation to IO::AppendElementGhostingInformation
     */
    void AppendElementGhostingInformation();

    /**
     * \brief Write the BACI internal node GIDs for each node
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
    void SetGeometryFromDiscretization();


   private:
    //! Discretization containing elements of which geometry and result data shall be visualized
    Teuchos::RCP<const DRT::Discretization> discretization_;

    //! The actual visualization writer object that additionally stores the geometry and result data
    Teuchos::RCP<VisualizationManager> visualization_manager_;

    //! Node row and col maps the geometry of visualization writer is based on
    Teuchos::RCP<Epetra_Map> noderowmap_last_geometry_set_;
    Teuchos::RCP<Epetra_Map> nodecolmap_last_geometry_set_;
  };

  /**
   * \brief Add a vector with the length of num_proc to each element, which contains a 1 for the
   * ranks that ghost the element.
   * @param discretization (in) Discretization
   * @param visualization_manager (in/out) Visualization writer object
   * @param is_beam (in) If beam or non-beam elements should be output
   */
  void AppendElementGhostingInformation(const DRT::Discretization& discretization,
      VisualizationManager& visualization_manager, bool is_beam = false);

}  // namespace IO
BACI_NAMESPACE_CLOSE

#endif
