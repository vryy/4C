/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Write nodal coordinates and values at the nodes to a visualization file

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/


#ifndef BACI_IO_DISCRETIZATION_VISUALIZATION_WRITER_NODES_H
#define BACI_IO_DISCRETIZATION_VISUALIZATION_WRITER_NODES_H


#include <Teuchos_RCP.hpp>


// forward declarations
class Epetra_Vector;
class Epetra_MultiVector;
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
  /**
   * @brief This object allows to write nodal visualization output for a discretization
   */
  class DiscretizationVisualizationWriterNodes
  {
   public:
    /// Constructor
    DiscretizationVisualizationWriterNodes(
        const Teuchos::RCP<const DRT::Discretization>& discretization,
        VisualizationParameters parameters);

    /// Destructor
    virtual ~DiscretizationVisualizationWriterNodes() = default;

    /** \brief Append result vector (flat vector) with num_dof values per node to output data
     */
    void AppendDofBasedResultDataVector(const Teuchos::RCP<Epetra_Vector>& result_data_dofbased,
        unsigned int result_num_dofs_per_node, const std::string& resultname);

    /** \brief Append result multi vector with num_components columns per node to output data
     */
    void AppendNodeBasedResultDataVector(
        const Teuchos::RCP<Epetra_MultiVector>& result_data_nodebased,
        unsigned int result_num_components_per_node, const std::string& resultname);

    /**
     * \brief Write the visualization files to disk
     */
    void WriteToDisk(const double visualziation_time, const int visualization_step);

    /** \brief Set geometry data from nodes based on reference configuration
     */
    void SetGeometryFromDiscretization();

   private:
    //! discretization containing nodes of which geometry and result data shall be visualized
    Teuchos::RCP<const DRT::Discretization> discretization_;

    //! the actual vtp writer object that additionally stores the geometry and result data
    Teuchos::RCP<VisualizationManager> visualization_manager_;
  };
}  // namespace IO
#endif
