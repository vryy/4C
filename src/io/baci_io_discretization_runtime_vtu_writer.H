/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Write visualization output for a discretization in vtk/vtu format at runtime

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef BACI_IO_DISCRETIZATION_RUNTIME_VTU_WRITER_H
#define BACI_IO_DISCRETIZATION_RUNTIME_VTU_WRITER_H

/*-----------------------------------------------------------------------------------------------*/
/* headers */

#include "baci_io_visualization_manager.H"

#include <Teuchos_RCP.hpp>

/*-----------------------------------------------------------------------------------------------*/
/* forward declarations */
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_Map;


namespace DRT
{
  class Discretization;
}

/*-----------------------------------------------------------------------------------------------*/
/* namespace */

/*!
 * \brief This object allows to write visualization output for a discretization
 *        - in vtk/vtu format
 *        - at runtime
 *        - in parallel
 *        - binary-encoded
 *
 * \author grill
 * \date 03/17
 */
class DiscretizationRuntimeVtuWriter
{
 public:
  /// Constructor
  DiscretizationRuntimeVtuWriter(const Teuchos::RCP<const DRT::Discretization>& discretization,
      IO::VisualizationParameters parameters);

  /// Destructor
  virtual ~DiscretizationRuntimeVtuWriter() = default;
  /**
   * @brief reset state depending if the maps changed or not
   */
  void Reset();

  /** \brief append result vector with num_dof values per node to output data
   *
   *  \note This method only accepts results based on col maps
   *
   *  \author grill
   *  \date 03/17 */
  void AppendDofBasedResultDataVector(const Teuchos::RCP<Epetra_Vector>& result_data_dofbased,
      unsigned int result_num_dofs_per_node, unsigned int read_result_data_from_dofindex,
      const std::string& resultname);

  /** \brief append result vector with num_components values per node to output data
   *
   *  \note This method only accepts results based on col maps
   *
   *  \author grill
   *  \date 03/17 */
  void AppendNodeBasedResultDataVector(
      const Teuchos::RCP<Epetra_MultiVector>& result_data_nodebased,
      unsigned int result_num_components_per_node, const std::string& resultname);

  /** \brief append result vector with num_components values per element to output data
   *
   *  \note This method only accepts results based on row maps
   *
   *  \author grill
   *  \date 03/17 */
  void AppendElementBasedResultDataVector(
      const Teuchos::RCP<Epetra_MultiVector>& result_data_elementbased,
      unsigned int result_num_components_per_element, const std::string& resultname);

  /**
   * \brief Write owner id for each element.
   *
   * In order to process all solid elements (beam elements are handled by a sepperate output
   * writer), it is sufficient that each processor writes all elements in his row map. This can
   * easily be achieved by appending the PID of the writing processor.
   *
   * @param resultname (in) Name of the owner field in the VTK file.
   */
  void AppendElementOwner(const std::string resultname);

  /**
   * \brief Write the BACI internal element GIDs for each element.
   *
   * @param resultname (in) Name of the field in the VTK file.
   */
  void AppendElementGID(const std::string& resultname);

  /**
   * \brief Write ghosting information to the elements.
   */
  void AppendElementGhostingInformation();

  /**
   * \brief Write the BACI internal node GIDs for each element.
   *
   * @param resultname (in) Name of the field in the VTK file.
   */
  void AppendNodeGID(const std::string& resultname);

  /**
   * \brief Write the visualization files to disk
   */
  void WriteToDisk(const double visualziation_time, const int visualization_step);

 private:
  /** \brief determine and set geometry data from elements based on reference configuration
   *
   *  \author grill
   *  \date 03/17 */
  void SetGeometryFromDiscretizationStandard();


 private:
  //! discretization containing elements of which geometry and result data shall be visualized
  Teuchos::RCP<const DRT::Discretization> discretization_;

  //! the actual vtu writer object that additionally stores the geometry and result data
  Teuchos::RCP<IO::VisualizationManager> visualization_manager_;

  //! node row and col maps the geometry of runtime vtu writer is based on
  Teuchos::RCP<Epetra_Map> noderowmap_last_geometry_set_;
  Teuchos::RCP<Epetra_Map> nodecolmap_last_geometry_set_;
};

namespace IO
{
  /**
   * \brief Add a vector with the length of num_proc to each element, which contains a 1 for the
   * ranks that ghost the element.
   * @param discretization (in) Discretization
   * @param runtime_vtuwriter (out) VTU writer
   * @param is_beam (in) If beam or solid elements should be output.
   */
  void AppendElementGhostingInformation(
      const Teuchos::RCP<const DRT::Discretization>& discretization,
      const Teuchos::RCP<IO::VisualizationManager>& visualization_manager, bool is_beam = false);
}  // namespace IO

#endif
