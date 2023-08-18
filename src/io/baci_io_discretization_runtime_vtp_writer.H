/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Write visualization output for a discretization in vtk/vtp format at runtime

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef BACI_IO_DISCRETIZATION_RUNTIME_VTP_WRITER_H
#define BACI_IO_DISCRETIZATION_RUNTIME_VTP_WRITER_H

/*-----------------------------------------------------------------------------------------------*/
/* headers */

#include "baci_io_visualization_manager.H"

#include <Teuchos_RCP.hpp>

/*-----------------------------------------------------------------------------------------------*/
/* forward declarations */
class Epetra_Vector;
class Epetra_MultiVector;


namespace DRT
{
  class Discretization;
}

/*-----------------------------------------------------------------------------------------------*/
/* namespace */

/*!
 * \brief This object allows to write visualization output for a discretization
 *        - in vtk/vtp format
 *        - at runtime
 *        - in parallel
 *        - binary-encoded
 *
 * \author grill
 * \date 04/17
 */
class DiscretizationRuntimeVtpWriter
{
 public:
  /// Constructor
  DiscretizationRuntimeVtpWriter(const Teuchos::RCP<const DRT::Discretization>& discretization,
      IO::VisualizationParameters parameters);

  /// Destructor
  virtual ~DiscretizationRuntimeVtpWriter() = default;

  /** \brief append result vector with num_dof values per node to output data
   *
   *  \author grill
   *  \date 04/17 */
  void AppendDofBasedResultDataVector(const Teuchos::RCP<Epetra_Vector>& result_data_dofbased,
      unsigned int result_num_dofs_per_node, const std::string& resultname);

  /** \brief append result vector with num_components values per node to output data
   *
   *  \author grill
   *  \date 04/17 */
  void AppendNodeBasedResultDataVector(
      const Teuchos::RCP<Epetra_MultiVector>& result_data_nodebased,
      unsigned int result_num_components_per_node, const std::string& resultname);

  /**
   * \brief Write the visualization files to disk
   */
  void WriteToDisk(const double visualziation_time, const int visualization_step);

  /** \brief determine and set geometry data from nodes based on reference configuration
   *
   *  \author grill
   *  \date 04/17 */
  void SetGeometryFromParticleDiscretization();

 private:
  //! discretization containing nodes of which geometry and result data shall be visualized
  Teuchos::RCP<const DRT::Discretization> discretization_;

  //! the actual vtp writer object that additionally stores the geometry and result data
  Teuchos::RCP<IO::VisualizationManager> visualization_manager_;
};

#endif
