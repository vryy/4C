/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Write visualization output for a discretization in vtk/vtu format at runtime

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef IO_DISCRETIZATION_RUNTIME_VTU_WRITER_H
#define IO_DISCRETIZATION_RUNTIME_VTU_WRITER_H

/*-----------------------------------------------------------------------------------------------*/
/* headers */

#include <Teuchos_RCP.hpp>

/*-----------------------------------------------------------------------------------------------*/
/* forward declarations */
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_Map;

class RuntimeVtuWriter;

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
  DiscretizationRuntimeVtuWriter();

  /// Destructor
  virtual ~DiscretizationRuntimeVtuWriter() {}

  /** \brief initialize object with all required data
   *
   *  \author grill
   *  \date 03/17 */
  void Initialize(Teuchos::RCP<const DRT::Discretization> discretization,
      unsigned int max_number_timesteps_to_be_written, double time, bool write_binary_output);

  /** \brief reset current simulation time and time step number
   *
   *  \author grill
   *  \date 03/17 */
  void ResetTimeAndTimeStep(double time, unsigned int timestep);

  /** \brief append result vector with num_dof values per node to output data
   *
   *  \author grill
   *  \date 03/17 */
  void AppendDofBasedResultDataVector(const Teuchos::RCP<Epetra_Vector>& result_data_dofbased,
      unsigned int result_num_dofs_per_node, unsigned int read_result_data_from_dofindex,
      const std::string& resultname);

  /** \brief append result vector with num_components values per node to output data
   *
   *  \author grill
   *  \date 03/17 */
  void AppendNodeBasedResultDataVector(
      const Teuchos::RCP<Epetra_MultiVector>& result_data_nodebased,
      unsigned int result_num_components_per_node, const std::string& resultname);

  /** \brief append result vector with num_components values per element to output data
   *
   *  \author grill
   *  \date 03/17 */
  void AppendElementBasedResultDataVector(
      const Teuchos::RCP<Epetra_MultiVector>& result_data_elementbased,
      unsigned int result_num_components_per_element, const std::string& resultname);

  /** \brief write all required VTU files to filesystem
   *
   *  \author grill
   *  \date 03/17 */
  void WriteFiles();

  /** \brief write a VTK collection file summarizing all written files of this object to filesystem
   *
   *  \author grill
   *  \date 03/17 */
  void WriteCollectionFileOfAllWrittenFiles();


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
  Teuchos::RCP<RuntimeVtuWriter> runtime_vtuwriter_;

  //! node row and col maps the geometry of runtime vtu writer is based on
  Teuchos::RCP<Epetra_Map> noderowmap_last_geometry_set_;
  Teuchos::RCP<Epetra_Map> nodecolmap_last_geometry_set_;
};

#endif
