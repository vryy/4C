/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Write visualization output for a discretization in vtk/vtp format at runtime

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef IO_DISCRETIZATION_RUNTIME_VTP_WRITER_H
#define IO_DISCRETIZATION_RUNTIME_VTP_WRITER_H

/*-----------------------------------------------------------------------------------------------*/
/* headers */

#include <Teuchos_RCP.hpp>

/*-----------------------------------------------------------------------------------------------*/
/* forward declarations */
class Epetra_Vector;
class Epetra_MultiVector;

class RuntimeVtpWriter;

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
  DiscretizationRuntimeVtpWriter();

  /// Destructor
  virtual ~DiscretizationRuntimeVtpWriter() {}

  /** \brief initialize object with all required data
   *
   *  \author grill
   *  \date 04/17 */
  void Initialize(Teuchos::RCP<const DRT::Discretization> discretization,
      const std::string& geometry_name, unsigned int max_number_timesteps_to_be_written,
      double time, bool write_binary_output);

  /** \brief reset current simulation time and time step number
   *
   *  \author grill
   *  \date 04/17 */
  void ResetTimeAndTimeStep(double time, unsigned int timestep);

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

  /** \brief write all required VTP files to filesystem
   *
   *  \author grill
   *  \date 04/17 */
  void WriteFiles();

  /** \brief write a VTP collection file summarizing all written files of this object to filesystem
   *
   *  \author grill
   *  \date 04/17 */
  void WriteCollectionFileOfAllWrittenFiles(const std::string& geometry_name);


  /** \brief determine and set geometry data from nodes based on reference configuration
   *
   *  \author grill
   *  \date 04/17 */
  void SetGeometryFromParticleDiscretization();

 private:
  //! discretization containing nodes of which geometry and result data shall be visualized
  Teuchos::RCP<const DRT::Discretization> discretization_;

  //! the actual vtp writer object that additionally stores the geometry and result data
  Teuchos::RCP<RuntimeVtpWriter> runtime_vtpwriter_;
};

#endif
