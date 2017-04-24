/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_discretization_runtime_vtu_writer.H

\brief Write visualization output for a beam discretization in vtk/vtu format at runtime

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef BEAM_DISCRETIZATION_RUNTIME_VTU_WRITER_H
#define BEAM_DISCRETIZATION_RUNTIME_VTU_WRITER_H

/*-----------------------------------------------------------------------------------------------*/
/* headers */

#include <Teuchos_RCP.hpp>

/*-----------------------------------------------------------------------------------------------*/
/* forward declarations */
namespace GEO{
namespace MESHFREE{
  class BoundingBox;
}
}
class Epetra_Vector;
class Epetra_MultiVector;

class RuntimeVtuWriter;

namespace DRT
{
  class Discretization;

namespace ELEMENTS
{
    class Beam3Base;
}
}

/*-----------------------------------------------------------------------------------------------*/
/* namespace */

/*!
 * \brief This object allows to write visualization output for a beam discretization
 *        - in vtk/vtu format (i.e. as an unstructured grid)
 *        - at runtime
 *        - in parallel
 *        - binary-encoded
 *
 *  Note: The special thing about beams is that they use 'non-standard' interpolation schemes,
 *        e.g. cubic Hermite polynomials for the interpolation of the centerline geometry.
 *        Thus, the geometry of one element cannot be represented by one simple vtk cell type.
 *        Moreover, the results that we want to visualize are special, e.g. the triad field,
 *        cross-section resultants (axial force, shear forces, bending moments, torque), ...
 *
 * \author grill
 * \date 03/17
 */
class BeamDiscretizationRuntimeVtuWriter
{

public:
  /// Constructor
  BeamDiscretizationRuntimeVtuWriter();

  /// Destructor
  virtual ~BeamDiscretizationRuntimeVtuWriter()
  {}

  /** \brief initialize object with all required data
   *
   *  \author grill
   *  \date 03/17 */
  void Initialize(
      Teuchos::RCP<DRT::Discretization> discretization,
      Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector,
      bool use_absolute_positions_for_point_coordinates,
      Teuchos::RCP<const GEO::MESHFREE::BoundingBox> const& periodic_boundingbox,
      unsigned int max_number_timesteps_to_be_written,
      double time,
      bool write_binary_output );


  /** \brief reset current simulation time and time step number
   *
   *  \author grill
   *  \date 03/17 */
  void ResetTimeAndTimeStep(
      double time,
      unsigned int timestep );


  /** \brief append triad field determined from given displacement state to output data
   *
   *  \author grill
   *  \date 03/17 */
  void AppendTriadField(
      Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector);

    /** \brief append discplacement state
   *
   *  \author grill
   *  \date 03/17 */
  void AppendDisplacementField(
      Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector);

  /** \brief append tangent vector field determined from given displacement state to output data
   *
   *  \author grill
   *  \date 03/17 */
//  void AppendTangentVectorField(
//      Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector);


  /** \brief append information about element owning processor to output data
   *
   *  \author grill
   *  \date 03/17 */
  void AppendElementOwningProcessor();

  /** \brief append circular cross-section radius of elements to output data
   *
   *  \author grill
   *  \date 03/17 */
  void AppendElementCircularCrossSectionRadius();

  /** \brief append a vector field defining orientation and radius of a circular cross-section to output data
   *
   *  \author grill
   *  \date 03/17 */
  void AppendPointCircularCrossSectionInformationVector(
      Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector);

  /** \brief append material cross-section strain resultant values at Gauss points to output data
   *
   *  \author grill
   *  \date 03/17 */
  void AppendGaussPointMaterialCrossSectionStrainResultants();

  /** \brief append material cross-section stress resultant values at Gauss points to output data
   *
   *  \author grill
   *  \date 03/17 */
  void AppendGaussPointMaterialCrossSectionStressResultants();

  /** \brief append current internal energy of the elements to output data
   *
   *  \author grill
   *  \date 03/17 */
  void AppendElementElasticEnergy();


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

  /** \brief determine and set geometry data from beam elements based on given displacement state
   *
   *  \author grill
   *  \date 03/17 */
  void SetGeometryFromBeamDiscretization(
      Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector);

private:

  /** \brief insert all values of a given std::vector at the end of another given std::vector
   *
   *  \author grill
   *  \date 03/17 */
  // Todo template <typename T>
  void InsertVectorValuesAtBackOfOtherVector(
      const std::vector<double> & vector_input,
      std::vector<double> & vector_output);


private:

  //! discretization containing beam elements of which geometry and result data shall be visualized
  Teuchos::RCP<const DRT::Discretization> discretization_;

  //! all local row indices of beam elements in the given discretization
  std::vector<unsigned int> local_row_indices_beam_elements_;

  //! periodic bounding box object
  Teuchos::RCP<const GEO::MESHFREE::BoundingBox> periodic_boundingbox_;

  //! number of points for each element ( in case of periodic boundary conditions
  //! not equal to BEAMSVTUVISUALSUBSEGMENTS + 1)
  std::vector<int> num_cells_per_element_;

  //! the actual vtu writer object that additionally stores the geometry and result data
  Teuchos::RCP<RuntimeVtuWriter> runtime_vtuwriter_;

  //! flag indicating whether to use absolute positions for point coordinates
  // false: use reference position
  bool use_absolute_positions_;

};

#endif
