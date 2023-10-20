/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Write visualization output for a beam discretization in vtk/vtu format at runtime

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef BACI_BEAM3_DISCRETIZATION_RUNTIME_VTU_WRITER_H
#define BACI_BEAM3_DISCRETIZATION_RUNTIME_VTU_WRITER_H

/*-----------------------------------------------------------------------------------------------*/
/* headers */

#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_io_visualization_manager.H"

#include <Teuchos_RCP.hpp>

/*-----------------------------------------------------------------------------------------------*/
/* forward declarations */
namespace CORE::GEO
{
  namespace MESHFREE
  {
    class BoundingBox;
  }
}  // namespace CORE::GEO
class Epetra_Comm;
class Epetra_Vector;
class Epetra_MultiVector;

namespace DRT
{
  class Discretization;

  namespace ELEMENTS
  {
    class Beam3Base;
  }
}  // namespace DRT

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
  BeamDiscretizationRuntimeVtuWriter(
      IO::VisualizationParameters parameters, const Epetra_Comm& comm);

  /// Destructor
  virtual ~BeamDiscretizationRuntimeVtuWriter() = default;
  /** \brief initialize object with all required data
   *
   *  \author grill
   *  \date 03/17 */
  void Initialize(Teuchos::RCP<DRT::Discretization> discretization,
      bool use_absolute_positions_for_point_coordinates, const unsigned int n_subsegments,
      Teuchos::RCP<const CORE::GEO::MESHFREE::BoundingBox> const& periodic_boundingbox);

  /** \brief append triad field determined from given displacement state to output data
   *
   *  \author grill
   *  \date 03/17 */
  void AppendTriadField(Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector);

  /** \brief append discplacement state
   *
   *  \author grill
   *  \date 03/17 */
  void AppendDisplacementField(Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector);

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

  /**
   * \brief Append the baci interal GIDs to all beam elements.
   */
  void AppendElementGID();

  /**
   * \brief Append the element ghosting information.
   */
  void AppendElementGhostingInformation();

  /** \brief append internal (elastic) energy of element
   *
   *  \author eichinger
   *  \date 01/18 */
  void AppendElementInternalEnergy();

  /** \brief append kinetic energy of element
   *
   *  \author eichinger
   *  \date 01/18 */
  void AppendElementKineticEnergy();

  /** \brief append information about to which filament an element belonging
   *
   *  \author eichinger
   *  \date 05/17 */
  void AppendElementFilamentIdAndType();

  /** \brief append circular cross-section radius of elements to output data
   *
   *  \author grill
   *  \date 03/17 */
  void AppendElementCircularCrossSectionRadius();

  /** \brief append a vector field defining orientation and radius of a circular cross-section to
   * output data
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

  /**
   * \brief Append interpolated GP values of the material cross-section strain resultants.
   */
  void AppendGaussPointMaterialCrossSectionStrainResultantsContinuous();

  /** \brief append material cross-section stress resultant values at Gauss points to output data
   *
   *  \author grill
   *  \date 03/17 */
  void AppendGaussPointMaterialCrossSectionStressResultants();

  /**
   * \brief Append interpolated GP values of the material cross-section stress resultants.
   */
  void AppendGaussPointMaterialCrossSectionStressResultantsContinuous();

  /** \brief append spatial cross-section stress resultant values at Gauss points to output data
   *
   *  \author grill
   *  \date 03/17 */
  void AppendGaussPointSpatialCrossSectionStressResultants();

  /**
   * \brief Append interpolated GP values of the spatial cross-section strain resultants.
   */
  void AppendGaussPointSpatialCrossSectionStressResultantsContinuous();

  /** \brief append element Orientation parameter with respect to x,y,z axis
   *
   *  \author eichinger
   *  \date 08/17 */
  void AppendElementOrientationParamater(
      Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector);

  /** \brief append sum of all element (node 0) internal energy cut in direction cut_dim
   *
   *  \author eichinger
   *  \date 08/17 */
  void AppendRVECrosssectionForces(
      Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector);

  /** \brief append current internal energy of the elements to output data
   *
   *  \author grill
   *  \date 03/17 */
  void AppendElementElasticEnergy();

  /**
   * \brief append the reference element length of the beam for the Hermitian interpolation.
   */
  void AppendRefLength();

  /**
   * \brief Write the visualization files to disk
   */
  void WriteToDisk(const double visualization_time, const int visualization_step);

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
      const std::vector<double>& vector_input, std::vector<double>& vector_output);

  /**
   * \brief Get the global (on all ranks) maximum value of Gauss Point resultants for stress /
   * strain output.
   **/
  int GetGlobalNumberOfGaussPointsPerBeam(unsigned int my_num_gp) const;

  /**
   * \brief Calculate the polynomial coefficients to interpolate the Gauss point values.
   * @param gauss_rule (in) Type of Gauss rule
   * @param gauss_point_values (in) Values at Gauss points
   * @param coefficients (out) Polynomial coefficients of the interpolated the Gauss point values.
   */
  void CalcInterpolationPolynomialCoefficients(const CORE::DRT::UTILS::GaussRule1D& gauss_rule,
      const std::vector<double>& gauss_point_values,
      std::vector<double>& polynomial_coefficients) const;

  /**
   * \brief Evaluate polynomial defined by its polynomial coefficients.
   * @param polynomial_coefficients (in)
   * @param xi (in) Point where the polynomial is evalated.
   * @return Interpolated value
   */
  double EvaluatePolynomialCoefficients(
      const std::vector<double>& polynomial_coefficients, const double& xi) const;

 private:
  /**
   * \brief Type of fields for continuous stress / strain output.
   */
  enum class StressStrainField
  {
    material_strain,
    material_stress,
  };

  /**
   * \brief Interpolate Gauss point values for stress / strain resultants along the beam, for each
   * visualization point along the beam centerline.
   * @param stress_strain_field (in) Type of stress / strain to append.
   */
  void AppendContinuousStressStrainResultants(const StressStrainField stress_strain_field);

 private:
  //! discretization containing beam elements of which geometry and result data shall be visualized
  Teuchos::RCP<const DRT::Discretization> discretization_;

  //! all local row indices of beam elements in the given discretization
  std::vector<unsigned int> local_row_indices_beam_elements_;

  //! periodic bounding box object
  Teuchos::RCP<const CORE::GEO::MESHFREE::BoundingBox> periodic_boundingbox_;

  //! number of points for each element ( in case of periodic boundary conditions
  //! not equal to 1)
  std::vector<int> num_cells_per_element_;

  //! the actual vtu writer object that additionally stores the geometry and result data
  Teuchos::RCP<IO::VisualizationManager> visualization_manager_;

  //! flag indicating whether to use absolute positions for point coordinates
  // false: use reference position
  bool use_absolute_positions_;

  //! number of visual subsegments for the beam polyline visialization.
  unsigned int n_subsegments_;
};

#endif
