/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Write visualization output for a beam discretization in vtk/vtu format at runtime

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef FOUR_C_BEAM3_DISCRETIZATION_RUNTIME_VTU_WRITER_HPP
#define FOUR_C_BEAM3_DISCRETIZATION_RUNTIME_VTU_WRITER_HPP

/*-----------------------------------------------------------------------------------------------*/
/* headers */

#include "4C_config.hpp"

#include "4C_fem_general_utils_integration.hpp"
#include "4C_io_visualization_manager.hpp"

#include <Teuchos_RCP.hpp>

class Epetra_Comm;
class Epetra_Vector;
class Epetra_MultiVector;

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*/
/* forward declarations */
namespace Core::Geo
{
  namespace MeshFree
  {
    class BoundingBox;
  }
}  // namespace Core::Geo

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace ELEMENTS
  {
    class Beam3Base;
  }
}  // namespace Discret

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
class BeamDiscretizationRuntimeOutputWriter
{
 public:
  /// Constructor
  BeamDiscretizationRuntimeOutputWriter(
      Core::IO::VisualizationParameters parameters, const Epetra_Comm& comm);

  /// Destructor
  virtual ~BeamDiscretizationRuntimeOutputWriter() = default;
  /** \brief initialize object with all required data
   *
   *  \author grill
   *  \date 03/17 */
  void Initialize(Teuchos::RCP<Core::FE::Discretization> discretization,
      bool use_absolute_positions_for_point_coordinates, const unsigned int n_subsegments,
      Teuchos::RCP<const Core::Geo::MeshFree::BoundingBox> const& periodic_boundingbox);

  /** \brief append triad field determined from given displacement state to output data
   *
   *  \author grill
   *  \date 03/17 */
  void AppendTriadField(Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector);

  /** \brief append discplacement state
   *
   *  \author grill
   *  \date 03/17 */
  void append_displacement_field(
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
  void append_element_owning_processor();

  /**
   * \brief Append the 4C interal GIDs to all beam elements.
   */
  void AppendElementGID();

  /**
   * \brief Append the element ghosting information.
   */
  void append_element_ghosting_information();

  /** \brief append internal (elastic) energy of element
   *
   *  \author eichinger
   *  \date 01/18 */
  void append_element_internal_energy();

  /** \brief append kinetic energy of element
   *
   *  \author eichinger
   *  \date 01/18 */
  void append_element_kinetic_energy();

  /** \brief append information about to which filament an element belonging
   *
   *  \author eichinger
   *  \date 05/17 */
  void append_element_filament_id_and_type();

  /** \brief append circular cross-section radius of elements to output data
   *
   *  \author grill
   *  \date 03/17 */
  void append_element_circular_cross_section_radius();

  /** \brief append a vector field defining orientation and radius of a circular cross-section to
   * output data
   *
   *  \author grill
   *  \date 03/17 */
  void append_point_circular_cross_section_information_vector(
      Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector);

  /** \brief append material cross-section strain resultant values at Gauss points to output data
   *
   *  \author grill
   *  \date 03/17 */
  void append_gauss_point_material_cross_section_strain_resultants();

  /**
   * \brief Append interpolated GP values of the material cross-section strain resultants.
   */
  void append_gauss_point_material_cross_section_strain_resultants_continuous();

  /** \brief append material cross-section stress resultant values at Gauss points to output data
   *
   *  \author grill
   *  \date 03/17 */
  void append_gauss_point_material_cross_section_stress_resultants();

  /**
   * \brief Append interpolated GP values of the material cross-section stress resultants.
   */
  void append_gauss_point_material_cross_section_stress_resultants_continuous();

  /** \brief append spatial cross-section stress resultant values at Gauss points to output data
   *
   *  \author grill
   *  \date 03/17 */
  void append_gauss_point_spatial_cross_section_stress_resultants();

  /**
   * \brief Append interpolated GP values of the spatial cross-section strain resultants.
   */
  void append_gauss_point_spatial_cross_section_stress_resultants_continuous();

  /** \brief append element Orientation parameter with respect to x,y,z axis
   *
   *  \author eichinger
   *  \date 08/17 */
  void append_element_orientation_paramater(
      Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector);

  /** \brief append sum of all element (node 0) internal energy cut in direction cut_dim
   *
   *  \author eichinger
   *  \date 08/17 */
  void append_rve_crosssection_forces(
      Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector);

  /** \brief append current internal energy of the elements to output data
   *
   *  \author grill
   *  \date 03/17 */
  void append_element_elastic_energy();

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
  void set_geometry_from_beam_discretization(
      Teuchos::RCP<const Epetra_Vector> const& displacement_state_vector);

 private:
  /** \brief insert all values of a given std::vector at the end of another given std::vector
   *
   *  \author grill
   *  \date 03/17 */
  // Todo template <typename T>
  void insert_vector_values_at_back_of_other_vector(
      const std::vector<double>& vector_input, std::vector<double>& vector_output);

  /**
   * \brief Get the global (on all ranks) maximum value of Gauss Point resultants for stress /
   * strain output.
   **/
  int get_global_number_of_gauss_points_per_beam(unsigned int my_num_gp) const;

  /**
   * \brief Calculate the polynomial coefficients to interpolate the Gauss point values.
   * @param gauss_rule (in) Type of Gauss rule
   * @param gauss_point_values (in) Values at Gauss points
   * @param coefficients (out) Polynomial coefficients of the interpolated the Gauss point values.
   */
  void calc_interpolation_polynomial_coefficients(const Core::FE::GaussRule1D& gauss_rule,
      const std::vector<double>& gauss_point_values,
      std::vector<double>& polynomial_coefficients) const;

  /**
   * \brief Evaluate polynomial defined by its polynomial coefficients.
   * @param polynomial_coefficients (in)
   * @param xi (in) Point where the polynomial is evalated.
   * @return Interpolated value
   */
  double evaluate_polynomial_coefficients(
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
  void append_continuous_stress_strain_resultants(const StressStrainField stress_strain_field);

 private:
  //! discretization containing beam elements of which geometry and result data shall be visualized
  Teuchos::RCP<const Core::FE::Discretization> discretization_;

  //! all local row indices of beam elements in the given discretization
  std::vector<unsigned int> local_row_indices_beam_elements_;

  //! periodic bounding box object
  Teuchos::RCP<const Core::Geo::MeshFree::BoundingBox> periodic_boundingbox_;

  //! number of points for each element ( in case of periodic boundary conditions
  //! not equal to 1)
  std::vector<int> num_cells_per_element_;

  //! the actual vtu writer object that additionally stores the geometry and result data
  Teuchos::RCP<Core::IO::VisualizationManager> visualization_manager_;

  //! flag indicating whether to use absolute positions for point coordinates
  // false: use reference position
  bool use_absolute_positions_;

  //! number of visual subsegments for the beam polyline visialization.
  unsigned int n_subsegments_;
};

FOUR_C_NAMESPACE_CLOSE

#endif
