#include "4C_beaminteraction_beam_to_solid_surface_contact_pair_base.hpp"

#include "4C_beam3_base.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_contact_params.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_visualization_output_params.hpp"
#include "4C_beaminteraction_beam_to_solid_utils.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_base.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_contact_params.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_element_faces.hpp"
#include "4C_geometry_pair_factory.hpp"
#include "4C_geometry_pair_line_to_surface.hpp"
#include "4C_geometry_pair_scalar_types.hpp"

#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
template <typename ScalarType, typename Beam, typename Surface>
BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<ScalarType, Beam,
    Surface>::BeamToSolidSurfaceContactPairBase()
    : base_class()
{
  // Empty constructor.
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Solid>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<ScalarType, Beam, Solid>::reset_state(
    const std::vector<double>& beam_centerline_dofvec,
    const std::vector<double>& solid_nodal_dofvec)
{
  // Clean the segments, as they will be re-evaluated in each iteration.
  this->line_to_3D_segments_.clear();

  // Set the current position of the beam element.
  const int n_patch_dof = face_element_->get_patch_gid().size();
  for (unsigned int i = 0; i < Beam::n_dof_; i++)
    this->ele1pos_.element_position_(i) = Core::FADUtils::HigherOrderFadValue<ScalarType>::apply(
        Beam::n_dof_ + n_patch_dof, i, beam_centerline_dofvec[i]);
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Surface>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<ScalarType, Beam, Surface>::pre_evaluate()
{
  // Call pre_evaluate on the geometry Pair.
  cast_geometry_pair()->pre_evaluate(
      this->ele1pos_, this->face_element_->get_face_element_data(), this->line_to_3D_segments_);
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Surface>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<ScalarType, Beam,
    Surface>::get_pair_visualization(Teuchos::RCP<BeamToSolidVisualizationOutputWriterBase>
                                         visualization_writer,
    Teuchos::ParameterList& visualization_params) const
{
  // Get visualization of base class.
  base_class::get_pair_visualization(visualization_writer, visualization_params);

  // Add segmentation and integration point data.
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_segmentation =
      visualization_writer->get_visualization_writer("btss-contact-segmentation");
  if (visualization_segmentation != Teuchos::null)
  {
    std::vector<GEOMETRYPAIR::ProjectionPoint1DTo3D<ScalarType>> points;
    for (const auto& segment : this->line_to_3D_segments_)
      for (const auto& segmentation_point : {segment.get_start_point(), segment.get_end_point()})
        points.push_back(segmentation_point);
    add_visualization_integration_points(*visualization_segmentation, points, visualization_params);
  }

  Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization>
      visualization_integration_points =
          visualization_writer->get_visualization_writer("btss-contact-integration-points");
  if (visualization_integration_points != Teuchos::null)
  {
    std::vector<GEOMETRYPAIR::ProjectionPoint1DTo3D<ScalarType>> points;
    for (const auto& segment : this->line_to_3D_segments_)
      for (const auto& segmentation_point : (segment.get_projection_points()))
        points.push_back(segmentation_point);
    add_visualization_integration_points(
        *visualization_integration_points, points, visualization_params);
  }
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Surface>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<ScalarType, Beam, Surface>::
    add_visualization_integration_points(
        BEAMINTERACTION::BeamToSolidOutputWriterVisualization& visualization_writer,
        const std::vector<GEOMETRYPAIR::ProjectionPoint1DTo3D<ScalarType>>& points,
        const Teuchos::ParameterList& visualization_params) const
{
  auto& visualization_data = visualization_writer.get_visualization_data();

  // Setup variables.
  Core::LinAlg::Matrix<3, 1, ScalarType> X_beam, u_beam;

  // Get beam cross-section diameter.
  auto beam_ptr = dynamic_cast<const Discret::ELEMENTS::Beam3Base*>(this->element1());
  const double beam_cross_section_radius =
      beam_ptr->get_circular_cross_section_radius_for_interactions();

  // Get the visualization vectors.
  std::vector<double>& point_coordinates = visualization_data.get_point_coordinates();
  std::vector<double>& displacement = visualization_data.get_point_data<double>("displacement");
  std::vector<double>& surface_normal_data =
      visualization_data.get_point_data<double>("surface_normal");
  std::vector<double>& gap_data = visualization_data.get_point_data<double>("gap");
  std::vector<double>& force_data = visualization_data.get_point_data<double>("force");

  const Teuchos::RCP<const BeamToSolidSurfaceVisualizationOutputParams>& output_params_ptr =
      visualization_params.get<Teuchos::RCP<const BeamToSolidSurfaceVisualizationOutputParams>>(
          "btss-output_params_ptr");
  const bool write_unique_ids = output_params_ptr->get_write_unique_i_ds_flag();
  std::vector<int>* pair_beam_id = nullptr;
  std::vector<int>* pair_solid_id = nullptr;
  if (write_unique_ids)
  {
    pair_beam_id = &(visualization_data.get_point_data<int>("uid_0_pair_beam_id"));
    pair_solid_id = &(visualization_data.get_point_data<int>("uid_1_pair_solid_id"));
  }

  for (const auto& point : points)
  {
    const auto [r_beam, r_surface, surface_normal, gap] =
        this->evaluate_contact_kinematics_at_projection_point(point, beam_cross_section_radius);
    GEOMETRYPAIR::evaluate_position<Beam>(point.get_eta(), this->ele1posref_, X_beam);

    u_beam = r_beam;
    u_beam -= X_beam;
    const auto force = penalty_force(gap, *this->params()->beam_to_solid_surface_contact_params());

    for (unsigned int dim = 0; dim < 3; dim++)
    {
      point_coordinates.push_back(Core::FADUtils::cast_to_double(X_beam(dim)));
      displacement.push_back(Core::FADUtils::cast_to_double(u_beam(dim)));
      surface_normal_data.push_back(Core::FADUtils::cast_to_double(surface_normal(dim)));
      force_data.push_back(Core::FADUtils::cast_to_double(force * surface_normal(dim)));
    }
    gap_data.push_back(Core::FADUtils::cast_to_double(gap));

    if (write_unique_ids)
    {
      pair_beam_id->push_back(this->element1()->id());
      pair_solid_id->push_back(this->element2()->id());
    }
  }
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Surface>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<ScalarType, Beam,
    Surface>::create_geometry_pair(const Core::Elements::Element* element1,
    const Core::Elements::Element* element2,
    const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>& geometry_evaluation_data_ptr)
{
  this->geometry_pair_ =
      GEOMETRYPAIR::geometry_pair_line_to_surface_factory_fad<ScalarType, Beam, Surface>(
          element1, element2, geometry_evaluation_data_ptr);
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Surface>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<ScalarType, Beam,
    Surface>::set_face_element(Teuchos::RCP<GEOMETRYPAIR::FaceElement>& face_element)
{
  face_element_ = Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::FaceElementTemplate<Surface, ScalarType>>(
      face_element, true);

  // Set the number of (centerline) degrees of freedom for the beam element in the face element
  face_element_->set_number_of_dof_other_element(
      Utils::get_number_of_element_centerline_dof(this->element1()));

  // If the solid surface is the surface of a 3D volume we set the face element here. Otherwise we
  // simply set the same element again.
  cast_geometry_pair()->set_element2(face_element_->get_element());
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Surface>
Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToSurface<ScalarType, Beam, Surface>>
BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<ScalarType, Beam, Surface>::cast_geometry_pair()
    const
{
  return Teuchos::rcp_dynamic_cast<
      GEOMETRYPAIR::GeometryPairLineToSurface<ScalarType, Beam, Surface>>(
      this->geometry_pair_, true);
}

/**
 *
 */
template <typename ScalarType, typename Beam, typename Surface>
std::tuple<Core::LinAlg::Matrix<3, 1, ScalarType>, Core::LinAlg::Matrix<3, 1, ScalarType>,
    Core::LinAlg::Matrix<3, 1, ScalarType>, ScalarType>
BEAMINTERACTION::BeamToSolidSurfaceContactPairBase<ScalarType, Beam, Surface>::
    evaluate_contact_kinematics_at_projection_point(
        const GEOMETRYPAIR::ProjectionPoint1DTo3D<ScalarType>& projection_point,
        const double beam_cross_section_radius) const
{
  // Get the projection coordinates
  const auto& xi = projection_point.get_xi();
  const auto& eta = projection_point.get_eta();

  // Get the surface normal vector
  Core::LinAlg::Matrix<3, 1, ScalarType> surface_normal;
  GEOMETRYPAIR::evaluate_surface_normal<Surface>(
      xi, this->face_element_->get_face_element_data(), surface_normal);

  // Evaluate the current position of beam and solid
  Core::LinAlg::Matrix<3, 1, ScalarType> r_beam;
  Core::LinAlg::Matrix<3, 1, ScalarType> r_surface;
  GEOMETRYPAIR::evaluate_position<Beam>(eta, this->ele1pos_, r_beam);
  GEOMETRYPAIR::evaluate_position<Surface>(
      xi, this->face_element_->get_face_element_data(), r_surface);

  // Evaluate the gap function
  Core::LinAlg::Matrix<3, 1, ScalarType> r_rel;
  r_rel = r_beam;
  r_rel -= r_surface;
  ScalarType gap = r_rel.dot(surface_normal) - beam_cross_section_radius;

  return {r_beam, r_surface, surface_normal, gap};
}

/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_tri3>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_tri6>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_quad4>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_quad8>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_quad9>;
  template class BeamToSolidSurfaceContactPairBase<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_line2, t_nurbs9>, t_line2, t_nurbs9>;

  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_line2,
      t_tri3>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_line2,
      t_tri6>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_line2,
      t_quad4>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_line2,
      t_quad8>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_line2,
      t_quad9>;
  template class BeamToSolidSurfaceContactPairBase<
      line_to_surface_patch_scalar_type_fixed_size<t_line2, t_nurbs9>, t_line2, t_nurbs9>;


  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_tri3>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_tri6>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad4>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad8>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad9>;
  template class BeamToSolidSurfaceContactPairBase<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_hermite, t_nurbs9>, t_hermite,
      t_nurbs9>;

  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_tri3>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_tri6>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_quad4>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_quad8>;
  template class BeamToSolidSurfaceContactPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_quad9>;
  template class BeamToSolidSurfaceContactPairBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;

}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE
