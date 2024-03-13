/*----------------------------------------------------------------------*/
/*! \file

\brief Base meshtying element for meshtying between a 3D beam and a surface element.

\level 3
*/


#include "baci_beaminteraction_beam_to_solid_surface_meshtying_pair_base.hpp"

#include "baci_beaminteraction_beam_to_solid_surface_meshtying_params.hpp"
#include "baci_beaminteraction_beam_to_solid_surface_visualization_output_params.hpp"
#include "baci_beaminteraction_beam_to_solid_visualization_output_writer_base.hpp"
#include "baci_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"
#include "baci_beaminteraction_calc_utils.hpp"
#include "baci_beaminteraction_contact_params.hpp"
#include "baci_geometry_pair_element_evaluation_functions.hpp"
#include "baci_geometry_pair_element_faces.hpp"
#include "baci_geometry_pair_factory.hpp"
#include "baci_geometry_pair_line_to_surface.hpp"
#include "baci_geometry_pair_scalar_types.hpp"

BACI_NAMESPACE_OPEN


/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam,
    surface>::BeamToSolidSurfaceMeshtyingPairBase()
    : base_class(), meshtying_is_evaluated_(false)
{
  // Empty constructor.
}

/**
 *
 */
template <typename scalar_type_fad, typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type_fad, beam, solid>::ResetState(
    const std::vector<double>& beam_centerline_dofvec,
    const std::vector<double>& solid_nodal_dofvec)
{
  // Beam element.
  const int n_patch_dof = face_element_->GetPatchGID().size();
  for (unsigned int i = 0; i < beam::n_dof_; i++)
    this->ele1pos_.element_position_(i) =
        CORE::FADUTILS::HigherOrderFadValue<scalar_type_fad>::apply(
            beam::n_dof_ + n_patch_dof, i, beam_centerline_dofvec[i]);
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam, surface>::PreEvaluate()
{
  // Call PreEvaluate on the geometry Pair.
  if (!meshtying_is_evaluated_)
  {
    CastGeometryPair()->PreEvaluate(this->ele1posref_,
        this->face_element_->GetFaceReferenceElementData(), this->line_to_3D_segments_);
  }
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam,
    surface>::GetPairVisualization(Teuchos::RCP<BeamToSolidVisualizationOutputWriterBase>
                                       visualization_writer,
    Teuchos::ParameterList& visualization_params) const
{
  // Get visualization of base class.
  base_class::GetPairVisualization(visualization_writer, visualization_params);

  // Add segmentation and integration point data.
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_segmentation =
      visualization_writer->GetVisualizationWriter("btssc-segmentation");
  if (visualization_segmentation != Teuchos::null)
  {
    std::vector<GEOMETRYPAIR::ProjectionPoint1DTo3D<double>> points;
    for (const auto& segment : this->line_to_3D_segments_)
      for (const auto& segmentation_point : {segment.GetStartPoint(), segment.GetEndPoint()})
        points.push_back(segmentation_point);
    AddVisualizationIntegrationPoints(visualization_segmentation, points, visualization_params);
  }

  Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization>
      visualization_integration_points =
          visualization_writer->GetVisualizationWriter("btssc-integration-points");
  if (visualization_integration_points != Teuchos::null)
  {
    std::vector<GEOMETRYPAIR::ProjectionPoint1DTo3D<double>> points;
    for (const auto& segment : this->line_to_3D_segments_)
      for (const auto& segmentation_point : (segment.GetProjectionPoints()))
        points.push_back(segmentation_point);
    AddVisualizationIntegrationPoints(
        visualization_integration_points, points, visualization_params);
  }
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam, surface>::
    AddVisualizationIntegrationPoints(
        const Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization>&
            visualization_writer,
        const std::vector<GEOMETRYPAIR::ProjectionPoint1DTo3D<double>>& points,
        const Teuchos::ParameterList& visualization_params) const
{
  auto& visualization_data = visualization_writer->GetVisualizationData();

  // Setup variables.
  CORE::LINALG::Matrix<3, 1, scalar_type> X_beam, u_beam, r_beam, r_solid, projection_dir;

  // Get the visualization vectors.
  std::vector<double>& point_coordinates = visualization_data.GetPointCoordinates();
  std::vector<double>& displacement = visualization_data.GetPointData<double>("displacement");
  std::vector<double>& projection_direction =
      visualization_data.GetPointData<double>("projection_direction");

  const Teuchos::RCP<const BeamToSolidSurfaceVisualizationOutputParams>& output_params_ptr =
      visualization_params.get<Teuchos::RCP<const BeamToSolidSurfaceVisualizationOutputParams>>(
          "btssc-output_params_ptr");
  const bool write_unique_ids = output_params_ptr->GetWriteUniqueIDsFlag();
  std::vector<double>* pair_beam_id = nullptr;
  std::vector<double>* pair_solid_id = nullptr;
  if (write_unique_ids)
  {
    pair_beam_id = &(visualization_data.GetPointData<double>("uid_0_pair_beam_id"));
    pair_solid_id = &(visualization_data.GetPointData<double>("uid_1_pair_solid_id"));
  }

  for (const auto& point : points)
  {
    GEOMETRYPAIR::EvaluatePosition<beam>(point.GetEta(), this->ele1posref_, X_beam);
    GEOMETRYPAIR::EvaluatePosition<beam>(point.GetEta(), this->ele1pos_, r_beam);
    u_beam = r_beam;
    u_beam -= X_beam;

    GEOMETRYPAIR::EvaluatePosition<surface>(
        point.GetXi(), this->face_element_->GetFaceElementData(), r_solid);
    projection_dir = r_solid;
    projection_dir -= r_beam;

    for (unsigned int dim = 0; dim < 3; dim++)
    {
      point_coordinates.push_back(CORE::FADUTILS::CastToDouble(X_beam(dim)));
      displacement.push_back(CORE::FADUTILS::CastToDouble(u_beam(dim)));
      projection_direction.push_back(CORE::FADUTILS::CastToDouble(projection_dir(dim)));
    }

    if (write_unique_ids)
    {
      pair_beam_id->push_back(this->Element1()->Id());
      pair_solid_id->push_back(this->Element2()->Id());
    }
  }
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam,
    surface>::CreateGeometryPair(const DRT::Element* element1, const DRT::Element* element2,
    const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>& geometry_evaluation_data_ptr)
{
  this->geometry_pair_ = GEOMETRYPAIR::GeometryPairLineToSurfaceFactory<double, beam, surface>(
      element1, element2, geometry_evaluation_data_ptr);
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam,
    surface>::SetFaceElement(Teuchos::RCP<GEOMETRYPAIR::FaceElement>& face_element)
{
  face_element_ =
      Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::FaceElementTemplate<surface, scalar_type>>(
          face_element, true);

  // Set the number of (centerline) degrees of freedom for the beam element in the face element
  face_element_->SetNumberOfDofOtherElement(
      UTILS::GetNumberOfElementCenterlineDof(this->Element1()));

  // The second element in the pair has to be the face element.
  CastGeometryPair()->SetElement2(face_element_->GetDrtFaceElement());
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToSurface<double, beam, surface>>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam, surface>::CastGeometryPair()
    const
{
  return Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::GeometryPairLineToSurface<double, beam, surface>>(
      this->geometry_pair_, true);
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
CORE::LINALG::Matrix<3, 1, scalar_type>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam, surface>::EvaluateCoupling(
    const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& evaluation_point) const
{
  using namespace INPAR::BEAMTOSOLID;

  CORE::LINALG::Matrix<3, 1, scalar_type> r_beam(true);
  CORE::LINALG::Matrix<3, 1, scalar_type> r_surface(true);

  const BeamToSolidSurfaceCoupling coupling_type =
      this->Params()->BeamToSolidSurfaceMeshtyingParams()->GetCouplingType();
  switch (coupling_type)
  {
    case BeamToSolidSurfaceCoupling::displacement:
    case BeamToSolidSurfaceCoupling::displacement_fad:
    {
      // In this case we have to substract the reference position from the DOF vectors.
      auto beam_dof = this->ele1pos_;
      auto surface_dof = this->face_element_->GetFaceElementData();

      for (unsigned int i_dof_beam = 0; i_dof_beam < beam::n_dof_; i_dof_beam++)
        beam_dof.element_position_(i_dof_beam) -= this->ele1posref_.element_position_(i_dof_beam);
      for (unsigned int i_dof_surface = 0; i_dof_surface < surface::n_dof_; i_dof_surface++)
        surface_dof.element_position_(i_dof_surface) -=
            this->face_element_->GetFaceReferenceElementData().element_position_(i_dof_surface);

      GEOMETRYPAIR::EvaluatePosition<beam>(evaluation_point.GetEta(), beam_dof, r_beam);
      GEOMETRYPAIR::EvaluatePosition<surface>(evaluation_point.GetXi(), surface_dof, r_surface);

      r_beam -= r_surface;
      return r_beam;
    }
    case BeamToSolidSurfaceCoupling::reference_configuration_forced_to_zero:
    case BeamToSolidSurfaceCoupling::reference_configuration_forced_to_zero_fad:
    {
      GEOMETRYPAIR::EvaluatePosition<beam>(evaluation_point.GetEta(), this->ele1pos_, r_beam);
      GEOMETRYPAIR::EvaluatePosition<surface>(
          evaluation_point.GetXi(), this->face_element_->GetFaceElementData(), r_surface);

      r_beam -= r_surface;
      return r_beam;
    }
    case BeamToSolidSurfaceCoupling::consistent_fad:
    {
      GEOMETRYPAIR::EvaluatePosition<beam>(evaluation_point.GetEta(), this->ele1pos_, r_beam);
      GEOMETRYPAIR::EvaluateSurfacePosition<surface>(
          evaluation_point.GetXi(), this->face_element_->GetFaceElementData(), r_surface);

      r_beam -= r_surface;
      return r_beam;
    }
    default:
    {
      dserror("Got unexpected coupling type.");
      return r_beam;
    }
  }
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
std::vector<int>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam, surface>::GetPairGID(
    const DRT::Discretization& discret) const
{
  // Get the beam centerline GIDs.
  CORE::LINALG::Matrix<beam::n_dof_, 1, int> beam_centerline_gid;
  UTILS::GetElementCenterlineGIDIndices(discret, this->Element1(), beam_centerline_gid);

  // Get the patch (in this case just the one face element) GIDs.
  const std::vector<int>& patch_gid = this->face_element_->GetPatchGID();
  std::vector<int> pair_gid;
  pair_gid.resize(beam::n_dof_ + patch_gid.size());

  // Combine beam and solid GIDs into one vector.
  for (unsigned int i_dof_beam = 0; i_dof_beam < beam::n_dof_; i_dof_beam++)
    pair_gid[i_dof_beam] = beam_centerline_gid(i_dof_beam);
  for (unsigned int i_dof_patch = 0; i_dof_patch < patch_gid.size(); i_dof_patch++)
    pair_gid[beam::n_dof_ + i_dof_patch] = patch_gid[i_dof_patch];

  return pair_gid;
}


/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  template class BeamToSolidSurfaceMeshtyingPairBase<
      line_to_surface_scalar_type<t_hermite, t_quad4>, t_hermite, t_quad4>;
  template class BeamToSolidSurfaceMeshtyingPairBase<
      line_to_surface_scalar_type<t_hermite, t_quad8>, t_hermite, t_quad8>;
  template class BeamToSolidSurfaceMeshtyingPairBase<
      line_to_surface_scalar_type<t_hermite, t_quad9>, t_hermite, t_quad9>;
  template class BeamToSolidSurfaceMeshtyingPairBase<line_to_surface_scalar_type<t_hermite, t_tri3>,
      t_hermite, t_tri3>;
  template class BeamToSolidSurfaceMeshtyingPairBase<line_to_surface_scalar_type<t_hermite, t_tri6>,
      t_hermite, t_tri6>;
  template class BeamToSolidSurfaceMeshtyingPairBase<
      line_to_surface_scalar_type<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;

  template class BeamToSolidSurfaceMeshtyingPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_quad4>;
  template class BeamToSolidSurfaceMeshtyingPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_quad8>;
  template class BeamToSolidSurfaceMeshtyingPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_quad9>;
  template class BeamToSolidSurfaceMeshtyingPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_tri3>;
  template class BeamToSolidSurfaceMeshtyingPairBase<line_to_surface_patch_scalar_type, t_hermite,
      t_tri6>;
  template class BeamToSolidSurfaceMeshtyingPairBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;
  template class BeamToSolidSurfaceMeshtyingPairBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex8>, t_hermite, t_quad4>;
  template class BeamToSolidSurfaceMeshtyingPairBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex20>, t_hermite, t_quad8>;
  template class BeamToSolidSurfaceMeshtyingPairBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex27>, t_hermite, t_quad9>;
}  // namespace BEAMINTERACTION

BACI_NAMESPACE_CLOSE
