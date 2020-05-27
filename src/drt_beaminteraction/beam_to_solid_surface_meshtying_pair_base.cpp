/*----------------------------------------------------------------------*/
/*! \file

\brief Base meshtying element for meshtying between a 3D beam and a surface element.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_surface_meshtying_pair_base.H"

#include "beam_to_solid_vtu_output_writer_base.H"
#include "beam_to_solid_vtu_output_writer_visualization.H"
#include "beam_contact_params.H"
#include "beam_to_solid_surface_meshtying_params.H"
#include "../drt_geometry_pair/geometry_pair_line_to_surface.H"
#include "../drt_geometry_pair/geometry_pair_element_functions.H"
#include "../drt_geometry_pair/geometry_pair_factory.H"
#include "../drt_geometry_pair/geometry_pair_element_faces.H"
#include "../drt_geometry_pair/geometry_pair_scalar_types.H"


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
    this->ele1pos_(i) = FADUTILS::HigherOrderFadValue<scalar_type_fad>::apply(
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
        this->face_element_->GetFaceReferencePosition(), this->line_to_3D_segments_,
        this->face_element_->GetReferenceNormals());
  }
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam,
    surface>::GetPairVisualization(Teuchos::RCP<BeamToSolidVtuOutputWriterBase>
                                       visualization_writer,
    const Teuchos::ParameterList& visualization_params) const
{
  // Get visualization of base class.
  base_class::GetPairVisualization(visualization_writer, visualization_params);

  // If a writer exists for segmentation point data, add the segmentation point data.
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization>
      visualization_segmentation = visualization_writer->GetVisualizationWriter("segmentation");
  if (visualization_segmentation != Teuchos::null)
  {
    // Setup variables.
    LINALG::Matrix<3, 1, scalar_type> X_beam, u_beam, r_beam, r_solid, projection_dir;

    // Get the visualization vectors.
    std::vector<double>& point_coordinates =
        visualization_segmentation->GetMutablePointCoordinateVector();
    std::vector<double>& displacement =
        visualization_segmentation->GetMutablePointDataVector("displacement");
    std::vector<double>& projection_direction =
        visualization_segmentation->GetMutablePointDataVector("projection_direction");

    // Loop over the segments on the beam.
    for (const auto& segment : this->line_to_3D_segments_)
    {
      // Add the left and right boundary point of the segment.
      for (const auto& segmentation_point : {segment.GetStartPoint(), segment.GetEndPoint()})
      {
        GEOMETRYPAIR::EvaluatePosition<beam>(
            segmentation_point.GetEta(), this->ele1posref_, X_beam, this->Element1());
        GEOMETRYPAIR::EvaluatePosition<beam>(
            segmentation_point.GetEta(), this->ele1pos_, r_beam, this->Element1());
        u_beam = r_beam;
        u_beam -= X_beam;

        GEOMETRYPAIR::EvaluatePosition<surface>(segmentation_point.GetXi(),
            this->face_element_->GetFacePosition(), r_solid,
            this->face_element_->GetDrtFaceElement());
        projection_dir = r_solid;
        projection_dir -= r_beam;

        for (unsigned int dim = 0; dim < 3; dim++)
        {
          point_coordinates.push_back(FADUTILS::CastToDouble(X_beam(dim)));
          displacement.push_back(FADUTILS::CastToDouble(u_beam(dim)));
          projection_direction.push_back(FADUTILS::CastToDouble(projection_dir(dim)));
        }
      }
    }
  }
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam,
    surface>::CreateGeometryPair(const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>&
        geometry_evaluation_data_ptr)
{
  // Call the method of the base class.
  BeamContactPair::CreateGeometryPair(geometry_evaluation_data_ptr);

  // Set up the geometry pair, it will be initialized in the Init call of the base class.
  this->geometry_pair_ = GEOMETRYPAIR::GeometryPairLineToSurfaceFactory<double, beam, surface>(
      geometry_evaluation_data_ptr);
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

  // The second element in the pair has to be the face element.
  CastGeometryPair()->SetElement2(face_element_->GetDrtFaceElement());
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
double BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam, surface>::GetEnergy()
    const
{
  return FADUTILS::CastToDouble(GetPenaltyPotential());
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
scalar_type BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam,
    surface>::GetPenaltyPotential() const
{
  using namespace INPAR::BEAMTOSOLID;

  // If there are no intersection segments, no penalty potential exists for this pair.
  if (this->line_to_3D_segments_.size() == 0) return 0.0;

  // Set the DOF vectors, depending on the desired coupling type.
  LINALG::Matrix<beam::n_dof_, 1, scalar_type> beam_dof_fad;
  LINALG::Matrix<surface::n_dof_, 1, scalar_type> surface_dof_fad;
  const INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling coupling_type =
      this->Params()->BeamToSolidSurfaceMeshtyingParams()->GetCouplingType();
  if (coupling_type == BeamToSolidSurfaceCoupling::reference_configuration_forced_to_zero or
      coupling_type == BeamToSolidSurfaceCoupling::reference_configuration_forced_to_zero_fad or
      coupling_type == BeamToSolidSurfaceCoupling::consistent_fad)
  {
    beam_dof_fad = this->ele1pos_;
    surface_dof_fad = this->face_element_->GetFacePosition();
  }
  else if (coupling_type == BeamToSolidSurfaceCoupling::displacement or
           coupling_type == BeamToSolidSurfaceCoupling::displacement_fad)
  {
    beam_dof_fad = this->ele1pos_;
    for (unsigned int i_dof_beam = 0; i_dof_beam < beam::n_dof_; i_dof_beam++)
      beam_dof_fad(i_dof_beam) -= this->ele1posref_(i_dof_beam);
    surface_dof_fad = this->face_element_->GetFacePosition();
    for (unsigned int i_dof_surface = 0; i_dof_surface < surface::n_dof_; i_dof_surface++)
      surface_dof_fad(i_dof_surface) -=
          this->face_element_->GetFaceReferencePosition()(i_dof_surface);
  }
  else
    dserror(
        "BeamToSolidSurfaceMeshtyingPairGaussPoint::EvaluateAndAssemble: Got unexpected "
        "surface coupling type.");


  // Initialize variables for position and potential.
  LINALG::Matrix<3, 1, double> dr_beam_ref;
  LINALG::Matrix<3, 1, scalar_type> coupling_vector_beam;
  LINALG::Matrix<3, 1, scalar_type> coupling_vector_surface;
  scalar_type potential = 0.0;

  // Initialize scalar variables.
  double segment_jacobian, beam_segmentation_factor;
  double penalty_parameter =
      this->Params()->BeamToSolidSurfaceMeshtyingParams()->GetPenaltyParameter();

  // Integrate over segments.
  for (unsigned int i_segment = 0; i_segment < this->line_to_3D_segments_.size(); i_segment++)
  {
    // Factor to account for a segment length not from -1 to 1.
    beam_segmentation_factor = 0.5 * this->line_to_3D_segments_[i_segment].GetSegmentLength();

    // Gauss point loop.
    for (unsigned int i_gp = 0;
         i_gp < this->line_to_3D_segments_[i_segment].GetProjectionPoints().size(); i_gp++)
    {
      // Get the current Gauss point.
      const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& projected_gauss_point =
          this->line_to_3D_segments_[i_segment].GetProjectionPoints()[i_gp];

      // Get the Jacobian in the reference configuration.
      GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(
          projected_gauss_point.GetEta(), this->ele1posref_, dr_beam_ref, this->Element1());

      // Jacobian including the segment length.
      segment_jacobian = dr_beam_ref.Norm2() * beam_segmentation_factor;

      // Evaluate the coupling position for the beam.
      GEOMETRYPAIR::EvaluatePosition<beam>(
          projected_gauss_point.GetEta(), beam_dof_fad, coupling_vector_beam, this->Element1());

      // Evaluate the coupling position for the surface.
      if (coupling_type == BeamToSolidSurfaceCoupling::reference_configuration_forced_to_zero or
          coupling_type == BeamToSolidSurfaceCoupling::reference_configuration_forced_to_zero_fad or
          coupling_type == BeamToSolidSurfaceCoupling::displacement or
          coupling_type == BeamToSolidSurfaceCoupling::displacement_fad)
      {
        GEOMETRYPAIR::EvaluatePosition<surface>(projected_gauss_point.GetXi(), surface_dof_fad,
            coupling_vector_surface, this->face_element_->GetDrtFaceElement());
      }
      else if (coupling_type == BeamToSolidSurfaceCoupling::consistent_fad)
      {
        GEOMETRYPAIR::EvaluateSurfacePosition<surface>(projected_gauss_point.GetXi(),
            surface_dof_fad, coupling_vector_surface, this->face_element_->GetDrtFaceElement(),
            this->face_element_->GetCurrentNormals());
      }
      else
        dserror(
            "BeamToSolidSurfaceMeshtyingPairGaussPoint::EvaluateAndAssemble: Got unexpected "
            "surface coupling type.");

      // Calculate the difference between the coupling vectors and add the corresponding term to the
      // potential.
      coupling_vector_surface -= coupling_vector_beam;
      potential += projected_gauss_point.GetGaussWeight() * segment_jacobian *
                   coupling_vector_surface.Dot(coupling_vector_surface) * penalty_parameter * 0.5;
    }
  }

  return potential;
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
};


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
      line_to_surface_patch_nurbs_scalar_type<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;
}  // namespace BEAMINTERACTION
