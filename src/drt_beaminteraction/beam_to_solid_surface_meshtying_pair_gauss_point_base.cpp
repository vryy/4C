/*----------------------------------------------------------------------*/
/*! \file

\brief Gauss point to segment mesh tying element for between a 3D beam and a surface element.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_surface_meshtying_pair_gauss_point_base.H"

#include "beam_contact_params.H"
#include "beam_to_solid_surface_meshtying_params.H"
#include "../drt_geometry_pair/geometry_pair_scalar_types.H"
#include "../drt_geometry_pair/geometry_pair_element_functions.H"
#include "../drt_geometry_pair/geometry_pair_element_faces.H"



/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointBase<scalar_type, beam,
    surface>::BeamToSolidSurfaceMeshtyingPairGaussPointBase()
    : base_class()
{
  // Empty constructor.
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
double BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointBase<scalar_type, beam,
    surface>::GetEnergy() const
{
  return FADUTILS::CastToDouble(GetPenaltyPotential());
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
scalar_type BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointBase<scalar_type, beam,
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
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  template class BeamToSolidSurfaceMeshtyingPairGaussPointBase<
      line_to_surface_scalar_type<t_hermite, t_tri3>, t_hermite, t_tri3>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointBase<
      line_to_surface_scalar_type<t_hermite, t_tri6>, t_hermite, t_tri6>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointBase<
      line_to_surface_scalar_type<t_hermite, t_quad4>, t_hermite, t_quad4>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointBase<
      line_to_surface_scalar_type<t_hermite, t_quad8>, t_hermite, t_quad8>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointBase<
      line_to_surface_scalar_type<t_hermite, t_quad9>, t_hermite, t_quad9>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointBase<
      line_to_surface_scalar_type<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;

  template class BeamToSolidSurfaceMeshtyingPairGaussPointBase<line_to_surface_patch_scalar_type,
      t_hermite, t_tri3>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointBase<line_to_surface_patch_scalar_type,
      t_hermite, t_tri6>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointBase<line_to_surface_patch_scalar_type,
      t_hermite, t_quad4>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointBase<line_to_surface_patch_scalar_type,
      t_hermite, t_quad8>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointBase<line_to_surface_patch_scalar_type,
      t_hermite, t_quad9>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointBase<
      line_to_surface_patch_nurbs_scalar_type<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;
}  // namespace BEAMINTERACTION
