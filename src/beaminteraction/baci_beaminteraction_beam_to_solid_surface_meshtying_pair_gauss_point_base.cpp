/*----------------------------------------------------------------------*/
/*! \file

\brief Gauss point to segment mesh tying element for between a 3D beam and a surface element.

\level 3
*/


#include "baci_beaminteraction_beam_to_solid_surface_meshtying_pair_gauss_point_base.hpp"

#include "baci_beaminteraction_beam_to_solid_surface_meshtying_params.hpp"
#include "baci_beaminteraction_contact_params.hpp"
#include "baci_geometry_pair_element_evaluation_functions.hpp"
#include "baci_geometry_pair_element_faces.hpp"
#include "baci_geometry_pair_scalar_types.hpp"

FOUR_C_NAMESPACE_OPEN



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
  return CORE::FADUTILS::CastToDouble(GetPenaltyPotential());
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

  // Initialize variables for position and potential.
  CORE::LINALG::Matrix<3, 1, double> dr_beam_ref;
  CORE::LINALG::Matrix<3, 1, scalar_type> coupling_vector;
  scalar_type potential = 0.0;

  // Initialize scalar variables.
  double segment_jacobian = 0.0;
  double beam_segmentation_factor = 0.0;
  double penalty_parameter =
      this->Params()->BeamToSolidSurfaceMeshtyingParams()->GetPenaltyParameter();

  // Integrate over segments.
  const unsigned int n_segments = this->line_to_3D_segments_.size();
  for (unsigned int i_segment = 0; i_segment < n_segments; i_segment++)
  {
    // Factor to account for the integration segment length.
    beam_segmentation_factor = 0.5 * this->line_to_3D_segments_[i_segment].GetSegmentLength();

    // Gauss point loop.
    const unsigned int n_gp = this->line_to_3D_segments_[i_segment].GetProjectionPoints().size();
    for (unsigned int i_gp = 0; i_gp < n_gp; i_gp++)
    {
      // Get the current Gauss point.
      const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& projected_gauss_point =
          this->line_to_3D_segments_[i_segment].GetProjectionPoints()[i_gp];

      // Get the Jacobian in the reference configuration.
      GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(
          projected_gauss_point.GetEta(), this->ele1posref_, dr_beam_ref);

      // Jacobian including the segment length.
      segment_jacobian = dr_beam_ref.Norm2() * beam_segmentation_factor;

      // Evaluate the coupling vector.
      coupling_vector = this->EvaluateCoupling(projected_gauss_point);

      // Calculate the difference between the coupling vectors and add the
      // corresponding term to the potential.
      potential += projected_gauss_point.GetGaussWeight() * segment_jacobian *
                   coupling_vector.Dot(coupling_vector) * penalty_parameter * 0.5;
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
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex8>, t_hermite, t_quad4>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex20>, t_hermite, t_quad8>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex27>, t_hermite, t_quad9>;
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE
