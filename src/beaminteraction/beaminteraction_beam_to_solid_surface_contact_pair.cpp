/*----------------------------------------------------------------------*/
/*! \file

\brief Contact element for contact between a 3D beam and a surface element.

\level 3
*/


#include "beaminteraction_beam_to_solid_surface_contact_pair.H"

#include "beaminteraction_contact_params.H"
#include "beaminteraction_beam_to_solid_surface_contact_params.H"
#include "beaminteraction_beam_to_solid_vtu_output_writer_base.H"
#include "beaminteraction_beam_to_solid_vtu_output_writer_visualization.H"
#include "beaminteraction_beam_to_solid_surface_vtk_output_params.H"
#include "beaminteraction_beam_to_solid_utils.H"
#include "beaminteraction_calc_utils.H"
#include "geometry_pair_line_to_surface.H"
#include "geometry_pair_element_functions.H"
#include "geometry_pair_factory.H"
#include "geometry_pair_element_faces.H"
#include "geometry_pair_scalar_types.H"

#include <Epetra_FEVector.h>


/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
BEAMINTERACTION::BeamToSolidSurfaceContactPairGapVariation<scalar_type, beam,
    surface>::BeamToSolidSurfaceContactPairGapVariation()
    : base_class()
{
  // Empty constructor.
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairGapVariation<scalar_type, beam,
    surface>::EvaluateAndAssemble(const Teuchos::RCP<const ::DRT::Discretization>& discret,
    const Teuchos::RCP<Epetra_FEVector>& force_vector,
    const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix,
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector)
{
  // Call Evaluate on the geometry Pair.
  this->CastGeometryPair()->Evaluate(this->ele1pos_, this->face_element_->GetFacePosition(),
      this->line_to_3D_segments_, this->face_element_->GetCurrentNormals());

  // If there are no intersection segments, no contact terms will be assembled.
  const unsigned int n_segments = this->line_to_3D_segments_.size();
  if (n_segments == 0) return;

  // Get beam cross-section diameter.
  auto beam_ptr = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(this->Element1());
  const double beam_cross_section_radius = beam_ptr->GetCircularCrossSectionRadiusForInteractions();

  // Initialize variables for contact kinematics.
  LINALG::Matrix<3, 1, scalar_type> dr_beam_ref;
  LINALG::Matrix<3, 1, scalar_type> surface_normal;
  LINALG::Matrix<3, 1, scalar_type> r_beam;
  LINALG::Matrix<3, 1, scalar_type> r_surface;
  LINALG::Matrix<3, 1, scalar_type> r_rel;
  LINALG::Matrix<1, beam::n_nodes_ * beam::n_val_, scalar_type> N_beam;
  LINALG::Matrix<1, surface::n_nodes_ * surface::n_val_, scalar_type> N_surface;
  LINALG::Matrix<beam::n_dof_ + surface::n_dof_, 1, scalar_type> gap_variation_times_normal;
  LINALG::Matrix<beam::n_dof_ + surface::n_dof_, 1, scalar_type> pair_force_vector;
  scalar_type gap = 0.0;
  scalar_type segment_jacobian = 0.0;
  scalar_type beam_segmentation_factor = 0.0;

  // GIDs of the pair and the force vector acting on the pair.
  const std::vector<int> pair_gid = this->GetPairGID(*discret);

  // Integrate over segments.
  for (unsigned int i_segment = 0; i_segment < n_segments; i_segment++)
  {
    // Factor to account for the integration segment length.
    beam_segmentation_factor = 0.5 * this->line_to_3D_segments_[i_segment].GetSegmentLength();

    // Gauss point loop.
    const unsigned int n_gp = this->line_to_3D_segments_[i_segment].GetProjectionPoints().size();
    for (unsigned int i_gp = 0; i_gp < n_gp; i_gp++)
    {
      // Get the current Gauss point.
      const auto& projected_gauss_point =
          this->line_to_3D_segments_[i_segment].GetProjectionPoints()[i_gp];

      // Get the fixed and variable projection coordinates.
      const auto& xi = projected_gauss_point.GetXi();
      const auto& eta = projected_gauss_point.GetEta();

      // Get the Jacobian in the reference configuration.
      GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(
          eta, this->ele1posref_, dr_beam_ref, this->Element1());

      // Jacobian including the segment length.
      segment_jacobian = FADUTILS::VectorNorm(dr_beam_ref) * beam_segmentation_factor;

      // Get the surface normal vector.
      GEOMETRYPAIR::EvaluateSurfaceNormal<surface>(xi, this->face_element_->GetFacePosition(),
          surface_normal, this->face_element_->GetDrtFaceElement(),
          this->face_element_->GetCurrentNormals());

      // Evaluate the current position of beam and solid.
      GEOMETRYPAIR::EvaluatePosition<beam>(eta, this->ele1pos_, r_beam, this->Element1());
      GEOMETRYPAIR::EvaluatePosition<surface>(xi, this->face_element_->GetFacePosition(), r_surface,
          this->face_element_->GetDrtFaceElement());

      // Evaluate the gap function.
      r_rel = r_beam;
      r_rel -= r_surface;
      gap = r_rel.Dot(surface_normal) - beam_cross_section_radius;

      // Get the shape function matrices.
      beam::EvaluateShapeFunction(
          N_beam, eta, std::integral_constant<unsigned int, beam::dim_>{}, this->Element1());
      surface::EvaluateShapeFunction(N_surface, xi,
          std::integral_constant<unsigned int, surface::dim_>{},
          this->face_element_->GetDrtFaceElement());

      // Calculate the variation of the gap function multiplied with the surface normal vector.
      for (unsigned int i_shape = 0; i_shape < N_beam.N(); i_shape++)
      {
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        {
          gap_variation_times_normal(i_shape * 3 + i_dim) = N_beam(i_shape) * surface_normal(i_dim);
        }
      }
      for (unsigned int i_shape = 0; i_shape < N_surface.N(); i_shape++)
      {
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        {
          gap_variation_times_normal(N_beam.N() * 3 + i_shape * 3 + i_dim) =
              -1.0 * N_surface(i_shape) * surface_normal(i_dim);
        }
      }

      // Get the contact force.
      scalar_type force = PenaltyForce(gap, this->Params()->BeamToSolidSurfaceContactParams());

      // Add the Gauss point contributions to the pair force vector.
      gap_variation_times_normal.Scale(
          force * projected_gauss_point.GetGaussWeight() * segment_jacobian);
      pair_force_vector -= gap_variation_times_normal;
    }
  }

  // If given, assemble force terms into the global vector.
  if (force_vector != Teuchos::null)
  {
    std::vector<double> force_pair_double(pair_gid.size(), 0.0);
    for (unsigned int j_dof = 0; j_dof < pair_force_vector.M(); j_dof++)
      force_pair_double[j_dof] = FADUTILS::CastToDouble(pair_force_vector(j_dof));
    force_vector->SumIntoGlobalValues(pair_gid.size(), pair_gid.data(), force_pair_double.data());
  }

  // If given, assemble force terms into the global stiffness matrix.
  if (stiffness_matrix != Teuchos::null)
    for (unsigned int i_dof = 0; i_dof < pair_force_vector.M(); i_dof++)
      for (unsigned int j_dof = 0; j_dof < pair_gid.size(); j_dof++)
        stiffness_matrix->FEAssemble(FADUTILS::CastToDouble(pair_force_vector(i_dof).dx(j_dof)),
            pair_gid[i_dof], pair_gid[j_dof]);
}


/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
BEAMINTERACTION::BeamToSolidSurfaceContactPairPotential<scalar_type, beam,
    surface>::BeamToSolidSurfaceContactPairPotential()
    : base_class()
{
  // Empty constructor.
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairPotential<scalar_type, beam,
    surface>::EvaluateAndAssemble(const Teuchos::RCP<const ::DRT::Discretization>& discret,
    const Teuchos::RCP<Epetra_FEVector>& force_vector,
    const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix,
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector)
{
  // Call Evaluate on the geometry Pair.
  this->CastGeometryPair()->Evaluate(this->ele1pos_, this->face_element_->GetFacePosition(),
      this->line_to_3D_segments_, this->face_element_->GetCurrentNormals());

  // If there are no intersection segments, no contact terms will be assembled.
  const unsigned int n_segments = this->line_to_3D_segments_.size();
  if (n_segments == 0) return;

  // Get beam cross-section diameter.
  auto beam_ptr = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(this->Element1());
  const double beam_cross_section_radius = beam_ptr->GetCircularCrossSectionRadiusForInteractions();

  // Initialize variables for contact kinematics.
  LINALG::Matrix<3, 1, scalar_type> dr_beam_ref;
  LINALG::Matrix<3, 1, scalar_type> surface_normal;
  LINALG::Matrix<3, 1, scalar_type> r_beam;
  LINALG::Matrix<3, 1, scalar_type> r_surface;
  LINALG::Matrix<3, 1, scalar_type> r_rel;
  scalar_type potential = 0.0;
  scalar_type gap = 0.0;
  scalar_type segment_jacobian = 0.0;
  scalar_type beam_segmentation_factor = 0.0;

  // GIDs of the pair and the force vector acting on the pair.
  const std::vector<int> pair_gid = this->GetPairGID(*discret);

  // Integrate over segments.
  for (unsigned int i_segment = 0; i_segment < n_segments; i_segment++)
  {
    // Factor to account for the integration segment length.
    beam_segmentation_factor = 0.5 * this->line_to_3D_segments_[i_segment].GetSegmentLength();

    // Gauss point loop.
    const unsigned int n_gp = this->line_to_3D_segments_[i_segment].GetProjectionPoints().size();
    for (unsigned int i_gp = 0; i_gp < n_gp; i_gp++)
    {
      // Get the current Gauss point.
      const auto& projected_gauss_point =
          this->line_to_3D_segments_[i_segment].GetProjectionPoints()[i_gp];

      // Get the fixed and variable projection coordinates.
      const auto& xi = projected_gauss_point.GetXi();
      const auto& eta = projected_gauss_point.GetEta();

      // Get the Jacobian in the reference configuration.
      GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(
          eta, this->ele1posref_, dr_beam_ref, this->Element1());

      // Jacobian including the segment length.
      segment_jacobian = FADUTILS::VectorNorm(dr_beam_ref) * beam_segmentation_factor;

      // Get the surface normal vector.
      GEOMETRYPAIR::EvaluateSurfaceNormal<surface>(xi, this->face_element_->GetFacePosition(),
          surface_normal, this->face_element_->GetDrtFaceElement(),
          this->face_element_->GetCurrentNormals());

      // Evaluate the current position of beam and solid.
      GEOMETRYPAIR::EvaluatePosition<beam>(eta, this->ele1pos_, r_beam, this->Element1());
      GEOMETRYPAIR::EvaluatePosition<surface>(xi, this->face_element_->GetFacePosition(), r_surface,
          this->face_element_->GetDrtFaceElement());

      // Evaluate the gap function.
      r_rel = r_beam;
      r_rel -= r_surface;
      gap = r_rel.Dot(surface_normal) - beam_cross_section_radius;

      // Get the contact force.
      potential += projected_gauss_point.GetGaussWeight() * segment_jacobian *
                   PenaltyPotential(gap, this->Params()->BeamToSolidSurfaceContactParams());
    }
  }

  // If given, assemble force terms into the global vector.
  if (force_vector != Teuchos::null)
  {
    std::vector<double> force_pair_double(pair_gid.size(), 0.0);
    for (unsigned int j_dof = 0; j_dof < pair_gid.size(); j_dof++)
      force_pair_double[j_dof] = FADUTILS::CastToDouble(potential.dx(j_dof));
    force_vector->SumIntoGlobalValues(pair_gid.size(), pair_gid.data(), force_pair_double.data());
  }

  // If given, assemble force terms into the global stiffness matrix.
  if (stiffness_matrix != Teuchos::null)
    for (unsigned int i_dof = 0; i_dof < pair_gid.size(); i_dof++)
      for (unsigned int j_dof = 0; j_dof < pair_gid.size(); j_dof++)
        stiffness_matrix->FEAssemble(FADUTILS::CastToDouble(potential.dx(i_dof).dx(j_dof)),
            pair_gid[i_dof], pair_gid[j_dof]);
}


/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  template class BeamToSolidSurfaceContactPairGapVariation<
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_tri3>;
  template class BeamToSolidSurfaceContactPairGapVariation<
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_tri6>;
  template class BeamToSolidSurfaceContactPairGapVariation<
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad4>;
  template class BeamToSolidSurfaceContactPairGapVariation<
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad8>;
  template class BeamToSolidSurfaceContactPairGapVariation<
      line_to_surface_patch_scalar_type_1st_order, t_hermite, t_quad9>;
  template class BeamToSolidSurfaceContactPairGapVariation<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_hermite, t_nurbs9>, t_hermite,
      t_nurbs9>;

  template class BeamToSolidSurfaceContactPairPotential<line_to_surface_patch_scalar_type,
      t_hermite, t_tri3>;
  template class BeamToSolidSurfaceContactPairPotential<line_to_surface_patch_scalar_type,
      t_hermite, t_tri6>;
  template class BeamToSolidSurfaceContactPairPotential<line_to_surface_patch_scalar_type,
      t_hermite, t_quad4>;
  template class BeamToSolidSurfaceContactPairPotential<line_to_surface_patch_scalar_type,
      t_hermite, t_quad8>;
  template class BeamToSolidSurfaceContactPairPotential<line_to_surface_patch_scalar_type,
      t_hermite, t_quad9>;
  template class BeamToSolidSurfaceContactPairPotential<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;

}  // namespace BEAMINTERACTION
