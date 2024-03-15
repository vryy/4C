/*----------------------------------------------------------------------*/
/*! \file

\brief Gauss point to segment mesh tying element for between a 3D beam and a surface element.

\level 3
*/


#include "baci_beaminteraction_beam_to_solid_surface_meshtying_pair_gauss_point.hpp"

#include "baci_beaminteraction_beam_to_solid_surface_meshtying_params.hpp"
#include "baci_beaminteraction_calc_utils.hpp"
#include "baci_beaminteraction_contact_params.hpp"
#include "baci_geometry_pair_element_evaluation_functions.hpp"
#include "baci_geometry_pair_element_faces.hpp"
#include "baci_geometry_pair_factory.hpp"
#include "baci_geometry_pair_line_to_surface.hpp"

#include <Epetra_FEVector.h>

BACI_NAMESPACE_OPEN


/**
 *
 */
template <typename beam, typename surface>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPoint<beam,
    surface>::BeamToSolidSurfaceMeshtyingPairGaussPoint()
    : base_class()
{
  // Empty constructor.
}

/**
 *
 */
template <typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPoint<beam, surface>::EvaluateAndAssemble(
    const Teuchos::RCP<const DRT::Discretization>& discret,
    const Teuchos::RCP<Epetra_FEVector>& force_vector,
    const Teuchos::RCP<CORE::LINALG::SparseMatrix>& stiffness_matrix,
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector)
{
  // Call Evaluate on the geometry Pair. Only do this once for mesh tying.
  if (!this->meshtying_is_evaluated_)
  {
    this->CastGeometryPair()->Evaluate(this->ele1posref_,
        this->face_element_->GetFaceReferenceElementData(), this->line_to_3D_segments_);
    this->meshtying_is_evaluated_ = true;
  }

  // If there are no intersection segments, no coupling terms will be assembled.
  if (this->line_to_3D_segments_.size() == 0) return;

  // Initialize variables for position and force vectors.
  CORE::LINALG::Matrix<3, 1, double> dr_beam_ref;
  CORE::LINALG::Matrix<3, 1, scalar_type> coupling_vector;
  CORE::LINALG::Matrix<3, 1, scalar_type> force;
  CORE::LINALG::Matrix<beam::n_dof_ + surface::n_dof_, 1, scalar_type> force_pair(true);

  // Initialize scalar variables.
  double segment_jacobian = 0.0;
  double beam_segmentation_factor = 0.0;
  double penalty_parameter =
      this->Params()->BeamToSolidSurfaceMeshtyingParams()->GetPenaltyParameter();

  // Calculate the mesh tying forces.
  // Loop over segments.
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

      // Calculate the force in this Gauss point. The sign of the force calculated here is the one
      // that acts on the beam.
      coupling_vector = this->EvaluateCoupling(projected_gauss_point);
      force = coupling_vector;
      force.Scale(penalty_parameter);

      // The force vector is in R3, we need to calculate the equivalent nodal forces on the element
      // dof. This is done with the virtual work equation $F \delta r = f \delta q$.
      for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
        for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
          force_pair(i_dof) += force(i_dir) * coupling_vector(i_dir).dx(i_dof) *
                               projected_gauss_point.GetGaussWeight() * segment_jacobian;
      for (unsigned int i_dof = 0; i_dof < surface::n_dof_; i_dof++)
        for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
          force_pair(i_dof + beam::n_dof_) +=
              force(i_dir) * coupling_vector(i_dir).dx(i_dof + beam::n_dof_) *
              projected_gauss_point.GetGaussWeight() * segment_jacobian;
    }
  }

  // Get the pair GIDs.
  CORE::LINALG::Matrix<beam::n_dof_ + surface::n_dof_, 1, int> pair_gid;
  {
    // Get the beam centerline GIDs.
    CORE::LINALG::Matrix<beam::n_dof_, 1, int> beam_centerline_gid;
    UTILS::GetElementCenterlineGIDIndices(*discret, this->Element1(), beam_centerline_gid);

    // Get the patch (in this case just the one face element) GIDs.
    const std::vector<int>& patch_gid = this->face_element_->GetPatchGID();

    // Combine beam and solid GIDs into one vector.
    for (unsigned int i_dof_beam = 0; i_dof_beam < beam::n_dof_; i_dof_beam++)
      pair_gid(i_dof_beam) = beam_centerline_gid(i_dof_beam);
    for (unsigned int i_dof_solid = 0; i_dof_solid < surface::n_dof_; i_dof_solid++)
      pair_gid(beam::n_dof_ + i_dof_solid) = patch_gid[i_dof_solid];
  }

  // If given, assemble force terms into the global vector.
  if (force_vector != Teuchos::null)
  {
    const auto force_pair_double = CORE::FADUTILS::CastToDouble(force_pair);
    force_vector->SumIntoGlobalValues(
        beam::n_dof_ + surface::n_dof_, pair_gid.A(), force_pair_double.A());
  }

  // If given, assemble force terms into the global stiffness matrix.
  if (stiffness_matrix != Teuchos::null)
    for (unsigned int i_dof = 0; i_dof < beam::n_dof_ + surface::n_dof_; i_dof++)
      for (unsigned int j_dof = 0; j_dof < beam::n_dof_ + surface::n_dof_; j_dof++)
        stiffness_matrix->FEAssemble(CORE::FADUTILS::CastToDouble(force_pair(i_dof).dx(j_dof)),
            pair_gid(i_dof), pair_gid(j_dof));
}


/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  template class BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, t_tri3>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, t_tri6>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, t_quad4>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, t_quad8>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, t_quad9>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, t_nurbs9>;
}  // namespace BEAMINTERACTION

BACI_NAMESPACE_CLOSE
