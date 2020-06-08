/*----------------------------------------------------------------------*/
/*! \file

\brief Gauss point to segment mesh tying element for between a 3D beam and a surface element.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_surface_meshtying_pair_gauss_point.H"

#include "beam_contact_params.H"
#include "beam_to_solid_surface_meshtying_params.H"
#include "beaminteraction_calc_utils.H"
#include "../drt_geometry_pair/geometry_pair_line_to_surface.H"
#include "../drt_geometry_pair/geometry_pair_element_functions.H"
#include "../drt_geometry_pair/geometry_pair_factory.H"
#include "../drt_geometry_pair/geometry_pair_element_faces.H"

#include "Epetra_FEVector.h"


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
    const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix)
{
  // Call Evaluate on the geometry Pair. Only do this once for mesh tying.
  if (!this->meshtying_is_evaluated_)
  {
    this->CastGeometryPair()->Evaluate(this->ele1posref_,
        this->face_element_->GetFaceReferencePosition(), this->line_to_3D_segments_,
        this->face_element_->GetReferenceNormals());
    this->meshtying_is_evaluated_ = true;
  }

  // If there are no intersection segments, no coupling terms will be assembled.
  if (this->line_to_3D_segments_.size() == 0) return;

  // Initialize variables for position and force vectors.
  LINALG::Matrix<3, 1, double> dr_beam_ref;
  LINALG::Matrix<3, 1, scalar_type> coupling_vector_beam;
  LINALG::Matrix<3, 1, scalar_type> coupling_vector_surface;
  LINALG::Matrix<3, 1, scalar_type> force;
  LINALG::Matrix<beam::n_dof_ + surface::n_dof_, 1, scalar_type> force_pair(true);

  // Set the DOF vectors, depending on the desired coupling type.
  LINALG::Matrix<beam::n_dof_, 1, scalar_type> beam_dof_fad;
  LINALG::Matrix<surface::n_dof_, 1, scalar_type> surface_dof_fad;
  const INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling coupling_type =
      this->Params()->BeamToSolidSurfaceMeshtyingParams()->GetCouplingType();
  if (coupling_type ==
      INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::reference_configuration_forced_to_zero)
  {
    // Couple the positions -> this will result in an initial stress of the system.
    beam_dof_fad = this->ele1pos_;
    surface_dof_fad = this->face_element_->GetFacePosition();
  }
  else if (coupling_type == INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::displacement)
  {
    // Couple the displacements -> this will result in a non-fulfilment of the conservation of
    // angular momentum.
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

  // Initialize scalar variables.
  double segment_jacobian, beam_segmentation_factor;
  double penalty_parameter =
      this->Params()->BeamToSolidSurfaceMeshtyingParams()->GetPenaltyParameter();

  // Calculate the mesh tying forces.
  // Loop over segments.
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

      // Evaluate the positions / displacements depending on the coupling type. In both cases, the
      // terms are exclusively formulated with 0 normal direction, i.e. directly on the surface.
      GEOMETRYPAIR::EvaluatePosition<beam>(
          projected_gauss_point.GetEta(), beam_dof_fad, coupling_vector_beam, this->Element1());
      GEOMETRYPAIR::EvaluatePosition<surface>(projected_gauss_point.GetXi(), surface_dof_fad,
          coupling_vector_surface, this->face_element_->GetDrtFaceElement());

      // Calculate the force in this Gauss point. The sign of the force calculated here is the one
      // that acts on the beam.
      force = coupling_vector_surface;
      force -= coupling_vector_beam;
      force.Scale(penalty_parameter);

      // The force vector is in R3, we need to calculate the equivalent nodal forces on the element
      // dof. This is done with the virtual work equation $F \delta r = f \delta q$.
      for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
        for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
          force_pair(i_dof) -= force(i_dir) * coupling_vector_beam(i_dir).dx(i_dof) *
                               projected_gauss_point.GetGaussWeight() * segment_jacobian;
      for (unsigned int i_dof = 0; i_dof < surface::n_dof_; i_dof++)
        for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
          force_pair(i_dof + beam::n_dof_) +=
              force(i_dir) * coupling_vector_surface(i_dir).dx(i_dof + beam::n_dof_) *
              projected_gauss_point.GetGaussWeight() * segment_jacobian;
    }
  }

  // Get the pair GIDs.
  LINALG::Matrix<beam::n_dof_ + surface::n_dof_, 1, int> pair_gid;
  {
    // Get the beam centerline GIDs.
    LINALG::Matrix<beam::n_dof_, 1, int> beam_centerline_gid;
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
    const auto force_pair_double = FADUTILS::CastToDouble(force_pair);
    force_vector->SumIntoGlobalValues(
        beam::n_dof_ + surface::n_dof_, pair_gid.A(), force_pair_double.A());
  }

  // If given, assemble force terms into the global stiffness matrix.
  if (stiffness_matrix != Teuchos::null)
    for (unsigned int i_dof = 0; i_dof < beam::n_dof_ + surface::n_dof_; i_dof++)
      for (unsigned int j_dof = 0; j_dof < beam::n_dof_ + surface::n_dof_; j_dof++)
        stiffness_matrix->FEAssemble(
            FADUTILS::CastToDouble(force_pair(i_dof).dx(j_dof)), pair_gid(i_dof), pair_gid(j_dof));
}


/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  template class BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, GEOMETRYPAIR::t_tri3>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, GEOMETRYPAIR::t_tri6>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, GEOMETRYPAIR::t_quad4>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, GEOMETRYPAIR::t_quad8>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, GEOMETRYPAIR::t_quad9>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPoint<t_hermite, GEOMETRYPAIR::t_nurbs9>;
}  // namespace BEAMINTERACTION
