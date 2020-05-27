/*----------------------------------------------------------------------*/
/*! \file

\brief Gauss point to segment mesh tying element for between a 3D beam and a surface element. The
coupling terms are evaluated using FAD.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_surface_meshtying_pair_gauss_point_FAD.H"

#include "beam_contact_params.H"
#include "beam_to_solid_surface_meshtying_params.H"
#include "beaminteraction_calc_utils.H"
#include "../drt_geometry_pair/geometry_pair_line_to_surface.H"
#include "../drt_geometry_pair/geometry_pair_element_functions.H"
#include "../drt_geometry_pair/geometry_pair_factory.H"
#include "../drt_geometry_pair/geometry_pair_element_faces.H"
#include "../drt_geometry_pair/geometry_pair_scalar_types.H"

#include "Epetra_FEVector.h"


/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointFAD<scalar_type, beam,
    surface>::BeamToSolidSurfaceMeshtyingPairGaussPointFAD()
    : base_class()
{
  // Empty constructor.
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPointFAD<scalar_type, beam,
    surface>::EvaluateAndAssemble(const Teuchos::RCP<const DRT::Discretization>& discret,
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

  // Set the DOF vectors, depending on the desired coupling type.
  LINALG::Matrix<beam::n_dof_, 1, scalar_type> beam_dof_fad;
  LINALG::Matrix<surface::n_dof_, 1, scalar_type> surface_dof_fad;
  const INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling coupling_type =
      this->Params()->BeamToSolidSurfaceMeshtyingParams()->GetCouplingType();
  if (coupling_type == INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::
                           reference_configuration_forced_to_zero_fad or
      coupling_type == INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::consistent_fad)
  {
    // Couple the positions -> this will result in an initial stress of the system.
    beam_dof_fad = this->ele1pos_;
    surface_dof_fad = this->face_element_->GetFacePosition();
  }
  else if (coupling_type == INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::displacement_fad)
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
      if (coupling_type == INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::
                               reference_configuration_forced_to_zero_fad or
          coupling_type == INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::displacement_fad)
      {
        GEOMETRYPAIR::EvaluatePosition<surface>(projected_gauss_point.GetXi(), surface_dof_fad,
            coupling_vector_surface, this->face_element_->GetDrtFaceElement());
      }
      else if (coupling_type == INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::consistent_fad)
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

  // Get the pair GIDs.
  std::vector<int> pair_gid;
  {
    // Get the beam centerline GIDs.
    LINALG::Matrix<beam::n_dof_, 1, int> beam_centerline_gid;
    UTILS::GetElementCenterlineGIDIndices(*discret, this->Element1(), beam_centerline_gid);

    // Get the patch (in this case just the one face element) GIDs.
    const std::vector<int>& patch_gid = this->face_element_->GetPatchGID();
    pair_gid.resize(beam::n_dof_ + patch_gid.size());

    // Combine beam and solid GIDs into one vector.
    for (unsigned int i_dof_beam = 0; i_dof_beam < beam::n_dof_; i_dof_beam++)
      pair_gid[i_dof_beam] = beam_centerline_gid(i_dof_beam);
    for (unsigned int i_dof_patch = 0; i_dof_patch < patch_gid.size(); i_dof_patch++)
      pair_gid[beam::n_dof_ + i_dof_patch] = patch_gid[i_dof_patch];
  }

  // If given, assemble force terms into the global vector.
  if (force_vector != Teuchos::null)
  {
    std::vector<double> force_pair_double(pair_gid.size());
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

  template class BeamToSolidSurfaceMeshtyingPairGaussPointFAD<line_to_surface_patch_scalar_type,
      t_hermite, t_quad4>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointFAD<line_to_surface_patch_scalar_type,
      t_hermite, t_quad8>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointFAD<line_to_surface_patch_scalar_type,
      t_hermite, t_quad9>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointFAD<line_to_surface_patch_scalar_type,
      t_hermite, t_tri3>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointFAD<line_to_surface_patch_scalar_type,
      t_hermite, t_tri6>;
  template class BeamToSolidSurfaceMeshtyingPairGaussPointFAD<
      line_to_surface_patch_nurbs_scalar_type<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;
}  // namespace BEAMINTERACTION
