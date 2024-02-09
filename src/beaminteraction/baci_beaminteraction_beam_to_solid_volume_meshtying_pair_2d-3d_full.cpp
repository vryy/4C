/*----------------------------------------------------------------------*/
/*! \file

\brief Class for full 2D-3D beam-to-solid volume mesh tying based on a Simo-Reissner beam element.

\level 3
*/


#include "baci_beaminteraction_beam_to_solid_volume_meshtying_pair_2d-3d_full.hpp"

#include "baci_beam3_reissner.hpp"
#include "baci_beaminteraction_beam_to_solid_utils.hpp"
#include "baci_beaminteraction_beam_to_solid_volume_meshtying_params.hpp"
#include "baci_beaminteraction_calc_utils.hpp"
#include "baci_beaminteraction_contact_params.hpp"
#include "baci_geometry_pair_element_functions.hpp"
#include "baci_geometry_pair_line_to_3D_evaluation_data.hpp"
#include "baci_geometry_pair_line_to_volume_gauss_point_projection_cross_section.hpp"
#include "baci_geometry_pair_utility_classes.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"
#include "baci_linalg_utils_densematrix_inverse.hpp"

#include <Epetra_FEVector.h>

BACI_NAMESPACE_OPEN


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair2D3DFull<beam, solid>::PreEvaluate()
{
  // Call PreEvaluate on the geometry Pair.
  if (!this->meshtying_is_evaluated_)
  {
    this->CastGeometryPair()->PreEvaluate(this->ele1posref_, this->ele2posref_,
        this->line_to_3D_segments_, &triad_interpolation_scheme_ref_);
  }
}

/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair2D3DFull<beam, solid>::EvaluateAndAssemble(
    const Teuchos::RCP<const DRT::Discretization>& discret,
    const Teuchos::RCP<Epetra_FEVector>& force_vector,
    const Teuchos::RCP<CORE::LINALG::SparseMatrix>& stiffness_matrix,
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector)
{
  // Call Evaluate on the geometry Pair. Only do this once for mesh tying.
  if (!this->meshtying_is_evaluated_)
  {
    this->CastGeometryPair()->Evaluate(
        this->ele1posref_, this->ele2posref_, this->line_to_3D_segments_);
    this->meshtying_is_evaluated_ = true;
  }

  // If there are no segments, this pair has no contribution. Also there can be no more than one
  // segment.
  if (this->line_to_3D_segments_.size() == 0)
    return;
  else if (this->line_to_3D_segments_.size() > 1)
    dserror("There can be a maximum of one segment!");

  // Check that the beam element is a Simo--Reissner beam.
  auto beam_ele = dynamic_cast<const DRT::ELEMENTS::Beam3r*>(this->Element1());
  if (beam_ele == nullptr)
    dserror("GetBeamTriadInterpolationScheme is only implemented for SR beams.");

  // Get the vector with the projection points for this pair.
  const std::vector<GEOMETRYPAIR::ProjectionPoint1DTo3D<double>>& projection_points =
      this->line_to_3D_segments_[0].GetProjectionPoints();

  // If there are no projection points, return no contact status.
  if (projection_points.size() == 0) return;

  // Set the FAD variables for the beam and solid DOFs.
  auto set_q_fad = [](const auto& q_original, auto& q_fad, unsigned int fad_offset = 0)
  {
    for (unsigned int i_dof = 0; i_dof < q_original.numRows(); i_dof++)
      q_fad(i_dof) = CORE::FADUTILS::HigherOrderFadValue<scalar_type_pair>::apply(
          n_dof_fad_, fad_offset + i_dof, CORE::FADUTILS::CastToDouble(q_original(i_dof)));
  };
  CORE::LINALG::Matrix<beam::n_dof_, 1, scalar_type_pair> q_beam;
  CORE::LINALG::Matrix<solid::n_dof_, 1, scalar_type_pair> q_solid;
  set_q_fad(this->ele1pos_, q_beam);
  set_q_fad(this->ele2pos_, q_solid, beam::n_dof_);

  // Initialize pair wise vectors and matrices.
  CORE::LINALG::Matrix<n_dof_pair_, 1, double> force_pair(true);
  CORE::LINALG::Matrix<n_dof_pair_, n_dof_pair_, double> stiff_pair(true);
  CORE::LINALG::Matrix<n_dof_pair_, 1, scalar_type_pair> force_pair_local;
  CORE::LINALG::Matrix<n_dof_pair_, 3, double> d_force_d_psi;
  CORE::LINALG::Matrix<n_dof_pair_, n_dof_rot_, double> local_stiffness_rot;

  // Shape function matrices.
  CORE::LINALG::Matrix<3, solid::n_dof_, scalar_type_pair> N;
  CORE::LINALG::Matrix<3, beam::n_dof_, scalar_type_pair> H;
  CORE::LINALG::Matrix<3, n_dof_rot_, scalar_type_pair> L;
  std::vector<CORE::LINALG::Matrix<3, 3, double>> I_tilde_vector;
  CORE::LINALG::Matrix<3, n_dof_rot_, double> I_tilde;

  // Initialize vector and matrix variables for the Gauss integration.
  CORE::LINALG::Matrix<3, 1, double> dr_beam_ref;
  CORE::LINALG::Matrix<3, 1, scalar_type_pair> cross_section_vector_ref;
  CORE::LINALG::Matrix<3, 1, scalar_type_pair> cross_section_vector_current;
  CORE::LINALG::Matrix<3, 1, scalar_type_pair> pos_beam;
  CORE::LINALG::Matrix<3, 1, scalar_type_pair> pos_solid;
  CORE::LINALG::Matrix<4, 1, double> quaternion_double;
  CORE::LINALG::Matrix<3, 1, double> rotation_vector_double;
  CORE::LINALG::Matrix<4, 1, scalar_type_pair> quaternion_fad;
  CORE::LINALG::Matrix<3, 1, scalar_type_pair> rotation_vector_fad;
  CORE::LINALG::Matrix<3, 3, scalar_type_pair> triad_fad;
  CORE::LINALG::Matrix<3, 3, double> T_beam_double;
  CORE::LINALG::Matrix<3, n_dof_rot_, double> T_times_I_tilde;
  CORE::LINALG::Matrix<solid::n_dof_, 1, scalar_type_pair> temp_solid_force;
  CORE::LINALG::Matrix<beam::n_dof_, 1, scalar_type_pair> temp_beam_force;
  CORE::LINALG::Matrix<n_dof_rot_, 1, scalar_type_pair> temp_beam_force_rot;

  // Initialize scalar variables.
  double eta = 1e10;
  double eta_last_gauss_point = 1e10;
  double beam_jacobian = 0.0;
  const double penalty_parameter =
      this->Params()->BeamToSolidVolumeMeshtyingParams()->GetPenaltyParameter();

  // Calculate the mesh tying forces.
  // Loop over segments.
  for (unsigned int i_integration_point = 0; i_integration_point < projection_points.size();
       i_integration_point++)
  {
    // Get the current Gauss point.
    const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& projected_gauss_point =
        projection_points[i_integration_point];
    eta = projected_gauss_point.GetEta();

    // Evaluate all beam specific terms. This only has to be done if the Gauss point position on the
    // beam has changed compared to the last Gauss point.
    if (std::abs(eta - eta_last_gauss_point) > 1e-10)
    {
      GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(
          eta, this->ele1posref_, dr_beam_ref, this->Element1());
      beam_jacobian = 0.5 * dr_beam_ref.Norm2();

      GEOMETRYPAIR::EvaluateShapeFunctionMatrix<beam>(eta, H, this->Element1());
      pos_beam.Multiply(H, q_beam);

      GEOMETRYPAIR::EvaluateShapeFunctionMatrix<GEOMETRYPAIR::t_line3>(eta, L);

      triad_interpolation_scheme_.GetNodalGeneralizedRotationInterpolationMatricesAtXi(
          I_tilde_vector, eta);
      for (unsigned int i_node = 0; i_node < n_nodes_rot_; i_node++)
        for (unsigned int i_dim_0 = 0; i_dim_0 < 3; i_dim_0++)
          for (unsigned int i_dim_1 = 0; i_dim_1 < 3; i_dim_1++)
            I_tilde(i_dim_0, i_node * 3 + i_dim_1) = I_tilde_vector[i_node](i_dim_0, i_dim_1);

      // Get the rotation vector at this Gauss point.
      triad_interpolation_scheme_.GetInterpolatedQuaternionAtXi(quaternion_double, eta);
      CORE::LARGEROTATIONS::quaterniontoangle(quaternion_double, rotation_vector_double);
      T_beam_double = CORE::LARGEROTATIONS::Tmatrix(rotation_vector_double);
      set_q_fad(rotation_vector_double, rotation_vector_fad, beam::n_dof_ + solid::n_dof_);
      CORE::LARGEROTATIONS::angletoquaternion(rotation_vector_fad, quaternion_fad);
      CORE::LARGEROTATIONS::quaterniontotriad(quaternion_fad, triad_fad);
    }

    // Get the shape function matrices.
    GEOMETRYPAIR::EvaluateShapeFunctionMatrix<solid>(
        projected_gauss_point.GetXi(), N, this->Element2());
    pos_solid.Multiply(N, q_solid);

    // Get the cross section vector.
    cross_section_vector_ref(0) = 0.0;
    cross_section_vector_ref(1) = projected_gauss_point.GetEtaCrossSection()(0);
    cross_section_vector_ref(2) = projected_gauss_point.GetEtaCrossSection()(1);
    cross_section_vector_current.Multiply(triad_fad, cross_section_vector_ref);

    // Reset the local force vector for this Gauss point.
    force_pair_local.PutScalar(0.0);

    // Numerical integration factor for this Gauss point.
    const double integration_factor =
        projected_gauss_point.GetGaussWeight() * beam_jacobian * penalty_parameter;

    // The following calculations are based on Steinbrecher, Popp, Meier: "Consistent coupling of
    // positions and rotations for embedding 1D Cosserat beams into 3D solid volumes", eq. 97-98. Be
    // aware, that there is a typo in eq. 98 where the derivative is taken with respect to the
    // rotation angle and not the rotational DOFs.
    auto r_diff = pos_beam;
    r_diff += cross_section_vector_current;
    r_diff -= pos_solid;

    // Evaluate the force on the solid DOFs.
    temp_solid_force.MultiplyTN(N, r_diff);
    temp_solid_force.Scale(-1.0);
    for (unsigned int i_dof = 0; i_dof < solid::n_dof_; i_dof++)
      force_pair_local(i_dof + beam::n_dof_) += temp_solid_force(i_dof);

    // Evaluate the force on the positional beam DOFs.
    temp_beam_force.MultiplyTN(H, r_diff);
    for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
      force_pair_local(i_dof) += temp_beam_force(i_dof);

    // Evaluate the force on the rotational beam DOFs.
    // In comparison to the mentioned paper, the relative cross section vector is also contained
    // here, but it cancels out in the cross product with itself.
    CORE::LINALG::Matrix<3, 1, scalar_type_pair> temp_beam_rot_cross;
    temp_beam_rot_cross.CrossProduct(cross_section_vector_current, r_diff);
    temp_beam_force_rot.MultiplyTN(L, temp_beam_rot_cross);
    for (unsigned int i_dof = 0; i_dof < n_dof_rot_; i_dof++)
      force_pair_local(i_dof + beam::n_dof_ + solid::n_dof_) += temp_beam_force_rot(i_dof);

    // Add to pair force contributions.
    force_pair_local.Scale(integration_factor);
    force_pair += CORE::FADUTILS::CastToDouble(force_pair_local);

    // The rotational stiffness contributions have to be handled separately due to the non-addidive
    // nature of the rotational DOFs.
    for (unsigned int i_dof = 0; i_dof < n_dof_pair_; i_dof++)
      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        d_force_d_psi(i_dof, i_dir) =
            force_pair_local(i_dof).dx(beam::n_dof_ + solid::n_dof_ + i_dir);
    T_times_I_tilde.Multiply(T_beam_double, I_tilde);
    local_stiffness_rot.MultiplyNN(d_force_d_psi, T_times_I_tilde);

    // Add the full stiffness contribution from this Gauss point.
    for (unsigned int i_dof = 0; i_dof < n_dof_pair_; i_dof++)
    {
      for (unsigned int j_dof = 0; j_dof < n_dof_pair_; j_dof++)
      {
        if (j_dof < beam::n_dof_ + solid::n_dof_)
          stiff_pair(i_dof, j_dof) += force_pair_local(i_dof).dx(j_dof);
        else
          stiff_pair(i_dof, j_dof) +=
              local_stiffness_rot(i_dof, j_dof - beam::n_dof_ - solid::n_dof_);
      }
    }

    // Set the eta value for this Gauss point.
    eta_last_gauss_point = eta;
  }

  // Get the GIDs of this pair.
  CORE::LINALG::Matrix<n_dof_pair_, 1, int> gid_pair;

  // Beam centerline GIDs.
  CORE::LINALG::Matrix<beam::n_dof_, 1, int> beam_centerline_gid;
  UTILS::GetElementCenterlineGIDIndices(*discret, this->Element1(), beam_centerline_gid);
  for (unsigned int i_dof_beam = 0; i_dof_beam < beam::n_dof_; i_dof_beam++)
    gid_pair(i_dof_beam) = beam_centerline_gid(i_dof_beam);

  // Solid GIDs.
  std::vector<int> lm, lmowner, lmstride;
  this->Element2()->LocationVector(*discret, lm, lmowner, lmstride);
  for (unsigned int i = 0; i < solid::n_dof_; i++) gid_pair(i + beam::n_dof_) = lm[i];

  // Beam rot GIDs.
  this->Element1()->LocationVector(*discret, lm, lmowner, lmstride);
  std::array<int, 9> rot_dof_indices = {3, 4, 5, 12, 13, 14, 18, 19, 20};
  for (unsigned int i = 0; i < n_dof_rot_; i++)
    gid_pair(i + beam::n_dof_ + solid::n_dof_) = lm[rot_dof_indices[i]];

  // If given, assemble force terms into the global force vector.
  if (force_vector != Teuchos::null)
    force_vector->SumIntoGlobalValues(gid_pair.numRows(), gid_pair.A(), force_pair.A());

  // If given, assemble force terms into the global stiffness matrix.
  if (stiffness_matrix != Teuchos::null)
    for (unsigned int i_dof = 0; i_dof < n_dof_pair_; i_dof++)
      for (unsigned int j_dof = 0; j_dof < n_dof_pair_; j_dof++)
        stiffness_matrix->FEAssemble(stiff_pair(i_dof, j_dof), gid_pair(i_dof), gid_pair(j_dof));
}

/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair2D3DFull<beam, solid>::ResetRotationState(
    const DRT::Discretization& discret, const Teuchos::RCP<const Epetra_Vector>& ia_discolnp)
{
  GetBeamTriadInterpolationScheme(discret, ia_discolnp, this->Element1(),
      triad_interpolation_scheme_, this->triad_interpolation_scheme_ref_);
}

/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair2D3DFull<beam, solid>::GetTriadAtXiDouble(
    const double xi, CORE::LINALG::Matrix<3, 3, double>& triad, const bool reference) const
{
  if (reference)
  {
    this->triad_interpolation_scheme_ref_.GetInterpolatedTriadAtXi(triad, xi);
  }
  else
  {
    this->triad_interpolation_scheme_.GetInterpolatedTriadAtXi(triad, xi);
  }
}


/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  template class BeamToSolidVolumeMeshtyingPair2D3DFull<t_hermite, t_hex8>;
  template class BeamToSolidVolumeMeshtyingPair2D3DFull<t_hermite, t_hex20>;
  template class BeamToSolidVolumeMeshtyingPair2D3DFull<t_hermite, t_hex27>;
  template class BeamToSolidVolumeMeshtyingPair2D3DFull<t_hermite, t_tet4>;
  template class BeamToSolidVolumeMeshtyingPair2D3DFull<t_hermite, t_tet10>;
}  // namespace BEAMINTERACTION

BACI_NAMESPACE_CLOSE
