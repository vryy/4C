/*----------------------------------------------------------------------*/
/*! \file

\brief Meshtying element for rotational meshtying between a 3D beam and a 3D solid element.

\level 3
*/


#include "4C_beaminteraction_beam_to_solid_volume_meshtying_pair_mortar_rotation.hpp"

#include "4C_beam3_reissner.hpp"
#include "4C_beam3_triad_interpolation_local_rotation_vectors.hpp"
#include "4C_beaminteraction_beam_to_solid_mortar_manager.hpp"
#include "4C_beaminteraction_beam_to_solid_utils.hpp"
#include "4C_beaminteraction_beam_to_solid_volume_meshtying_params.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_contact_params.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_line_to_volume.hpp"

#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN



/**
 *
 */
template <typename beam, typename solid, typename mortar, typename mortar_rot>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortarRotation<beam, solid, mortar,
    mortar_rot>::BeamToSolidVolumeMeshtyingPairMortarRotation()
    : BeamToSolidVolumeMeshtyingPairMortar<beam, solid, mortar>()
{
  // Set the number of rotational mortar DOFs.
  this->n_mortar_rot_ = mortar_rot::n_dof_;
}

/**
 *
 */
template <typename beam, typename solid, typename mortar, typename mortar_rot>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortarRotation<beam, solid, mortar,
    mortar_rot>::evaluate_and_assemble_mortar_contributions(const Core::FE::Discretization& discret,
    const BeamToSolidMortarManager* mortar_manager,
    Core::LinAlg::SparseMatrix& global_constraint_lin_beam,
    Core::LinAlg::SparseMatrix& global_constraint_lin_solid,
    Core::LinAlg::SparseMatrix& global_force_beam_lin_lambda,
    Core::LinAlg::SparseMatrix& global_force_solid_lin_lambda, Epetra_FEVector& global_constraint,
    Epetra_FEVector& global_kappa, Core::LinAlg::SparseMatrix& global_kappa_lin_beam,
    Core::LinAlg::SparseMatrix& global_kappa_lin_solid, Epetra_FEVector& global_lambda_active,
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector)
{
  // Call the base method.
  base_class::evaluate_and_assemble_mortar_contributions(discret, mortar_manager,
      global_constraint_lin_beam, global_constraint_lin_solid, global_force_beam_lin_lambda,
      global_force_solid_lin_lambda, global_constraint, global_kappa, global_kappa_lin_beam,
      global_kappa_lin_solid, global_lambda_active, displacement_vector);

  // If there are no intersection segments, return as no contact can occur.
  if (this->line_to_3D_segments_.size() == 0) return;

  // Get the beam triad interpolation schemes.
  LargeRotations::TriadInterpolationLocalRotationVectors<3, double> triad_interpolation_scheme;
  LargeRotations::TriadInterpolationLocalRotationVectors<3, double> ref_triad_interpolation_scheme;
  GetBeamTriadInterpolationScheme(discret, displacement_vector, this->Element1(),
      triad_interpolation_scheme, ref_triad_interpolation_scheme);

  // Set the FAD variables for the solid DOFs. For the terms calculated here we only need first
  // order derivatives.
  auto q_solid =
      GEOMETRYPAIR::InitializeElementData<solid, scalar_type_rot_1st>::Initialize(this->Element2());
  for (unsigned int i_solid = 0; i_solid < solid::n_dof_; i_solid++)
    q_solid.element_position_(i_solid) =
        Core::FADUtils::HigherOrderFadValue<scalar_type_rot_1st>::apply(3 + solid::n_dof_,
            3 + i_solid, Core::FADUtils::CastToDouble(this->ele2pos_.element_position_(i_solid)));

  // Initialize local matrices.
  Core::LinAlg::Matrix<mortar_rot::n_dof_, 1, double> local_g(true);
  Core::LinAlg::Matrix<mortar_rot::n_dof_, n_dof_rot_, double> local_G_B(true);
  Core::LinAlg::Matrix<mortar_rot::n_dof_, solid::n_dof_, double> local_G_S(true);
  Core::LinAlg::Matrix<n_dof_rot_, mortar_rot::n_dof_, double> local_FB_L(true);
  Core::LinAlg::Matrix<solid::n_dof_, mortar_rot::n_dof_, double> local_FS_L(true);
  Core::LinAlg::Matrix<mortar_rot::n_dof_, 1, double> local_kappa(true);

  const auto rot_coupling_type =
      this->Params()->beam_to_solid_volume_meshtying_params()->get_rotational_coupling_type();
  if (rot_coupling_type == Inpar::BeamToSolid::BeamToSolidRotationCoupling::fix_triad_2d)
  {
    // In the case of "fix_triad_2d" we couple both, the ey and ez direction to the beam. Therefore,
    // we have to evaluate the coupling terms w.r.t both of those coupling types.
    evaluate_rotational_coupling_terms(
        Inpar::BeamToSolid::BeamToSolidRotationCoupling::deformation_gradient_y_2d, q_solid,
        triad_interpolation_scheme, ref_triad_interpolation_scheme, local_g, local_G_B, local_G_S,
        local_FB_L, local_FS_L, local_kappa);
    evaluate_rotational_coupling_terms(
        Inpar::BeamToSolid::BeamToSolidRotationCoupling::deformation_gradient_z_2d, q_solid,
        triad_interpolation_scheme, ref_triad_interpolation_scheme, local_g, local_G_B, local_G_S,
        local_FB_L, local_FS_L, local_kappa);
  }
  else
    evaluate_rotational_coupling_terms(rot_coupling_type, q_solid, triad_interpolation_scheme,
        ref_triad_interpolation_scheme, local_g, local_G_B, local_G_S, local_FB_L, local_FS_L,
        local_kappa);

  // Get the GIDs of the solid and beam.
  std::vector<int> lm_beam, gid_solid, lmowner, lmstride;
  this->Element1()->LocationVector(discret, lm_beam, lmowner, lmstride);
  this->Element2()->LocationVector(discret, gid_solid, lmowner, lmstride);
  // Local indices of rotational DOFs for the Simo--Reissner beam element.
  const int rot_dof_indices[] = {3, 4, 5, 12, 13, 14, 18, 19, 20};
  Core::LinAlg::Matrix<n_dof_rot_, 1, int> gid_rot;
  for (unsigned int i = 0; i < n_dof_rot_; i++) gid_rot(i) = lm_beam[rot_dof_indices[i]];

  // Get the Lagrange multiplier GIDs.
  const auto& [_, lambda_gid_rot] = mortar_manager->LocationVector(*this);

  // Assemble into the global vectors
  global_constraint.SumIntoGlobalValues(lambda_gid_rot.size(), lambda_gid_rot.data(), local_g.A());
  global_kappa.SumIntoGlobalValues(lambda_gid_rot.size(), lambda_gid_rot.data(), local_kappa.A());
  local_kappa.PutScalar(1.0);
  global_lambda_active.SumIntoGlobalValues(
      lambda_gid_rot.size(), lambda_gid_rot.data(), local_kappa.A());

  // Assemble into global matrices.
  for (unsigned int i_dof_lambda = 0; i_dof_lambda < mortar_rot::n_dof_; i_dof_lambda++)
  {
    for (unsigned int i_dof_rot = 0; i_dof_rot < n_dof_rot_; i_dof_rot++)
    {
      global_constraint_lin_beam.FEAssemble(
          local_G_B(i_dof_lambda, i_dof_rot), lambda_gid_rot[i_dof_lambda], gid_rot(i_dof_rot));
      global_force_beam_lin_lambda.FEAssemble(
          local_FB_L(i_dof_rot, i_dof_lambda), gid_rot(i_dof_rot), lambda_gid_rot[i_dof_lambda]);
    }
    for (unsigned int i_dof_solid = 0; i_dof_solid < solid::n_dof_; i_dof_solid++)
    {
      global_constraint_lin_solid.FEAssemble(local_G_S(i_dof_lambda, i_dof_solid),
          lambda_gid_rot[i_dof_lambda], gid_solid[i_dof_solid]);
      global_force_solid_lin_lambda.FEAssemble(local_FS_L(i_dof_solid, i_dof_lambda),
          gid_solid[i_dof_solid], lambda_gid_rot[i_dof_lambda]);
    }
  }
}

/**
 *
 */
template <typename beam, typename solid, typename mortar, typename mortar_rot>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortarRotation<beam, solid, mortar,
    mortar_rot>::evaluate_rotational_coupling_terms(  //
    const Inpar::BeamToSolid::BeamToSolidRotationCoupling& rot_coupling_type,
    const GEOMETRYPAIR::ElementData<solid, scalar_type_rot_1st>& q_solid,
    const LargeRotations::TriadInterpolationLocalRotationVectors<3, double>&
        triad_interpolation_scheme,
    const LargeRotations::TriadInterpolationLocalRotationVectors<3, double>&
        ref_triad_interpolation_scheme,
    Core::LinAlg::Matrix<mortar_rot::n_dof_, 1, double>& local_g,
    Core::LinAlg::Matrix<mortar_rot::n_dof_, n_dof_rot_, double>& local_G_B,
    Core::LinAlg::Matrix<mortar_rot::n_dof_, solid::n_dof_, double>& local_G_S,
    Core::LinAlg::Matrix<n_dof_rot_, mortar_rot::n_dof_, double>& local_FB_L,
    Core::LinAlg::Matrix<solid::n_dof_, mortar_rot::n_dof_, double>& local_FS_L,
    Core::LinAlg::Matrix<mortar_rot::n_dof_, 1, double>& local_kappa) const
{
  // Initialize variables.
  Core::LinAlg::Matrix<3, 1, double> dr_beam_ref;
  Core::LinAlg::Matrix<4, 1, double> quaternion_beam_double;
  Core::LinAlg::Matrix<3, 1, double> psi_beam_double;
  Core::LinAlg::Matrix<3, 1, scalar_type_rot_1st> psi_beam;
  Core::LinAlg::Matrix<3, 1, scalar_type_rot_1st> psi_solid;
  Core::LinAlg::Matrix<3, 1, scalar_type_rot_1st> psi_rel;
  Core::LinAlg::Matrix<4, 1, scalar_type_rot_1st> quaternion_beam;
  Core::LinAlg::Matrix<4, 1, scalar_type_rot_1st> quaternion_beam_inv;
  Core::LinAlg::Matrix<4, 1, double> quaternion_beam_ref;
  Core::LinAlg::Matrix<4, 1, scalar_type_rot_1st> quaternion_solid;
  Core::LinAlg::Matrix<4, 1, scalar_type_rot_1st> quaternion_rel;
  Core::LinAlg::Matrix<3, 3, double> T_beam;
  Core::LinAlg::Matrix<3, 3, double> T_solid;
  Core::LinAlg::Matrix<3, 3, double> T_solid_inv;
  Core::LinAlg::Matrix<3, 3, double> T_rel;

  Core::LinAlg::Matrix<mortar_rot::n_nodes_, 1, double> lambda_shape_functions;
  Core::LinAlg::Matrix<3, mortar_rot::n_dof_, double> lambda_shape_functions_full(true);
  Core::LinAlg::SerialDenseVector L_i(3);
  Core::LinAlg::Matrix<3, n_dof_rot_, double> L_full(true);
  std::vector<Core::LinAlg::Matrix<3, 3, double>> I_beam_tilde;
  Core::LinAlg::Matrix<3, n_dof_rot_, double> I_beam_tilde_full;
  Core::LinAlg::Matrix<3, n_dof_rot_, double> T_beam_times_I_beam_tilde_full;
  Core::LinAlg::Matrix<3, mortar_rot::n_dof_, double> T_rel_tr_times_lambda_shape;
  Core::LinAlg::Matrix<3, mortar_rot::n_dof_, double> T_solid_mtr_times_T_rel_tr_times_lambda_shape;
  Core::LinAlg::Matrix<n_dof_rot_, mortar_rot::n_dof_, double> d_fb_d_lambda_gp;
  Core::LinAlg::Matrix<solid::n_dof_, mortar_rot::n_dof_, double> d_fs_d_lambda_gp;
  Core::LinAlg::Matrix<mortar_rot::n_dof_, 1, scalar_type_rot_1st> g_gp;
  Core::LinAlg::Matrix<3, solid::n_dof_, double> d_psi_solid_d_q_solid;
  Core::LinAlg::Matrix<mortar_rot::n_dof_, 3, double> d_g_d_psi_beam;
  Core::LinAlg::Matrix<mortar_rot::n_dof_, n_dof_rot_, double> d_g_d_psi_beam_times_T_beam_I;
  Core::LinAlg::Matrix<mortar_rot::n_dof_, solid::n_dof_, double> d_g_d_q_solid;

  // Initialize scalar variables.
  double segment_jacobian = 0.0;
  double beam_segmentation_factor = 0.0;

  // Calculate the meshtying forces.
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

      // Get the jacobian in the reference configuration.
      GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(
          projected_gauss_point.GetEta(), this->ele1posref_, dr_beam_ref);

      // Jacobian including the segment length.
      segment_jacobian = dr_beam_ref.Norm2() * beam_segmentation_factor;

      // Calculate the rotation vector of this cross section.
      triad_interpolation_scheme.get_interpolated_quaternion_at_xi(
          quaternion_beam_double, projected_gauss_point.GetEta());
      Core::LargeRotations::quaterniontoangle(quaternion_beam_double, psi_beam_double);
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        psi_beam(i_dim) = Core::FADUtils::HigherOrderFadValue<scalar_type_rot_1st>::apply(
            3 + solid::n_dof_, i_dim, psi_beam_double(i_dim));
      Core::LargeRotations::angletoquaternion(psi_beam, quaternion_beam);
      quaternion_beam_inv = Core::LargeRotations::inversequaternion(quaternion_beam);

      // Get the solid rotation vector.
      ref_triad_interpolation_scheme.get_interpolated_quaternion_at_xi(
          quaternion_beam_ref, projected_gauss_point.GetEta());
      GetSolidRotationVector<solid>(rot_coupling_type, projected_gauss_point.GetXi(),
          this->ele2posref_, q_solid, quaternion_beam_ref, psi_solid);
      Core::LargeRotations::angletoquaternion(psi_solid, quaternion_solid);

      // Calculate the relative rotation vector.
      Core::LargeRotations::quaternionproduct(
          quaternion_beam_inv, quaternion_solid, quaternion_rel);
      Core::LargeRotations::quaterniontoangle(quaternion_rel, psi_rel);

      // Calculate the transformation matrices.
      T_rel = Core::LargeRotations::Tmatrix(Core::FADUtils::CastToDouble(psi_rel));
      T_beam = Core::LargeRotations::Tmatrix(Core::FADUtils::CastToDouble(psi_beam));
      T_solid = Core::LargeRotations::Tmatrix(Core::FADUtils::CastToDouble(psi_solid));
      T_solid_inv = T_solid;
      Core::LinAlg::Inverse(T_solid_inv);

      // Evaluate shape functions.
      GEOMETRYPAIR::EvaluateShapeFunction<mortar_rot>::Evaluate(
          lambda_shape_functions, projected_gauss_point.GetEta());
      for (unsigned int i_node = 0; i_node < mortar_rot::n_nodes_; i_node++)
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          lambda_shape_functions_full(i_dim, 3 * i_node + i_dim) = lambda_shape_functions(i_node);

      Core::FE::shape_function_1D(L_i, projected_gauss_point.GetEta(), Core::FE::CellType::line3);
      for (unsigned int i_node = 0; i_node < 3; i_node++)
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          L_full(i_dim, 3 * i_node + i_dim) = L_i(i_node);

      triad_interpolation_scheme.get_nodal_generalized_rotation_interpolation_matrices_at_xi(
          I_beam_tilde, projected_gauss_point.GetEta());
      for (unsigned int i_node = 0; i_node < 3; i_node++)
        for (unsigned int i_dim_0 = 0; i_dim_0 < 3; i_dim_0++)
          for (unsigned int i_dim_1 = 0; i_dim_1 < 3; i_dim_1++)
            I_beam_tilde_full(i_dim_0, i_node * 3 + i_dim_1) =
                I_beam_tilde[i_node](i_dim_0, i_dim_1);

      // Solid angle derived w.r.t. the solid DOFs.
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        for (unsigned int i_solid = 0; i_solid < solid::n_dof_; i_solid++)
          d_psi_solid_d_q_solid(i_dim, i_solid) = psi_solid(i_dim).dx(3 + i_solid);

      // Calculate the force terms derived w.r.t. the Lagrange multipliers.
      T_rel_tr_times_lambda_shape.MultiplyTN(T_rel, lambda_shape_functions_full);
      d_fb_d_lambda_gp.MultiplyTN(L_full, T_rel_tr_times_lambda_shape);
      d_fb_d_lambda_gp.Scale(-1.0 * projected_gauss_point.GetGaussWeight() * segment_jacobian);

      T_solid_mtr_times_T_rel_tr_times_lambda_shape.MultiplyTN(
          T_solid_inv, T_rel_tr_times_lambda_shape);
      d_fs_d_lambda_gp.MultiplyTN(
          d_psi_solid_d_q_solid, T_solid_mtr_times_T_rel_tr_times_lambda_shape);
      d_fs_d_lambda_gp.Scale(projected_gauss_point.GetGaussWeight() * segment_jacobian);

      // Constraint vector.
      g_gp.PutScalar(0.0);
      for (unsigned int i_row = 0; i_row < mortar_rot::n_dof_; i_row++)
        for (unsigned int i_col = 0; i_col < 3; i_col++)
          g_gp(i_row) += lambda_shape_functions_full(i_col, i_row) * psi_rel(i_col);
      g_gp.Scale(projected_gauss_point.GetGaussWeight() * segment_jacobian);

      // Derivatives of constraint vector.
      T_beam_times_I_beam_tilde_full.Multiply(T_beam, I_beam_tilde_full);

      for (unsigned int i_lambda = 0; i_lambda < mortar_rot::n_dof_; i_lambda++)
        for (unsigned int i_psi = 0; i_psi < 3; i_psi++)
          d_g_d_psi_beam(i_lambda, i_psi) = g_gp(i_lambda).dx(i_psi);
      d_g_d_psi_beam_times_T_beam_I.Multiply(d_g_d_psi_beam, T_beam_times_I_beam_tilde_full);

      for (unsigned int i_lambda = 0; i_lambda < mortar_rot::n_dof_; i_lambda++)
        for (unsigned int i_solid = 0; i_solid < solid::n_dof_; i_solid++)
          d_g_d_q_solid(i_lambda, i_solid) = g_gp(i_lambda).dx(3 + i_solid);

      // Add to output matrices and vector.
      local_g += Core::FADUtils::CastToDouble(g_gp);
      local_G_B += d_g_d_psi_beam_times_T_beam_I;
      local_G_S += d_g_d_q_solid;
      local_FB_L += d_fb_d_lambda_gp;
      local_FS_L += d_fs_d_lambda_gp;

      // Calculate the scaling entries.
      for (unsigned int i_mortar_node = 0; i_mortar_node < mortar_rot::n_nodes_; i_mortar_node++)
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          local_kappa(i_mortar_node * 3 + i_dim) += lambda_shape_functions(i_mortar_node) *
                                                    projected_gauss_point.GetGaussWeight() *
                                                    segment_jacobian;
    }
  }
}

/**
 *
 */
template <typename beam, typename solid, typename mortar, typename mortar_rot>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortarRotation<beam, solid, mortar,
    mortar_rot>::EvaluateAndAssemble(const Core::FE::Discretization& discret,
    const BeamToSolidMortarManager* mortar_manager,
    const Teuchos::RCP<Epetra_FEVector>& force_vector,
    const Teuchos::RCP<Core::LinAlg::SparseMatrix>& stiffness_matrix,
    const Epetra_Vector& global_lambda, const Epetra_Vector& displacement_vector)
{
  // Call the base method.
  base_class::EvaluateAndAssemble(
      discret, mortar_manager, force_vector, stiffness_matrix, global_lambda, displacement_vector);

  // If there are no intersection segments, return as no contact can occur.
  if (this->line_to_3D_segments_.size() == 0) return;

  // This pair only gives contributions to the stiffness matrix.
  if (stiffness_matrix == Teuchos::null) return;

  // Get the beam triad interpolation schemes.
  LargeRotations::TriadInterpolationLocalRotationVectors<3, double> triad_interpolation_scheme;
  LargeRotations::TriadInterpolationLocalRotationVectors<3, double> ref_triad_interpolation_scheme;
  GetBeamTriadInterpolationScheme(discret, Teuchos::rcpFromRef(displacement_vector),
      this->Element1(), triad_interpolation_scheme, ref_triad_interpolation_scheme);

  // Set the FAD variables for the solid DOFs. For the terms calculated here we only need first
  // order derivatives.
  auto q_solid =
      GEOMETRYPAIR::InitializeElementData<solid, scalar_type_rot_2nd>::Initialize(this->Element2());
  for (unsigned int i_solid = 0; i_solid < solid::n_dof_; i_solid++)
    q_solid.element_position_(i_solid) =
        Core::FADUtils::HigherOrderFadValue<scalar_type_rot_2nd>::apply(3 + solid::n_dof_,
            3 + i_solid, Core::FADUtils::CastToDouble(this->ele2pos_.element_position_(i_solid)));

  // Get the rotational Lagrange multipliers for this pair.
  const auto& [_, lambda_gid_rot] = mortar_manager->LocationVector(*this);

  std::vector<double> lambda_rot_double;
  Core::FE::ExtractMyValues(global_lambda, lambda_rot_double, lambda_gid_rot);
  Core::LinAlg::Matrix<mortar_rot::n_dof_, 1, double> lambda_rot;
  for (unsigned int i_dof = 0; i_dof < mortar_rot::n_dof_; i_dof++)
    lambda_rot(i_dof) = lambda_rot_double[i_dof];

  // Initialize local matrices.
  Core::LinAlg::Matrix<n_dof_rot_, n_dof_rot_, double> local_stiff_BB(true);
  Core::LinAlg::Matrix<n_dof_rot_, solid::n_dof_, double> local_stiff_BS(true);
  Core::LinAlg::Matrix<solid::n_dof_, n_dof_rot_, double> local_stiff_SB(true);
  Core::LinAlg::Matrix<solid::n_dof_, solid::n_dof_, double> local_stiff_SS(true);

  const auto rot_coupling_type =
      this->Params()->beam_to_solid_volume_meshtying_params()->get_rotational_coupling_type();
  if (rot_coupling_type == Inpar::BeamToSolid::BeamToSolidRotationCoupling::fix_triad_2d)
  {
    // In the case of "fix_triad_2d" we couple both, the ey and ez direction to the beam. Therefore,
    // we have to evaluate the coupling terms w.r.t both of those coupling types.
    evaluate_rotational_coupling_stiff_terms(
        Inpar::BeamToSolid::BeamToSolidRotationCoupling::deformation_gradient_y_2d, q_solid,
        lambda_rot, triad_interpolation_scheme, ref_triad_interpolation_scheme, local_stiff_BB,
        local_stiff_BS, local_stiff_SB, local_stiff_SS);
    evaluate_rotational_coupling_stiff_terms(
        Inpar::BeamToSolid::BeamToSolidRotationCoupling::deformation_gradient_z_2d, q_solid,
        lambda_rot, triad_interpolation_scheme, ref_triad_interpolation_scheme, local_stiff_BB,
        local_stiff_BS, local_stiff_SB, local_stiff_SS);
  }
  else
    evaluate_rotational_coupling_stiff_terms(rot_coupling_type, q_solid, lambda_rot,
        triad_interpolation_scheme, ref_triad_interpolation_scheme, local_stiff_BB, local_stiff_BS,
        local_stiff_SB, local_stiff_SS);

  // Get the GIDs of the solid and beam.
  std::vector<int> lm_beam, gid_solid, lmowner, lmstride;
  this->Element1()->LocationVector(discret, lm_beam, lmowner, lmstride);
  this->Element2()->LocationVector(discret, gid_solid, lmowner, lmstride);
  int rot_dof_indices[] = {3, 4, 5, 12, 13, 14, 18, 19, 20};
  Core::LinAlg::Matrix<n_dof_rot_, 1, int> gid_rot;
  for (unsigned int i = 0; i < n_dof_rot_; i++) gid_rot(i) = lm_beam[rot_dof_indices[i]];

  // Assemble into global matrix.
  for (unsigned int i_dof_beam = 0; i_dof_beam < n_dof_rot_; i_dof_beam++)
  {
    for (unsigned int j_dof_beam = 0; j_dof_beam < n_dof_rot_; j_dof_beam++)
      stiffness_matrix->FEAssemble(
          local_stiff_BB(i_dof_beam, j_dof_beam), gid_rot(i_dof_beam), gid_rot(j_dof_beam));
    for (unsigned int j_dof_solid = 0; j_dof_solid < solid::n_dof_; j_dof_solid++)
      stiffness_matrix->FEAssemble(
          local_stiff_BS(i_dof_beam, j_dof_solid), gid_rot(i_dof_beam), gid_solid[j_dof_solid]);
  }
  for (unsigned int i_dof_solid = 0; i_dof_solid < solid::n_dof_; i_dof_solid++)
  {
    for (unsigned int j_dof_beam = 0; j_dof_beam < n_dof_rot_; j_dof_beam++)
      stiffness_matrix->FEAssemble(
          local_stiff_SB(i_dof_solid, j_dof_beam), gid_solid[i_dof_solid], gid_rot(j_dof_beam));
    for (unsigned int j_dof_solid = 0; j_dof_solid < solid::n_dof_; j_dof_solid++)
      stiffness_matrix->FEAssemble(
          local_stiff_SS(i_dof_solid, j_dof_solid), gid_solid[i_dof_solid], gid_solid[j_dof_solid]);
  }
}

/**
 *
 */
template <typename beam, typename solid, typename mortar, typename mortar_rot>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortarRotation<beam, solid, mortar,
    mortar_rot>::
    evaluate_rotational_coupling_stiff_terms(
        const Inpar::BeamToSolid::BeamToSolidRotationCoupling& rot_coupling_type,
        const GEOMETRYPAIR::ElementData<solid, scalar_type_rot_2nd>& q_solid,
        Core::LinAlg::Matrix<mortar_rot::n_dof_, 1, double>& lambda_rot,
        const LargeRotations::TriadInterpolationLocalRotationVectors<3, double>&
            triad_interpolation_scheme,
        const LargeRotations::TriadInterpolationLocalRotationVectors<3, double>&
            ref_triad_interpolation_scheme,
        Core::LinAlg::Matrix<n_dof_rot_, n_dof_rot_, double>& local_stiff_BB,
        Core::LinAlg::Matrix<n_dof_rot_, solid::n_dof_, double>& local_stiff_BS,
        Core::LinAlg::Matrix<solid::n_dof_, n_dof_rot_, double>& local_stiff_SB,
        Core::LinAlg::Matrix<solid::n_dof_, solid::n_dof_, double>& local_stiff_SS) const
{
  // Initialize variables.
  Core::LinAlg::Matrix<3, 1, double> dr_beam_ref;
  Core::LinAlg::Matrix<4, 1, double> quaternion_beam_double;
  Core::LinAlg::Matrix<3, 1, double> psi_beam_double;
  Core::LinAlg::Matrix<3, 1, scalar_type_rot_1st> psi_beam;
  Core::LinAlg::Matrix<3, 1, scalar_type_rot_2nd> psi_solid;
  Core::LinAlg::Matrix<3, 1, scalar_type_rot_1st> psi_solid_val;
  Core::LinAlg::Matrix<3, 1, scalar_type_rot_1st> psi_rel;
  Core::LinAlg::Matrix<4, 1, scalar_type_rot_1st> quaternion_beam;
  Core::LinAlg::Matrix<4, 1, scalar_type_rot_1st> quaternion_beam_inv;
  Core::LinAlg::Matrix<4, 1, double> quaternion_beam_ref;
  Core::LinAlg::Matrix<4, 1, scalar_type_rot_1st> quaternion_solid;
  Core::LinAlg::Matrix<4, 1, scalar_type_rot_1st> quaternion_rel;
  Core::LinAlg::Matrix<3, 3, double> T_beam;
  Core::LinAlg::Matrix<3, 3, scalar_type_rot_1st> T_solid;
  Core::LinAlg::Matrix<3, 3, scalar_type_rot_1st> T_solid_inv;
  Core::LinAlg::Matrix<3, 3, scalar_type_rot_1st> T_rel;

  Core::LinAlg::Matrix<mortar_rot::n_nodes_, 1, double> lambda_shape_functions;
  Core::LinAlg::Matrix<3, mortar_rot::n_dof_, scalar_type_rot_1st> lambda_shape_functions_full(
      true);
  Core::LinAlg::SerialDenseVector L_i(3);
  Core::LinAlg::Matrix<3, n_dof_rot_, scalar_type_rot_1st> L_full(true);
  std::vector<Core::LinAlg::Matrix<3, 3, double>> I_beam_tilde;
  Core::LinAlg::Matrix<3, n_dof_rot_, double> I_beam_tilde_full;
  Core::LinAlg::Matrix<3, n_dof_rot_, double> T_beam_times_I_beam_tilde_full;
  Core::LinAlg::Matrix<3, mortar_rot::n_dof_, scalar_type_rot_1st> T_rel_tr_times_lambda_shape;
  Core::LinAlg::Matrix<3, mortar_rot::n_dof_, scalar_type_rot_1st>
      T_solid_mtr_times_T_rel_tr_times_lambda_shape;
  Core::LinAlg::Matrix<n_dof_rot_, mortar_rot::n_dof_, scalar_type_rot_1st> d_fb_d_lambda_gp;
  Core::LinAlg::Matrix<solid::n_dof_, mortar_rot::n_dof_, scalar_type_rot_1st> d_fs_d_lambda_gp;
  Core::LinAlg::Matrix<3, solid::n_dof_, scalar_type_rot_1st> d_psi_solid_d_q_solid;
  Core::LinAlg::Matrix<mortar_rot::n_dof_, 3, double> d_g_d_psi_beam;
  Core::LinAlg::Matrix<mortar_rot::n_dof_, solid::n_dof_, double> d_g_d_q_solid;
  Core::LinAlg::Matrix<n_dof_rot_, 1, scalar_type_rot_1st> f_beam;
  Core::LinAlg::Matrix<solid::n_dof_, 1, scalar_type_rot_1st> f_solid;
  Core::LinAlg::Matrix<n_dof_rot_, 3, double> d_f_beam_d_phi;
  Core::LinAlg::Matrix<solid::n_dof_, 3, double> d_f_solid_d_phi;
  Core::LinAlg::Matrix<n_dof_rot_, n_dof_rot_, double>
      d_f_beam_d_phi_times_T_beam_times_I_beam_tilde_full;
  Core::LinAlg::Matrix<solid::n_dof_, n_dof_rot_, double>
      d_f_solid_d_phi_times_T_beam_times_I_beam_tilde_full;

  // Initialize scalar variables.
  double segment_jacobian, beam_segmentation_factor;

  // Calculate the meshtying forces.
  // Loop over segments.
  const unsigned int n_segments = this->line_to_3D_segments_.size();
  for (unsigned int i_segment = 0; i_segment < n_segments; i_segment++)
  {
    // Factor to account for the integration segment length.
    beam_segmentation_factor = 0.5 * this->line_to_3D_segments_[i_segment].GetSegmentLength();

    // Gauss point loop.
    for (unsigned int i_gp = 0;
         i_gp < this->line_to_3D_segments_[i_segment].GetProjectionPoints().size(); i_gp++)
    {
      // Get the current Gauss point.
      const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& projected_gauss_point =
          this->line_to_3D_segments_[i_segment].GetProjectionPoints()[i_gp];

      // Get the jacobian in the reference configuration.
      GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(
          projected_gauss_point.GetEta(), this->ele1posref_, dr_beam_ref);

      // Jacobian including the segment length.
      segment_jacobian = dr_beam_ref.Norm2() * beam_segmentation_factor;

      // Calculate the rotation vector of this cross section.
      triad_interpolation_scheme.get_interpolated_quaternion_at_xi(
          quaternion_beam_double, projected_gauss_point.GetEta());
      Core::LargeRotations::quaterniontoangle(quaternion_beam_double, psi_beam_double);
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        psi_beam(i_dim) = Core::FADUtils::HigherOrderFadValue<scalar_type_rot_1st>::apply(
            3 + solid::n_dof_, i_dim, psi_beam_double(i_dim));
      Core::LargeRotations::angletoquaternion(psi_beam, quaternion_beam);
      quaternion_beam_inv = Core::LargeRotations::inversequaternion(quaternion_beam);

      // Get the solid rotation vector.
      ref_triad_interpolation_scheme.get_interpolated_quaternion_at_xi(
          quaternion_beam_ref, projected_gauss_point.GetEta());
      GetSolidRotationVector<solid>(rot_coupling_type, projected_gauss_point.GetXi(),
          this->ele2posref_, q_solid, quaternion_beam_ref, psi_solid);
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        psi_solid_val(i_dim) = psi_solid(i_dim).val();
      Core::LargeRotations::angletoquaternion(psi_solid_val, quaternion_solid);

      // Calculate the relative rotation vector.
      Core::LargeRotations::quaternionproduct(
          quaternion_beam_inv, quaternion_solid, quaternion_rel);
      Core::LargeRotations::quaterniontoangle(quaternion_rel, psi_rel);

      // Calculate the transformation matrices.
      T_rel = Core::LargeRotations::Tmatrix(psi_rel);
      T_beam = Core::LargeRotations::Tmatrix(Core::FADUtils::CastToDouble(psi_beam));
      T_solid = Core::LargeRotations::Tmatrix(psi_solid_val);
      T_solid_inv = T_solid;
      Core::LinAlg::Inverse(T_solid_inv);

      // Evaluate shape functions.
      GEOMETRYPAIR::EvaluateShapeFunction<mortar_rot>::Evaluate(
          lambda_shape_functions, projected_gauss_point.GetEta());
      for (unsigned int i_node = 0; i_node < mortar_rot::n_nodes_; i_node++)
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          lambda_shape_functions_full(i_dim, 3 * i_node + i_dim) = lambda_shape_functions(i_node);

      Core::FE::shape_function_1D(L_i, projected_gauss_point.GetEta(), Core::FE::CellType::line3);
      for (unsigned int i_node = 0; i_node < 3; i_node++)
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          L_full(i_dim, 3 * i_node + i_dim) = L_i(i_node);

      triad_interpolation_scheme.get_nodal_generalized_rotation_interpolation_matrices_at_xi(
          I_beam_tilde, projected_gauss_point.GetEta());
      for (unsigned int i_node = 0; i_node < 3; i_node++)
        for (unsigned int i_dim_0 = 0; i_dim_0 < 3; i_dim_0++)
          for (unsigned int i_dim_1 = 0; i_dim_1 < 3; i_dim_1++)
            I_beam_tilde_full(i_dim_0, i_node * 3 + i_dim_1) =
                I_beam_tilde[i_node](i_dim_0, i_dim_1);

      // Solid angle derived w.r.t. the solid DOFs.
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        for (unsigned int i_solid = 0; i_solid < solid::n_dof_; i_solid++)
          d_psi_solid_d_q_solid(i_dim, i_solid) = psi_solid(i_dim).dx(3 + i_solid);

      // Calculate the force terms derived w.r.t. the Lagrange multipliers.
      T_rel_tr_times_lambda_shape.MultiplyTN(T_rel, lambda_shape_functions_full);
      d_fb_d_lambda_gp.MultiplyTN(L_full, T_rel_tr_times_lambda_shape);
      d_fb_d_lambda_gp.Scale(-1.0 * projected_gauss_point.GetGaussWeight() * segment_jacobian);

      T_solid_mtr_times_T_rel_tr_times_lambda_shape.MultiplyTN(
          T_solid_inv, T_rel_tr_times_lambda_shape);
      d_fs_d_lambda_gp.MultiplyTN(
          d_psi_solid_d_q_solid, T_solid_mtr_times_T_rel_tr_times_lambda_shape);
      d_fs_d_lambda_gp.Scale(projected_gauss_point.GetGaussWeight() * segment_jacobian);

      // Calculate the force vectors.
      f_beam.PutScalar(0.0);
      for (unsigned int i_row = 0; i_row < n_dof_rot_; i_row++)
        for (unsigned int i_col = 0; i_col < mortar_rot::n_dof_; i_col++)
          f_beam(i_row) += d_fb_d_lambda_gp(i_row, i_col) * lambda_rot(i_col);
      f_solid.PutScalar(0.0);
      for (unsigned int i_row = 0; i_row < solid::n_dof_; i_row++)
        for (unsigned int i_col = 0; i_col < mortar_rot::n_dof_; i_col++)
          f_solid(i_row) += d_fs_d_lambda_gp(i_row, i_col) * lambda_rot(i_col);

      // Derivatives of the force vectors.
      for (unsigned int i_row = 0; i_row < n_dof_rot_; i_row++)
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          d_f_beam_d_phi(i_row, i_dim) = f_beam(i_row).dx(i_dim);
      for (unsigned int i_row = 0; i_row < solid::n_dof_; i_row++)
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          d_f_solid_d_phi(i_row, i_dim) = f_solid(i_row).dx(i_dim);

      T_beam_times_I_beam_tilde_full.Multiply(T_beam, I_beam_tilde_full);
      d_f_beam_d_phi_times_T_beam_times_I_beam_tilde_full.Multiply(
          d_f_beam_d_phi, T_beam_times_I_beam_tilde_full);
      d_f_solid_d_phi_times_T_beam_times_I_beam_tilde_full.Multiply(
          d_f_solid_d_phi, T_beam_times_I_beam_tilde_full);

      // Add to output matrices and vector.
      local_stiff_BB += d_f_beam_d_phi_times_T_beam_times_I_beam_tilde_full;
      for (unsigned int i_beam = 0; i_beam < n_dof_rot_; i_beam++)
        for (unsigned int j_solid = 0; j_solid < solid::n_dof_; j_solid++)
          local_stiff_BS(i_beam, j_solid) += f_beam(i_beam).dx(3 + j_solid);
      local_stiff_SB += d_f_solid_d_phi_times_T_beam_times_I_beam_tilde_full;
      for (unsigned int i_solid = 0; i_solid < solid::n_dof_; i_solid++)
        for (unsigned int j_solid = 0; j_solid < solid::n_dof_; j_solid++)
          local_stiff_SS(i_solid, j_solid) += f_solid(i_solid).dx(3 + j_solid);
    }
  }
}


/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

#define initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(            \
    mortar, mortar_rot)                                                                     \
  template class BeamToSolidVolumeMeshtyingPairMortarRotation<t_hermite, t_hex8, mortar,    \
      mortar_rot>;                                                                          \
  template class BeamToSolidVolumeMeshtyingPairMortarRotation<t_hermite, t_hex20, mortar,   \
      mortar_rot>;                                                                          \
  template class BeamToSolidVolumeMeshtyingPairMortarRotation<t_hermite, t_hex27, mortar,   \
      mortar_rot>;                                                                          \
  template class BeamToSolidVolumeMeshtyingPairMortarRotation<t_hermite, t_tet4, mortar,    \
      mortar_rot>;                                                                          \
  template class BeamToSolidVolumeMeshtyingPairMortarRotation<t_hermite, t_tet10, mortar,   \
      mortar_rot>;                                                                          \
  template class BeamToSolidVolumeMeshtyingPairMortarRotation<t_hermite, t_nurbs27, mortar, \
      mortar_rot>;

  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line2, t_line2);
  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line2, t_line3);
  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line2, t_line4);

  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line3, t_line2);
  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line3, t_line3);
  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line3, t_line4);

  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line4, t_line2);
  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line4, t_line3);
  initialize_template_beam_to_solid_volume_meshtying_pair_mortar_rotation(t_line4, t_line4);
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE
