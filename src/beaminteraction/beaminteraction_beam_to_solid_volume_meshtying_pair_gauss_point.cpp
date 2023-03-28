/*----------------------------------------------------------------------*/
/*! \file

\brief Meshtying element for meshtying between a 3D beam and a 3D solid element.

\level 3
*/


#include "beaminteraction_beam_to_solid_volume_meshtying_pair_gauss_point.H"

#include "linalg_utils_densematrix_inverse.H"
#include "linalg_serialdensematrix.H"
#include "linalg_serialdensevector.H"

#include "beam3_reissner.H"
#include "beam3_kirchhoff.H"
#include "beam3_euler_bernoulli.H"

#include "beaminteraction_contact_params.H"
#include "beaminteraction_beam_to_solid_volume_meshtying_params.H"
#include "beaminteraction_beam_to_solid_utils.H"
#include "geometry_pair_element_functions.H"
#include "geometry_pair_line_to_volume.H"
#include "beam3_reissner.H"
#include "beam3_triad_interpolation_local_rotation_vectors.H"

#include <Epetra_FEVector.h>


/**
 *
 */
template <typename beam, typename solid>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<beam,
    solid>::BeamToSolidVolumeMeshtyingPairGaussPoint()
    : BeamToSolidVolumeMeshtyingPairBase<beam, solid>()
{
  // Empty constructor.
}


/**
 *
 */
template <typename beam, typename solid>
bool BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<beam, solid>::Evaluate(
    LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
    LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
    LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22)
{
  // Call Evaluate on the geometry Pair. Only do this once for meshtying.
  if (!this->meshtying_is_evaluated_)
  {
    LINALG::Matrix<beam::n_dof_, 1, double> beam_coupling_ref;
    LINALG::Matrix<solid::n_dof_, 1, double> solid_coupling_ref;
    this->GetCouplingReferencePosition(beam_coupling_ref, solid_coupling_ref);
    this->CastGeometryPair()->Evaluate(
        beam_coupling_ref, solid_coupling_ref, this->line_to_3D_segments_);
    this->meshtying_is_evaluated_ = true;
  }

  // If there are no intersection segments, return no contact status.
  if (this->line_to_3D_segments_.size() == 0) return false;

  // Initialize variables for position and force vectors.
  LINALG::Matrix<3, 1, double> dr_beam_ref;
  LINALG::Matrix<3, 1, scalar_type> r_beam;
  LINALG::Matrix<3, 1, scalar_type> r_solid;
  LINALG::Matrix<3, 1, scalar_type> force;
  LINALG::Matrix<beam::n_dof_, 1, scalar_type> force_element_1(true);
  LINALG::Matrix<solid::n_dof_, 1, scalar_type> force_element_2(true);

  // Initialize scalar variables.
  double segment_jacobian = 0.0;
  double beam_segmentation_factor = 0.0;
  double penalty_parameter =
      this->Params()->BeamToSolidVolumeMeshtyingParams()->GetPenaltyParameter();

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
          projected_gauss_point.GetEta(), this->ele1posref_, dr_beam_ref, this->Element1());

      // Jacobian including the segment length.
      segment_jacobian = dr_beam_ref.Norm2() * beam_segmentation_factor;

      // Get the current positions on beam and solid.
      GEOMETRYPAIR::EvaluatePosition<beam>(
          projected_gauss_point.GetEta(), this->ele1pos_, r_beam, this->Element1());
      GEOMETRYPAIR::EvaluatePosition<solid>(
          projected_gauss_point.GetXi(), this->ele2pos_, r_solid, this->Element2());

      // Calculate the force in this Gauss point. The sign of the force calculated here is the one
      // that acts on the beam.
      force = r_solid;
      force -= r_beam;
      force.Scale(penalty_parameter);

      // The force vector is in R3, we need to calculate the equivalent nodal forces on the element
      // dof. This is done with the virtual work equation $F \delta r = f \delta q$.
      for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
        for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
          force_element_1(i_dof) += force(i_dir) * r_beam(i_dir).dx(i_dof) *
                                    projected_gauss_point.GetGaussWeight() * segment_jacobian;
      for (unsigned int i_dof = 0; i_dof < solid::n_dof_; i_dof++)
        for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
          force_element_2(i_dof) -= force(i_dir) * r_solid(i_dir).dx(i_dof + beam::n_dof_) *
                                    projected_gauss_point.GetGaussWeight() * segment_jacobian;
    }
  }


  // Fill in the entries for the local matrices and vectors.
  {
    // Resize and initialize the return variables.
    if (forcevec1 != nullptr) forcevec1->Size(beam::n_dof_);
    if (forcevec2 != nullptr) forcevec2->Size(solid::n_dof_);
    if (stiffmat11 != nullptr) stiffmat11->Shape(beam::n_dof_, beam::n_dof_);
    if (stiffmat12 != nullptr) stiffmat12->Shape(beam::n_dof_, solid::n_dof_);
    if (stiffmat21 != nullptr) stiffmat21->Shape(solid::n_dof_, beam::n_dof_);
    if (stiffmat22 != nullptr) stiffmat22->Shape(solid::n_dof_, solid::n_dof_);

    if (forcevec1 != nullptr && forcevec2 != nullptr)
    {
      // $f_1$
      for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
        (*forcevec1)(i_dof) = FADUTILS::CastToDouble(force_element_1(i_dof));
      // $f_2$
      for (unsigned int i_dof = 0; i_dof < solid::n_dof_; i_dof++)
        (*forcevec2)(i_dof) = FADUTILS::CastToDouble(force_element_2(i_dof));
    }

    if (stiffmat11 != nullptr && stiffmat12 != nullptr && stiffmat21 != nullptr &&
        stiffmat22 != nullptr)
    {
      // $k_{11}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < beam::n_dof_; i_dof_1++)
        for (unsigned int i_dof_2 = 0; i_dof_2 < beam::n_dof_; i_dof_2++)
          (*stiffmat11)(i_dof_1, i_dof_2) =
              -FADUTILS::CastToDouble(force_element_1(i_dof_1).dx(i_dof_2));

      // $k_{12}, k_{21}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < beam::n_dof_; i_dof_1++)
      {
        for (unsigned int i_dof_2 = 0; i_dof_2 < solid::n_dof_; i_dof_2++)
        {
          (*stiffmat12)(i_dof_1, i_dof_2) =
              -FADUTILS::CastToDouble(force_element_1(i_dof_1).dx(beam::n_dof_ + i_dof_2));
          (*stiffmat21)(i_dof_2, i_dof_1) =
              -FADUTILS::CastToDouble(force_element_2(i_dof_2).dx(i_dof_1));
        }
      }

      // $k_{22}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < solid::n_dof_; i_dof_1++)
        for (unsigned int i_dof_2 = 0; i_dof_2 < solid::n_dof_; i_dof_2++)
          (*stiffmat22)(i_dof_1, i_dof_2) =
              -FADUTILS::CastToDouble(force_element_2(i_dof_1).dx(beam::n_dof_ + i_dof_2));
    }
  }

  // Return true as there are meshtying contributions.
  return true;
}


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<beam, solid>::EvaluateAndAssemble(
    const Teuchos::RCP<const DRT::Discretization>& discret,
    const Teuchos::RCP<Epetra_FEVector>& force_vector,
    const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix,
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector)
{
  // This function only gives contributions for rotational coupling.
  auto rot_coupling_type =
      this->Params()->BeamToSolidVolumeMeshtyingParams()->GetRotationalCouplingType();
  if (rot_coupling_type == INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::none) return;

  // Call Evaluate on the geometry Pair. Only do this once for meshtying.
  if (!this->meshtying_is_evaluated_)
  {
    LINALG::Matrix<beam::n_dof_, 1, double> beam_coupling_ref;
    LINALG::Matrix<solid::n_dof_, 1, double> solid_coupling_ref;
    this->GetCouplingReferencePosition(beam_coupling_ref, solid_coupling_ref);
    this->CastGeometryPair()->Evaluate(
        beam_coupling_ref, solid_coupling_ref, this->line_to_3D_segments_);
    this->meshtying_is_evaluated_ = true;
  }

  // If there are no intersection segments, return no contact status.
  if (this->line_to_3D_segments_.size() == 0) return;

  // Get the beam triad interpolation schemes.
  LARGEROTATIONS::TriadInterpolationLocalRotationVectors<3, double> triad_interpolation_scheme;
  LARGEROTATIONS::TriadInterpolationLocalRotationVectors<3, double> ref_triad_interpolation_scheme;
  GetBeamTriadInterpolationScheme(*discret, displacement_vector, this->Element1(),
      triad_interpolation_scheme, ref_triad_interpolation_scheme);

  // Set the FAD variables for the solid DOFs.
  LINALG::Matrix<solid::n_dof_, 1, scalar_type_rot_2nd> q_solid(true);
  for (unsigned int i_solid = 0; i_solid < solid::n_dof_; i_solid++)
    q_solid(i_solid) = FADUTILS::HigherOrderFadValue<scalar_type_rot_2nd>::apply(
        3 + solid::n_dof_, 3 + i_solid, FADUTILS::CastToDouble(this->ele2pos_(i_solid)));


  // Initialize local matrices.
  LINALG::Matrix<n_dof_pair_, 1, double> local_force(true);
  LINALG::Matrix<n_dof_pair_, n_dof_pair_, double> local_stiff(true);


  if (rot_coupling_type == INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::fix_triad_2d)
  {
    // In the case of "fix_triad_2d" we couple both, the ey and ez direction to the beam. Therefore,
    // we have to evaluate the coupling terms w.r.t both of those coupling types.
    EvaluateRotationalCouplingTerms(
        INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::deformation_gradient_y_2d, q_solid,
        triad_interpolation_scheme, ref_triad_interpolation_scheme, local_force, local_stiff);
    EvaluateRotationalCouplingTerms(
        INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::deformation_gradient_z_2d, q_solid,
        triad_interpolation_scheme, ref_triad_interpolation_scheme, local_force, local_stiff);
  }
  else
    EvaluateRotationalCouplingTerms(rot_coupling_type, q_solid, triad_interpolation_scheme,
        ref_triad_interpolation_scheme, local_force, local_stiff);


  // Get the GIDs of this pair.
  // The first 9 entries in the vector will be the rotational DOFs of the beam, the other entries
  // are the solid DOFs.
  std::vector<int> lm_beam, lm_solid, lmowner, lmstride;
  this->Element1()->LocationVector(*discret, lm_beam, lmowner, lmstride);
  this->Element2()->LocationVector(*discret, lm_solid, lmowner, lmstride);
  std::array<int, 9> rot_dof_indices = {3, 4, 5, 12, 13, 14, 18, 19, 20};
  LINALG::Matrix<n_dof_pair_, 1, int> gid_pair;
  for (unsigned int i = 0; i < n_dof_rot_; i++) gid_pair(i) = lm_beam[rot_dof_indices[i]];
  for (unsigned int i = 0; i < solid::n_dof_; i++) gid_pair(i + n_dof_rot_) = lm_solid[i];


  // If given, assemble force terms into the global force vector.
  if (force_vector != Teuchos::null)
    force_vector->SumIntoGlobalValues(gid_pair.M(), gid_pair.A(), local_force.A());

  // If given, assemble force terms into the global stiffness matrix.
  if (stiffness_matrix != Teuchos::null)
    for (unsigned int i_dof = 0; i_dof < n_dof_pair_; i_dof++)
      for (unsigned int j_dof = 0; j_dof < n_dof_pair_; j_dof++)
        stiffness_matrix->FEAssemble(
            FADUTILS::CastToDouble(local_stiff(i_dof, j_dof)), gid_pair(i_dof), gid_pair(j_dof));
}

/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<beam,
    solid>::EvaluateRotationalCouplingTerms(  //
    const INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling& rot_coupling_type,
    const LINALG::Matrix<solid::n_dof_, 1, scalar_type_rot_2nd>& q_solid,
    const LARGEROTATIONS::TriadInterpolationLocalRotationVectors<3, double>&
        triad_interpolation_scheme,
    const LARGEROTATIONS::TriadInterpolationLocalRotationVectors<3, double>&
        ref_triad_interpolation_scheme,
    LINALG::Matrix<n_dof_pair_, 1, double>& local_force,
    LINALG::Matrix<n_dof_pair_, n_dof_pair_, double>& local_stiff) const
{
  // Initialize variables.
  LINALG::Matrix<3, 1, double> dr_beam_ref;
  LINALG::Matrix<4, 1, double> quaternion_beam_double;
  LINALG::Matrix<3, 1, double> psi_beam_double;
  LINALG::Matrix<3, 1, scalar_type_rot_1st> psi_beam;
  LINALG::Matrix<3, 1, scalar_type_rot_2nd> psi_solid;
  LINALG::Matrix<3, 1, scalar_type_rot_1st> psi_solid_val;
  LINALG::Matrix<3, 1, scalar_type_rot_1st> psi_rel;
  LINALG::Matrix<4, 1, scalar_type_rot_1st> quaternion_beam;
  LINALG::Matrix<4, 1, scalar_type_rot_1st> quaternion_beam_inv;
  LINALG::Matrix<4, 1, double> quaternion_beam_ref;
  LINALG::Matrix<4, 1, scalar_type_rot_1st> quaternion_solid;
  LINALG::Matrix<4, 1, scalar_type_rot_1st> quaternion_rel;
  LINALG::Matrix<3, 3, double> T_beam;
  LINALG::Matrix<3, 3, scalar_type_rot_1st> T_solid;
  LINALG::Matrix<3, 3, scalar_type_rot_1st> T_solid_inv;
  LINALG::Matrix<3, 1, scalar_type_rot_1st> potential_variation;
  LINALG::Matrix<n_dof_rot_, 1, scalar_type_rot_1st> fc_beam_gp;
  LINALG::Matrix<3, solid::n_dof_, scalar_type_rot_1st> d_psi_solid_d_q_solid;
  LINALG::Matrix<3, 1, scalar_type_rot_1st> Tinv_solid_times_potential_variation;
  LINALG::Matrix<solid::n_dof_, 1, scalar_type_rot_1st> fc_solid_gp;
  Epetra_SerialDenseVector L_i(3);
  LINALG::Matrix<n_dof_rot_, 3, double> d_fc_beam_d_psi_beam;
  LINALG::Matrix<solid::n_dof_, 3, double> d_fc_solid_d_psi_beam;
  std::vector<LINALG::Matrix<3, 3, double>> I_beam_tilde;
  LINALG::Matrix<3, n_dof_rot_, double> I_beam_tilde_full;
  LINALG::Matrix<3, n_dof_rot_, double> T_beam_times_I_beam_tilde_full;
  LINALG::Matrix<n_dof_rot_, n_dof_rot_, double> stiff_beam_beam_gp;
  LINALG::Matrix<solid::n_dof_, n_dof_rot_, double> stiff_solid_beam_gp;
  LINALG::Matrix<n_dof_rot_, solid::n_dof_, double> stiff_beam_solid_gp;

  // Initialize scalar variables.
  double segment_jacobian = 0.0;
  double beam_segmentation_factor = 0.0;
  double rotational_penalty_parameter =
      this->Params()->BeamToSolidVolumeMeshtyingParams()->GetRotationalCouplingPenaltyParameter();

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
          projected_gauss_point.GetEta(), this->ele1posref_, dr_beam_ref, this->Element1());

      // Jacobian including the segment length.
      segment_jacobian = dr_beam_ref.Norm2() * beam_segmentation_factor;

      // Calculate the rotation vector of this cross section.
      triad_interpolation_scheme.GetInterpolatedQuaternionAtXi(
          quaternion_beam_double, projected_gauss_point.GetEta());
      LARGEROTATIONS::quaterniontoangle(quaternion_beam_double, psi_beam_double);
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        psi_beam(i_dim) = FADUTILS::HigherOrderFadValue<scalar_type_rot_1st>::apply(
            3 + solid::n_dof_, i_dim, psi_beam_double(i_dim));
      LARGEROTATIONS::angletoquaternion(psi_beam, quaternion_beam);
      quaternion_beam_inv = LARGEROTATIONS::inversequaternion(quaternion_beam);

      // Get the solid rotation vector.
      ref_triad_interpolation_scheme.GetInterpolatedQuaternionAtXi(
          quaternion_beam_ref, projected_gauss_point.GetEta());
      GetSolidRotationVector<solid>(rot_coupling_type, projected_gauss_point.GetXi(),
          this->ele2posref_, q_solid, quaternion_beam_ref, psi_solid, this->Element2());
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        psi_solid_val(i_dim) = psi_solid(i_dim).val();
      LARGEROTATIONS::angletoquaternion(psi_solid_val, quaternion_solid);

      // Calculate the relative rotation vector.
      LARGEROTATIONS::quaternionproduct(quaternion_beam_inv, quaternion_solid, quaternion_rel);
      LARGEROTATIONS::quaterniontoangle(quaternion_rel, psi_rel);

      // Calculate the transformation matrices.
      T_beam = LARGEROTATIONS::Tmatrix(FADUTILS::CastToDouble(psi_beam));
      T_solid = LARGEROTATIONS::Tmatrix(psi_solid_val);

      // Force terms.
      DRT::UTILS::shape_function_1D(L_i, projected_gauss_point.GetEta(), DRT::Element::line3);
      potential_variation = psi_rel;
      potential_variation.Scale(rotational_penalty_parameter);
      for (unsigned int i_node = 0; i_node < 3; i_node++)
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          fc_beam_gp(3 * i_node + i_dim) = -1.0 * L_i(i_node) * potential_variation(i_dim) *
                                           projected_gauss_point.GetGaussWeight() *
                                           segment_jacobian;
      for (unsigned int i_dof = 0; i_dof < n_dof_rot_; i_dof++)
        local_force(i_dof) += FADUTILS::CastToDouble(fc_beam_gp(i_dof));

      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        for (unsigned int i_solid = 0; i_solid < solid::n_dof_; i_solid++)
          d_psi_solid_d_q_solid(i_dim, i_solid) = psi_solid(i_dim).dx(3 + i_solid);
      T_solid_inv = T_solid;
      LINALG::Inverse(T_solid_inv);
      Tinv_solid_times_potential_variation.MultiplyTN(T_solid_inv, potential_variation);
      fc_solid_gp.MultiplyTN(d_psi_solid_d_q_solid, Tinv_solid_times_potential_variation);
      fc_solid_gp.Scale(projected_gauss_point.GetGaussWeight() * segment_jacobian);
      for (unsigned int i_dof = 0; i_dof < solid::n_dof_; i_dof++)
        local_force(n_dof_rot_ + i_dof) += FADUTILS::CastToDouble(fc_solid_gp(i_dof));


      // Stiffness terms.
      for (unsigned int i_beam_dof = 0; i_beam_dof < n_dof_rot_; i_beam_dof++)
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          d_fc_beam_d_psi_beam(i_beam_dof, i_dim) = fc_beam_gp(i_beam_dof).dx(i_dim);
      for (unsigned int i_solid_dof = 0; i_solid_dof < solid::n_dof_; i_solid_dof++)
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          d_fc_solid_d_psi_beam(i_solid_dof, i_dim) = fc_solid_gp(i_solid_dof).dx(i_dim);

      triad_interpolation_scheme.GetNodalGeneralizedRotationInterpolationMatricesAtXi(
          I_beam_tilde, projected_gauss_point.GetEta());
      for (unsigned int i_node = 0; i_node < 3; i_node++)
        for (unsigned int i_dim_0 = 0; i_dim_0 < 3; i_dim_0++)
          for (unsigned int i_dim_1 = 0; i_dim_1 < 3; i_dim_1++)
            I_beam_tilde_full(i_dim_0, i_node * 3 + i_dim_1) =
                I_beam_tilde[i_node](i_dim_0, i_dim_1);

      T_beam_times_I_beam_tilde_full.Multiply(FADUTILS::CastToDouble(T_beam), I_beam_tilde_full);
      stiff_beam_beam_gp.Multiply(d_fc_beam_d_psi_beam, T_beam_times_I_beam_tilde_full);
      for (unsigned int i_dof = 0; i_dof < n_dof_rot_; i_dof++)
        for (unsigned int j_dof = 0; j_dof < n_dof_rot_; j_dof++)
          local_stiff(i_dof, j_dof) += stiff_beam_beam_gp(i_dof, j_dof);

      stiff_solid_beam_gp.Multiply(d_fc_solid_d_psi_beam, T_beam_times_I_beam_tilde_full);
      for (unsigned int i_dof = 0; i_dof < solid::n_dof_; i_dof++)
        for (unsigned int j_dof = 0; j_dof < n_dof_rot_; j_dof++)
          local_stiff(i_dof + n_dof_rot_, j_dof) += stiff_solid_beam_gp(i_dof, j_dof);

      for (unsigned int i_dof = 0; i_dof < n_dof_rot_; i_dof++)
        for (unsigned int j_dof = 0; j_dof < solid::n_dof_; j_dof++)
          local_stiff(i_dof, j_dof + n_dof_rot_) += fc_beam_gp(i_dof).dx(3 + j_dof);

      for (unsigned int i_dof = 0; i_dof < solid::n_dof_; i_dof++)
        for (unsigned int j_dof = 0; j_dof < solid::n_dof_; j_dof++)
          local_stiff(i_dof + n_dof_rot_, j_dof + n_dof_rot_) += fc_solid_gp(i_dof).dx(3 + j_dof);
    }
  }
}


/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  template class BeamToSolidVolumeMeshtyingPairGaussPoint<t_hermite, t_hex8>;
  template class BeamToSolidVolumeMeshtyingPairGaussPoint<t_hermite, t_hex20>;
  template class BeamToSolidVolumeMeshtyingPairGaussPoint<t_hermite, t_hex27>;
  template class BeamToSolidVolumeMeshtyingPairGaussPoint<t_hermite, t_tet4>;
  template class BeamToSolidVolumeMeshtyingPairGaussPoint<t_hermite, t_tet10>;
  template class BeamToSolidVolumeMeshtyingPairGaussPoint<t_hermite, t_nurbs27>;
}  // namespace BEAMINTERACTION
