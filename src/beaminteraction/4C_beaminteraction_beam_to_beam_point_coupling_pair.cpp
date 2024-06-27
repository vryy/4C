/*----------------------------------------------------------------------*/
/*! \file

\brief Mesh tying element to couple points of two 3D beam elements together.

\level 3
*/


#include "4C_beaminteraction_beam_to_beam_point_coupling_pair.hpp"

#include "4C_beam3_reissner.hpp"
#include "4C_beam3_triad_interpolation_local_rotation_vectors.hpp"
#include "4C_beaminteraction_beam_to_solid_utils.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_geometry_pair_access_traits.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_sparsematrix.hpp"

#include <Epetra_FEVector.h>

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
template <typename beam>
BEAMINTERACTION::BeamToBeamPointCouplingPair<beam>::BeamToBeamPointCouplingPair(
    double penalty_parameter_rot, double penalty_parameter_pos,
    std::array<double, 2> pos_in_parameterspace)
    : BeamContactPair(),
      penalty_parameter_pos_(penalty_parameter_pos),
      penalty_parameter_rot_(penalty_parameter_rot),
      position_in_parameterspace_(pos_in_parameterspace)
{
  // Empty constructor.
}

/**
 *
 */
template <typename beam>
void BEAMINTERACTION::BeamToBeamPointCouplingPair<beam>::setup()
{
  // This pair only works for Simo Reissner beam elements.
  const auto check_simo_reissner_beam = [](auto element)
  {
    const bool is_sr_beam = dynamic_cast<const Discret::ELEMENTS::Beam3r*>(element) != nullptr;
    if (!is_sr_beam)
      FOUR_C_THROW("The BeamToBeamPointCouplingPair only works for Simo Reissner beams");
  };
  check_simo_reissner_beam(this->Element1());
  check_simo_reissner_beam(this->Element2());

  this->issetup_ = true;
}

/**
 *
 */
template <typename beam>
void BEAMINTERACTION::BeamToBeamPointCouplingPair<beam>::EvaluateAndAssemble(
    const Teuchos::RCP<const Core::FE::Discretization>& discret,
    const Teuchos::RCP<Epetra_FEVector>& force_vector,
    const Teuchos::RCP<Core::LinAlg::SparseMatrix>& stiffness_matrix,
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector)
{
  evaluate_and_assemble_positional_coupling(
      discret, force_vector, stiffness_matrix, displacement_vector);
  evaluate_and_assemble_rotational_coupling(
      discret, force_vector, stiffness_matrix, displacement_vector);
}

/**
 *
 */
template <typename beam>
void BEAMINTERACTION::BeamToBeamPointCouplingPair<beam>::evaluate_and_assemble_positional_coupling(
    const Teuchos::RCP<const Core::FE::Discretization>& discret,
    const Teuchos::RCP<Epetra_FEVector>& force_vector,
    const Teuchos::RCP<Core::LinAlg::SparseMatrix>& stiffness_matrix,
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector) const
{
  const std::array<const Core::Elements::Element*, 2> beam_ele = {
      this->Element1(), this->Element2()};

  // Initialize variables for evaluation of the positional coupling terms.
  std::array<Core::LinAlg::Matrix<beam::n_dof_, 1, int>, 2> gid_pos;
  std::array<GEOMETRYPAIR::ElementData<beam, scalar_type_pos>, 2> beam_pos = {
      GEOMETRYPAIR::InitializeElementData<beam, scalar_type_pos>::initialize(beam_ele[0]),
      GEOMETRYPAIR::InitializeElementData<beam, scalar_type_pos>::initialize(beam_ele[1])};
  std::array<Core::LinAlg::Matrix<3, 1, scalar_type_pos>, 2> r;
  Core::LinAlg::Matrix<3, 1, scalar_type_pos> force;
  std::array<Core::LinAlg::Matrix<beam::n_dof_, 1, scalar_type_pos>, 2> force_element;

  // Evaluate individual positions in the two beams.
  for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
  {
    // Get GIDs of the beams positional DOF.
    std::vector<int> lm_beam, lm_solid, lmowner, lmstride;
    beam_ele[i_beam]->LocationVector(*discret, lm_beam, lmowner, lmstride);
    const std::array<int, 12> pos_dof_indices = {0, 1, 2, 6, 7, 8, 9, 10, 11, 15, 16, 17};
    for (unsigned int i = 0; i < beam::n_dof_; i++)
      gid_pos[i_beam](i) = lm_beam[pos_dof_indices[i]];

    // Set current nodal positions (and tangents) for beam element
    std::vector<double> element_posdofvec_absolutevalues(beam::n_dof_, 0.0);
    std::vector<int> lm(beam::n_dof_);
    for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++) lm[i_dof] = gid_pos[i_beam](i_dof);
    BEAMINTERACTION::UTILS::ExtractPosDofVecAbsoluteValues(
        *discret, beam_ele[i_beam], displacement_vector, element_posdofvec_absolutevalues);
    for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
      beam_pos[i_beam].element_position_(i_dof) =
          Core::FADUtils::HigherOrderFadValue<scalar_type_pos>::apply(2 * beam::n_dof_,
              i_beam * beam::n_dof_ + i_dof, element_posdofvec_absolutevalues[i_dof]);

    // Evaluate the position of the coupling point.
    GEOMETRYPAIR::EvaluatePosition<beam>(
        position_in_parameterspace_[i_beam], beam_pos[i_beam], r[i_beam]);
  }

  // Calculate the force between the two beams.
  force = r[1];
  force -= r[0];
  force.scale(penalty_parameter_pos_);

  for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
  {
    double factor;
    if (i_beam == 0)
      factor = -1.0;
    else
      factor = 1.0;

    // The force vector is in R3, we need to calculate the equivalent nodal forces on the
    // element dof. This is done with the virtual work equation $F \delta r = f \delta q$.
    force_element[i_beam].put_scalar(0.0);
    for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        force_element[i_beam](i_dof) +=
            factor * force(i_dir) * r[i_beam](i_dir).dx(i_beam * beam::n_dof_ + i_dof);

    // Add the coupling force to the global force vector.
    if (force_vector != Teuchos::null)
      force_vector->SumIntoGlobalValues(gid_pos[i_beam].numRows(), gid_pos[i_beam].data(),
          Core::FADUtils::CastToDouble(force_element[i_beam]).data());
  }

  // Evaluate and assemble the coupling stiffness terms.
  if (stiffness_matrix != Teuchos::null)
  {
    for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
    {
      for (unsigned int j_beam = 0; j_beam < 2; j_beam++)
      {
        for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
          for (unsigned int j_dof = 0; j_dof < beam::n_dof_; j_dof++)
            stiffness_matrix->FEAssemble(
                force_element[i_beam](i_dof).dx(j_beam * beam::n_dof_ + j_dof),
                gid_pos[i_beam](i_dof), gid_pos[j_beam](j_dof));
      }
    }
  }
}

/**
 *
 */
template <typename beam>
void BEAMINTERACTION::BeamToBeamPointCouplingPair<beam>::evaluate_and_assemble_rotational_coupling(
    const Teuchos::RCP<const Core::FE::Discretization>& discret,
    const Teuchos::RCP<Epetra_FEVector>& force_vector,
    const Teuchos::RCP<Core::LinAlg::SparseMatrix>& stiffness_matrix,
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector) const
{
  const std::array<const Core::Elements::Element*, 2> beam_ele = {
      this->Element1(), this->Element2()};

  // Declare variables for evaluation of the rotational coupling terms.
  std::array<Core::LinAlg::Matrix<n_dof_rot_, 1, int>, 2> gid_rot;
  std::array<Core::LinAlg::Matrix<4, 1, double>, 2> quaternion_ref;
  std::array<Core::LinAlg::Matrix<4, 1, scalar_type_rot>, 2> quaternion;
  std::array<Core::LinAlg::SerialDenseVector, 2> L_i = {
      Core::LinAlg::SerialDenseVector(3), Core::LinAlg::SerialDenseVector(3)};
  std::array<Core::LinAlg::Matrix<3, n_dof_rot_, double>, 2> T_times_I_tilde_full;
  std::array<Core::LinAlg::Matrix<n_dof_rot_, 1, scalar_type_rot>, 2> moment_nodal_load;
  std::array<std::array<Core::LinAlg::Matrix<n_dof_rot_, 3, double>, 2>, 2>
      d_moment_nodal_load_d_psi;

  // Evaluate individual rotational fields in the two beams.
  for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
  {
    // Get GIDs of the beams rotational DOF.
    std::vector<int> lm_beam, lm_solid, lmowner, lmstride;
    beam_ele[i_beam]->LocationVector(*discret, lm_beam, lmowner, lmstride);
    const std::array<int, 9> rot_dof_indices = {3, 4, 5, 12, 13, 14, 18, 19, 20};
    for (unsigned int i = 0; i < n_dof_rot_; i++) gid_rot[i_beam](i) = lm_beam[rot_dof_indices[i]];

    // Get the triad interpolation schemes for the two beams.
    LargeRotations::TriadInterpolationLocalRotationVectors<3, double> triad_interpolation_scheme;
    LargeRotations::TriadInterpolationLocalRotationVectors<3, double>
        ref_triad_interpolation_scheme;
    BEAMINTERACTION::GetBeamTriadInterpolationScheme(*discret, displacement_vector,
        beam_ele[i_beam], triad_interpolation_scheme, ref_triad_interpolation_scheme);

    // Calculate the rotation vector of the beam cross sections and its FAD representation.
    Core::LinAlg::Matrix<4, 1, double> quaternion_double;
    Core::LinAlg::Matrix<3, 1, double> psi_double;
    Core::LinAlg::Matrix<3, 1, scalar_type_rot> psi;
    triad_interpolation_scheme.get_interpolated_quaternion_at_xi(
        quaternion_double, position_in_parameterspace_[i_beam]);
    Core::LargeRotations::quaterniontoangle(quaternion_double, psi_double);
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      psi(i_dim) = Core::FADUtils::HigherOrderFadValue<scalar_type_rot>::apply(
          6, i_beam * rot_dim_ + i_dim, psi_double(i_dim));
    Core::LargeRotations::angletoquaternion(psi, quaternion[i_beam]);

    ref_triad_interpolation_scheme.get_interpolated_quaternion_at_xi(
        quaternion_ref[i_beam], position_in_parameterspace_[i_beam]);

    // Transformation matrix.
    Core::LinAlg::Matrix<3, 3, double> T = Core::LargeRotations::Tmatrix(psi_double);

    // Interpolation matrices.
    std::vector<Core::LinAlg::Matrix<3, 3, double>> I_tilde;
    Core::LinAlg::Matrix<3, n_dof_rot_, double> I_tilde_full;
    Core::FE::shape_function_1D(
        L_i[i_beam], position_in_parameterspace_[i_beam], Core::FE::CellType::line3);
    triad_interpolation_scheme.get_nodal_generalized_rotation_interpolation_matrices_at_xi(
        I_tilde, position_in_parameterspace_[i_beam]);
    for (unsigned int i_node = 0; i_node < 3; i_node++)
      for (unsigned int i_dim_0 = 0; i_dim_0 < 3; i_dim_0++)
        for (unsigned int i_dim_1 = 0; i_dim_1 < 3; i_dim_1++)
          I_tilde_full(i_dim_0, i_node * 3 + i_dim_1) = I_tilde[i_node](i_dim_0, i_dim_1);
    T_times_I_tilde_full[i_beam].multiply(T, I_tilde_full);
  }

  // Get the relative rotation vector between the two cross sections.
  Core::LinAlg::Matrix<4, 1, scalar_type_rot> temp_quaternion_1, temp_quaternion_2, quaternion_rel;
  Core::LinAlg::Matrix<3, 1, scalar_type_rot> psi_rel;
  Core::LinAlg::Matrix<4, 1, scalar_type_rot> quaternion_0_inv =
      Core::LargeRotations::inversequaternion(quaternion[0]);
  Core::LinAlg::Matrix<4, 1, double> quaternion_1_ref_inv =
      Core::LargeRotations::inversequaternion(quaternion_ref[1]);
  Core::LargeRotations::quaternionproduct(quaternion_0_inv, quaternion_ref[0], temp_quaternion_1);
  Core::LargeRotations::quaternionproduct(
      temp_quaternion_1, quaternion_1_ref_inv, temp_quaternion_2);
  Core::LargeRotations::quaternionproduct(temp_quaternion_2, quaternion[1], quaternion_rel);
  Core::LargeRotations::quaterniontoangle(quaternion_rel, psi_rel);

  // Evaluate and assemble the moment coupling terms.
  for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
  {
    for (unsigned int i_node = 0; i_node < 3; i_node++)
    {
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      {
        double factor;
        if (i_beam == 0)
          factor = -1.0;
        else
          factor = 1.0;
        moment_nodal_load[i_beam](3 * i_node + i_dim) =
            factor * penalty_parameter_rot_ * L_i[i_beam](i_node) * psi_rel(i_dim);
      }
    }

    if (force_vector != Teuchos::null)
      force_vector->SumIntoGlobalValues(gid_rot[i_beam].numRows(), gid_rot[i_beam].data(),
          Core::FADUtils::CastToDouble(moment_nodal_load[i_beam]).data());
  }

  // Evaluate and assemble the coupling stiffness terms.
  if (stiffness_matrix != Teuchos::null)
  {
    for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
    {
      for (unsigned int j_beam = 0; j_beam < 2; j_beam++)
        for (unsigned int i_beam_dof = 0; i_beam_dof < n_dof_rot_; i_beam_dof++)
          for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
            d_moment_nodal_load_d_psi[i_beam][j_beam](i_beam_dof, i_dim) =
                moment_nodal_load[i_beam](i_beam_dof).dx(j_beam * 3 + i_dim);

      for (unsigned int j_beam = 0; j_beam < 2; j_beam++)
      {
        Core::LinAlg::Matrix<n_dof_rot_, n_dof_rot_, double> moment_stiff_temp;
        moment_stiff_temp.multiply(
            d_moment_nodal_load_d_psi[i_beam][j_beam], T_times_I_tilde_full[j_beam]);

        for (unsigned int i_dof = 0; i_dof < n_dof_rot_; i_dof++)
          for (unsigned int j_dof = 0; j_dof < n_dof_rot_; j_dof++)
            stiffness_matrix->FEAssemble(
                moment_stiff_temp(i_dof, j_dof), gid_rot[i_beam](i_dof), gid_rot[j_beam](j_dof));
      }
    }
  }
}

/**
 *
 */
template <typename beam>
void BEAMINTERACTION::BeamToBeamPointCouplingPair<beam>::print(std::ostream& out) const
{
  check_init_setup();

  // Print some general information: Element IDs and dofvecs.
  out << "\n------------------------------------------------------------------------";
  out << "\nInstance of BeamToBeamPenaltyPointCouplingPair"
      << "\nBeam1 EleGID:  " << Element1()->Id() << "\nBeam2 EleGID: " << Element2()->Id();
  out << "------------------------------------------------------------------------\n";
}

/**
 *
 */
template <typename beam>
void BEAMINTERACTION::BeamToBeamPointCouplingPair<
    beam>::print_summary_one_line_per_active_segment_pair(std::ostream& out) const
{
  check_init_setup();

  out << "Beam-to-beam point coupling pair, beam1 gid: " << Element1()->Id()
      << " beam2 gid: " << Element2()->Id() << ", position in parameter space: ["
      << position_in_parameterspace_[0] << ", " << position_in_parameterspace_[1] << "]\n";
}

/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  template class BeamToBeamPointCouplingPair<t_hermite>;
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE
