// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_contact_beam_to_beam_point_coupling_pair.hpp"

#include "4C_beam3_base.hpp"
#include "4C_beam3_kirchhoff.hpp"
#include "4C_beam3_reissner.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_contact_beam_to_solid_mortar_manager.hpp"
#include "4C_beaminteraction_contact_beam_to_solid_utils.hpp"
#include "4C_beaminteraction_geometry_pair_access_traits.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_line_to_line.hpp"
#include "4C_linalg_fevector.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_utils_exceptions.hpp"

#include <array>
#include <cstddef>

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
void BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2, n_dof_beam_2>::setup()
{
  beam_elements_ = {this->element1(), this->element2()};

  if (parameters_.use_closest_point_projection)
  {
    evaluate_closest_point_projection();
  }

  this->issetup_ = true;
}

/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
void BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2,
    n_dof_beam_2>::evaluate_and_assemble(const std::shared_ptr<const Core::FE::Discretization>&
                                             discret,
    const std::shared_ptr<Core::LinAlg::FEVector<double>>& force_vector,
    const std::shared_ptr<Core::LinAlg::SparseMatrix>& stiffness_matrix,
    const std::shared_ptr<const Core::LinAlg::Vector<double>>& displacement_vector)
{
  check_init_setup();

  if (!parameters_.evaluate_pair) return;

  // Evaluate the kinematics of the coupling pair.
  const auto pair_kinematic = evaluate_kinematics(*discret, *displacement_vector);

  // Positional coupling terms
  const auto coupling_terms_position = evaluate_positional_coupling(pair_kinematic);

  // Rotational coupling terms
  const auto coupling_terms_rotation = evaluate_rotational_coupling(pair_kinematic);

  // Penalty regularization of Lagrange multipliers
  Core::LinAlg::Matrix<3, 1> lambda_position = coupling_terms_position.constraint;
  lambda_position.scale(parameters_.penalty_parameter_pos);
  Core::LinAlg::Matrix<3, 1> lambda_rotation = coupling_terms_rotation.constraint;
  lambda_rotation.scale(parameters_.penalty_parameter_rot);

  // Penalty regularization for stiffness
  Core::LinAlg::Matrix<12, 12> stiffness(Core::LinAlg::Initialization::zero);

  auto penalty_stiffness =
      [&](const auto& coupling_terms, const auto& lambda, const double penalty_parameter)
  {
    // Stiffness contribution
    Core::LinAlg::Matrix<12, 12> stiffness_temp;
    stiffness_temp.multiply_nn(
        coupling_terms.residuum_lin_lambda, coupling_terms.constraint_lin_kinematic);
    stiffness_temp.scale(penalty_parameter);
    stiffness += stiffness_temp;
  };

  penalty_stiffness(coupling_terms_position, lambda_position, parameters_.penalty_parameter_pos);
  penalty_stiffness(coupling_terms_rotation, lambda_rotation, parameters_.penalty_parameter_rot);

  // Add stiffness contributions due to coupling formulation
  add_coupling_stiffness(stiffness, pair_kinematic, coupling_terms_position, lambda_position,
      coupling_terms_rotation, lambda_rotation);

  // Map residuum and stiffness to element DOFs
  const auto [residuum_pair, stiffness_pair] =
      map_residuum_and_stiffness_to_pair_dof(stiffness, pair_kinematic, coupling_terms_position,
          lambda_position, coupling_terms_rotation, lambda_rotation);

  // Add the coupling terms into the global vector and matrix.
  if (force_vector != nullptr)
    force_vector->sum_into_global_values(
        pair_kinematic.pair_gid.size(), pair_kinematic.pair_gid.data(), residuum_pair.data());
  if (stiffness_matrix != nullptr)
  {
    for (unsigned int i_dof = 0; i_dof < pair_kinematic.pair_gid.size(); i_dof++)
    {
      for (unsigned int j_dof = 0; j_dof < pair_kinematic.pair_gid.size(); j_dof++)
      {
        stiffness_matrix->fe_assemble(stiffness_pair(i_dof, j_dof), pair_kinematic.pair_gid[i_dof],
            pair_kinematic.pair_gid[j_dof]);
      }
    }
  }
}

/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
void BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2, n_dof_beam_2>::print(
    std::ostream& out) const
{
  check_init_setup();

  // Print some general information: Element IDs and dofvecs.
  out << "\n------------------------------------------------------------------------";
  out << "\nInstance of BeamToBeamPenaltyPointCouplingPair"
      << "\nBeam1 EleGID:  " << element1()->id() << "\nBeam2 EleGID: " << element2()->id();
  out << "------------------------------------------------------------------------\n";
}

/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
void BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2,
    n_dof_beam_2>::print_summary_one_line_per_active_segment_pair(std::ostream& out) const
{
  check_init_setup();

  out << "Beam-to-beam point coupling pair, beam1 gid: " << element1()->id()
      << " beam2 gid: " << element2()->id() << ", position in parameter space: ["
      << parameters_.position_in_parameterspace[0] << ", "
      << parameters_.position_in_parameterspace[1] << "]\n";
}


/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
void BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2,
    n_dof_beam_2>::evaluate_and_assemble_mortar_contributions(const Core::FE::Discretization&
                                                                  discret,
    const BeamToSolidMortarManager* mortar_manager,
    Core::LinAlg::SparseMatrix& global_constraint_lin_beam,
    Core::LinAlg::SparseMatrix& global_constraint_lin_solid,
    Core::LinAlg::SparseMatrix& global_force_beam_lin_lambda,
    Core::LinAlg::SparseMatrix& global_force_solid_lin_lambda,
    Core::LinAlg::FEVector<double>& global_constraint, Core::LinAlg::FEVector<double>& global_kappa,
    Core::LinAlg::SparseMatrix& global_kappa_lin_beam,
    Core::LinAlg::SparseMatrix& global_kappa_lin_solid,
    Core::LinAlg::FEVector<double>& global_lambda_active,
    const std::shared_ptr<const Core::LinAlg::Vector<double>>& displacement_vector)
{
  check_init_setup();

  if (!parameters_.evaluate_pair) return;

  // Get the Lagrange multiplier GIDs.
  const auto& [lambda_gid_pos, lambda_gid_rot] = mortar_manager->location_vector(*this);


  auto get_pair_lambda_gid = [&](const std::vector<int>& lambda_gid_total)
  {
    // In this function we extract the pair GIDs from the given total lambda GID vector.
    // We also check that the size is consistent, i.e., a multiple of the spatial dimension.
    constexpr size_t n_dim = 3;

    FOUR_C_ASSERT_ALWAYS(lambda_gid_total.size() >= n_dim && lambda_gid_total.size() % n_dim == 0,
        "Expected a non-zero multiple of {} Lagrange multiplier GIDs, got {}!", n_dim,
        lambda_gid_total.size());

    // Get the GIDs that are not taken by another coupling pair.
    std::vector<int> available_gid;
    for (size_t i = 0; i < lambda_gid_total.size(); i++)
    {
      const int gid = lambda_gid_total[i];
      const int lid = mortar_manager->get_lambda_dof_row_map()->lid(gid);
      if (lid < 0) FOUR_C_THROW("Could not find the GID {} in the lambda dof row map!", gid);
      if (global_lambda_active.get_values()[lid] < 1e-5)
      {
        if (available_gid.size() == 0 && i % n_dim != 0)
          FOUR_C_THROW(
              "Invalid starting Lagrange multiplier index {} for coupling pair. Expected multiple "
              "of {}.",
              gid, n_dim);
        available_gid.push_back(gid);
      }
      if (available_gid.size() == n_dim) break;
    }

    // Perform some sanity checks.
    if (available_gid.size() != n_dim)
      FOUR_C_THROW(
          "Could not find enough available Lagrange multiplier GIDs for the coupling pair! "
          "Increase the value of MAX_NUMBER_OF_PAIRS_PER_ELEMENT");
    // Ensure a consecutive numbering
    for (size_t i = 1; i < available_gid.size(); i++)
    {
      if (available_gid[i - 1] + 1 != available_gid[i])
        FOUR_C_THROW(
            "Invalid Lagrange multiplier GIDs for coupling pair. Expected consecutive numbering. "
            "Got GIDs {} and {} for the entries {} and {}.",
            available_gid[0], available_gid[1], i - 1, i);
    }

    return available_gid;
  };

  const auto lambda_gid_pos_selected = get_pair_lambda_gid(lambda_gid_pos);
  const auto lambda_gid_rot_selected = get_pair_lambda_gid(lambda_gid_rot);
  for (unsigned int i = 0; i < 3; i++)
  {
    lambda_gid_[i] = lambda_gid_pos_selected[i];
    lambda_gid_[3 + i] = lambda_gid_rot_selected[i];
  }

  // Evaluate the kinematics of the coupling pair.
  const auto pair_kinematic = evaluate_kinematics(discret, *displacement_vector);

  // Positional coupling terms
  const auto coupling_terms_position = evaluate_positional_coupling(pair_kinematic);

  // Rotational coupling terms
  const auto coupling_terms_rotation = evaluate_rotational_coupling(pair_kinematic);

  // Assemble into the global constraint vectors
  for (unsigned int i = 0; i < 3; i++)
  {
    global_constraint.sum_into_global_value(
        lambda_gid_[i], 0, coupling_terms_position.constraint(i));
    global_constraint.sum_into_global_value(
        lambda_gid_[i + 3], 0, coupling_terms_rotation.constraint(i));
  }

  // Set the global active vector
  for (unsigned int i = 0; i < 6; i++)
  {
    global_kappa.sum_into_global_value(lambda_gid_[i], 0, 1.0);
    global_lambda_active.sum_into_global_value(lambda_gid_[i], 0, 1.0);
  }

  // Map the matrices to the pair DOFs
  Core::LinAlg::Matrix<3, n_dof_total> constraint_position_lin_kinematic_pair(
      Core::LinAlg::Initialization::zero);
  constraint_position_lin_kinematic_pair.multiply(
      coupling_terms_position.constraint_lin_kinematic, pair_kinematic.right_transformation_matrix);

  Core::LinAlg::Matrix<3, n_dof_total> constraint_rotation_lin_kinematic_pair(
      Core::LinAlg::Initialization::zero);
  constraint_rotation_lin_kinematic_pair.multiply(
      coupling_terms_rotation.constraint_lin_kinematic, pair_kinematic.right_transformation_matrix);

  Core::LinAlg::Matrix<n_dof_total, 3> residuum_position_lin_lambda_pair(
      Core::LinAlg::Initialization::zero);
  residuum_position_lin_lambda_pair.multiply(
      pair_kinematic.left_transformation_matrix, coupling_terms_position.residuum_lin_lambda);

  Core::LinAlg::Matrix<n_dof_total, 3> residuum_rotation_lin_lambda_pair(
      Core::LinAlg::Initialization::zero);
  residuum_rotation_lin_lambda_pair.multiply(
      pair_kinematic.left_transformation_matrix, coupling_terms_rotation.residuum_lin_lambda);

  // Assemble into global matrices.
  for (unsigned int i_dof_lambda = 0; i_dof_lambda < 3; i_dof_lambda++)
  {
    for (unsigned int i_dof_beam_pos = 0; i_dof_beam_pos < n_dof_total; i_dof_beam_pos++)
    {
      global_constraint_lin_beam.fe_assemble(
          constraint_position_lin_kinematic_pair(i_dof_lambda, i_dof_beam_pos),
          lambda_gid_[i_dof_lambda], pair_kinematic.pair_gid[i_dof_beam_pos]);
      global_constraint_lin_beam.fe_assemble(
          constraint_rotation_lin_kinematic_pair(i_dof_lambda, i_dof_beam_pos),
          lambda_gid_[3 + i_dof_lambda], pair_kinematic.pair_gid[i_dof_beam_pos]);

      global_force_beam_lin_lambda.fe_assemble(
          residuum_position_lin_lambda_pair(i_dof_beam_pos, i_dof_lambda),
          pair_kinematic.pair_gid[i_dof_beam_pos], lambda_gid_[i_dof_lambda]);
      global_force_beam_lin_lambda.fe_assemble(
          residuum_rotation_lin_lambda_pair(i_dof_beam_pos, i_dof_lambda),
          pair_kinematic.pair_gid[i_dof_beam_pos], lambda_gid_[3 + i_dof_lambda]);
    }
  }
}

/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
void BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2,
    n_dof_beam_2>::evaluate_and_assemble(const Core::FE::Discretization& discret,
    const BeamToSolidMortarManager* mortar_manager,
    const std::shared_ptr<Core::LinAlg::FEVector<double>>& force_vector,
    const std::shared_ptr<Core::LinAlg::SparseMatrix>& stiffness_matrix,
    const Core::LinAlg::Vector<double>& global_lambda,
    const Core::LinAlg::Vector<double>& displacement_vector)
{
  check_init_setup();

  if (!parameters_.evaluate_pair) return;

  // Get the Lagrange multiplier GIDs.
  auto pair_lambda = Core::FE::extract_values(global_lambda, lambda_gid_);

  // Evaluate the kinematics of the coupling pair.
  const auto pair_kinematic = evaluate_kinematics(discret, displacement_vector);

  // Positional coupling terms
  const auto coupling_terms_position = evaluate_positional_coupling(pair_kinematic);

  // Rotational coupling terms
  const auto coupling_terms_rotation = evaluate_rotational_coupling(pair_kinematic);

  // Store lambda for pair
  Core::LinAlg::Matrix<3, 1> lambda_position{Core::LinAlg::Initialization::zero};
  Core::LinAlg::Matrix<3, 1> lambda_rotation{Core::LinAlg::Initialization::zero};
  for (unsigned int i = 0; i < 3; i++)
  {
    lambda_position(i) = pair_lambda[i];
    lambda_rotation(i) = pair_lambda[3 + i];
  }

  // Add stiffness contributions due to coupling formulation
  Core::LinAlg::Matrix<12, 12> stiffness(Core::LinAlg::Initialization::zero);
  add_coupling_stiffness(stiffness, pair_kinematic, coupling_terms_position, lambda_position,
      coupling_terms_rotation, lambda_rotation);

  // Map residuum and stiffness to element DOFs
  const auto [residuum_pair, stiffness_pair] =
      map_residuum_and_stiffness_to_pair_dof(stiffness, pair_kinematic, coupling_terms_position,
          lambda_position, coupling_terms_rotation, lambda_rotation);

  // Add the coupling terms into the global vector and matrix.
  if (stiffness_matrix != nullptr)
  {
    for (unsigned int i_dof = 0; i_dof < pair_kinematic.pair_gid.size(); i_dof++)
    {
      for (unsigned int j_dof = 0; j_dof < pair_kinematic.pair_gid.size(); j_dof++)
      {
        stiffness_matrix->fe_assemble(stiffness_pair(i_dof, j_dof), pair_kinematic.pair_gid[i_dof],
            pair_kinematic.pair_gid[j_dof]);
      }
    }
  }
}


/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
void BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2,
    n_dof_beam_2>::evaluate_closest_point_projection()
{
  const std::array<const Core::Elements::Element*, 2> beam_elements = {
      this->element1(), this->element2()};

  auto get_beam_ref_data = [&]<typename Element, unsigned int i_beam>()
  {
    auto beam_pos_ref =
        GeometryPair::InitializeElementData<Element, double>::initialize(beam_elements[i_beam]);

    const auto* beam_element =
        dynamic_cast<const Discret::Elements::Beam3Base*>(beam_elements[i_beam]);
    std::vector<double> zero_reference_displacement(
        std::get<i_beam>(std::make_tuple(n_dof_beam_1, n_dof_beam_2)), 0.0);
    std::vector<double> reference_dof_centerline(Element::n_dof_, 0.0);
    beam_element->extract_centerline_dof_values_from_element_state_vector(
        zero_reference_displacement, reference_dof_centerline, true);

    // Get the current and reference position.
    for (unsigned int i_dof = 0; i_dof < Element::n_dof_; i_dof++)
    {
      beam_pos_ref.element_position_(i_dof) = reference_dof_centerline[i_dof];
    }

    return beam_pos_ref;
  };

  // Closest point projection between the two curves
  const auto beam_pos_ref_1 = get_beam_ref_data.template operator()<Beam1, 0>();
  const auto beam_pos_ref_2 = get_beam_ref_data.template operator()<Beam2, 1>();
  const auto projection_result =
      GeometryPair::line_to_line_closest_point_projection(beam_pos_ref_1, beam_pos_ref_2,
          parameters_.position_in_parameterspace[0], parameters_.position_in_parameterspace[1]);

  if (projection_result != GeometryPair::ProjectionResult::projection_found_valid)
  {
    // No projection was found
    parameters_.evaluate_pair = false;
    return;
  }

  // Check the projection distance
  Core::LinAlg::Matrix<3, 1> diff{Core::LinAlg::Initialization::zero};
  Core::LinAlg::Matrix<3, 1> r;
  GeometryPair::evaluate_position<Beam1>(
      parameters_.position_in_parameterspace[0], beam_pos_ref_1, r);
  diff -= r;
  GeometryPair::evaluate_position<Beam2>(
      parameters_.position_in_parameterspace[1], beam_pos_ref_2, r);
  diff += r;

  double beam_radii = 0.0;
  for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
  {
    const auto* beam_ptr = dynamic_cast<const Discret::Elements::Beam3Base*>(beam_elements[i_beam]);
    beam_radii += beam_ptr->get_circular_cross_section_radius_for_interactions();
  }
  if (parameters_.projection_valid_factor * beam_radii < diff.norm2())
  {
    parameters_.evaluate_pair = false;
    return;
  }

  // Make sure that we have unique pairs for projections directly on nodes
  if (not line_to_line_evaluation_data_->evaluate_projection_coordinates(
          beam_elements, parameters_.position_in_parameterspace))
  {
    parameters_.evaluate_pair = false;
    return;
  }

  // All checks passed, this pair shall be evaluated.
  parameters_.evaluate_pair = true;
}


/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
BeamInteraction::BeamToBeamKinematic<n_dof_beam_1, n_dof_beam_2>
BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2,
    n_dof_beam_2>::evaluate_kinematics(const Core::FE::Discretization& discret,
    const Core::LinAlg::Vector<double>& displacement_vector) const
{
  BeamInteraction::BeamToBeamKinematic<n_dof_beam_1, n_dof_beam_2> pair_kinematic;

  // Evaluate all quantities at the current state for both beams.
  {
    for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
    {
      const auto* beam_element =
          dynamic_cast<const Discret::Elements::Beam3Base*>(beam_elements_[i_beam]);

      std::vector<int> lm, lmowner, lmstride;
      beam_element->location_vector(discret, lm, lmowner, lmstride);
      auto eledisp = Core::FE::extract_values(displacement_vector, lm);
      auto eledisp_ref = std::vector<double>(eledisp.size(), 0.0);
      pair_kinematic.element_displacement[i_beam] = eledisp;

      for (unsigned int i_dof = 0; i_dof < n_dof_beam[i_beam]; i_dof++)
        pair_kinematic.pair_gid[i_dof + i_beam * n_dof_beam_1] = lm[i_dof];

      beam_element->get_pos_at_xi(
          pair_kinematic.r[i_beam], parameters_.position_in_parameterspace[i_beam], eledisp);
      beam_element->get_pos_at_xi(pair_kinematic.r_ref[i_beam],
          parameters_.position_in_parameterspace[i_beam], eledisp_ref);

      Core::LinAlg::Matrix<3, 3> triad;
      Core::LinAlg::Matrix<3, 1> psi_double;
      Core::LinAlg::Matrix<3, 1, scalar_type_rot> psi;
      beam_element->get_triad_at_xi(triad, parameters_.position_in_parameterspace[i_beam], eledisp);
      Core::LinAlg::Matrix<4, 1, double> quaternion_double;
      Core::LargeRotations::triadtoquaternion(triad, quaternion_double);
      Core::LargeRotations::quaterniontoangle(quaternion_double, psi_double);
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        psi(i_dim) = Core::FADUtils::HigherOrderFadValue<scalar_type_rot>::apply(
            6, i_beam * 3 + i_dim, psi_double(i_dim));
      Core::LargeRotations::angletoquaternion(psi, pair_kinematic.cross_section_quaternion[i_beam]);

      Core::LinAlg::Matrix<3, 3> triad_ref;
      beam_element->get_triad_at_xi(
          triad_ref, parameters_.position_in_parameterspace[i_beam], eledisp_ref);
      Core::LargeRotations::triadtoquaternion(
          triad_ref, pair_kinematic.cross_section_quaternion_ref[i_beam]);

      Core::LinAlg::SerialDenseMatrix trafomatrix_left;
      trafomatrix_left.shape(6, n_dof_beam[i_beam]);
      beam_element->get_generalized_interpolation_matrix_variations_at_xi(
          trafomatrix_left, parameters_.position_in_parameterspace[i_beam], eledisp);

      Core::LinAlg::SerialDenseMatrix trafomatrix_right;
      trafomatrix_right.shape(6, n_dof_beam[i_beam]);
      beam_element->get_generalized_interpolation_matrix_increments_at_xi(
          trafomatrix_right, parameters_.position_in_parameterspace[i_beam], eledisp);

      for (unsigned int i_dof = 0; i_dof < n_dof_beam[i_beam]; i_dof++)
      {
        for (unsigned int i_dir = 0; i_dir < 6; i_dir++)
        {
          const unsigned int i_row = i_dir + i_beam * 6;
          const unsigned int i_col = i_dof + i_beam * n_dof_beam_1;
          pair_kinematic.left_transformation_matrix(i_col, i_row) = trafomatrix_left(i_dir, i_dof);
          pair_kinematic.right_transformation_matrix(i_row, i_col) =
              trafomatrix_right(i_dir, i_dof);
        }
      }
    }
  }

  return pair_kinematic;
}

/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
BeamInteraction::BeamToBeamCouplingTerms<n_dof_beam_1, n_dof_beam_2>
BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2,
    n_dof_beam_2>::evaluate_positional_coupling(const BeamToBeamKinematic<n_dof_beam_1,
    n_dof_beam_2>& pair_kinematic) const
{
  // Note: If the two points are not at the same position in the reference configuration, then even
  // for the positional coupling terms there will be contributions to the rotational equilibrium,
  // i.e., moments. To keep things simple, we first evaluate the terms that will affect the
  // positional DOFs and then we compute the terms that will affect the rotational DOFs.

  BeamToBeamCouplingTerms<n_dof_beam_1, n_dof_beam_2> coupling_terms;

  // Constraint (pure positional terms)
  coupling_terms.constraint = pair_kinematic.r[1];
  coupling_terms.constraint -= pair_kinematic.r[0];

  // Linearizations (pure positional terms)
  for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
  {
    coupling_terms.constraint_lin_kinematic(i_dir, i_dir) = -1.0;
    coupling_terms.constraint_lin_kinematic(i_dir, 6 + i_dir) = 1.0;

    coupling_terms.residuum_lin_lambda(i_dir, i_dir) = -1.0;
    coupling_terms.residuum_lin_lambda(6 + i_dir, i_dir) = 1.0;
  }

  // Terms arising due to non matching points in the reference configuration.
  {
    // Evaluate reference offset in the material and spatial frames.
    Core::LinAlg::Matrix<3, 1, double> director_21_ref_spatial = pair_kinematic.r_ref[1];
    director_21_ref_spatial -= pair_kinematic.r_ref[0];
    std::array<Core::LinAlg::Matrix<3, 1, double>, 2> director_21_ref_material;
    for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
    {
      Core::LinAlg::Matrix<3, 3, double> rotation_matrix;
      Core::LargeRotations::quaterniontotriad(
          pair_kinematic.cross_section_quaternion_ref[i_beam], rotation_matrix);
      director_21_ref_material[i_beam].multiply_tn(rotation_matrix, director_21_ref_spatial);
    }

    // Add to constraint vector.
    // We add a half of the reference offset vector in each frame of the two beams. By doing it this
    // way, we obtain a coupling formulation invariant to the ordering of the cross-sections in this
    // pair.
    std::array<Core::LinAlg::Matrix<3, 1, double>, 2> director_21_desired_spatial_half;
    for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
    {
      Core::LinAlg::Matrix<3, 3, double> rotation_matrix;
      Core::LargeRotations::quaterniontotriad(
          Core::FADUtils::cast_to_double(pair_kinematic.cross_section_quaternion[i_beam]),
          rotation_matrix);
      director_21_desired_spatial_half[i_beam].multiply_nn(
          rotation_matrix, director_21_ref_material[i_beam]);
      director_21_desired_spatial_half[i_beam].scale(0.5);
      coupling_terms.constraint -= director_21_desired_spatial_half[i_beam];
    }

    // Current director from point 1 to 2.
    Core::LinAlg::Matrix<3, 1, double> director_21_spatial_half = pair_kinematic.r[1];
    director_21_spatial_half -= pair_kinematic.r[0];
    director_21_spatial_half.scale(0.5);
    Core::LinAlg::Matrix<3, 3, double> skew_director_21_spatial_half;
    Core::LargeRotations::computespin(skew_director_21_spatial_half, director_21_spatial_half);

    // Add the linearizations of the previous constraint terms.
    for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
    {
      Core::LinAlg::Matrix<3, 3, double> skew_director_21_desired_spatial_half;
      Core::LargeRotations::computespin(
          skew_director_21_desired_spatial_half, director_21_desired_spatial_half[i_beam]);

      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
      {
        for (unsigned int j_dir = 0; j_dir < 3; j_dir++)
        {
          coupling_terms.constraint_lin_kinematic(i_dir, j_dir + 3 + 6 * i_beam) +=
              skew_director_21_desired_spatial_half(i_dir, j_dir);

          coupling_terms.residuum_lin_lambda(j_dir + 3 + 6 * i_beam, i_dir) -=
              skew_director_21_spatial_half(j_dir, i_dir);
        }
      }
    }
  }

  return coupling_terms;
}

/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
BeamInteraction::BeamToBeamCouplingTerms<n_dof_beam_1, n_dof_beam_2>
BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2,
    n_dof_beam_2>::evaluate_rotational_coupling(const BeamToBeamKinematic<n_dof_beam_1,
    n_dof_beam_2>& pair_kinematic) const
{
  BeamToBeamCouplingTerms<n_dof_beam_1, n_dof_beam_2> coupling_terms;

  // Get the relative rotation vector between the two cross sections.
  Core::LinAlg::Matrix<4, 1, scalar_type_rot> temp_quaternion_1, temp_quaternion_2, quaternion_rel;
  Core::LinAlg::Matrix<3, 1, scalar_type_rot> psi_rel;
  Core::LinAlg::Matrix<4, 1, scalar_type_rot> quaternion_0_inv =
      Core::LargeRotations::inversequaternion(pair_kinematic.cross_section_quaternion[0]);
  Core::LinAlg::Matrix<4, 1, double> quaternion_1_ref_inv =
      Core::LargeRotations::inversequaternion(pair_kinematic.cross_section_quaternion_ref[1]);
  Core::LargeRotations::quaternionproduct(
      quaternion_0_inv, pair_kinematic.cross_section_quaternion_ref[0], temp_quaternion_1);
  Core::LargeRotations::quaternionproduct(
      temp_quaternion_1, quaternion_1_ref_inv, temp_quaternion_2);
  Core::LargeRotations::quaternionproduct(
      temp_quaternion_2, pair_kinematic.cross_section_quaternion[1], quaternion_rel);
  Core::LargeRotations::quaterniontoangle(quaternion_rel, psi_rel);

  // Transformation matrix
  const Core::LinAlg::Matrix<3, 3, scalar_type_rot> T_psi_rel =
      Core::LargeRotations::tmatrix(psi_rel);

  coupling_terms.constraint = Core::FADUtils::cast_to_double(psi_rel);
  for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
  {
    const double load_factor = i_beam == 0 ? -1.0 : 1.0;

    // Get the rotation angle
    Core::LinAlg::Matrix<3, 1, double> psi_beam;
    Core::LargeRotations::quaterniontoangle(
        Core::FADUtils::cast_to_double(pair_kinematic.cross_section_quaternion[i_beam]), psi_beam);

    const Core::LinAlg::Matrix<3, 3, double> T_psi_beam = Core::LargeRotations::tmatrix(psi_beam);
    Core::LinAlg::Matrix<3, 3, double> d_psi_rel_d_psi_beam;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
    {
      for (unsigned int j_dim = 0; j_dim < 3; j_dim++)
      {
        d_psi_rel_d_psi_beam(i_dim, j_dim) = psi_rel(i_dim).dx(j_dim + i_beam * 3);
      }
    }

    Core::LinAlg::Matrix<3, 3, double> constraint_lin_psi_beam;
    constraint_lin_psi_beam.multiply_nn(d_psi_rel_d_psi_beam, T_psi_beam);
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
    {
      for (unsigned int j_dim = 0; j_dim < 3; j_dim++)
      {
        coupling_terms.constraint_lin_kinematic(i_dim, j_dim + 3 + i_beam * 6) =
            constraint_lin_psi_beam(i_dim, j_dim);
        // Note: Here we insert the transposed, thus the switched ordering in the target but not
        // the source.
        coupling_terms.residuum_lin_lambda(j_dim + 3 + i_beam * 6, i_dim) =
            load_factor * Core::FADUtils::cast_to_double(T_psi_rel(i_dim, j_dim));
      }
    }
    for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < 3; j++)
      {
        for (unsigned int l = 0; l < 3; l++)
        {
          for (unsigned int k = 0; k < 3; k++)
          {
            coupling_terms.evaluation_data_rotation[i_beam][i][j][l] +=
                T_psi_rel(j, i).dx(k + 3 * i_beam) * T_psi_beam(k, l);
          }
        }
      }
    }
  }

  return coupling_terms;
}


/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
void BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2,
    n_dof_beam_2>::add_coupling_stiffness(Core::LinAlg::Matrix<12, 12>& stiffness,
    const BeamToBeamKinematic<n_dof_beam_1, n_dof_beam_2>& pair_kinematic,
    const BeamToBeamCouplingTerms<n_dof_beam_1, n_dof_beam_2>& coupling_terms_position,
    const Core::LinAlg::Matrix<3, 1>& lambda_position,
    const BeamToBeamCouplingTerms<n_dof_beam_1, n_dof_beam_2>& coupling_terms_rotation,
    const Core::LinAlg::Matrix<3, 1>& lambda_rotation)
{
  Core::LinAlg::Matrix<3, 3, double> skew_lambda;
  Core::LargeRotations::computespin(skew_lambda, lambda_position);
  for (unsigned int i_beam_rot = 0; i_beam_rot < 2; i_beam_rot++)
  {
    for (unsigned int i_beam_pos = 0; i_beam_pos < 2; i_beam_pos++)
    {
      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
      {
        for (unsigned int j_dir = 0; j_dir < 3; j_dir++)
        {
          const double factor = (i_beam_pos == 0) ? -0.5 : 0.5;
          stiffness(i_dir + 3 + 6 * i_beam_rot, j_dir + 6 * i_beam_pos) +=
              factor * skew_lambda(i_dir, j_dir);
        }
      }
    }
  }

  for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
  {
    for (unsigned int i = 0; i < 3; i++)
    {
      for (unsigned int j = 0; j < 3; j++)
      {
        for (unsigned int l = 0; l < 3; l++)
        {
          stiffness(i + 3, l + 3 + 6 * i_beam) -=
              coupling_terms_rotation.evaluation_data_rotation[i_beam][i][j][l] *
              lambda_rotation(j);
          stiffness(i + 9, l + 3 + 6 * i_beam) +=
              coupling_terms_rotation.evaluation_data_rotation[i_beam][i][j][l] *
              lambda_rotation(j);
        }
      }
    }
  }
}

/**
 *
 */
template <typename Beam1, unsigned int n_dof_beam_1, typename Beam2, unsigned int n_dof_beam_2>
std::pair<Core::LinAlg::Matrix<n_dof_beam_1 + n_dof_beam_2, 1>,
    Core::LinAlg::Matrix<n_dof_beam_1 + n_dof_beam_2, n_dof_beam_1 + n_dof_beam_2>>
BeamInteraction::BeamToBeamPointCouplingPair<Beam1, n_dof_beam_1, Beam2,
    n_dof_beam_2>::map_residuum_and_stiffness_to_pair_dof(Core::LinAlg::Matrix<12, 12>& stiffness,
    const BeamToBeamKinematic<n_dof_beam_1, n_dof_beam_2>& pair_kinematic,
    const BeamToBeamCouplingTerms<n_dof_beam_1, n_dof_beam_2>& coupling_terms_position,
    const Core::LinAlg::Matrix<3, 1>& lambda_position,
    const BeamToBeamCouplingTerms<n_dof_beam_1, n_dof_beam_2>& coupling_terms_rotation,
    const Core::LinAlg::Matrix<3, 1>& lambda_rotation) const
{
  // Map residuum and stiffness to element DOFs
  Core::LinAlg::Matrix<12, 1> residuum(Core::LinAlg::Initialization::zero);

  // Coupling residuum contribution
  Core::LinAlg::Matrix<12, 1> residuum_temp;
  residuum_temp.multiply_nn(coupling_terms_position.residuum_lin_lambda, lambda_position);
  residuum += residuum_temp;
  residuum_temp.multiply_nn(coupling_terms_rotation.residuum_lin_lambda, lambda_rotation);
  residuum += residuum_temp;

  // Pair residuum and stiffness
  Core::LinAlg::Matrix<n_dof_total, 1> residuum_pair{Core::LinAlg::Initialization::zero};
  Core::LinAlg::Matrix<n_dof_total, n_dof_total> stiffness_pair{Core::LinAlg::Initialization::zero};

  // Map the residuum and stiffness to the pair DOFs
  residuum_pair.multiply(pair_kinematic.left_transformation_matrix, residuum);
  Core::LinAlg::Matrix<n_dof_total, 12> temp_matrix;
  temp_matrix.multiply(pair_kinematic.left_transformation_matrix, stiffness);
  stiffness_pair.multiply(temp_matrix, pair_kinematic.right_transformation_matrix);

  // Add stiffness contributions due to the beam formulation
  for (unsigned int i_beam = 0; i_beam < 2; i_beam++)
  {
    const auto* beam_element =
        dynamic_cast<const Discret::Elements::Beam3Base*>(beam_elements_[i_beam]);

    Core::LinAlg::SerialDenseVector force(6, true);
    for (unsigned int i_dim = 0; i_dim < 6; i_dim++) force(i_dim) = residuum(i_dim + 6 * i_beam);

    Core::LinAlg::SerialDenseMatrix stiffness_beam;
    stiffness_beam.shape(n_dof_beam[i_beam], n_dof_beam[i_beam]);
    beam_element->get_stiffmat_resulting_from_generalized_interpolation_matrix_at_xi(stiffness_beam,
        parameters_.position_in_parameterspace[i_beam], pair_kinematic.element_displacement[i_beam],
        force);

    for (unsigned int i_dof = 0; i_dof < n_dof_beam[i_beam]; i_dof++)
    {
      for (unsigned int j_dof = 0; j_dof < n_dof_beam[i_beam]; j_dof++)
      {
        stiffness_pair(i_dof + i_beam * n_dof_beam_1, j_dof + i_beam * n_dof_beam_1) +=
            stiffness_beam(i_dof, j_dof);
      }
    }
  }

  return {residuum_pair, stiffness_pair};
}

/**
 *
 */
std::unique_ptr<BeamInteraction::BeamContactPair>
BeamInteraction::beam_to_beam_point_coupling_pair_factory(
    const std::array<const Core::Elements::Element*, 2>& element_ptrs,
    const BeamToBeamPointCouplingPairParameters& parameters,
    const std::shared_ptr<GeometryPair::LineToLineEvaluationData>& line_to_line_evaluation_data)
{
  const auto get_element_info =
      [&](size_t i_beam) -> std::tuple<bool, Core::FE::CellType, unsigned int>
  {
    const auto* reissner_beam =
        dynamic_cast<const Discret::Elements::Beam3r*>(element_ptrs[i_beam]);
    const auto* kirchhoff_beam =
        dynamic_cast<const Discret::Elements::Beam3k*>(element_ptrs[i_beam]);

    unsigned int n_dof = 0;
    bool is_hermite = false;
    if (reissner_beam != nullptr)
    {
      is_hermite = reissner_beam->hermite_centerline_interpolation();
      if (is_hermite)
      {
        if (reissner_beam->shape() != Core::FE::CellType::line3)
        {
          FOUR_C_THROW(
              "The beam_to_beam_point_coupling_pair_factory only works for Simo Reissner beam "
              "elements with Hermite shape functions of base type line3");
        }
        n_dof = 21;
      }
      else
      {
        n_dof = reissner_beam->num_centerline_nodes() * 6;
      }
    }
    else if (kirchhoff_beam != nullptr)
    {
      if (kirchhoff_beam->shape() != Core::FE::CellType::line3 or
          not kirchhoff_beam->hermite_centerline_interpolation())
      {
        FOUR_C_THROW(
            "The beam_to_beam_point_coupling_pair_factory only works for Kirchhoff beam elements "
            "with Hermite shape functions");
      }
      is_hermite = true;
      n_dof = 15;
    }
    else
    {
      FOUR_C_THROW(
          "The beam_to_beam_point_coupling_pair_factory only works for Simo Reissner or "
          "Kirchhoff beam elements");
    }
    return {is_hermite, element_ptrs[i_beam]->shape(), n_dof};
  };

  const auto element_info_1 = get_element_info(0);
  const auto element_info_2 = get_element_info(1);

  auto dispatch_beam_type = []<typename F>(bool is_hermite, unsigned int n_dof,
                                F&& f) -> std::unique_ptr<BeamInteraction::BeamContactPair>
  {
    if (is_hermite)
    {
      if (n_dof == 21) return f.template operator()<GeometryPair::t_hermite, 21>();
      if (n_dof == 15) return f.template operator()<GeometryPair::t_hermite, 15>();
    }
    else
    {
      if (n_dof == 12) return f.template operator()<GeometryPair::t_line2, 12>();
      if (n_dof == 18) return f.template operator()<GeometryPair::t_line3, 18>();
      if (n_dof == 24) return f.template operator()<GeometryPair::t_line4, 24>();
      if (n_dof == 30) return f.template operator()<GeometryPair::t_line5, 30>();
    }

    FOUR_C_THROW("Unsupported beam type with ndof {}", n_dof);
  };

  return dispatch_beam_type(std::get<0>(element_info_1), std::get<2>(element_info_1),
      [&parameters, &line_to_line_evaluation_data, &element_info_2,
          &dispatch_beam_type]<typename Beam1, unsigned int n_dof_1_temp>()
          -> std::unique_ptr<BeamInteraction::BeamContactPair>
      {
        return dispatch_beam_type(std::get<0>(element_info_2), std::get<2>(element_info_2),
            [&parameters,
                &line_to_line_evaluation_data]<typename Beam2, unsigned int n_dof_2_temp>()
                -> std::unique_ptr<BeamInteraction::BeamContactPair>
            {
              using Pair = BeamToBeamPointCouplingPair<Beam1, n_dof_1_temp, Beam2, n_dof_2_temp>;
              return std::make_unique<Pair>(parameters, line_to_line_evaluation_data);
            });
      });
}


FOUR_C_NAMESPACE_CLOSE
