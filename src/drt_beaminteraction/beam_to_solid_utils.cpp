/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for beam-to-solid interactions

\level 3

*/
// End doxygen header.


#include "beam_to_solid_utils.H"

#include "beam_contact_pair.H"
#include "beaminteraction_calc_utils.H"
#include "beam_to_solid_mortar_manager.H"
#include "beam_to_solid_surface_contact_params.H"
#include "../drt_geometry_pair/geometry_pair.H"
#include "../drt_geometry_pair/geometry_pair_element.H"
#include "../drt_geometry_pair/geometry_pair_scalar_types.H"
#include "../headers/FAD_utils.H"
#include "../drt_geometry_pair/geometry_pair_element_functions.H"
#include "../drt_fem_general/largerotations.H"
#include "../drt_inpar/inpar_beam_to_solid.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/triad_interpolation_local_rotation_vectors.H"

#include <Epetra_FEVector.h>


/**
 *
 */
template <typename scalar_type>
scalar_type BEAMINTERACTION::PenaltyForce(const scalar_type& gap,
    const Teuchos::RCP<const BeamToSolidSurfaceContactParams>& contact_params)
{
  const INPAR::BEAMTOSOLID::BeamToSolidSurfaceContactPenaltyLaw penalty_law =
      contact_params->GetPenaltyLaw();
  const double penalty_parameter = contact_params->GetPenaltyParameter();

  scalar_type penalty_force = 0.0;

  switch (penalty_law)
  {
    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceContactPenaltyLaw::linear:
    {
      if (gap < 0.0)
      {
        penalty_force = -gap * penalty_parameter;
      }
      break;
    }
    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceContactPenaltyLaw::linear_quadratic:
    {
      const double penalty_parameter_g0 = contact_params->GetPenaltyParameterG0();

      if (gap < 0.0)
      {
        penalty_force = 0.5 * (-2.0 * gap + penalty_parameter_g0) * penalty_parameter;
      }
      else if (gap < penalty_parameter_g0)
      {
        penalty_force = 0.5 * std::pow(gap - penalty_parameter_g0, 2) * penalty_parameter /
                        penalty_parameter_g0;
      }
      break;
    }
    default:
      dserror("Got unexpected penalty law.");
      break;
  }

  return penalty_force;
}

/**
 *
 */
template <typename scalar_type>
scalar_type BEAMINTERACTION::PenaltyPotential(const scalar_type& gap,
    const Teuchos::RCP<const BeamToSolidSurfaceContactParams>& contact_params)
{
  const INPAR::BEAMTOSOLID::BeamToSolidSurfaceContactPenaltyLaw penalty_law =
      contact_params->GetPenaltyLaw();
  const double penalty_parameter = contact_params->GetPenaltyParameter();

  scalar_type penalty_potential = 0.0;

  switch (penalty_law)
  {
    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceContactPenaltyLaw::linear:
    {
      if (gap < 0.0)
      {
        penalty_potential = 0.5 * std::pow(gap, 2) * penalty_parameter;
      }
      break;
    }
    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceContactPenaltyLaw::linear_quadratic:
    {
      const double penalty_parameter_g0 = contact_params->GetPenaltyParameterG0();

      if (gap < 0.0)
      {
        penalty_potential = penalty_parameter / 6.0 *
                            (3.0 * std::pow(gap, 2) - 3.0 * gap * penalty_parameter_g0 +
                                std::pow(penalty_parameter_g0, 2));
      }
      else if (gap < penalty_parameter_g0)
      {
        penalty_potential = penalty_parameter / (6.0 * penalty_parameter_g0) *
                            std::pow(penalty_parameter_g0 - gap, 3);
      }
      break;
    }
    default:
      dserror("Got unexpected penalty law.");
      break;
  }

  return penalty_potential;
}

/**
 *
 */
void BEAMINTERACTION::MortarShapeFunctionsToNumberOfLagrangeValues(
    const INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions shape_function,
    unsigned int& n_lambda_node, unsigned int& n_lambda_element)
{
  switch (shape_function)
  {
    case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::none:
    {
      n_lambda_node = 0;
      n_lambda_element = 0;
      return;
    }
    case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line2:
    {
      n_lambda_node = 1 * 3;
      n_lambda_element = 0 * 3;
      return;
    }
    case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line3:
    {
      n_lambda_node = 1 * 3;
      n_lambda_element = 1 * 3;
      return;
    }
    case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line4:
    {
      n_lambda_node = 1 * 3;
      n_lambda_element = 2 * 3;
      return;
    }
    default:
      dserror("Mortar shape function not implemented!");
  }
}

/**
 *
 */
void BEAMINTERACTION::GetMortarGID(const BeamToSolidMortarManager* mortar_manager,
    const BEAMINTERACTION::BeamContactPair* contact_pair, const unsigned int n_mortar_pos,
    const unsigned int n_mortar_rot, std::vector<int>* lambda_gid_pos,
    std::vector<int>* lambda_gid_rot)
{
  std::vector<int> lambda_total;
  mortar_manager->LocationVector(contact_pair, lambda_total);

#if DEBUG
  if (lambda_total.size() != n_mortar_pos + n_mortar_rot)
    dserror("BEAMINTERACTION::GetMortarGID the local and global GID sizes do not match.");
#endif

  unsigned int n_nodal_dof = 3;
  if (n_mortar_rot > 0) n_nodal_dof = 6;

  if (lambda_gid_pos != nullptr)
  {
    // Get the Lagrange multiplier GIDs for positional coupling.
    lambda_gid_pos->clear();
    for (unsigned int i_node = 0; i_node < 2; i_node++)
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        lambda_gid_pos->push_back(lambda_total[n_nodal_dof * i_node + i_dim]);
    for (unsigned int i_dof = 6; i_dof < n_mortar_pos; i_dof++)
      lambda_gid_pos->push_back(lambda_total[n_nodal_dof * 2 - 6 + i_dof]);
  }

  if (lambda_gid_rot != nullptr)
  {
    // Get the Lagrange multiplier GIDs for rotational coupling.
    lambda_gid_rot->clear();
    for (unsigned int i_node = 0; i_node < 2; i_node++)
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        lambda_gid_rot->push_back(lambda_total[3 + n_nodal_dof * i_node + i_dim]);
    const unsigned int offset = n_mortar_pos - 2 * n_nodal_dof;
    for (unsigned int i_dof = 6; i_dof < n_mortar_rot; i_dof++)
      lambda_gid_rot->push_back(lambda_total[offset + n_nodal_dof * 2 + i_dof]);
  }
}

/**
 *
 */
void BEAMINTERACTION::GetBeamTriadInterpolationScheme(const DRT::Discretization& discret,
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector, const DRT::Element* ele,
    LARGEROTATIONS::TriadInterpolationLocalRotationVectors<3, double>& triad_interpolation_scheme,
    LARGEROTATIONS::TriadInterpolationLocalRotationVectors<3, double>&
        ref_triad_interpolation_scheme)
{
  // Check that the beam element is a SR beam.
  auto beam_ele = dynamic_cast<const DRT::ELEMENTS::Beam3r*>(ele);
  if (beam_ele == nullptr)
    dserror("GetBeamTriadInterpolationScheme is only implemented for SR beams.");

  // Get the rotations of the beam rotation nodes.
  std::vector<double> beam_displacement_vector_full_double;
  BEAMINTERACTION::UTILS::GetCurrentElementDis(
      discret, beam_ele, displacement_vector, beam_displacement_vector_full_double);

  // Create object for triad interpolation schemes.
  std::vector<LINALG::Matrix<4, 1, double>> nodal_quaternions(3);
  beam_ele->GetNodalTriadsFromFullDispVecOrFromDispTheta<3, double>(
      beam_displacement_vector_full_double, nodal_quaternions);
  triad_interpolation_scheme.Reset(nodal_quaternions);

  std::vector<double> beam_displacement_vector_full_ref(
      beam_displacement_vector_full_double.size(), 0.0);
  beam_ele->GetNodalTriadsFromFullDispVecOrFromDispTheta<3, double>(
      beam_displacement_vector_full_ref, nodal_quaternions);
  ref_triad_interpolation_scheme.Reset(nodal_quaternions);
}

/**
 *
 */
template <typename solid, typename scalar_type>
void BEAMINTERACTION::GetSolidRotationVector(
    const INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling& rot_coupling_type,
    const LINALG::Matrix<3, 1, double>& xi,
    const LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
    const LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
    const LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
    LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element)
{
  switch (rot_coupling_type)
  {
    case INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::polar_decomposition_2d:
      GetSolidRotationVectorPolarDecomposition2D<solid>(
          xi, q_solid_ref, q_solid, quaternion_beam_ref, psi_solid, element);
      return;
    case INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::deformation_gradient_average_2d:
    case INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::deformation_gradient_y_2d:
    case INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::deformation_gradient_z_2d:
      GetSolidRotationVectorDeformationGradient2D<solid>(
          rot_coupling_type, xi, q_solid_ref, q_solid, quaternion_beam_ref, psi_solid, element);
      return;
    case INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::deformation_gradient_3d_local_1:
    case INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::deformation_gradient_3d_local_2:
    case INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::deformation_gradient_3d_local_3:
      GetSolidRotationVectorDeformationGradient3D<solid>(
          rot_coupling_type, xi, q_solid_ref, q_solid, quaternion_beam_ref, psi_solid, element);
      return;
    case INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::deformation_gradient_3d_general:
      GetSolidRotationVectorDeformationGradient3DGeneral<solid>(
          xi, q_solid_ref, q_solid, quaternion_beam_ref, psi_solid, element);
      return;
    case INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::
        deformation_gradient_3d_general_in_cross_section_plane:
      GetSolidRotationVectorDeformationGradient3DGeneralInCrossSectionPlane<solid>(
          xi, q_solid_ref, q_solid, quaternion_beam_ref, psi_solid, element);
      return;
    case INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::deformation_gradient_3d_base_1:
      GetSolidRotationVectorDeformationGradient3DBase1<solid>(
          xi, q_solid_ref, q_solid, quaternion_beam_ref, psi_solid, element);
      return;
    default:
      dserror("Got unexpected rotational coupling type");
      break;
  }
}

/**
 *
 */
template <typename solid, typename scalar_type>
void BEAMINTERACTION::GetSolidRotationVectorDeformationGradient3DGeneral(
    const LINALG::Matrix<3, 1, double>& xi,
    const LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
    const LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
    const LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
    LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element)
{
  // Get basis vectors of reference triad in the current configuration.
  LINALG::Matrix<4, 1, scalar_type> quaternion_beam_ref_fad;
  for (unsigned int i = 0; i < 4; i++) quaternion_beam_ref_fad(i) = quaternion_beam_ref(i);
  LINALG::Matrix<3, 3, scalar_type> ref_triad;
  LARGEROTATIONS::quaterniontotriad(quaternion_beam_ref_fad, ref_triad);
  LINALG::Matrix<3, 3, scalar_type> deformation_gradient;
  GEOMETRYPAIR::EvaluateDeformationGradient<solid>(
      xi, q_solid_ref, q_solid, deformation_gradient, element);
  LINALG::Matrix<3, 3, scalar_type> deformed_basis;
  deformed_basis.Multiply(deformation_gradient, ref_triad);

  // Average of deformed basis vectors.
  LINALG::Matrix<3, 1, scalar_type> average_vector(true);
  for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
    for (unsigned int j_vec = 0; j_vec < 3; j_vec++)
      average_vector(i_dim) += deformed_basis(i_dim, j_vec);
  average_vector.Scale(1.0 / FADUTILS::Norm(average_vector));

  // Project the deformed basis vectors on the plane.
  LINALG::Matrix<3, 1, scalar_type> projected_basis[3];
  LINALG::Matrix<3, 1, scalar_type> temp_vec;
  scalar_type projection;
  for (unsigned int i_basis = 0; i_basis < 3; i_basis++)
  {
    projection = 0.0;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      projection += average_vector(i_dim) * deformed_basis(i_dim, i_basis);

    temp_vec = average_vector;
    temp_vec.Scale(projection);

    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      projected_basis[i_basis](i_dim) = deformed_basis(i_dim, i_basis) - temp_vec(i_dim);

    projected_basis[i_basis].Scale(1.0 / FADUTILS::Norm(projected_basis[i_basis]));
  }

  // Calculate angles between the projected basis vectors.
  scalar_type alpha_21 = std::acos(projected_basis[0].Dot(projected_basis[1]));
  scalar_type alpha_31 = 2.0 * M_PI - acos(projected_basis[0].Dot(projected_basis[2]));

  // Minimum relative angle.
  scalar_type alpha = 1.0 / 3.0 * (alpha_21 + alpha_31 - 2.0 * M_PI);

  // Rotate up the first basis vector.
  LINALG::Matrix<3, 1, scalar_type> rot_vec;
  LINALG::Matrix<4, 1, scalar_type> rot_quat;
  LINALG::Matrix<3, 3, scalar_type> rot_mat;
  LINALG::Matrix<3, 1, scalar_type> start_vec;

  rot_vec.CrossProduct(projected_basis[0], average_vector);
  rot_vec.Scale(0.5 * (M_PI - 2.0 * acos(1.0 / sqrt(3.0))));  // No need to normalize before, should
                                                              // already be length one.
  LARGEROTATIONS::angletoquaternion(rot_vec, rot_quat);
  LARGEROTATIONS::quaterniontotriad(rot_quat, rot_mat);
  start_vec.Multiply(rot_mat, projected_basis[0]);

  // Rotate to the new basis vectors.
  LINALG::Matrix<3, 3, scalar_type> new_basis;
  for (unsigned int i_basis = 0; i_basis < 3; i_basis++)
  {
    rot_vec = average_vector;
    rot_vec.Scale(alpha + i_basis * M_PI * 2.0 / 3.0);
    LARGEROTATIONS::angletoquaternion(rot_vec, rot_quat);
    LARGEROTATIONS::quaterniontotriad(rot_quat, rot_mat);

    temp_vec.Multiply(rot_mat, start_vec);
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++) new_basis(i_dim, i_basis) = temp_vec(i_dim);
  }

  // Get the rotation angle.
  LARGEROTATIONS::triadtoquaternion(new_basis, rot_quat);
  LARGEROTATIONS::quaterniontoangle(rot_quat, psi_solid);
}


/**
 *
 */
template <typename solid, typename scalar_type>
void BEAMINTERACTION::GetSolidRotationVectorDeformationGradient3DGeneralInCrossSectionPlane(
    const LINALG::Matrix<3, 1, double>& xi,
    const LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
    const LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
    const LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
    LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element)
{
  // Get basis vectors of reference triad in the current configuration.
  LINALG::Matrix<3, 3, double> ref_triad;
  LARGEROTATIONS::quaterniontotriad(quaternion_beam_ref, ref_triad);
  LINALG::Matrix<3, 3, scalar_type> deformation_gradient;
  GEOMETRYPAIR::EvaluateDeformationGradient<solid>(
      xi, q_solid_ref, q_solid, deformation_gradient, element);

  // Get the rotation angle.
  GetSolidRotationVectorDeformationGradient3DGeneralInCrossSectionPlane(
      deformation_gradient, ref_triad, psi_solid);
}

/**
 *
 */
template <typename scalar_type>
void BEAMINTERACTION::GetSolidRotationVectorDeformationGradient3DGeneralInCrossSectionPlane(
    const LINALG::Matrix<3, 3, scalar_type>& F, const LINALG::Matrix<3, 3, double>& beam_ref_triad,
    LINALG::Matrix<3, 1, scalar_type>& psi_solid)
{
  // Get the current cross section basis vectors.
  std::array<LINALG::Matrix<3, 1, scalar_type>, 2> cross_section_basis_vector;
  for (unsigned int i_basis = 0; i_basis < 2; i_basis++)
  {
    cross_section_basis_vector[i_basis].PutScalar(0.0);
    for (unsigned int i_row = 0; i_row < 3; i_row++)
      for (unsigned int i_col = 0; i_col < 3; i_col++)
        cross_section_basis_vector[i_basis](i_row) +=
            F(i_row, i_col) * beam_ref_triad(i_col, i_basis + 1);
  }

  // Get the normal on the cross section.
  LINALG::Matrix<3, 1, scalar_type> cross_section_normal_vector;
  cross_section_normal_vector.CrossProduct(
      cross_section_basis_vector[0], cross_section_basis_vector[1]);
  cross_section_normal_vector.Scale(1.0 / FADUTILS::VectorNorm(cross_section_normal_vector));

  // Average the current cross section basis vectors.
  LINALG::Matrix<3, 1, scalar_type> cross_section_average_vector(true);
  for (unsigned int i_basis = 0; i_basis < 2; i_basis++)
  {
    cross_section_basis_vector[i_basis].Scale(
        1.0 / FADUTILS::VectorNorm(cross_section_basis_vector[i_basis]));
    cross_section_average_vector += cross_section_basis_vector[i_basis];
  }
  cross_section_average_vector.Scale(1.0 / FADUTILS::VectorNorm(cross_section_average_vector));

  // Calculate the current solid triad.
  LINALG::Matrix<3, 1, scalar_type> cross_section_third_vector;
  cross_section_third_vector.CrossProduct(
      cross_section_normal_vector, cross_section_average_vector);
  LINALG::Matrix<3, 3, scalar_type> solid_triad_ref;
  for (unsigned int i_row = 0; i_row < 3; i_row++)
  {
    solid_triad_ref(i_row, 0) = cross_section_normal_vector(i_row);
    solid_triad_ref(i_row, 1) = cross_section_average_vector(i_row);
    solid_triad_ref(i_row, 2) = cross_section_third_vector(i_row);
  }

  // Calculate the relative matrix.
  LINALG::Matrix<3, 3, scalar_type> beam_ref_triad_fad;
  for (unsigned int i_row = 0; i_row < 3; i_row++)
    for (unsigned int i_col = 0; i_col < 3; i_col++)
      beam_ref_triad_fad(i_row, i_col) = beam_ref_triad(i_row, i_col);

  LINALG::Matrix<3, 1, scalar_type> rot_vec(true);
  LINALG::Matrix<4, 1, scalar_type> rot_quat(true);
  LINALG::Matrix<3, 3, scalar_type> solid_triad_rel(true);

  rot_vec(0) = -M_PI_4;
  LARGEROTATIONS::angletoquaternion(rot_vec, rot_quat);
  LARGEROTATIONS::quaterniontotriad(rot_quat, solid_triad_rel);

  // Get the rotation angle.
  LINALG::Matrix<3, 3, scalar_type> solid_triad;
  solid_triad.Multiply(solid_triad_ref, solid_triad_rel);
  LARGEROTATIONS::triadtoquaternion(solid_triad, rot_quat);
  LARGEROTATIONS::quaterniontoangle(rot_quat, psi_solid);
}

/**
 *
 */
template <typename solid, typename scalar_type>
void BEAMINTERACTION::GetSolidRotationVectorDeformationGradient3DBase1(
    const LINALG::Matrix<3, 1, double>& xi,
    const LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
    const LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
    const LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
    LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element)
{
  // Get basis vectors of reference triad in the current configuration.
  LINALG::Matrix<4, 1, scalar_type> quaternion_beam_ref_fad;
  for (unsigned int i = 0; i < 4; i++) quaternion_beam_ref_fad(i) = quaternion_beam_ref(i);
  LINALG::Matrix<3, 3, scalar_type> ref_triad;
  LARGEROTATIONS::quaterniontotriad(quaternion_beam_ref_fad, ref_triad);
  LINALG::Matrix<3, 3, scalar_type> deformation_gradient;
  GEOMETRYPAIR::EvaluateDeformationGradient<solid>(
      xi, q_solid_ref, q_solid, deformation_gradient, element);
  LINALG::Matrix<3, 3, scalar_type> deformed_basis;
  deformed_basis.Multiply(deformation_gradient, ref_triad);

  // Average of deformed basis vectors.
  LINALG::Matrix<3, 1, scalar_type> normalized_base_1(true);
  for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
    normalized_base_1(i_dim) = deformed_basis(i_dim, 0);
  normalized_base_1.Scale(1.0 / FADUTILS::Norm(normalized_base_1));

  // Project the deformed basis vectors on the plane.
  LINALG::Matrix<3, 1, scalar_type> projected_basis[2];
  LINALG::Matrix<3, 1, scalar_type> temp_vec;
  scalar_type projection;
  for (unsigned int i_basis = 1; i_basis < 3; i_basis++)
  {
    projection = 0.0;
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      projection += normalized_base_1(i_dim) * deformed_basis(i_dim, i_basis);

    temp_vec = normalized_base_1;
    temp_vec.Scale(projection);

    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      projected_basis[i_basis - 1](i_dim) = deformed_basis(i_dim, i_basis) - temp_vec(i_dim);

    projected_basis[i_basis - 1].Scale(1.0 / FADUTILS::Norm(projected_basis[i_basis - 1]));
  }

  // Calculate angles between the projected basis vectors.
  scalar_type alpha_32 = acos(projected_basis[0].Dot(projected_basis[1]));

  // Rotation angle for base 2.
  scalar_type alpha = 0.5 * (alpha_32 - 0.5 * M_PI);

  // Construct the new basis.
  // Rotate to the new basis vectors.
  LINALG::Matrix<3, 3, scalar_type> new_basis;
  for (unsigned int i_dim = 0; i_dim < 3; i_dim++) new_basis(i_dim, 0) = normalized_base_1(i_dim);

  LINALG::Matrix<3, 1, scalar_type> rot_vec;
  LINALG::Matrix<4, 1, scalar_type> rot_quat;
  LINALG::Matrix<3, 3, scalar_type> rot_mat;
  for (unsigned int i_basis = 0; i_basis < 2; i_basis++)
  {
    rot_vec = normalized_base_1;
    rot_vec.Scale(alpha + i_basis * M_PI * 0.5);
    LARGEROTATIONS::angletoquaternion(rot_vec, rot_quat);
    LARGEROTATIONS::quaterniontotriad(rot_quat, rot_mat);

    temp_vec.Multiply(rot_mat, projected_basis[0]);
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      new_basis(i_dim, i_basis + 1) = temp_vec(i_dim);
  }

  // Get the rotation angle.
  LARGEROTATIONS::triadtoquaternion(new_basis, rot_quat);
  LARGEROTATIONS::quaterniontoangle(rot_quat, psi_solid);
}

/**
 *
 */
template <typename solid, typename scalar_type>
void BEAMINTERACTION::GetSolidRotationVectorDeformationGradient3D(
    const INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling& rot_coupling_type,
    const LINALG::Matrix<3, 1, double>& xi,
    const LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
    const LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
    const LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
    LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element)
{
  // Get basis vectors of reference triad in the current configuration.
  LINALG::Matrix<4, 1, scalar_type> quaternion_beam_ref_fad;
  for (unsigned int i = 0; i < 4; i++) quaternion_beam_ref_fad(i) = quaternion_beam_ref(i);
  LINALG::Matrix<3, 3, scalar_type> ref_triad;
  LARGEROTATIONS::quaterniontotriad(quaternion_beam_ref_fad, ref_triad);
  LINALG::Matrix<3, 3, scalar_type> deformation_gradient;
  GEOMETRYPAIR::EvaluateDeformationGradient<solid>(
      xi, q_solid_ref, q_solid, deformation_gradient, element);
  LINALG::Matrix<3, 3, scalar_type> deformed_basis;
  deformed_basis.Multiply(deformation_gradient, ref_triad);

  // Get the order of the basis vector to be used for the construction.
  unsigned int local_basis_vector_construction_order[2];
  unsigned int constructed_basis_vector_to_triad_order[3];
  switch (rot_coupling_type)
  {
    case INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::deformation_gradient_3d_local_1:
      local_basis_vector_construction_order[0] = 0;
      local_basis_vector_construction_order[1] = 2;
      constructed_basis_vector_to_triad_order[0] = 0;
      constructed_basis_vector_to_triad_order[1] = 1;
      constructed_basis_vector_to_triad_order[2] = 2;
      break;
    case INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::deformation_gradient_3d_local_2:
      local_basis_vector_construction_order[0] = 1;
      local_basis_vector_construction_order[1] = 0;
      constructed_basis_vector_to_triad_order[0] = 1;
      constructed_basis_vector_to_triad_order[1] = 2;
      constructed_basis_vector_to_triad_order[2] = 0;
      break;
    case INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::deformation_gradient_3d_local_3:
      local_basis_vector_construction_order[0] = 2;
      local_basis_vector_construction_order[1] = 1;
      constructed_basis_vector_to_triad_order[0] = 2;
      constructed_basis_vector_to_triad_order[1] = 0;
      constructed_basis_vector_to_triad_order[2] = 1;
      break;
    default:
      dserror("Unexpected coupling type for GetSolidRotationVectorDeformationGradient3D");
  }

  // Basis vectors and the triad of the solid.
  LINALG::Matrix<3, 1, scalar_type> construction_vector;
  LINALG::Matrix<3, 1, scalar_type> triad_basis_vectors[3];
  LINALG::Matrix<3, 3, scalar_type> solid_triad;

  // Set the first basis vector.
  for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
    triad_basis_vectors[0](i_dim) = deformed_basis(i_dim, local_basis_vector_construction_order[0]);

  // Construct the second basis vector.
  for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
    construction_vector(i_dim) = deformed_basis(i_dim, local_basis_vector_construction_order[1]);
  triad_basis_vectors[1].CrossProduct(construction_vector, triad_basis_vectors[0]);

  // Construct the third basis vector.
  triad_basis_vectors[2].CrossProduct(triad_basis_vectors[0], triad_basis_vectors[1]);

  // Norm all vectors and add them to the solid triad.
  scalar_type norm;
  for (unsigned int i_basis = 0; i_basis < 3; i_basis++)
  {
    triad_basis_vectors[i_basis].Scale(1.0 / FADUTILS::Norm(triad_basis_vectors[i_basis]));
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      solid_triad(i_dim, constructed_basis_vector_to_triad_order[i_basis]) =
          triad_basis_vectors[i_basis](i_dim);
  }

  // Convert the triad into a rotation vector.
  LINALG::Matrix<4, 1, scalar_type> solid_quaternion;
  LARGEROTATIONS::triadtoquaternion(solid_triad, solid_quaternion);
  LARGEROTATIONS::quaterniontoangle(solid_quaternion, psi_solid);
}

/**
 *
 */
template <typename solid, typename scalar_type>
void BEAMINTERACTION::GetSolidRotationVectorPolarDecomposition2D(
    const LINALG::Matrix<3, 1, double>& xi,
    const LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
    const LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
    const LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
    LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element)
{
  // Get the deformation gradient in the solid.
  LINALG::Matrix<3, 3, scalar_type> deformation_gradient;
  GEOMETRYPAIR::EvaluateDeformationGradient<solid>(
      xi, q_solid_ref, q_solid, deformation_gradient, element);

  // Check that the assumption of plane rotations is full filled.
  CheckPlaneRotations(deformation_gradient, quaternion_beam_ref);

  // Reference rotation of beam cross-section in plane.
  LINALG::Matrix<3, 1, double> beam_ref_psi;
  LARGEROTATIONS::quaterniontoangle(quaternion_beam_ref, beam_ref_psi);
  double reference_rotation_beam = FADUTILS::VectorNorm(beam_ref_psi);

  // Perform a polar decomposition of the 2D deformation gradient.
  scalar_type solid_angle = 0.0;
  {
    LINALG::Matrix<2, 2, scalar_type> F, R, U, U_times_U;
    for (unsigned int dim_y = 0; dim_y < 2; dim_y++)
      for (unsigned int dim_z = 0; dim_z < 2; dim_z++)
        F(dim_y, dim_z) = deformation_gradient(dim_y + 1, dim_z + 1);

    // Compute U*U.
    U_times_U.MultiplyTN(F, F);

    // We have to calculate the square root of the matrix U*U here.
    U.Clear();
    if (abs(FADUTILS::CastToDouble(U_times_U(0, 0) - U_times_U(1, 1))) < 1e-10 and
        abs(FADUTILS::CastToDouble(U_times_U(0, 1))) < 1e-10)
    {
      U(0, 0) = FADUTILS::sqrt(U_times_U(0, 0));
      U(1, 1) = FADUTILS::sqrt(U_times_U(1, 1));
    }
    else
    {
      scalar_type v[16];
      v[15] = U_times_U(0, 0) - U_times_U(1, 1);
      v[5] = 4.0 * std::pow(U_times_U(0, 1), 2.0) + (v[15] * v[15]);
      v[12] = 1.0 / FADUTILS::sqrt(v[5]);
      v[6] = FADUTILS::sqrt(v[5]);
      v[7] = 0.3535533905932738e0 * v[12];
      v[8] = FADUTILS::sqrt<scalar_type>(U_times_U(0, 0) + U_times_U(1, 1) - v[6]);
      v[9] = -U_times_U(0, 0) + U_times_U(1, 1) + v[6];
      v[10] = v[15] + v[6];
      v[11] = FADUTILS::sqrt<scalar_type>(2.0 * U_times_U(0, 0) + v[9]);
      v[13] = -0.7071067811865476e0 * U_times_U(0, 1) * v[12] * (-v[11] + v[8]);
      U(0, 0) = v[7] * (v[10] * v[11] + v[8] * v[9]);
      U(0, 1) = v[13];
      U(1, 0) = v[13];
      U(1, 1) = v[7] * (v[10] * v[8] + v[11] * v[9]);
    }

    // Compute R.
    Inverse(U);
    R.Multiply(F, U);
    solid_angle = atan2(-R(0, 1), R(0, 0));
  }

  // Return the solid rotation vector.
  psi_solid.PutScalar(0.0);
  psi_solid(0) = solid_angle + reference_rotation_beam;
}

/**
 *
 */
template <typename solid, typename scalar_type>
void BEAMINTERACTION::GetSolidRotationVectorDeformationGradient2D(
    const INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling& rot_coupling_type,
    const LINALG::Matrix<3, 1, double>& xi,
    const LINALG::Matrix<solid::n_dof_, 1, double>& q_solid_ref,
    const LINALG::Matrix<solid::n_dof_, 1, scalar_type>& q_solid,
    const LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
    LINALG::Matrix<3, 1, scalar_type>& psi_solid, const DRT::Element* element)
{
  // Get the deformation gradient in the solid.
  LINALG::Matrix<3, 3, scalar_type> deformation_gradient;
  GEOMETRYPAIR::EvaluateDeformationGradient<solid>(
      xi, q_solid_ref, q_solid, deformation_gradient, element);

  // Check that the assumption of plane rotations is full filled.
  CheckPlaneRotations(deformation_gradient, quaternion_beam_ref);

  // Reference rotation of beam cross-section in plane.
  LINALG::Matrix<3, 1, double> beam_ref_psi;
  LARGEROTATIONS::quaterniontoangle(quaternion_beam_ref, beam_ref_psi);
  double reference_rotation_beam = FADUTILS::VectorNorm(beam_ref_psi);

  // Get the rotation of the solid.
  scalar_type angle;
  switch (rot_coupling_type)
  {
    case INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::deformation_gradient_y_2d:
      angle = atan2(deformation_gradient(2, 1), deformation_gradient(1, 1));
      break;
    case INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::deformation_gradient_z_2d:
      angle = atan2(deformation_gradient(2, 2), deformation_gradient(1, 2)) - M_PI * 0.5;
      break;
    case INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling::deformation_gradient_average_2d:
      angle = 0.5 * (atan2(deformation_gradient(2, 1), deformation_gradient(1, 1)) +
                        atan2(deformation_gradient(2, 2), deformation_gradient(1, 2))) -
              M_PI * 0.25;
      break;
    default:
      dserror("Unexpected coupling type for GetSolidRotationVectorDeformationGradient2D");
      break;
  }
  psi_solid.PutScalar(0.0);
  psi_solid(0) = reference_rotation_beam + angle;
}

/**
 *
 */
template <typename scalar_type>
void BEAMINTERACTION::CheckPlaneRotations(
    const LINALG::Matrix<3, 3, scalar_type> deformation_gradient,
    const LINALG::Matrix<4, 1, double>& quaternion_beam_ref)
{
  // Check that the solid as well as the reference beam rotations are plane.
  const double tol = 1e-10;
  double out_of_plane_values = 0.0;
  for (unsigned int i = 0; i < 2; i++)
  {
    out_of_plane_values += pow(FADUTILS::CastToDouble(deformation_gradient)(i + 1, 0), 2.0);
    out_of_plane_values += pow(FADUTILS::CastToDouble(deformation_gradient)(0, i + 1), 2.0);
  }
  if (FADUTILS::sqrt(out_of_plane_values) > tol)
    dserror("The solid deformation is not just in plane. Out of plane value: %f, tolerance: %f",
        FADUTILS::sqrt(out_of_plane_values), tol);
  LINALG::Matrix<2, 1, double> projection_on_x;
  LINALG::Matrix<3, 1, double> beam_ref_psi;
  LARGEROTATIONS::quaterniontoangle(quaternion_beam_ref, beam_ref_psi);
  for (unsigned int i = 0; i < 2; i++) projection_on_x(i) = beam_ref_psi(i + 1);
  if (FADUTILS::VectorNorm(projection_on_x) > tol)
    dserror("The beam reference rotation is not just in plane. Projection value: %f, tolerance: %f",
        FADUTILS::VectorNorm(projection_on_x), tol);
}

/**
 *
 */
template <typename beam, typename other, typename mortar>
void BEAMINTERACTION::AssembleLocalMortarContributions(const BEAMINTERACTION::BeamContactPair* pair,
    const DRT::Discretization& discret, const BeamToSolidMortarManager* mortar_manager,
    LINALG::SparseMatrix& global_G_B, LINALG::SparseMatrix& global_G_S,
    LINALG::SparseMatrix& global_FB_L, LINALG::SparseMatrix& global_FS_L,
    Epetra_FEVector& global_constraint, Epetra_FEVector& global_kappa,
    Epetra_FEVector& global_lambda_active,
    const LINALG::Matrix<mortar::n_dof_, beam::n_dof_, double>& local_D,
    const LINALG::Matrix<mortar::n_dof_, other::n_dof_, double>& local_M,
    const LINALG::Matrix<mortar::n_dof_, 1, double>& local_kappa,
    const LINALG::Matrix<mortar::n_dof_, 1, double>& local_constraint,
    const unsigned int n_mortar_rot)
{
  // Get the GIDs of the Lagrange multipliers.
  std::vector<int> lambda_row;
  GetMortarGID(mortar_manager, pair, mortar::n_dof_, n_mortar_rot, &lambda_row, nullptr);

  // Get the beam centerline GIDs.
  LINALG::Matrix<beam::n_dof_, 1, int> beam_centerline_gid;
  UTILS::GetElementCenterlineGIDIndices(discret, pair->Element1(), beam_centerline_gid);

  // Get the other GIDs.
  // We call this function on the element pointer of the geometry pair, since for face elements,
  // the element pointer of the beam contact pair is to the volume element and only the element
  // pointer of the geometry pair is to the face element.
  std::vector<int> other_row;
  std::vector<int> dummy_1;
  std::vector<int> dummy_2;
  pair->GeometryPair()->Element2()->LocationVector(discret, other_row, dummy_1, dummy_2);

  // Assemble into the global matrices. All contributions here are assumed to be symmetric.
  for (unsigned int i_lambda = 0; i_lambda < mortar::n_dof_; ++i_lambda)
  {
    for (unsigned int i_beam = 0; i_beam < beam::n_dof_; ++i_beam)
    {
      global_G_B.FEAssemble(
          local_D(i_lambda, i_beam), lambda_row[i_lambda], beam_centerline_gid(i_beam));
      global_FB_L.FEAssemble(
          local_D(i_lambda, i_beam), beam_centerline_gid(i_beam), lambda_row[i_lambda]);
    }
    for (unsigned int i_other = 0; i_other < other::n_dof_; ++i_other)
    {
      global_G_S.FEAssemble(-local_M(i_lambda, i_other), lambda_row[i_lambda], other_row[i_other]);
      global_FS_L.FEAssemble(-local_M(i_lambda, i_other), other_row[i_other], lambda_row[i_lambda]);
    }
  }
  global_kappa.SumIntoGlobalValues(mortar::n_dof_, &lambda_row[0], local_kappa.A());
  global_constraint.SumIntoGlobalValues(mortar::n_dof_, &lambda_row[0], local_constraint.A());

  // Set all entries in the local kappa vector to 1 and add them to the active vector.
  LINALG::Matrix<mortar::n_dof_, 1, double> local_kappa_active;
  local_kappa_active.PutScalar(1.0);
  global_lambda_active.SumIntoGlobalValues(mortar::n_dof_, &lambda_row[0], local_kappa_active.A());
}


/**
 * Explicit template initialization of template functions.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  // Helper types for the macro initialization. The compiler has troubles inserting the templated
  // typenames into the macros.
  using line_to_surface_patch_scalar_type_fixed_size_1st_order_nurbs_9 =
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_hermite, t_nurbs9>;
  using line_to_surface_patch_scalar_type_fixed_size_nurbs_9 =
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>;

#define initialize_template_penalty(scalar_type)                                       \
  template scalar_type PenaltyForce<scalar_type>(                                      \
      const scalar_type&, const Teuchos::RCP<const BeamToSolidSurfaceContactParams>&); \
  template scalar_type PenaltyPotential<scalar_type>(                                  \
      const scalar_type&, const Teuchos::RCP<const BeamToSolidSurfaceContactParams>&);

  initialize_template_penalty(double);
  initialize_template_penalty(line_to_surface_patch_scalar_type_1st_order);
  initialize_template_penalty(line_to_surface_patch_scalar_type_fixed_size_1st_order_nurbs_9);
  initialize_template_penalty(line_to_surface_patch_scalar_type);
  initialize_template_penalty(line_to_surface_patch_scalar_type_fixed_size_nurbs_9);


#define initialize_template_get_solid_rotation_vector(a, fad_order)                                \
  template void GetSolidRotationVector<a,                                                          \
      FADUTILS::HigherOrderFadType<fad_order, Sacado::Fad::SLFad<double, 3 + a::n_dof_>>::type>(   \
      const INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling&, const LINALG::Matrix<3, 1, double>&, \
      const LINALG::Matrix<a::n_dof_, 1, double>&,                                                 \
      const LINALG::Matrix<a::n_dof_, 1,                                                           \
          FADUTILS::HigherOrderFadType<fad_order,                                                  \
              Sacado::Fad::SLFad<double, 3 + a::n_dof_>>::type>&,                                  \
      const LINALG::Matrix<4, 1, double>&,                                                         \
      LINALG::Matrix<3, 1,                                                                         \
          FADUTILS::HigherOrderFadType<fad_order,                                                  \
              Sacado::Fad::SLFad<double, 3 + a::n_dof_>>::type>&,                                  \
      const DRT::Element* element);

  initialize_template_get_solid_rotation_vector(t_hex8, 1);
  initialize_template_get_solid_rotation_vector(t_hex20, 1);
  initialize_template_get_solid_rotation_vector(t_hex27, 1);
  initialize_template_get_solid_rotation_vector(t_tet4, 1);
  initialize_template_get_solid_rotation_vector(t_tet10, 1);
  initialize_template_get_solid_rotation_vector(t_nurbs27, 1);

  initialize_template_get_solid_rotation_vector(t_hex8, 2);
  initialize_template_get_solid_rotation_vector(t_hex20, 2);
  initialize_template_get_solid_rotation_vector(t_hex27, 2);
  initialize_template_get_solid_rotation_vector(t_tet4, 2);
  initialize_template_get_solid_rotation_vector(t_tet10, 2);
  initialize_template_get_solid_rotation_vector(t_nurbs27, 2);

#define initialize_template_get_surface_rotation_vector(a, fad_order)                            \
  template void GetSolidRotationVectorDeformationGradient3DGeneralInCrossSectionPlane<           \
      FADUTILS::HigherOrderFadType<fad_order, Sacado::Fad::SLFad<double, 3 + a::n_dof_>>::type>( \
      const LINALG::Matrix<3, 3,                                                                 \
          FADUTILS::HigherOrderFadType<fad_order,                                                \
              Sacado::Fad::SLFad<double, 3 + a::n_dof_>>::type>&,                                \
      const LINALG::Matrix<3, 3, double>&,                                                       \
      LINALG::Matrix<3, 1,                                                                       \
          FADUTILS::HigherOrderFadType<fad_order,                                                \
              Sacado::Fad::SLFad<double, 3 + a::n_dof_>>::type>&);

  initialize_template_get_surface_rotation_vector(t_quad4, 1);
  initialize_template_get_surface_rotation_vector(t_quad8, 1);
  initialize_template_get_surface_rotation_vector(t_quad9, 1);  // Also initializes nurbs9
  initialize_template_get_surface_rotation_vector(t_tri3, 1);
  initialize_template_get_surface_rotation_vector(t_tri6, 1);

  initialize_template_get_surface_rotation_vector(t_quad4, 2);
  initialize_template_get_surface_rotation_vector(t_quad8, 2);
  initialize_template_get_surface_rotation_vector(t_quad9, 2);  // Also initializes nurbs9
  initialize_template_get_surface_rotation_vector(t_tri3, 2);
  initialize_template_get_surface_rotation_vector(t_tri6, 2);

#define initialize_template_assemble_local_mortar_contributions(beam, other, mortar)    \
  template void AssembleLocalMortarContributions<beam, other, mortar>(                  \
      const BEAMINTERACTION::BeamContactPair*, const DRT::Discretization&,              \
      const BeamToSolidMortarManager*, LINALG::SparseMatrix&, LINALG::SparseMatrix&,    \
      LINALG::SparseMatrix&, LINALG::SparseMatrix&, Epetra_FEVector&, Epetra_FEVector&, \
      Epetra_FEVector&, const LINALG::Matrix<mortar::n_dof_, beam::n_dof_, double>&,    \
      const LINALG::Matrix<mortar::n_dof_, other::n_dof_, double>&,                     \
      const LINALG::Matrix<mortar::n_dof_, 1, double>&,                                 \
      const LINALG::Matrix<mortar::n_dof_, 1, double>&, const unsigned int);

  initialize_template_assemble_local_mortar_contributions(t_hermite, t_hex8, t_line2);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_hex8, t_line3);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_hex8, t_line4);

  initialize_template_assemble_local_mortar_contributions(t_hermite, t_hex20, t_line2);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_hex20, t_line3);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_hex20, t_line4);

  initialize_template_assemble_local_mortar_contributions(t_hermite, t_hex27, t_line2);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_hex27, t_line3);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_hex27, t_line4);

  initialize_template_assemble_local_mortar_contributions(t_hermite, t_tet4, t_line2);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_tet4, t_line3);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_tet4, t_line4);

  initialize_template_assemble_local_mortar_contributions(t_hermite, t_tet10, t_line2);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_tet10, t_line3);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_tet10, t_line4);

  initialize_template_assemble_local_mortar_contributions(t_hermite, t_nurbs27, t_line2);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_nurbs27, t_line3);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_nurbs27, t_line4);

  initialize_template_assemble_local_mortar_contributions(t_hermite, t_quad4, t_line2);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_quad4, t_line3);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_quad4, t_line4);

  initialize_template_assemble_local_mortar_contributions(t_hermite, t_quad8, t_line2);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_quad8, t_line3);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_quad8, t_line4);

  initialize_template_assemble_local_mortar_contributions(t_hermite, t_quad9, t_line2);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_quad9, t_line3);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_quad9, t_line4);

  initialize_template_assemble_local_mortar_contributions(t_hermite, t_tri3, t_line2);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_tri3, t_line3);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_tri3, t_line4);

  initialize_template_assemble_local_mortar_contributions(t_hermite, t_tri6, t_line2);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_tri6, t_line3);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_tri6, t_line4);

  initialize_template_assemble_local_mortar_contributions(t_hermite, t_nurbs9, t_line2);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_nurbs9, t_line3);
  initialize_template_assemble_local_mortar_contributions(t_hermite, t_nurbs9, t_line4);
}  // namespace BEAMINTERACTION
