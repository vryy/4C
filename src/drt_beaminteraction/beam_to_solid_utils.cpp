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
#include "../drt_geometry_pair/geometry_pair.H"
#include "../drt_geometry_pair/geometry_pair_element.H"
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
    default:
      dserror("Got unexpected rotational coupling type");
      break;
  }
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
  LINALG::Matrix<3, 1, scalar_type> t[3];
  LINALG::Matrix<3, 3, scalar_type> solid_triad;

  // Set the first basis vector.
  for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
    t[0](i_dim) = deformed_basis(i_dim, local_basis_vector_construction_order[0]);

  // Construct the second basis vector.
  for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
    construction_vector(i_dim) = deformed_basis(i_dim, local_basis_vector_construction_order[1]);
  t[1].CrossProduct(construction_vector, t[0]);

  // Construct the third basis vector.
  t[2].CrossProduct(t[0], t[1]);

  // Norm all vectors and add them to the solid triad.
  scalar_type norm;
  for (unsigned int i_basis = 0; i_basis < 3; i_basis++)
  {
    t[i_basis].Scale(1.0 / FADUTILS::Norm(t[i_basis]));
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
      solid_triad(i_dim, constructed_basis_vector_to_triad_order[i_basis]) = t[i_basis](i_dim);
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
      v[5] = 4e0 * pow(U_times_U(0, 1), 2.0) + (v[15] * v[15]);
      v[12] = 1e0 / FADUTILS::sqrt(v[5]);
      v[6] = FADUTILS::sqrt(v[5]);
      v[7] = 0.3535533905932738e0 * v[12];
      v[8] = FADUTILS::sqrt<scalar_type>(U_times_U(0, 0) + U_times_U(1, 1) - v[6]);
      v[9] = -U_times_U(0, 0) + U_times_U(1, 1) + v[6];
      v[10] = v[15] + v[6];
      v[11] = FADUTILS::sqrt<scalar_type>(2e0 * U_times_U(0, 0) + v[9]);
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
    LINALG::SparseMatrix& global_GB, LINALG::SparseMatrix& global_GS,
    LINALG::SparseMatrix& global_FB, LINALG::SparseMatrix& global_FS,
    Epetra_FEVector& global_constraint, Epetra_FEVector& global_kappa,
    Epetra_FEVector& global_lambda_active,
    const LINALG::Matrix<mortar::n_dof_, beam::n_dof_, double>& local_D,
    const LINALG::Matrix<mortar::n_dof_, other::n_dof_, double>& local_M,
    const LINALG::Matrix<mortar::n_dof_, 1, double>& local_kappa,
    const LINALG::Matrix<mortar::n_dof_, 1, double>& local_constraint)
{
  // Get the GIDs of the Lagrange multipliers.
  std::vector<int> lambda_row;
  mortar_manager->LocationVector(pair, lambda_row);

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
      global_GB.FEAssemble(
          local_D(i_lambda, i_beam), lambda_row[i_lambda], beam_centerline_gid(i_beam));
      global_FB.FEAssemble(
          local_D(i_lambda, i_beam), beam_centerline_gid(i_beam), lambda_row[i_lambda]);
    }
    for (unsigned int i_other = 0; i_other < other::n_dof_; ++i_other)
    {
      global_GS.FEAssemble(-local_M(i_lambda, i_other), lambda_row[i_lambda], other_row[i_other]);
      global_FS.FEAssemble(-local_M(i_lambda, i_other), other_row[i_other], lambda_row[i_lambda]);
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

#define initialize_template_get_solid_rotation_vector(a)                                           \
  template void GetSolidRotationVector<a,                                                          \
      FADUTILS::HigherOrderFadType<2, Sacado::Fad::SLFad<double, 3 + a::n_dof_>>::type>(           \
      const INPAR::BEAMTOSOLID::BeamToSolidRotationCoupling&, const LINALG::Matrix<3, 1, double>&, \
      const LINALG::Matrix<a::n_dof_, 1, double>&,                                                 \
      const LINALG::Matrix<a::n_dof_, 1,                                                           \
          FADUTILS::HigherOrderFadType<2, Sacado::Fad::SLFad<double, 3 + a::n_dof_>>::type>&,      \
      const LINALG::Matrix<4, 1, double>&,                                                         \
      LINALG::Matrix<3, 1,                                                                         \
          FADUTILS::HigherOrderFadType<2, Sacado::Fad::SLFad<double, 3 + a::n_dof_>>::type>&,      \
      const DRT::Element* element);

  initialize_template_get_solid_rotation_vector(t_hex8);
  initialize_template_get_solid_rotation_vector(t_hex20);
  initialize_template_get_solid_rotation_vector(t_hex27);
  initialize_template_get_solid_rotation_vector(t_tet4);
  initialize_template_get_solid_rotation_vector(t_tet10);
  initialize_template_get_solid_rotation_vector(t_nurbs27);

#define initialize_template_assemble_local_mortar_contributions(beam, other, mortar)    \
  template void AssembleLocalMortarContributions<beam, other, mortar>(                  \
      const BEAMINTERACTION::BeamContactPair*, const DRT::Discretization&,              \
      const BeamToSolidMortarManager*, LINALG::SparseMatrix&, LINALG::SparseMatrix&,    \
      LINALG::SparseMatrix&, LINALG::SparseMatrix&, Epetra_FEVector&, Epetra_FEVector&, \
      Epetra_FEVector&, const LINALG::Matrix<mortar::n_dof_, beam::n_dof_, double>&,    \
      const LINALG::Matrix<mortar::n_dof_, other::n_dof_, double>&,                     \
      const LINALG::Matrix<mortar::n_dof_, 1, double>&,                                 \
      const LINALG::Matrix<mortar::n_dof_, 1, double>&);

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
