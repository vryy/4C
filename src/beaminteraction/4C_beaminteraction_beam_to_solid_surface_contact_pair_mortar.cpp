/*----------------------------------------------------------------------*/
/*! \file

\brief Contact element for contact between a 3D beam and a surface element.

\level 3
*/


#include "4C_beaminteraction_beam_to_solid_surface_contact_pair_mortar.hpp"

#include "4C_beam3_base.hpp"
#include "4C_beaminteraction_beam_to_solid_mortar_manager_contact.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_contact_params.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_visualization_output_params.hpp"
#include "4C_beaminteraction_beam_to_solid_utils.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_base.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"
#include "4C_beaminteraction_calc_utils.hpp"
#include "4C_beaminteraction_contact_params.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_element_faces.hpp"
#include "4C_geometry_pair_element_shape_functions.hpp"
#include "4C_geometry_pair_factory.hpp"
#include "4C_geometry_pair_line_to_surface.hpp"
#include "4C_geometry_pair_scalar_types.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_fad.hpp"

#include <Epetra_FEVector.h>
#include <Teuchos_RCPDecl.hpp>

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
BEAMINTERACTION::BeamToSolidSurfaceContactPairMortar<scalar_type, beam, surface,
    mortar>::BeamToSolidSurfaceContactPairMortar()
    : base_class()
{
  // Empty constructor.
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairMortar<scalar_type, beam, surface,
    mortar>::evaluate_and_assemble_mortar_contributions(const Core::FE::Discretization& discret,
    const BeamToSolidMortarManager* mortar_manager,
    Core::LinAlg::SparseMatrix& global_constraint_lin_beam,
    Core::LinAlg::SparseMatrix& global_constraint_lin_solid,
    Core::LinAlg::SparseMatrix& global_force_beam_lin_lambda,
    Core::LinAlg::SparseMatrix& global_force_solid_lin_lambda, Epetra_FEVector& global_constraint,
    Epetra_FEVector& global_kappa, Core::LinAlg::SparseMatrix& global_kappa_lin_beam,
    Core::LinAlg::SparseMatrix& global_kappa_lin_solid, Epetra_FEVector& global_lambda_active,
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector)
{
  // Call Evaluate on the geometry Pair
  this->cast_geometry_pair()->evaluate(
      this->ele1pos_, this->face_element_->GetFaceElementData(), this->line_to_3D_segments_);

  // If there are no intersection segments, no contact terms will be assembled
  const unsigned int n_segments = this->line_to_3D_segments_.size();
  if (n_segments == 0) return;

  // Pointer to the contact parameters and input parameters
  const auto contact_parameters = this->Params()->beam_to_solid_surface_contact_params();
  const auto contact_defined_on =
      contact_parameters->get_beam_to_solid_surface_contact_mortar_defined_in();

  // Get beam cross-section diameter
  const auto beam_ptr = dynamic_cast<const Discret::ELEMENTS::Beam3Base*>(this->Element1());
  const double beam_cross_section_radius =
      beam_ptr->get_circular_cross_section_radius_for_interactions();

  // Initialize variables for contact kinematics
  Core::LinAlg::Matrix<3, 1, scalar_type> surface_normal;
  Core::LinAlg::Matrix<3, 1, scalar_type> r_beam;
  Core::LinAlg::Matrix<3, 1, scalar_type> r_surface;
  Core::LinAlg::Matrix<3, 1, scalar_type> r_rel;
  scalar_type gap = 0.0;

  // Initialize variables for the shape function matrices
  Core::LinAlg::Matrix<mortar_trial::spatial_dim_, mortar_trial::n_dof_, scalar_type>
      N_lambda_trial;
  Core::LinAlg::Matrix<mortar::spatial_dim_, mortar::n_dof_, scalar_type> N_lambda;
  Core::LinAlg::Matrix<beam::spatial_dim_, beam::n_dof_, scalar_type> N_beam;
  Core::LinAlg::Matrix<surface::spatial_dim_, surface::n_dof_, scalar_type> N_surface;

  // Initialize variables for local vectors
  Core::LinAlg::Matrix<mortar_trial::n_dof_, 1, scalar_type> constraint_vector(true);
  Core::LinAlg::Matrix<mortar_trial::n_dof_, 1, scalar_type> kappa(true);
  Core::LinAlg::Matrix<3, mortar::n_dof_, scalar_type> normal_times_lambda_shape(true);
  Core::LinAlg::Matrix<beam::n_dof_, mortar::n_dof_, scalar_type>
      beam_shape_times_normal_times_lambda_shape_gp(true);
  Core::LinAlg::Matrix<surface::n_dof_, mortar::n_dof_, scalar_type>
      surface_shape_times_normal_times_lambda_shape_gp(true);
  beam_shape_times_normal_times_lambda_shape_.put_scalar(0.0);
  surface_shape_times_normal_times_lambda_shape_.put_scalar(0.0);

  // Integrate over segments
  for (const auto& segment : this->line_to_3D_segments_)
  {
    // Gauss point loop
    for (const auto& projected_gauss_point : segment.GetProjectionPoints())
    {
      // Get the projection coordinates
      const auto& xi = projected_gauss_point.GetXi();
      const auto& eta = projected_gauss_point.GetEta();

      // Get the current Gauss integration factor. This includes everything, e.g., Gauss weight,
      // segment Jacobian and beam Jacobian.
      const scalar_type gauss_factor = projected_gauss_point.GetGaussWeight() *  //
                                       0.5 * segment.GetSegmentLength() *
                                       get_jacobian_for_configuration(eta, contact_defined_on);

      // Get the surface normal vector
      GEOMETRYPAIR::EvaluateSurfaceNormal<surface>(
          xi, this->face_element_->GetFaceElementData(), surface_normal);

      // Evaluate the current position of beam and solid
      GEOMETRYPAIR::EvaluatePosition<beam>(eta, this->ele1pos_, r_beam);
      GEOMETRYPAIR::EvaluatePosition<surface>(
          xi, this->face_element_->GetFaceElementData(), r_surface);

      // Evaluate the gap function
      r_rel = r_beam;
      r_rel -= r_surface;
      gap = r_rel.dot(surface_normal) - beam_cross_section_radius;

      // Get the shape function matrices
      GEOMETRYPAIR::EvaluateShapeFunctionMatrix<mortar_trial>(N_lambda_trial, eta);
      GEOMETRYPAIR::EvaluateShapeFunctionMatrix<mortar>(N_lambda, eta);
      GEOMETRYPAIR::EvaluateShapeFunctionMatrix<beam>(
          N_beam, eta, this->ele1pos_.shape_function_data_);
      GEOMETRYPAIR::EvaluateShapeFunctionMatrix<surface>(
          N_surface, xi, this->face_element_->GetFaceElementData().shape_function_data_);

      // Weighted gap
      constraint_vector.update_t(gauss_factor * gap, N_lambda_trial, 1.0);

      // Force contributions
      normal_times_lambda_shape.multiply(surface_normal, N_lambda);

      beam_shape_times_normal_times_lambda_shape_gp.multiply_tn(N_beam, normal_times_lambda_shape);
      beam_shape_times_normal_times_lambda_shape_gp.scale(1.0 * gauss_factor);
      beam_shape_times_normal_times_lambda_shape_ += beam_shape_times_normal_times_lambda_shape_gp;

      surface_shape_times_normal_times_lambda_shape_gp.multiply_tn(
          N_surface, normal_times_lambda_shape);
      surface_shape_times_normal_times_lambda_shape_gp.scale(-1.0 * gauss_factor);
      surface_shape_times_normal_times_lambda_shape_ +=
          surface_shape_times_normal_times_lambda_shape_gp;

      // Scaling vector
      Core::LinAlg::Matrix<mortar_trial::spatial_dim_, 1, double> ones(true);
      ones.put_scalar(1.0);
      Core::LinAlg::Matrix<mortar_trial::n_dof_, 1, scalar_type> N_lambda_trial_flat(true);
      N_lambda_trial_flat.multiply_tn(N_lambda_trial, ones);
      N_lambda_trial_flat.scale(gauss_factor);
      kappa += N_lambda_trial_flat;
    }
  }

  // Get the beam centerline GIDs.
  Core::LinAlg::Matrix<beam::n_dof_, 1, int> beam_centerline_gid;
  BEAMINTERACTION::UTILS::GetElementCenterlineGIDIndices(
      discret, this->Element1(), beam_centerline_gid);

  // Get the patch GIDs.
  const std::vector<int>& patch_gid = this->face_element_->GetPatchGID();

  // Get the Lagrange multiplier GIDs.
  const auto& [lambda_gid_pos, _] = mortar_manager->LocationVector(*this);

  // Assemble into the matrix in the beam row and lambda column
  for (unsigned int i_beam = 0; i_beam < beam::n_dof_; i_beam++)
  {
    for (unsigned int i_lambda = 0; i_lambda < mortar::n_dof_; i_lambda++)
    {
      const double value = Core::FADUtils::CastToDouble(
          beam_shape_times_normal_times_lambda_shape_(i_beam, i_lambda));
      global_force_beam_lin_lambda.FEAssemble(
          value, beam_centerline_gid(i_beam), lambda_gid_pos[i_lambda]);
    }
  }

  // Assemble into the matrix in the surface row and lambda column
  for (unsigned int i_surface = 0; i_surface < surface::n_dof_; i_surface++)
  {
    for (unsigned int i_lambda = 0; i_lambda < mortar::n_dof_; i_lambda++)
    {
      const double value = Core::FADUtils::CastToDouble(
          surface_shape_times_normal_times_lambda_shape_(i_surface, i_lambda));
      global_force_solid_lin_lambda.FEAssemble(
          value, patch_gid[i_surface], lambda_gid_pos[i_lambda]);
    }
  }

  // Assemble into the in the lambda row
  for (unsigned int i_lambda = 0; i_lambda < mortar::n_dof_; i_lambda++)
  {
    // Assemble into the beam column
    for (unsigned int i_beam = 0; i_beam < beam::n_dof_; i_beam++)
    {
      const double value = constraint_vector(i_lambda).dx(i_beam);
      global_constraint_lin_beam.FEAssemble(
          value, lambda_gid_pos[i_lambda], beam_centerline_gid(i_beam));

      const double value_kappa_linearization = kappa(i_lambda).dx(i_beam);
      global_kappa_lin_beam.FEAssemble(
          value_kappa_linearization, lambda_gid_pos[i_lambda], beam_centerline_gid(i_beam));
    }

    // Assemble into the solid column
    for (unsigned int i_patch = 0; i_patch < patch_gid.size(); i_patch++)
    {
      const double value = constraint_vector(i_lambda).dx(beam::n_dof_ + i_patch);
      global_constraint_lin_solid.FEAssemble(value, lambda_gid_pos[i_lambda], patch_gid[i_patch]);

      const double value_kappa_linearization = kappa(i_lambda).dx(beam::n_dof_ + i_patch);
      global_kappa_lin_solid.FEAssemble(
          value_kappa_linearization, lambda_gid_pos[i_lambda], patch_gid[i_patch]);
    }
  }

  // Assemble into global coupling vector
  const auto constraint_vector_double = Core::FADUtils::CastToDouble(constraint_vector);
  global_constraint.SumIntoGlobalValues(
      lambda_gid_pos.size(), lambda_gid_pos.data(), constraint_vector_double.data());

  // Assemble into global kappa vector
  auto kappa_double = Core::FADUtils::CastToDouble(kappa);
  global_kappa.SumIntoGlobalValues(
      lambda_gid_pos.size(), lambda_gid_pos.data(), kappa_double.data());
  kappa_double.put_scalar(1.0);
  global_lambda_active.SumIntoGlobalValues(
      lambda_gid_pos.size(), lambda_gid_pos.data(), kappa_double.data());
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
void BEAMINTERACTION::BeamToSolidSurfaceContactPairMortar<scalar_type, beam, surface,
    mortar>::EvaluateAndAssemble(const Core::FE::Discretization& discret,
    const BeamToSolidMortarManager* mortar_manager,
    const Teuchos::RCP<Epetra_FEVector>& force_vector,
    const Teuchos::RCP<Core::LinAlg::SparseMatrix>& stiffness_matrix,
    const Epetra_Vector& global_lambda, const Epetra_Vector& displacement_vector)
{
  // At this point the pair is already evaluated in the current deformation state, so we don't have
  // to perform the projections or integration again, we can simply take the values previously
  // computed and multiply them with the Lagrange multipliers.

  // If there are no intersection segments, no contact terms will be assembled.
  const unsigned int n_segments = this->line_to_3D_segments_.size();
  if (n_segments == 0) return;

  // This pair only gives contributions to the stiffness matrix
  if (stiffness_matrix == Teuchos::null) return;

  // Get the Lagrange multipliers DOF vector for this pair
  const auto& [lambda_gid_pos, _] = mortar_manager->LocationVector(*this);
  std::vector<double> lambda_pos_vector;
  Core::FE::ExtractMyValues(global_lambda, lambda_pos_vector, lambda_gid_pos);
  const auto lambda_pos = Core::LinAlg::Matrix<mortar::n_dof_, 1, double>(lambda_pos_vector.data());

  // Multiply with the matrices evaluated in evaluate_and_assemble_mortar_contributions
  auto force_beam = Core::LinAlg::Matrix<beam::n_dof_, 1, scalar_type>(true);
  force_beam.multiply(beam_shape_times_normal_times_lambda_shape_, lambda_pos);
  auto force_surface = Core::LinAlg::Matrix<surface::n_dof_, 1, scalar_type>(true);
  force_surface.multiply(surface_shape_times_normal_times_lambda_shape_, lambda_pos);

  // Assemble the terms to the global stiffness matrix
  Core::LinAlg::Matrix<beam::n_dof_, 1, int> beam_centerline_gid;
  BEAMINTERACTION::UTILS::GetElementCenterlineGIDIndices(
      discret, this->Element1(), beam_centerline_gid);
  const std::vector<int>& patch_gid = this->face_element_->GetPatchGID();

  for (unsigned int i_beam = 0; i_beam < beam::n_dof_; i_beam++)
  {
    for (unsigned int j_beam = 0; j_beam < beam::n_dof_; j_beam++)
    {
      const double value = force_beam(i_beam).dx(j_beam);
      stiffness_matrix->FEAssemble(value, beam_centerline_gid(i_beam), beam_centerline_gid(j_beam));
    }
    for (unsigned int j_patch = 0; j_patch < patch_gid.size(); j_patch++)
    {
      const double value = force_beam(i_beam).dx(beam::n_dof_ + j_patch);
      stiffness_matrix->FEAssemble(value, beam_centerline_gid(i_beam), patch_gid[j_patch]);
    }
  }
  for (unsigned int i_surface = 0; i_surface < surface::n_dof_; i_surface++)
  {
    for (unsigned int j_beam = 0; j_beam < beam::n_dof_; j_beam++)
    {
      const double value = force_surface(i_surface).dx(j_beam);
      stiffness_matrix->FEAssemble(value, patch_gid[i_surface], beam_centerline_gid(j_beam));
    }
    for (unsigned int j_patch = 0; j_patch < patch_gid.size(); j_patch++)
    {
      const double value = force_surface(i_surface).dx(beam::n_dof_ + j_patch);
      stiffness_matrix->FEAssemble(value, patch_gid[i_surface], patch_gid[j_patch]);
    }
  }
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
scalar_type BEAMINTERACTION::BeamToSolidSurfaceContactPairMortar<scalar_type, beam, surface,
    mortar>::get_jacobian_for_configuration(const scalar_type& eta,
    const Inpar::BeamToSolid::BeamToSolidSurfaceContactMortarDefinedIn mortar_configuration) const
{
  Core::LinAlg::Matrix<3, 1, scalar_type> dr_beam;
  switch (mortar_configuration)
  {
    case Inpar::BeamToSolid::BeamToSolidSurfaceContactMortarDefinedIn::reference_configuration:
    {
      GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(eta, this->ele1posref_, dr_beam);
      return Core::FADUtils::VectorNorm(dr_beam);
    }
    case Inpar::BeamToSolid::BeamToSolidSurfaceContactMortarDefinedIn::current_configuration:
    {
      GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(eta, this->ele1pos_, dr_beam);
      return Core::FADUtils::VectorNorm(dr_beam);
    }
    default:
      FOUR_C_THROW("Got unexpected mortar configuration");
  }
}

/**
 * @brief Factory function templated on the type of beam element and the surface shape
 */
template <typename beam, typename surface, typename scalar_type>
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BeamToSolidSurfaceContactPairMortarFactoryTemplateBeamSurface(
    const Teuchos::RCP<const BEAMINTERACTION::BeamToSolidSurfaceContactParams>
        beam_to_surface_contact_params)
{
  using namespace GEOMETRYPAIR;

  switch (beam_to_surface_contact_params->get_mortar_shape_function_type())
  {
    case Inpar::BeamToSolid::BeamToSolidMortarShapefunctions::line2:
      return Teuchos::rcp(new BEAMINTERACTION::BeamToSolidSurfaceContactPairMortar<scalar_type,
          beam, surface, t_line2_scalar>);
    default:
      FOUR_C_THROW("Got unexpected mortar shape function");
  }
}

/**
 * @brief Factory function templated on the type of beam element
 */
template <typename beam>
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BeamToSolidSurfaceContactPairMortarFactoryTemplateBeam(
    const Teuchos::RCP<const BEAMINTERACTION::BeamToSolidSurfaceContactParams>
        beam_to_surface_contact_params,
    const Core::FE::CellType& surface_type)
{
  using namespace GEOMETRYPAIR;

  switch (surface_type)
  {
    case Core::FE::CellType::quad4:
      return BeamToSolidSurfaceContactPairMortarFactoryTemplateBeamSurface<beam, t_quad4,
          line_to_surface_patch_scalar_type_1st_order>(beam_to_surface_contact_params);
    default:
      FOUR_C_THROW("Got unexpected surface shape");
  }
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidSurfaceContactPairMortarFactory(
    const Teuchos::RCP<const BeamToSolidSurfaceContactParams> beam_to_surface_contact_params,
    const Core::FE::CellType& surface_type, const bool beam_is_hermite)
{
  using namespace GEOMETRYPAIR;

  if (beam_is_hermite)
  {
    return BeamToSolidSurfaceContactPairMortarFactoryTemplateBeam<t_hermite>(
        beam_to_surface_contact_params, surface_type);
  }
  else
  {
    FOUR_C_THROW("Beam-to-solid contact with mortar is not implemented for linear beam elements");
  }
}

FOUR_C_NAMESPACE_CLOSE