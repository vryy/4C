/*----------------------------------------------------------------------*/
/*! \file

\brief Mortar mesh tying element for between a 3D beam and a surface element, coupling terms are
evaluated with FAD.

\level 3
*/

#include "beaminteraction_beam_to_solid_surface_meshtying_pair_mortar_FAD.H"

#include "beam3_reissner.H"
#include "beaminteraction_contact_params.H"
#include "beaminteraction_beam_to_solid_surface_meshtying_params.H"
#include "beaminteraction_calc_utils.H"
#include "beaminteraction_beam_to_solid_mortar_manager.H"
#include "geometry_pair_line_to_surface.H"
#include "geometry_pair_element_faces.H"
#include "inpar_beam_to_solid.H"
#include "inpar_geometry_pair.H"
#include "beaminteraction_beam_to_solid_utils.H"
#include "beam3_triad_interpolation_local_rotation_vectors.H"

#include <Epetra_FEVector.h>


/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarFAD<scalar_type, beam, surface,
    mortar>::BeamToSolidSurfaceMeshtyingPairMortarFAD()
    : base_class()
{
  // Empty constructor.
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarFAD<scalar_type, beam, surface,
    mortar>::EvaluateAndAssemble(const ::DRT::Discretization& discret,
    const BeamToSolidMortarManager* mortar_manager,
    const Teuchos::RCP<Epetra_FEVector>& force_vector,
    const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix, const Epetra_Vector& global_lambda,
    const Epetra_Vector& displacement_vector)
{
  // Call Evaluate on the geometry Pair. Only do this once for meshtying.
  if (!this->meshtying_is_evaluated_)
  {
    this->CastGeometryPair()->Evaluate(this->ele1posref_,
        this->face_element_->GetFaceReferencePosition(), this->line_to_3D_segments_,
        this->face_element_->GetReferenceNormals());
    this->meshtying_is_evaluated_ = true;
  }

  // If there are no intersection segments, return no contact status.
  if (this->line_to_3D_segments_.size() == 0) return;

  // Get the positional Lagrange multipliers for this pair.
  std::vector<int> lambda_gid_pos;
  GetMortarGID(mortar_manager, this, mortar::n_dof_, this->n_mortar_rot_, &lambda_gid_pos, nullptr);
  std::vector<double> local_lambda_pos;
  DRT::UTILS::ExtractMyValues(global_lambda, local_lambda_pos, lambda_gid_pos);
  LINALG::Matrix<mortar::n_dof_, 1, double> q_lambda(local_lambda_pos.data());

  // Initialize variables for local values.
  LINALG::Matrix<3, 1, scalar_type> coupling_vector(true);
  LINALG::Matrix<3, 1, double> lambda(true);
  LINALG::Matrix<3, 1, double> dr_beam_ref(true);
  scalar_type potential = 0.0;

  // Initialize scalar variables.
  double segment_jacobian = 0.0;
  double beam_segmentation_factor = 0.0;

  // Loop over segments to evaluate the coupling potential.
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

      // Get the Gauss point contribution to the coupling potential.
      coupling_vector = this->EvaluateCoupling(projected_gauss_point);
      GEOMETRYPAIR::EvaluatePosition<mortar>(projected_gauss_point.GetEta(), q_lambda, lambda);
      for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
        potential += coupling_vector(i_dim) * lambda(i_dim) *
                     projected_gauss_point.GetGaussWeight() * segment_jacobian;
    }
  }

  // Get the pair GIDs.
  std::vector<int> pair_gid = this->GetPairGID(discret);

  // Add the terms to the global stiffness matrix.
  if (stiffness_matrix != Teuchos::null)
    for (unsigned int i_dof = 0; i_dof < pair_gid.size(); i_dof++)
      for (unsigned int j_dof = 0; j_dof < pair_gid.size(); j_dof++)
        stiffness_matrix->FEAssemble(CORE::FADUTILS::CastToDouble(potential.dx(i_dof).dx(j_dof)),
            pair_gid[i_dof], pair_gid[j_dof]);
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarFAD<scalar_type, beam, surface,
    mortar>::EvaluateAndAssembleMortarContributions(const ::DRT::Discretization& discret,
    const BeamToSolidMortarManager* mortar_manager, LINALG::SparseMatrix& global_G_B,
    LINALG::SparseMatrix& global_G_S, LINALG::SparseMatrix& global_FB_L,
    LINALG::SparseMatrix& global_FS_L, Epetra_FEVector& global_constraint,
    Epetra_FEVector& global_kappa, Epetra_FEVector& global_lambda_active,
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector)
{
  // Call Evaluate on the geometry Pair. Only do this once for meshtying.
  if (!this->meshtying_is_evaluated_)
  {
    this->CastGeometryPair()->Evaluate(this->ele1posref_,
        this->face_element_->GetFaceReferencePosition(), this->line_to_3D_segments_,
        this->face_element_->GetReferenceNormals());
    this->meshtying_is_evaluated_ = true;
  }

  // If there are no intersection segments, return no contact status.
  if (this->line_to_3D_segments_.size() == 0) return;

  // Initialize variables for local values.
  LINALG::Matrix<3, 1, scalar_type> coupling_vector(true);
  LINALG::Matrix<mortar::n_dof_, 1, scalar_type> constraint_vector(true);
  LINALG::Matrix<mortar::n_dof_, 1, double> local_kappa(true);
  LINALG::Matrix<3, 1, double> dr_beam_ref(true);
  LINALG::Matrix<1, mortar::n_nodes_ * mortar::n_val_, double> N_mortar(true);

  // Initialize scalar variables.
  double segment_jacobian = 0.0;
  double beam_segmentation_factor = 0.0;

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
          projected_gauss_point.GetEta(), this->ele1posref_, dr_beam_ref, this->Element1());

      // Jacobian including the segment length.
      segment_jacobian = dr_beam_ref.Norm2() * beam_segmentation_factor;

      // Get the mortar shape functions.
      N_mortar.Clear();
      mortar::EvaluateShapeFunction(
          N_mortar, projected_gauss_point.GetEta(), std::integral_constant<unsigned int, 1>{});

      // Fill in the local mortar scaling vector kappa.
      for (unsigned int i_mortar_node = 0; i_mortar_node < mortar::n_nodes_; i_mortar_node++)
        for (unsigned int i_mortar_val = 0; i_mortar_val < mortar::n_val_; i_mortar_val++)
          for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
            local_kappa(i_mortar_node * mortar::n_val_ * 3 + i_mortar_val * 3 + i_dim) +=
                N_mortar(i_mortar_node * mortar::n_val_ + i_mortar_val) *
                projected_gauss_point.GetGaussWeight() * segment_jacobian;

      // Get the constraint vector. This is the coupling potentials variation w.r.t the discrete
      // Lagrange multiplier DOFs.
      coupling_vector = this->EvaluateCoupling(projected_gauss_point);
      for (unsigned int i_mortar_node = 0; i_mortar_node < mortar::n_nodes_; i_mortar_node++)
        for (unsigned int i_mortar_val = 0; i_mortar_val < mortar::n_val_; i_mortar_val++)
          for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
            constraint_vector(i_mortar_node * mortar::n_val_ * 3 + i_mortar_val * 3 + i_dim) +=
                N_mortar(i_mortar_node * mortar::n_val_ + i_mortar_val) * coupling_vector(i_dim) *
                projected_gauss_point.GetGaussWeight() * segment_jacobian;
    }
  }


  // Get the beam centerline GIDs.
  LINALG::Matrix<beam::n_dof_, 1, int> beam_centerline_gid;
  BEAMINTERACTION::UTILS::GetElementCenterlineGIDIndices(
      discret, this->Element1(), beam_centerline_gid);

  // Get the patch GIDs.
  const std::vector<int>& patch_gid = this->face_element_->GetPatchGID();

  // Get the Lagrange multiplier GIDs.
  std::vector<int> lambda_gid_pos;
  GetMortarGID(mortar_manager, this, mortar::n_dof_, this->n_mortar_rot_, &lambda_gid_pos, nullptr);

  // Assemble into the matrices related to beam DOFs.
  for (unsigned int i_lambda = 0; i_lambda < mortar::n_dof_; i_lambda++)
    for (unsigned int i_beam = 0; i_beam < beam::n_dof_; i_beam++)
    {
      const double val = CORE::FADUTILS::CastToDouble(constraint_vector(i_lambda).dx(i_beam));
      global_G_B.FEAssemble(val, lambda_gid_pos[i_lambda], beam_centerline_gid(i_beam));
      global_FB_L.FEAssemble(val, beam_centerline_gid(i_beam), lambda_gid_pos[i_lambda]);
    }

  // Assemble into the matrices related to solid DOFs.
  for (unsigned int i_lambda = 0; i_lambda < mortar::n_dof_; i_lambda++)
    for (unsigned int i_patch = 0; i_patch < patch_gid.size(); i_patch++)
    {
      const double val =
          CORE::FADUTILS::CastToDouble(constraint_vector(i_lambda).dx(beam::n_dof_ + i_patch));
      global_G_S.FEAssemble(val, lambda_gid_pos[i_lambda], patch_gid[i_patch]);
      global_FS_L.FEAssemble(val, patch_gid[i_patch], lambda_gid_pos[i_lambda]);
    }

  // Assemble into global coupling vector.
  LINALG::Matrix<mortar::n_dof_, 1, double> constraint_vector_double =
      CORE::FADUTILS::CastToDouble(constraint_vector);
  global_constraint.SumIntoGlobalValues(
      lambda_gid_pos.size(), lambda_gid_pos.data(), constraint_vector_double.A());

  // Assemble into global kappa vector.
  global_kappa.SumIntoGlobalValues(lambda_gid_pos.size(), lambda_gid_pos.data(), local_kappa.A());

  // Assemble into global lambda active vector.
  local_kappa.PutScalar(1.0);
  global_lambda_active.SumIntoGlobalValues(
      lambda_gid_pos.size(), lambda_gid_pos.data(), local_kappa.A());
}


/**
 *
 */
template <typename surface, typename scalar_type_basis>
void GetSurfaceBasis(const LINALG::Matrix<3, 1, double>& xi,
    const LINALG::Matrix<surface::n_dof_, 1, scalar_type_basis>& q_solid,
    LINALG::Matrix<3, 3, scalar_type_basis>& surface_basis, const DRT::Element* element)
{
  // Calculate surface basis vectors in the surface plane.
  LINALG::Matrix<3, 2, scalar_type_basis> dr_surf(true);
  GEOMETRYPAIR::EvaluatePositionDerivative1<surface>(xi, q_solid, dr_surf, element);

  LINALG::Matrix<3, 1, scalar_type_basis> dr_surf_0;
  LINALG::Matrix<3, 1, scalar_type_basis> dr_surf_1;
  for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
  {
    dr_surf_0(i_dir) = dr_surf(i_dir, 0);
    dr_surf_1(i_dir) = dr_surf(i_dir, 1);
  }

  // Calculate normal on the basis vectors.
  LINALG::Matrix<3, 1, scalar_type_basis> element_surface_normal;
  element_surface_normal.CrossProduct(dr_surf_0, dr_surf_1);
  element_surface_normal.Scale(1.0 / CORE::FADUTILS::VectorNorm(element_surface_normal));

  // Put the new basis vectors in a matrix.
  for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
  {
    surface_basis(i_dir, 0) = dr_surf(i_dir, 0);
    surface_basis(i_dir, 1) = dr_surf(i_dir, 1);
    surface_basis(i_dir, 2) = element_surface_normal(i_dir);
  }
}

/**
 *
 */
template <typename surface, typename scalar_type_rot_vec>
void GetSurfaceRotationVectorAveraged(const LINALG::Matrix<3, 1, double>& xi,
    const LINALG::Matrix<surface::n_dof_, 1, double>& q_solid_ref,
    const LINALG::Matrix<surface::n_dof_, 1, scalar_type_rot_vec>& q_solid,
    const LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
    LINALG::Matrix<3, 1, scalar_type_rot_vec>& psi_solid, const DRT::Element* element)
{
  // Get beam basis vectors in reference configuration.
  LINALG::Matrix<3, 3, double> triad_beam_ref(true);
  CORE::LARGEROTATIONS::quaterniontotriad(quaternion_beam_ref, triad_beam_ref);

  // Calculate surface basis coordinate transformation matrix.
  LINALG::Matrix<3, 3, double> surface_basis_ref_inverse;
  GetSurfaceBasis<surface>(xi, q_solid_ref, surface_basis_ref_inverse, element);
  LINALG::Inverse(surface_basis_ref_inverse);

  // Calculate the current surface basis vectors.
  LINALG::Matrix<3, 3, scalar_type_rot_vec> surface_basis_current;
  GetSurfaceBasis<surface>(xi, q_solid, surface_basis_current, element);

  // Calculate the in plane surface deformation gradient.
  LINALG::Matrix<3, 3, scalar_type_rot_vec> surface_basis_ref_inverse_scalar_type;
  for (unsigned int i_row = 0; i_row < 3; i_row++)
    for (unsigned int i_col = 0; i_col < 3; i_col++)
      surface_basis_ref_inverse_scalar_type(i_row, i_col) = surface_basis_ref_inverse(i_row, i_col);
  LINALG::Matrix<3, 3, scalar_type_rot_vec> surface_F;
  surface_F.Multiply(surface_basis_current, surface_basis_ref_inverse_scalar_type);

  // Get the solid rotation vector from the deformation gradient via construction in the cross
  // section plane.
  BEAMINTERACTION::GetSolidRotationVectorDeformationGradient3DGeneralInCrossSectionPlane(
      surface_F, triad_beam_ref, psi_solid);
}

/**
 *
 */
template <typename surface, typename scalar_type_rot_vec>
void GetSurfaceRotationVectorCrossSectionDirector(const LINALG::Matrix<3, 1, double>& xi,
    const LINALG::Matrix<surface::n_dof_, 1, double>& q_solid_ref,
    const LINALG::Matrix<surface::n_dof_, 1, scalar_type_rot_vec>& q_solid,
    const LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
    LINALG::Matrix<3, 1, scalar_type_rot_vec>& psi_solid, const DRT::Element* element)
{
  // Get beam basis vectors in reference configuration.
  LINALG::Matrix<3, 3, double> triad_beam_ref(true);
  CORE::LARGEROTATIONS::quaterniontotriad(quaternion_beam_ref, triad_beam_ref);

  // Get the surface basis vectors in the reference configuration.
  LINALG::Matrix<3, 3, double> surface_basis_ref;
  GetSurfaceBasis<surface>(xi, q_solid_ref, surface_basis_ref, element);

  // Get the surface basis vectors in the current configuration.
  LINALG::Matrix<3, 3, scalar_type_rot_vec> surface_basis_current;
  GetSurfaceBasis<surface>(xi, q_solid, surface_basis_current, element);

  // Get the surface material director (the intersection between the beam cross-section and the
  // surface tangent plane).
  LINALG::Matrix<3, 1, double> surface_normal_ref;
  LINALG::Matrix<3, 1, double> beam_cross_section_normal_ref;
  for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
  {
    surface_normal_ref(i_dim) = surface_basis_ref(i_dim, 2);
    beam_cross_section_normal_ref(i_dim) = triad_beam_ref(i_dim, 0);
  }
  LINALG::Matrix<3, 1, double> surface_material_director_ref;
  surface_material_director_ref.CrossProduct(surface_normal_ref, beam_cross_section_normal_ref);
  surface_material_director_ref.Scale(
      1.0 / CORE::FADUTILS::VectorNorm(surface_material_director_ref));

  // Get the reference triad of the surface.
  LINALG::Matrix<3, 1, double> surface_material_director_perpendicular_ref;
  surface_material_director_perpendicular_ref.CrossProduct(
      surface_material_director_ref, surface_normal_ref);
  LINALG::Matrix<3, 3, double> surface_triad_ref;
  for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
  {
    surface_triad_ref(i_dim, 0) = surface_material_director_ref(i_dim);
    surface_triad_ref(i_dim, 1) = surface_normal_ref(i_dim);
    surface_triad_ref(i_dim, 2) = surface_material_director_perpendicular_ref(i_dim);
  }

  // Get the offset of the reference triad, so it matches the beam reference triad.
  LINALG::Matrix<3, 3, double> surface_triad_offset;
  surface_triad_offset.MultiplyTN(surface_triad_ref, triad_beam_ref);

  // Calculate the in plane surface deformation gradient.
  LINALG::Matrix<3, 3, double> surface_basis_ref_inverse;
  surface_basis_ref_inverse = surface_basis_ref;
  LINALG::Inverse(surface_basis_ref_inverse);
  LINALG::Matrix<3, 3, scalar_type_rot_vec> surface_basis_ref_inverse_scalar_type;
  for (unsigned int i_row = 0; i_row < 3; i_row++)
    for (unsigned int i_col = 0; i_col < 3; i_col++)
      surface_basis_ref_inverse_scalar_type(i_row, i_col) = surface_basis_ref_inverse(i_row, i_col);
  LINALG::Matrix<3, 3, scalar_type_rot_vec> surface_F;
  surface_F.Multiply(surface_basis_current, surface_basis_ref_inverse_scalar_type);

  // Get the current material director.
  LINALG::Matrix<3, 1, scalar_type_rot_vec> surface_material_director_ref_fad;
  for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
    surface_material_director_ref_fad(i_dim) = surface_material_director_ref(i_dim);
  LINALG::Matrix<3, 1, scalar_type_rot_vec> surface_material_director_current;
  surface_material_director_current.Multiply(surface_F, surface_material_director_ref_fad);
  surface_material_director_current.Scale(
      1.0 / CORE::FADUTILS::VectorNorm(surface_material_director_current));

  // Get the current triad of the surface.
  LINALG::Matrix<3, 1, scalar_type_rot_vec> surface_normal_current;
  for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
    surface_normal_current(i_dim) = surface_basis_current(i_dim, 2);
  LINALG::Matrix<3, 1, scalar_type_rot_vec> surface_material_director_perpendicular_current;
  surface_material_director_perpendicular_current.CrossProduct(
      surface_material_director_current, surface_normal_current);
  LINALG::Matrix<3, 3, scalar_type_rot_vec> surface_triad_current;
  for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
  {
    surface_triad_current(i_dim, 0) = surface_material_director_current(i_dim);
    surface_triad_current(i_dim, 1) = surface_normal_current(i_dim);
    surface_triad_current(i_dim, 2) = surface_material_director_perpendicular_current(i_dim);
  }

  // Add the offset to the surface triad.
  LINALG::Matrix<3, 3, scalar_type_rot_vec> surface_triad_offset_fad;
  for (unsigned int i_row = 0; i_row < 3; i_row++)
    for (unsigned int i_col = 0; i_col < 3; i_col++)
      surface_triad_offset_fad(i_row, i_col) = surface_triad_offset(i_row, i_col);
  LINALG::Matrix<3, 3, scalar_type_rot_vec> surface_triad_current_with_offset;
  surface_triad_current_with_offset.Multiply(surface_triad_current, surface_triad_offset_fad);

  // Get the rotation angle.
  LINALG::Matrix<4, 1, scalar_type_rot_vec> rot_quat;
  CORE::LARGEROTATIONS::triadtoquaternion(surface_triad_current_with_offset, rot_quat);
  CORE::LARGEROTATIONS::quaterniontoangle(rot_quat, psi_solid);

#ifdef DEBUG
  LINALG::Matrix<3, 1, scalar_type_rot_vec> current_normal;
  for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
    current_normal(i_dim) = surface_basis_current(i_dim, 2);
  if (abs(surface_material_director_current.Dot(current_normal)) > 1e-10)
    dserror("The current material director has to lie within the surface tangent plane.");
#endif
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
template <typename scalar_type_rot_vec>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarRotationFAD<scalar_type, beam, surface,
    mortar>::GetSurfaceRotationVector(const LINALG::Matrix<3, 1, double>& xi,
    const LINALG::Matrix<surface::n_dof_, 1, double>& q_solid_ref,
    const LINALG::Matrix<surface::n_dof_, 1, scalar_type_rot_vec>& q_solid,
    const LINALG::Matrix<4, 1, double>& quaternion_beam_ref,
    const INPAR::BEAMTOSOLID::BeamToSolidSurfaceRotationCoupling surface_triad_type,
    LINALG::Matrix<3, 1, scalar_type_rot_vec>& psi_solid, const DRT::Element* element) const
{
  switch (surface_triad_type)
  {
    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceRotationCoupling::averaged:
      GetSurfaceRotationVectorAveraged<surface>(
          xi, q_solid_ref, q_solid, quaternion_beam_ref, psi_solid, element);
      break;
    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceRotationCoupling::surface_cross_section_director:
      GetSurfaceRotationVectorCrossSectionDirector<surface>(
          xi, q_solid_ref, q_solid, quaternion_beam_ref, psi_solid, element);
      break;
    default:
      dserror("Please supply a suitable solid triad construction.");
      break;
  }
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarRotationFAD<scalar_type, beam, surface,
    mortar>::EvaluateAndAssemble(const ::DRT::Discretization& discret,
    const BeamToSolidMortarManager* mortar_manager,
    const Teuchos::RCP<Epetra_FEVector>& force_vector,
    const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix, const Epetra_Vector& global_lambda,
    const Epetra_Vector& displacement_vector)
{
  base_class::EvaluateAndAssemble(
      discret, mortar_manager, force_vector, stiffness_matrix, global_lambda, displacement_vector);

  // If there are no intersection segments, return as no contact can occur.
  if (this->line_to_3D_segments_.size() == 0) return;

  // This pair only gives contributions to the stiffness matrix.
  if (stiffness_matrix == Teuchos::null) return;

  // Get the beam triad interpolation schemes.
  LARGEROTATIONS::TriadInterpolationLocalRotationVectors<3, double> triad_interpolation_scheme;
  LARGEROTATIONS::TriadInterpolationLocalRotationVectors<3, double> ref_triad_interpolation_scheme;
  GetBeamTriadInterpolationScheme(discret, Teuchos::rcpFromRef(displacement_vector),
      this->Element1(), triad_interpolation_scheme, ref_triad_interpolation_scheme);

  // Set the FAD variables for the solid DOFs. For the terms calculated here we need second
  // order derivatives.
  LINALG::Matrix<surface::n_dof_, 1, scalar_type_rot_2nd> q_surface(true);
  for (unsigned int i_surface = 0; i_surface < surface::n_dof_; i_surface++)
    q_surface(i_surface) = CORE::FADUTILS::HigherOrderFadValue<scalar_type_rot_2nd>::apply(
        3 + surface::n_dof_, 3 + i_surface,
        CORE::FADUTILS::CastToDouble(this->face_element_->GetFacePosition()(i_surface)));

  // Get the rotational Lagrange multipliers for this pair.
  std::vector<int> lambda_gid_rot;
  GetMortarGID(mortar_manager, this, mortar::n_dof_, mortar::n_dof_, nullptr, &lambda_gid_rot);
  std::vector<double> lambda_rot_double;
  DRT::UTILS::ExtractMyValues(global_lambda, lambda_rot_double, lambda_gid_rot);
  LINALG::Matrix<mortar::n_dof_, 1, double> lambda_rot;
  for (unsigned int i_dof = 0; i_dof < mortar::n_dof_; i_dof++)
    lambda_rot(i_dof) = lambda_rot_double[i_dof];

  // Get the type of surface triad construction.
  const auto surface_triad_type =
      this->Params()->BeamToSolidSurfaceMeshtyingParams()->GetSurfaceTriadConstruction();

  // Initialize local matrices.
  LINALG::Matrix<n_dof_rot_, n_dof_rot_, double> local_stiff_BB(true);
  LINALG::Matrix<n_dof_rot_, surface::n_dof_, double> local_stiff_BS(true);
  LINALG::Matrix<surface::n_dof_, n_dof_rot_, double> local_stiff_SB(true);
  LINALG::Matrix<surface::n_dof_, surface::n_dof_, double> local_stiff_SS(true);

  // Evaluate the pair wise terms.
  {
    // Initialize variables.
    LINALG::Matrix<3, 1, double> dr_beam_ref;
    LINALG::Matrix<4, 1, double> quaternion_beam_double;
    LINALG::Matrix<3, 1, double> psi_beam_double;
    LINALG::Matrix<3, 1, scalar_type_rot_1st> psi_beam;
    LINALG::Matrix<3, 1, scalar_type_rot_2nd> psi_surface;
    LINALG::Matrix<3, 1, scalar_type_rot_1st> psi_surface_val;
    LINALG::Matrix<3, 1, scalar_type_rot_1st> psi_rel;
    LINALG::Matrix<4, 1, scalar_type_rot_1st> quaternion_beam;
    LINALG::Matrix<4, 1, scalar_type_rot_1st> quaternion_beam_inv;
    LINALG::Matrix<4, 1, double> quaternion_beam_ref;
    LINALG::Matrix<4, 1, scalar_type_rot_1st> quaternion_surface;
    LINALG::Matrix<4, 1, scalar_type_rot_1st> quaternion_rel;
    LINALG::Matrix<3, 3, double> T_beam;
    LINALG::Matrix<3, 3, scalar_type_rot_1st> T_surface;
    LINALG::Matrix<3, 3, scalar_type_rot_1st> T_surface_inv;
    LINALG::Matrix<3, 3, scalar_type_rot_1st> T_rel;

    LINALG::Matrix<mortar::n_nodes_, 1, double> lambda_shape_functions;
    LINALG::Matrix<3, mortar::n_dof_, scalar_type_rot_1st> lambda_shape_functions_full(true);
    Epetra_SerialDenseVector L_i(3);
    LINALG::Matrix<3, n_dof_rot_, scalar_type_rot_1st> L_full(true);
    std::vector<LINALG::Matrix<3, 3, double>> I_beam_tilde;
    LINALG::Matrix<3, n_dof_rot_, double> I_beam_tilde_full;
    LINALG::Matrix<3, n_dof_rot_, double> T_beam_times_I_beam_tilde_full;
    LINALG::Matrix<3, mortar::n_dof_, scalar_type_rot_1st> T_rel_tr_times_lambda_shape;
    LINALG::Matrix<3, mortar::n_dof_, scalar_type_rot_1st>
        T_surface_mtr_times_T_rel_tr_times_lambda_shape;
    LINALG::Matrix<n_dof_rot_, mortar::n_dof_, scalar_type_rot_1st> d_fb_d_lambda_gp;
    LINALG::Matrix<surface::n_dof_, mortar::n_dof_, scalar_type_rot_1st> d_fs_d_lambda_gp;
    LINALG::Matrix<3, surface::n_dof_, scalar_type_rot_1st> d_psi_surface_d_q_surface;
    LINALG::Matrix<mortar::n_dof_, 3, double> d_g_d_psi_beam;
    LINALG::Matrix<mortar::n_dof_, surface::n_dof_, double> d_g_d_q_surface;
    LINALG::Matrix<n_dof_rot_, 1, scalar_type_rot_1st> f_beam;
    LINALG::Matrix<surface::n_dof_, 1, scalar_type_rot_1st> f_surface;
    LINALG::Matrix<n_dof_rot_, 3, double> d_f_beam_d_phi;
    LINALG::Matrix<surface::n_dof_, 3, double> d_f_surface_d_phi;
    LINALG::Matrix<n_dof_rot_, n_dof_rot_, double>
        d_f_beam_d_phi_times_T_beam_times_I_beam_tilde_full;
    LINALG::Matrix<surface::n_dof_, n_dof_rot_, double>
        d_f_surface_d_phi_times_T_beam_times_I_beam_tilde_full;

    // Initialize scalar variables.
    double segment_jacobian = 0.0;
    double beam_segmentation_factor = 0.0;

    // Calculate the meshtying forces.
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

        // Get the jacobian in the reference configuration.
        GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(
            projected_gauss_point.GetEta(), this->ele1posref_, dr_beam_ref, this->Element1());

        // Jacobian including the segment length.
        segment_jacobian = dr_beam_ref.Norm2() * beam_segmentation_factor;

        // Calculate the rotation vector of this cross section.
        triad_interpolation_scheme.GetInterpolatedQuaternionAtXi(
            quaternion_beam_double, projected_gauss_point.GetEta());
        CORE::LARGEROTATIONS::quaterniontoangle(quaternion_beam_double, psi_beam_double);
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          psi_beam(i_dim) = CORE::FADUTILS::HigherOrderFadValue<scalar_type_rot_1st>::apply(
              3 + surface::n_dof_, i_dim, psi_beam_double(i_dim));
        CORE::LARGEROTATIONS::angletoquaternion(psi_beam, quaternion_beam);
        quaternion_beam_inv = CORE::LARGEROTATIONS::inversequaternion(quaternion_beam);

        // Get the surface rotation vector.
        ref_triad_interpolation_scheme.GetInterpolatedQuaternionAtXi(
            quaternion_beam_ref, projected_gauss_point.GetEta());
        GetSurfaceRotationVector(projected_gauss_point.GetXi(),
            this->face_element_->GetFaceReferencePosition(), q_surface, quaternion_beam_ref,
            surface_triad_type, psi_surface, this->face_element_->GetDrtFaceElement());
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          psi_surface_val(i_dim) = psi_surface(i_dim).val();
        CORE::LARGEROTATIONS::angletoquaternion(psi_surface_val, quaternion_surface);

        // Calculate the relative rotation vector.
        CORE::LARGEROTATIONS::quaternionproduct(
            quaternion_beam_inv, quaternion_surface, quaternion_rel);
        CORE::LARGEROTATIONS::quaterniontoangle(quaternion_rel, psi_rel);

        // Calculate the transformation matrices.
        T_rel = CORE::LARGEROTATIONS::Tmatrix(psi_rel);
        T_beam = CORE::LARGEROTATIONS::Tmatrix(CORE::FADUTILS::CastToDouble(psi_beam));
        T_surface = CORE::LARGEROTATIONS::Tmatrix(psi_surface_val);
        T_surface_inv = T_surface;
        LINALG::Inverse(T_surface_inv);

        // Evaluate mortar shape functions.
        mortar::EvaluateShapeFunction(lambda_shape_functions, projected_gauss_point.GetEta(),
            std::integral_constant<unsigned int, mortar::dim_>{});
        for (unsigned int i_node = 0; i_node < mortar::n_nodes_; i_node++)
          for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
            lambda_shape_functions_full(i_dim, 3 * i_node + i_dim) = lambda_shape_functions(i_node);

        // Get the shape functions for the interpolation of the beam rotations. This is currently
        // only implemented for 2nd order Lagrange interpolation (Beam3rHerm2Line3).
        const unsigned int n_nodes_rot = 3;
        CORE::DRT::UTILS::shape_function_1D(
            L_i, projected_gauss_point.GetEta(), DRT::Element::line3);
        for (unsigned int i_node = 0; i_node < n_nodes_rot; i_node++)
          for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
            L_full(i_dim, 3 * i_node + i_dim) = L_i(i_node);

        triad_interpolation_scheme.GetNodalGeneralizedRotationInterpolationMatricesAtXi(
            I_beam_tilde, projected_gauss_point.GetEta());
        for (unsigned int i_node = 0; i_node < n_nodes_rot; i_node++)
          for (unsigned int i_dim_0 = 0; i_dim_0 < 3; i_dim_0++)
            for (unsigned int i_dim_1 = 0; i_dim_1 < 3; i_dim_1++)
              I_beam_tilde_full(i_dim_0, i_node * 3 + i_dim_1) =
                  I_beam_tilde[i_node](i_dim_0, i_dim_1);

        // Solid angle derived w.r.t. the surface DOFs.
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          for (unsigned int i_surface = 0; i_surface < surface::n_dof_; i_surface++)
            d_psi_surface_d_q_surface(i_dim, i_surface) = psi_surface(i_dim).dx(3 + i_surface);

        // Calculate the force terms derived w.r.t. the Lagrange multipliers.
        T_rel_tr_times_lambda_shape.MultiplyTN(T_rel, lambda_shape_functions_full);
        d_fb_d_lambda_gp.MultiplyTN(L_full, T_rel_tr_times_lambda_shape);
        d_fb_d_lambda_gp.Scale(-1.0 * projected_gauss_point.GetGaussWeight() * segment_jacobian);

        T_surface_mtr_times_T_rel_tr_times_lambda_shape.MultiplyTN(
            T_surface_inv, T_rel_tr_times_lambda_shape);
        d_fs_d_lambda_gp.MultiplyTN(
            d_psi_surface_d_q_surface, T_surface_mtr_times_T_rel_tr_times_lambda_shape);
        d_fs_d_lambda_gp.Scale(projected_gauss_point.GetGaussWeight() * segment_jacobian);

        // Calculate the force vectors.
        f_beam.PutScalar(0.0);
        for (unsigned int i_row = 0; i_row < n_dof_rot_; i_row++)
          for (unsigned int i_col = 0; i_col < mortar::n_dof_; i_col++)
            f_beam(i_row) += d_fb_d_lambda_gp(i_row, i_col) * lambda_rot(i_col);
        f_surface.PutScalar(0.0);
        for (unsigned int i_row = 0; i_row < surface::n_dof_; i_row++)
          for (unsigned int i_col = 0; i_col < mortar::n_dof_; i_col++)
            f_surface(i_row) += d_fs_d_lambda_gp(i_row, i_col) * lambda_rot(i_col);

        // Derivatives of the force vectors.
        for (unsigned int i_row = 0; i_row < n_dof_rot_; i_row++)
          for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
            d_f_beam_d_phi(i_row, i_dim) = f_beam(i_row).dx(i_dim);
        for (unsigned int i_row = 0; i_row < surface::n_dof_; i_row++)
          for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
            d_f_surface_d_phi(i_row, i_dim) = f_surface(i_row).dx(i_dim);

        T_beam_times_I_beam_tilde_full.Multiply(T_beam, I_beam_tilde_full);
        d_f_beam_d_phi_times_T_beam_times_I_beam_tilde_full.Multiply(
            d_f_beam_d_phi, T_beam_times_I_beam_tilde_full);
        d_f_surface_d_phi_times_T_beam_times_I_beam_tilde_full.Multiply(
            d_f_surface_d_phi, T_beam_times_I_beam_tilde_full);

        // Add to output matrices and vector.
        local_stiff_BB += d_f_beam_d_phi_times_T_beam_times_I_beam_tilde_full;
        for (unsigned int i_beam = 0; i_beam < n_dof_rot_; i_beam++)
          for (unsigned int j_surface = 0; j_surface < surface::n_dof_; j_surface++)
            local_stiff_BS(i_beam, j_surface) += f_beam(i_beam).dx(3 + j_surface);
        local_stiff_SB += d_f_surface_d_phi_times_T_beam_times_I_beam_tilde_full;
        for (unsigned int i_surface = 0; i_surface < surface::n_dof_; i_surface++)
          for (unsigned int j_surface = 0; j_surface < surface::n_dof_; j_surface++)
            local_stiff_SS(i_surface, j_surface) += f_surface(i_surface).dx(3 + j_surface);
      }
    }
  }

  // Get the rotational GIDs of the surface and beam.
  std::vector<int> gid_surface;
  LINALG::Matrix<n_dof_rot_, 1, int> gid_rot;
  GetPairRotationalGIDs(discret, gid_surface, gid_rot);

  // Assemble into global matrix.
  for (unsigned int i_dof_beam = 0; i_dof_beam < n_dof_rot_; i_dof_beam++)
  {
    for (unsigned int j_dof_beam = 0; j_dof_beam < n_dof_rot_; j_dof_beam++)
      stiffness_matrix->FEAssemble(
          local_stiff_BB(i_dof_beam, j_dof_beam), gid_rot(i_dof_beam), gid_rot(j_dof_beam));
    for (unsigned int j_dof_surface = 0; j_dof_surface < surface::n_dof_; j_dof_surface++)
      stiffness_matrix->FEAssemble(local_stiff_BS(i_dof_beam, j_dof_surface), gid_rot(i_dof_beam),
          gid_surface[j_dof_surface]);
  }
  for (unsigned int i_dof_surface = 0; i_dof_surface < surface::n_dof_; i_dof_surface++)
  {
    for (unsigned int j_dof_beam = 0; j_dof_beam < n_dof_rot_; j_dof_beam++)
      stiffness_matrix->FEAssemble(local_stiff_SB(i_dof_surface, j_dof_beam),
          gid_surface[i_dof_surface], gid_rot(j_dof_beam));
    for (unsigned int j_dof_surface = 0; j_dof_surface < surface::n_dof_; j_dof_surface++)
      stiffness_matrix->FEAssemble(local_stiff_SS(i_dof_surface, j_dof_surface),
          gid_surface[i_dof_surface], gid_surface[j_dof_surface]);
  }
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarRotationFAD<scalar_type, beam, surface,
    mortar>::EvaluateAndAssembleMortarContributions(const ::DRT::Discretization& discret,
    const BeamToSolidMortarManager* mortar_manager, LINALG::SparseMatrix& global_GB,
    LINALG::SparseMatrix& global_GS, LINALG::SparseMatrix& global_FB,
    LINALG::SparseMatrix& global_FS, Epetra_FEVector& global_constraint,
    Epetra_FEVector& global_kappa, Epetra_FEVector& global_lambda_active,
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector)
{
  base_class::EvaluateAndAssembleMortarContributions(discret, mortar_manager, global_GB, global_GS,
      global_FB, global_FS, global_constraint, global_kappa, global_lambda_active,
      displacement_vector);

  // If there are no intersection segments, return as no contact can occur.
  if (this->line_to_3D_segments_.size() == 0) return;

  // Get the beam triad interpolation schemes.
  LARGEROTATIONS::TriadInterpolationLocalRotationVectors<3, double> triad_interpolation_scheme;
  LARGEROTATIONS::TriadInterpolationLocalRotationVectors<3, double> ref_triad_interpolation_scheme;
  GetBeamTriadInterpolationScheme(discret, displacement_vector, this->Element1(),
      triad_interpolation_scheme, ref_triad_interpolation_scheme);

  // Set the FAD variables for the surface DOFs. For the terms calculated here we only need first
  // order derivatives.
  LINALG::Matrix<surface::n_dof_, 1, scalar_type_rot_1st> q_surface(true);
  for (unsigned int i_surface = 0; i_surface < surface::n_dof_; i_surface++)
    q_surface(i_surface) = CORE::FADUTILS::HigherOrderFadValue<scalar_type_rot_1st>::apply(
        3 + surface::n_dof_, 3 + i_surface,
        CORE::FADUTILS::CastToDouble(this->face_element_->GetFacePosition()(i_surface)));

  // Initialize local matrices.
  LINALG::Matrix<mortar::n_dof_, 1, double> local_g(true);
  LINALG::Matrix<mortar::n_dof_, n_dof_rot_, double> local_GB(true);
  LINALG::Matrix<mortar::n_dof_, surface::n_dof_, double> local_GS(true);
  LINALG::Matrix<n_dof_rot_, mortar::n_dof_, double> local_FB(true);
  LINALG::Matrix<surface::n_dof_, mortar::n_dof_, double> local_FS(true);
  LINALG::Matrix<mortar::n_dof_, 1, double> local_kappa(true);

  // Get the type of surface triad construction.
  const auto surface_triad_type =
      this->Params()->BeamToSolidSurfaceMeshtyingParams()->GetSurfaceTriadConstruction();

  // Evaluate the mortar terms for this pair.
  {
    // Initialize variables.
    LINALG::Matrix<3, 1, double> dr_beam_ref;
    LINALG::Matrix<4, 1, double> quaternion_beam_double;
    LINALG::Matrix<3, 1, double> psi_beam_double;
    LINALG::Matrix<3, 1, scalar_type_rot_1st> psi_beam;
    LINALG::Matrix<3, 1, scalar_type_rot_1st> psi_surface;
    LINALG::Matrix<3, 1, scalar_type_rot_1st> psi_rel;
    LINALG::Matrix<4, 1, scalar_type_rot_1st> quaternion_beam;
    LINALG::Matrix<4, 1, scalar_type_rot_1st> quaternion_beam_inv;
    LINALG::Matrix<4, 1, double> quaternion_beam_ref;
    LINALG::Matrix<4, 1, scalar_type_rot_1st> quaternion_surface;
    LINALG::Matrix<4, 1, scalar_type_rot_1st> quaternion_rel;
    LINALG::Matrix<3, 3, double> T_beam;
    LINALG::Matrix<3, 3, double> T_surface;
    LINALG::Matrix<3, 3, double> T_surface_inv;
    LINALG::Matrix<3, 3, double> T_rel;

    LINALG::Matrix<mortar::n_nodes_, 1, double> lambda_shape_functions;
    LINALG::Matrix<3, mortar::n_dof_, double> lambda_shape_functions_full(true);
    Epetra_SerialDenseVector L_i(3);
    LINALG::Matrix<3, n_dof_rot_, double> L_full(true);
    std::vector<LINALG::Matrix<3, 3, double>> I_beam_tilde;
    LINALG::Matrix<3, n_dof_rot_, double> I_beam_tilde_full;
    LINALG::Matrix<3, n_dof_rot_, double> T_beam_times_I_beam_tilde_full;
    LINALG::Matrix<3, mortar::n_dof_, double> T_rel_tr_times_lambda_shape;
    LINALG::Matrix<3, mortar::n_dof_, double> T_surface_mtr_times_T_rel_tr_times_lambda_shape;
    LINALG::Matrix<n_dof_rot_, mortar::n_dof_, double> d_fb_d_lambda_gp;
    LINALG::Matrix<surface::n_dof_, mortar::n_dof_, double> d_fs_d_lambda_gp;
    LINALG::Matrix<mortar::n_dof_, 1, scalar_type_rot_1st> g_gp;
    LINALG::Matrix<3, surface::n_dof_, double> d_psi_surface_d_q_surface;
    LINALG::Matrix<mortar::n_dof_, 3, double> d_g_d_psi_beam;
    LINALG::Matrix<mortar::n_dof_, n_dof_rot_, double> d_g_d_psi_beam_times_T_beam_I;
    LINALG::Matrix<mortar::n_dof_, surface::n_dof_, double> d_g_d_q_surface;

    // Initialize scalar variables.
    double segment_jacobian = 0.0;
    double beam_segmentation_factor = 0.0;

    // Calculate the meshtying forces.
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

        // Get the jacobian in the reference configuration.
        GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(
            projected_gauss_point.GetEta(), this->ele1posref_, dr_beam_ref, this->Element1());

        // Jacobian including the segment length.
        segment_jacobian = dr_beam_ref.Norm2() * beam_segmentation_factor;

        // Calculate the rotation vector of this cross section.
        triad_interpolation_scheme.GetInterpolatedQuaternionAtXi(
            quaternion_beam_double, projected_gauss_point.GetEta());
        CORE::LARGEROTATIONS::quaterniontoangle(quaternion_beam_double, psi_beam_double);
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          psi_beam(i_dim) = CORE::FADUTILS::HigherOrderFadValue<scalar_type_rot_1st>::apply(
              3 + surface::n_dof_, i_dim, psi_beam_double(i_dim));
        CORE::LARGEROTATIONS::angletoquaternion(psi_beam, quaternion_beam);
        quaternion_beam_inv = CORE::LARGEROTATIONS::inversequaternion(quaternion_beam);

        // Get the surface rotation vector.
        ref_triad_interpolation_scheme.GetInterpolatedQuaternionAtXi(
            quaternion_beam_ref, projected_gauss_point.GetEta());
        GetSurfaceRotationVector(projected_gauss_point.GetXi(),
            this->face_element_->GetFaceReferencePosition(), q_surface, quaternion_beam_ref,
            surface_triad_type, psi_surface, this->face_element_->GetDrtFaceElement());
        CORE::LARGEROTATIONS::angletoquaternion(psi_surface, quaternion_surface);

        // Calculate the relative rotation vector.
        CORE::LARGEROTATIONS::quaternionproduct(
            quaternion_beam_inv, quaternion_surface, quaternion_rel);
        CORE::LARGEROTATIONS::quaterniontoangle(quaternion_rel, psi_rel);

        // Calculate the transformation matrices.
        T_rel = CORE::LARGEROTATIONS::Tmatrix(CORE::FADUTILS::CastToDouble(psi_rel));
        T_beam = CORE::LARGEROTATIONS::Tmatrix(CORE::FADUTILS::CastToDouble(psi_beam));
        T_surface = CORE::LARGEROTATIONS::Tmatrix(CORE::FADUTILS::CastToDouble(psi_surface));
        T_surface_inv = T_surface;
        LINALG::Inverse(T_surface_inv);

        // Evaluate mortar shape functions.
        mortar::EvaluateShapeFunction(lambda_shape_functions, projected_gauss_point.GetEta(),
            std::integral_constant<unsigned int, mortar::dim_>{});
        for (unsigned int i_node = 0; i_node < mortar::n_nodes_; i_node++)
          for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
            lambda_shape_functions_full(i_dim, 3 * i_node + i_dim) = lambda_shape_functions(i_node);

        // Get the shape functions for the interpolation of the beam rotations. This is currently
        // only implemented for 2nd order Lagrange interpolation (Beam3rHerm2Line3).
        const unsigned int n_nodes_rot = 3;
        CORE::DRT::UTILS::shape_function_1D(
            L_i, projected_gauss_point.GetEta(), DRT::Element::line3);
        for (unsigned int i_node = 0; i_node < n_nodes_rot; i_node++)
          for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
            L_full(i_dim, 3 * i_node + i_dim) = L_i(i_node);

        triad_interpolation_scheme.GetNodalGeneralizedRotationInterpolationMatricesAtXi(
            I_beam_tilde, projected_gauss_point.GetEta());
        for (unsigned int i_node = 0; i_node < n_nodes_rot; i_node++)
          for (unsigned int i_dim_0 = 0; i_dim_0 < 3; i_dim_0++)
            for (unsigned int i_dim_1 = 0; i_dim_1 < 3; i_dim_1++)
              I_beam_tilde_full(i_dim_0, i_node * 3 + i_dim_1) =
                  I_beam_tilde[i_node](i_dim_0, i_dim_1);

        // Solid angle derived w.r.t. the surface DOFs.
        for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
          for (unsigned int i_surface = 0; i_surface < surface::n_dof_; i_surface++)
            d_psi_surface_d_q_surface(i_dim, i_surface) = psi_surface(i_dim).dx(3 + i_surface);

        // Calculate the force terms derived w.r.t. the Lagrange multipliers.
        T_rel_tr_times_lambda_shape.MultiplyTN(T_rel, lambda_shape_functions_full);
        d_fb_d_lambda_gp.MultiplyTN(L_full, T_rel_tr_times_lambda_shape);
        d_fb_d_lambda_gp.Scale(-1.0 * projected_gauss_point.GetGaussWeight() * segment_jacobian);

        T_surface_mtr_times_T_rel_tr_times_lambda_shape.MultiplyTN(
            T_surface_inv, T_rel_tr_times_lambda_shape);
        d_fs_d_lambda_gp.MultiplyTN(
            d_psi_surface_d_q_surface, T_surface_mtr_times_T_rel_tr_times_lambda_shape);
        d_fs_d_lambda_gp.Scale(projected_gauss_point.GetGaussWeight() * segment_jacobian);

        // Constraint vector.
        g_gp.PutScalar(0.0);
        for (unsigned int i_row = 0; i_row < mortar::n_dof_; i_row++)
          for (unsigned int i_col = 0; i_col < 3; i_col++)
            g_gp(i_row) += lambda_shape_functions_full(i_col, i_row) * psi_rel(i_col);
        g_gp.Scale(projected_gauss_point.GetGaussWeight() * segment_jacobian);

        // Derivatives of constraint vector.
        T_beam_times_I_beam_tilde_full.Multiply(T_beam, I_beam_tilde_full);

        for (unsigned int i_lambda = 0; i_lambda < mortar::n_dof_; i_lambda++)
          for (unsigned int i_psi = 0; i_psi < 3; i_psi++)
            d_g_d_psi_beam(i_lambda, i_psi) = g_gp(i_lambda).dx(i_psi);
        d_g_d_psi_beam_times_T_beam_I.Multiply(d_g_d_psi_beam, T_beam_times_I_beam_tilde_full);

        for (unsigned int i_lambda = 0; i_lambda < mortar::n_dof_; i_lambda++)
          for (unsigned int i_surface = 0; i_surface < surface::n_dof_; i_surface++)
            d_g_d_q_surface(i_lambda, i_surface) = g_gp(i_lambda).dx(3 + i_surface);

        // Add to output matrices and vector.
        local_g += CORE::FADUTILS::CastToDouble(g_gp);
        local_GB += d_g_d_psi_beam_times_T_beam_I;
        local_GS += d_g_d_q_surface;
        local_FB += d_fb_d_lambda_gp;
        local_FS += d_fs_d_lambda_gp;

        // Calculate the scaling entries.
        for (unsigned int i_mortar_node = 0; i_mortar_node < mortar::n_nodes_; i_mortar_node++)
          for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
            local_kappa(i_mortar_node * 3 + i_dim) += lambda_shape_functions(i_mortar_node) *
                                                      projected_gauss_point.GetGaussWeight() *
                                                      segment_jacobian;
      }
    }
  }

  // Get the rotational GIDs of the surface and beam.
  std::vector<int> gid_surface;
  LINALG::Matrix<n_dof_rot_, 1, int> gid_rot;
  GetPairRotationalGIDs(discret, gid_surface, gid_rot);

  // Get the Lagrange multiplier GIDs.
  std::vector<int> lambda_gid_rot;
  GetMortarGID(mortar_manager, this, mortar::n_dof_, mortar::n_dof_, nullptr, &lambda_gid_rot);

  // Assemble into the global vectors
  global_constraint.SumIntoGlobalValues(lambda_gid_rot.size(), lambda_gid_rot.data(), local_g.A());
  global_kappa.SumIntoGlobalValues(lambda_gid_rot.size(), lambda_gid_rot.data(), local_kappa.A());
  local_kappa.PutScalar(1.0);
  global_lambda_active.SumIntoGlobalValues(
      lambda_gid_rot.size(), lambda_gid_rot.data(), local_kappa.A());

  // Assemble into global matrices.
  for (unsigned int i_dof_lambda = 0; i_dof_lambda < mortar::n_dof_; i_dof_lambda++)
  {
    for (unsigned int i_dof_rot = 0; i_dof_rot < n_dof_rot_; i_dof_rot++)
    {
      global_GB.FEAssemble(
          local_GB(i_dof_lambda, i_dof_rot), lambda_gid_rot[i_dof_lambda], gid_rot(i_dof_rot));
      global_FB.FEAssemble(
          local_FB(i_dof_rot, i_dof_lambda), gid_rot(i_dof_rot), lambda_gid_rot[i_dof_lambda]);
    }
    for (unsigned int i_dof_surface = 0; i_dof_surface < surface::n_dof_; i_dof_surface++)
    {
      global_GS.FEAssemble(local_GS(i_dof_lambda, i_dof_surface), lambda_gid_rot[i_dof_lambda],
          gid_surface[i_dof_surface]);
      global_FS.FEAssemble(local_FS(i_dof_surface, i_dof_lambda), gid_surface[i_dof_surface],
          lambda_gid_rot[i_dof_lambda]);
    }
  }
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarRotationFAD<scalar_type, beam, surface,
    mortar>::GetPairRotationalGIDs(const ::DRT::Discretization& discret,
    std::vector<int>& gid_surface, LINALG::Matrix<n_dof_rot_, 1, int>& gid_rot) const
{
  // Get the GIDs of the surface and beam.
  std::vector<int> lm_beam, lmowner, lmstride;
  this->Element1()->LocationVector(discret, lm_beam, lmowner, lmstride);
  this->face_element_->GetDrtFaceElement()->LocationVector(discret, gid_surface, lmowner, lmstride);
  std::array<int, n_dof_rot_> rot_dof_indices = {3, 4, 5, 12, 13, 14, 18, 19, 20};
  for (unsigned int i = 0; i < n_dof_rot_; i++) gid_rot(i) = lm_beam[rot_dof_indices[i]];
}


/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortarRotation(const bool rotational_coupling)
{
  using namespace BEAMINTERACTION;
  using namespace GEOMETRYPAIR;

  if (!rotational_coupling)
    return Teuchos::rcp(
        new BeamToSolidSurfaceMeshtyingPairMortarFAD<scalar_type, beam, surface, mortar>());
  else
    return Teuchos::rcp(
        new BeamToSolidSurfaceMeshtyingPairMortarRotationFAD<scalar_type, beam, surface, mortar>());
}

/**
 *
 */
template <typename mortar>
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortar(
    const DRT::Element::DiscretizationType surface_shape, const bool rotational_coupling)
{
  using namespace BEAMINTERACTION;
  using namespace GEOMETRYPAIR;

  switch (surface_shape)
  {
    case DRT::Element::tri3:
      return BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortarRotation<
          line_to_surface_patch_scalar_type, t_hermite, t_tri3, mortar>(rotational_coupling);
    case DRT::Element::tri6:
      return BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortarRotation<
          line_to_surface_patch_scalar_type, t_hermite, t_tri6, mortar>(rotational_coupling);
    case DRT::Element::quad4:
      return BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortarRotation<
          line_to_surface_patch_scalar_type, t_hermite, t_quad4, mortar>(rotational_coupling);
    case DRT::Element::quad8:
      return BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortarRotation<
          line_to_surface_patch_scalar_type, t_hermite, t_quad8, mortar>(rotational_coupling);
    case DRT::Element::quad9:
      return BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortarRotation<
          line_to_surface_patch_scalar_type, t_hermite, t_quad9, mortar>(rotational_coupling);
    case DRT::Element::nurbs9:
      return BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortarRotation<
          line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite, t_nurbs9,
          mortar>(rotational_coupling);
    default:
      dserror("Wrong element type for surface element.");
      return Teuchos::null;
  }
}

/**
 *
 */
template <typename mortar>
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortarXVolume(
    const DRT::Element::DiscretizationType surface_shape, const bool rotational_coupling)
{
  using namespace BEAMINTERACTION;
  using namespace GEOMETRYPAIR;

  switch (surface_shape)
  {
    case DRT::Element::quad4:
      return BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortarRotation<
          line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex8>, t_hermite, t_quad4,
          mortar>(rotational_coupling);
    case DRT::Element::quad8:
      return BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortarRotation<
          line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex20>, t_hermite, t_quad8,
          mortar>(rotational_coupling);
    case DRT::Element::quad9:
      return BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortarRotation<
          line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex27>, t_hermite, t_quad9,
          mortar>(rotational_coupling);
    default:
      dserror("Wrong element type for surface element.");
      return Teuchos::null;
  }
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarFADFactory(
    const DRT::Element::DiscretizationType surface_shape,
    const INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions mortar_shapefunction,
    const bool rotational_coupling,
    const INPAR::GEOMETRYPAIR::SurfaceNormals surface_normal_strategy)
{
  using namespace GEOMETRYPAIR;

  if (surface_normal_strategy == INPAR::GEOMETRYPAIR::SurfaceNormals::standard)
  {
    switch (mortar_shapefunction)
    {
      case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line2:
      {
        return BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortar<t_line2>(
            surface_shape, rotational_coupling);
      }
      case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line3:
      {
        return BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortar<t_line3>(
            surface_shape, rotational_coupling);
      }
      case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line4:
      {
        return BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortar<t_line4>(
            surface_shape, rotational_coupling);
      }
      default:
        dserror("Wrong mortar shape function.");
    }
  }
  else if (surface_normal_strategy == INPAR::GEOMETRYPAIR::SurfaceNormals::extended_volume)
  {
    switch (mortar_shapefunction)
    {
      case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line2:
        return BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortarXVolume<t_line2>(
            surface_shape, rotational_coupling);
      case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line3:
        return BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortarXVolume<t_line3>(
            surface_shape, rotational_coupling);
      case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line4:
        return BeamToSolidSurfaceMeshtyingPairMortarFADFactoryMortarXVolume<t_line4>(
            surface_shape, rotational_coupling);
      default:
        dserror("Wrong mortar shape function.");
    }
  }
  else
    dserror("Surface normal strategy not recognized.");

  return Teuchos::null;
}
