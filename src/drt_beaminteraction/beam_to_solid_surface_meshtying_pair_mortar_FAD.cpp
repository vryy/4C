/*----------------------------------------------------------------------*/
/*! \file

\brief Mortar mesh tying element for between a 3D beam and a surface element, coupling terms are
evaluated with FAD.

\level 3
*/


#include "beam_to_solid_surface_meshtying_pair_mortar_FAD.H"

#include "beaminteraction_calc_utils.H"
#include "beam_to_solid_mortar_manager.H"
#include "../drt_geometry_pair/geometry_pair_line_to_surface.H"
#include "../drt_geometry_pair/geometry_pair_element_faces.H"
#include "../drt_inpar/inpar_beam_to_solid.H"

#include "Epetra_FEVector.h"


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
    mortar>::EvaluateAndAssemble(const DRT::Discretization& discret,
    const BeamToSolidMortarManager* mortar_manager,
    const Teuchos::RCP<Epetra_FEVector>& force_vector,
    const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix, const Epetra_Vector& global_lambda)
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

  // Get the local Lagrange multiplier vector.
  std::vector<int> lambda_gid;
  mortar_manager->LocationVector(this, lambda_gid);
  std::vector<double> local_lambda;
  DRT::UTILS::ExtractMyValues(global_lambda, local_lambda, lambda_gid);
  LINALG::Matrix<mortar::n_dof_, 1, double> q_lambda(local_lambda.data());


  // Initialize variables for local values.
  LINALG::Matrix<3, 1, scalar_type> coupling_vector(true);
  LINALG::Matrix<3, 1, double> lambda(true);
  LINALG::Matrix<3, 1, double> dr_beam_ref(true);
  scalar_type potential = 0.0;

  // Initialize scalar variables.
  double segment_jacobian = 0.0;
  double beam_segmentation_factor = 0.0;

  // Loop over segments to evaluate the coupling potential.
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
        stiffness_matrix->FEAssemble(FADUTILS::CastToDouble(potential.dx(i_dof).dx(j_dof)),
            pair_gid[i_dof], pair_gid[j_dof]);
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarFAD<scalar_type, beam, surface,
    mortar>::EvaluateAndAssembleDM(const DRT::Discretization& discret,
    const BeamToSolidMortarManager* mortar_manager, LINALG::SparseMatrix& global_D,
    LINALG::SparseMatrix& global_M, Epetra_FEVector& global_constraint,
    Epetra_FEVector& global_kappa, Epetra_FEVector& global_lambda_active)
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
  std::vector<int> lambda_gid;
  mortar_manager->LocationVector(this, lambda_gid);

  // Assemble into the global D matrix.
  for (unsigned int i_lambda = 0; i_lambda < mortar::n_dof_; i_lambda++)
    for (unsigned int i_beam = 0; i_beam < beam::n_dof_; i_beam++)
      global_D.FEAssemble(FADUTILS::CastToDouble(constraint_vector(i_lambda).dx(i_beam)),
          lambda_gid[i_lambda], beam_centerline_gid(i_beam));

  // Assemble into the global M matrix.
  for (unsigned int i_lambda = 0; i_lambda < mortar::n_dof_; i_lambda++)
    for (unsigned int i_patch = 0; i_patch < patch_gid.size(); i_patch++)
      global_M.FEAssemble(
          -1.0 * FADUTILS::CastToDouble(constraint_vector(i_lambda).dx(beam::n_dof_ + i_patch)),
          lambda_gid[i_lambda], patch_gid[i_patch]);

  // Assemble into global coupling vector.
  LINALG::Matrix<mortar::n_dof_, 1, double> constraint_vector_double =
      FADUTILS::CastToDouble(constraint_vector);
  global_constraint.SumIntoGlobalValues(
      lambda_gid.size(), lambda_gid.data(), constraint_vector_double.A());

  // Assemble into global kappa vector.
  global_kappa.SumIntoGlobalValues(lambda_gid.size(), lambda_gid.data(), local_kappa.A());

  // Assemble into global lambda active vector.
  local_kappa.PutScalar(1.0);
  global_lambda_active.SumIntoGlobalValues(lambda_gid.size(), lambda_gid.data(), local_kappa.A());
}


/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarFADFactory(
    const DRT::Element::DiscretizationType surface_shape,
    const INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions mortar_shapefunction)
{
  using namespace GEOMETRYPAIR;

  switch (mortar_shapefunction)
  {
    case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line2:
    {
      switch (surface_shape)
      {
        case DRT::Element::tri3:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_tri3, t_line2>());
        case DRT::Element::tri6:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_tri6, t_line2>());
        case DRT::Element::quad4:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad4, t_line2>());
        case DRT::Element::quad8:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad8, t_line2>());
        case DRT::Element::quad9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad9, t_line2>());
        case DRT::Element::nurbs9:
          return Teuchos::rcp(new BeamToSolidSurfaceMeshtyingPairMortarFAD<
              line_to_surface_patch_nurbs_scalar_type<t_hermite, t_nurbs9>, t_hermite, t_nurbs9,
              t_line2>());
        default:
          dserror("Wrong element type for surface element.");
      }
      break;
    }
    case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line3:
    {
      switch (surface_shape)
      {
        case DRT::Element::tri3:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_tri3, t_line3>());
        case DRT::Element::tri6:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_tri6, t_line3>());
        case DRT::Element::quad4:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad4, t_line3>());
        case DRT::Element::quad8:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad8, t_line3>());
        case DRT::Element::quad9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad9, t_line3>());
        case DRT::Element::nurbs9:
          return Teuchos::rcp(new BeamToSolidSurfaceMeshtyingPairMortarFAD<
              line_to_surface_patch_nurbs_scalar_type<t_hermite, t_nurbs9>, t_hermite, t_nurbs9,
              t_line3>());
        default:
          dserror("Wrong element type for surface element.");
      }
      break;
    }
    case INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::line4:
    {
      switch (surface_shape)
      {
        case DRT::Element::tri3:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_tri3, t_line4>());
        case DRT::Element::tri6:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_tri6, t_line4>());
        case DRT::Element::quad4:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad4, t_line4>());
        case DRT::Element::quad8:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad8, t_line4>());
        case DRT::Element::quad9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortarFAD<line_to_surface_patch_scalar_type,
                  t_hermite, t_quad9, t_line4>());
        case DRT::Element::nurbs9:
          return Teuchos::rcp(new BeamToSolidSurfaceMeshtyingPairMortarFAD<
              line_to_surface_patch_nurbs_scalar_type<t_hermite, t_nurbs9>, t_hermite, t_nurbs9,
              t_line4>());
        default:
          dserror("Wrong element type for surface element.");
      }
      break;
    }
    default:
      dserror("Wrong mortar shape function.");
  }

  return Teuchos::null;
}
