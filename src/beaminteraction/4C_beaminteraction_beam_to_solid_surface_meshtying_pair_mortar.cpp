/*----------------------------------------------------------------------*/
/*! \file

\brief Mortar mesh tying element for between a 3D beam and a surface element.

\level 3
*/


#include "4C_beaminteraction_beam_to_solid_surface_meshtying_pair_mortar.hpp"

#include "4C_beaminteraction_beam_to_solid_surface_meshtying_params.hpp"
#include "4C_beaminteraction_beam_to_solid_utils.hpp"
#include "4C_beaminteraction_contact_params.hpp"
#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_element_faces.hpp"
#include "4C_geometry_pair_line_to_surface.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
template <typename beam, typename surface, typename mortar>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortar<beam, surface,
    mortar>::BeamToSolidSurfaceMeshtyingPairMortar()
    : base_class()
{
  // Empty constructor.
}

/**
 *
 */
template <typename beam, typename surface, typename mortar>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortar<beam, surface,
    mortar>::evaluate_and_assemble_mortar_contributions(const Discret::Discretization& discret,
    const BeamToSolidMortarManager* mortar_manager, Core::LinAlg::SparseMatrix& global_G_B,
    Core::LinAlg::SparseMatrix& global_G_S, Core::LinAlg::SparseMatrix& global_FB_L,
    Core::LinAlg::SparseMatrix& global_FS_L, Epetra_FEVector& global_constraint,
    Epetra_FEVector& global_kappa, Epetra_FEVector& global_lambda_active,
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector)
{
  // Call Evaluate on the geometry Pair. Only do this once for meshtying.
  if (!this->meshtying_is_evaluated_)
  {
    this->cast_geometry_pair()->Evaluate(this->ele1posref_,
        this->face_element_->get_face_reference_element_data(), this->line_to_3D_segments_);
    this->meshtying_is_evaluated_ = true;
  }

  // If there are no intersection segments, return no contact status.
  if (this->line_to_3D_segments_.size() == 0) return;

  // Initialize variables for local mortar matrices.
  Core::LinAlg::Matrix<mortar::n_dof_, beam::n_dof_, double> local_D(false);
  Core::LinAlg::Matrix<mortar::n_dof_, surface::n_dof_, double> local_M(false);
  Core::LinAlg::Matrix<mortar::n_dof_, 1, double> local_kappa(false);
  Core::LinAlg::Matrix<mortar::n_dof_, 1, double> local_constraint(false);

  // Evaluate the local mortar contributions.
  evaluate_dm(local_D, local_M, local_kappa, local_constraint);

  // Assemble into global matrices.
  AssembleLocalMortarContributions<beam, surface, mortar>(this, discret, mortar_manager, global_G_B,
      global_G_S, global_FB_L, global_FS_L, global_constraint, global_kappa, global_lambda_active,
      local_D, local_M, local_kappa, local_constraint);
}

/**
 *
 */
template <typename beam, typename surface, typename mortar>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortar<beam, surface, mortar>::evaluate_dm(
    Core::LinAlg::Matrix<mortar::n_dof_, beam::n_dof_, double>& local_D,
    Core::LinAlg::Matrix<mortar::n_dof_, surface::n_dof_, double>& local_M,
    Core::LinAlg::Matrix<mortar::n_dof_, 1, double>& local_kappa,
    Core::LinAlg::Matrix<mortar::n_dof_, 1, double>& local_constraint) const
{
  // Initialize the local mortar matrices.
  local_D.PutScalar(0.0);
  local_M.PutScalar(0.0);
  local_kappa.PutScalar(0.0);
  local_constraint.PutScalar(0.0);

  // Initialize variables for shape function values.
  Core::LinAlg::Matrix<1, mortar::n_nodes_ * mortar::n_val_, double> N_mortar(true);
  Core::LinAlg::Matrix<1, beam::n_nodes_ * beam::n_val_, double> N_beam(true);
  Core::LinAlg::Matrix<1, surface::n_nodes_ * surface::n_val_, double> N_surface(true);

  // Initialize variable for beam position derivative.
  Core::LinAlg::Matrix<3, 1, double> dr_beam_ref(true);

  // Initialize scalar variables.Clear
  double segment_jacobian = 0.0;
  double beam_segmentation_factor = 0.0;

  // Calculate the mortar matrices.
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

      // Get the shape function matrices.
      N_mortar.Clear();
      N_beam.Clear();
      N_surface.Clear();
      GEOMETRYPAIR::EvaluateShapeFunction<mortar>::Evaluate(
          N_mortar, projected_gauss_point.GetEta());
      GEOMETRYPAIR::EvaluateShapeFunction<beam>::Evaluate(
          N_beam, projected_gauss_point.GetEta(), this->ele1pos_.shape_function_data_);
      GEOMETRYPAIR::EvaluateShapeFunction<surface>::Evaluate(N_surface,
          projected_gauss_point.GetXi(),
          this->face_element_->GetFaceElementData().shape_function_data_);

      // Fill in the local templated mortar matrix D.
      for (unsigned int i_mortar_node = 0; i_mortar_node < mortar::n_nodes_; i_mortar_node++)
        for (unsigned int i_mortar_val = 0; i_mortar_val < mortar::n_val_; i_mortar_val++)
          for (unsigned int i_beam_node = 0; i_beam_node < beam::n_nodes_; i_beam_node++)
            for (unsigned int i_beam_val = 0; i_beam_val < beam::n_val_; i_beam_val++)
              for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
                local_D(i_mortar_node * mortar::n_val_ * 3 + i_mortar_val * 3 + i_dim,
                    i_beam_node * beam::n_val_ * 3 + i_beam_val * 3 + i_dim) +=
                    N_mortar(i_mortar_node * mortar::n_val_ + i_mortar_val) *
                    N_beam(i_beam_node * beam::n_val_ + i_beam_val) *
                    projected_gauss_point.GetGaussWeight() * segment_jacobian;

      // Fill in the local templated mortar matrix M.
      for (unsigned int i_mortar_node = 0; i_mortar_node < mortar::n_nodes_; i_mortar_node++)
        for (unsigned int i_mortar_val = 0; i_mortar_val < mortar::n_val_; i_mortar_val++)
          for (unsigned int i_surface_node = 0; i_surface_node < surface::n_nodes_;
               i_surface_node++)
            for (unsigned int i_surface_val = 0; i_surface_val < surface::n_val_; i_surface_val++)
              for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
                local_M(i_mortar_node * mortar::n_val_ * 3 + i_mortar_val * 3 + i_dim,
                    i_surface_node * surface::n_val_ * 3 + i_surface_val * 3 + i_dim) +=
                    N_mortar(i_mortar_node * mortar::n_val_ + i_mortar_val) *
                    N_surface(i_surface_node * surface::n_val_ + i_surface_val) *
                    projected_gauss_point.GetGaussWeight() * segment_jacobian;

      // Fill in the local templated mortar scaling vector kappa.
      for (unsigned int i_mortar_node = 0; i_mortar_node < mortar::n_nodes_; i_mortar_node++)
        for (unsigned int i_mortar_val = 0; i_mortar_val < mortar::n_val_; i_mortar_val++)
          for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
            local_kappa(i_mortar_node * mortar::n_val_ * 3 + i_mortar_val * 3 + i_dim) +=
                N_mortar(i_mortar_node * mortar::n_val_ + i_mortar_val) *
                projected_gauss_point.GetGaussWeight() * segment_jacobian;
    }
  }

  // Add the local constraint contributions. For this we multiply the local mortar matrices with the
  // positions / displacements to get the actual constraint terms for this pair.
  Core::LinAlg::Matrix<beam::n_dof_, 1, double> beam_coupling_dof(true);
  Core::LinAlg::Matrix<surface::n_dof_, 1, double> surface_coupling_dof(true);
  switch (this->Params()->beam_to_solid_surface_meshtying_params()->GetCouplingType())
  {
    case Inpar::BeamToSolid::BeamToSolidSurfaceCoupling::reference_configuration_forced_to_zero:
    {
      beam_coupling_dof = Core::FADUtils::CastToDouble(this->ele1pos_.element_position_);
      surface_coupling_dof =
          Core::FADUtils::CastToDouble(this->face_element_->GetFaceElementData().element_position_);
      break;
    }
    case Inpar::BeamToSolid::BeamToSolidSurfaceCoupling::displacement:
    {
      beam_coupling_dof = Core::FADUtils::CastToDouble(this->ele1pos_.element_position_);
      beam_coupling_dof -= this->ele1posref_.element_position_;
      surface_coupling_dof =
          Core::FADUtils::CastToDouble(this->face_element_->GetFaceElementData().element_position_);
      surface_coupling_dof -=
          this->face_element_->get_face_reference_element_data().element_position_;
      break;
    }
    default:
      FOUR_C_THROW("Wrong coupling type.");
  }
  for (unsigned int i_lambda = 0; i_lambda < mortar::n_dof_; i_lambda++)
  {
    for (unsigned int i_beam = 0; i_beam < beam::n_dof_; i_beam++)
      local_constraint(i_lambda) += local_D(i_lambda, i_beam) * beam_coupling_dof(i_beam);
    for (unsigned int i_surface = 0; i_surface < surface::n_dof_; i_surface++)
      local_constraint(i_lambda) -= local_M(i_lambda, i_surface) * surface_coupling_dof(i_surface);
  }
}


/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarFactory(
    const Core::FE::CellType surface_shape,
    const Inpar::BeamToSolid::BeamToSolidMortarShapefunctions mortar_shapefunction)
{
  using namespace GEOMETRYPAIR;

  switch (mortar_shapefunction)
  {
    case Inpar::BeamToSolid::BeamToSolidMortarShapefunctions::line2:
    {
      switch (surface_shape)
      {
        case Core::FE::CellType::tri3:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri3, t_line2>());
        case Core::FE::CellType::tri6:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri6, t_line2>());
        case Core::FE::CellType::quad4:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad4, t_line2>());
        case Core::FE::CellType::quad8:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad8, t_line2>());
        case Core::FE::CellType::quad9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad9, t_line2>());
        case Core::FE::CellType::nurbs9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_nurbs9, t_line2>());
        default:
          FOUR_C_THROW("Wrong element type for surface element.");
      }
      break;
    }
    case Inpar::BeamToSolid::BeamToSolidMortarShapefunctions::line3:
    {
      switch (surface_shape)
      {
        case Core::FE::CellType::tri3:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri3, t_line3>());
        case Core::FE::CellType::tri6:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri6, t_line3>());
        case Core::FE::CellType::quad4:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad4, t_line3>());
        case Core::FE::CellType::quad8:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad8, t_line3>());
        case Core::FE::CellType::quad9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad9, t_line3>());
        case Core::FE::CellType::nurbs9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_nurbs9, t_line3>());
        default:
          FOUR_C_THROW("Wrong element type for surface element.");
      }
      break;
    }
    case Inpar::BeamToSolid::BeamToSolidMortarShapefunctions::line4:
    {
      switch (surface_shape)
      {
        case Core::FE::CellType::tri3:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri3, t_line4>());
        case Core::FE::CellType::tri6:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri6, t_line4>());
        case Core::FE::CellType::quad4:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad4, t_line4>());
        case Core::FE::CellType::quad8:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad8, t_line4>());
        case Core::FE::CellType::quad9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad9, t_line4>());
        case Core::FE::CellType::nurbs9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_nurbs9, t_line4>());
        default:
          FOUR_C_THROW("Wrong element type for surface element.");
      }
      break;
    }
    default:
      FOUR_C_THROW("Wrong mortar shape function.");
  }

  return Teuchos::null;
}

FOUR_C_NAMESPACE_CLOSE
