/*----------------------------------------------------------------------*/
/*! \file

\brief Mortar mesh tying element for between a 3D beam and a surface element.

\level 3
*/


#include "beaminteraction_beam_to_solid_surface_meshtying_pair_mortar.H"

#include "beaminteraction_contact_params.H"
#include "beaminteraction_beam_to_solid_surface_meshtying_params.H"
#include "beaminteraction_beam_to_solid_utils.H"
#include "geometry_pair_line_to_surface.H"
#include "geometry_pair_element_functions.H"
#include "geometry_pair_element_faces.H"


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
    mortar>::EvaluateAndAssembleMortarContributions(const DRT::Discretization& discret,
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

  // Initialize variables for local mortar matrices.
  LINALG::Matrix<mortar::n_dof_, beam::n_dof_, double> local_D(false);
  LINALG::Matrix<mortar::n_dof_, surface::n_dof_, double> local_M(false);
  LINALG::Matrix<mortar::n_dof_, 1, double> local_kappa(false);
  LINALG::Matrix<mortar::n_dof_, 1, double> local_constraint(false);

  // Evaluate the local mortar contributions.
  EvaluateDM(local_D, local_M, local_kappa, local_constraint);

  // Assemble into global matrices.
  AssembleLocalMortarContributions<beam, surface, mortar>(this, discret, mortar_manager, global_G_B,
      global_G_S, global_FB_L, global_FS_L, global_constraint, global_kappa, global_lambda_active,
      local_D, local_M, local_kappa, local_constraint);
}

/**
 *
 */
template <typename beam, typename surface, typename mortar>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortar<beam, surface, mortar>::EvaluateDM(
    LINALG::Matrix<mortar::n_dof_, beam::n_dof_, double>& local_D,
    LINALG::Matrix<mortar::n_dof_, surface::n_dof_, double>& local_M,
    LINALG::Matrix<mortar::n_dof_, 1, double>& local_kappa,
    LINALG::Matrix<mortar::n_dof_, 1, double>& local_constraint) const
{
  // Initialize the local mortar matrices.
  local_D.PutScalar(0.0);
  local_M.PutScalar(0.0);
  local_kappa.PutScalar(0.0);
  local_constraint.PutScalar(0.0);

  // Initialize variables for shape function values.
  LINALG::Matrix<1, mortar::n_nodes_ * mortar::n_val_, double> N_mortar(true);
  LINALG::Matrix<1, beam::n_nodes_ * beam::n_val_, double> N_beam(true);
  LINALG::Matrix<1, surface::n_nodes_ * surface::n_val_, double> N_surface(true);

  // Initialize variable for beam position derivative.
  LINALG::Matrix<3, 1, double> dr_beam_ref(true);

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
          projected_gauss_point.GetEta(), this->ele1posref_, dr_beam_ref, this->Element1());

      // Jacobian including the segment length.
      segment_jacobian = dr_beam_ref.Norm2() * beam_segmentation_factor;

      // Get the shape function matrices.
      N_mortar.Clear();
      N_beam.Clear();
      N_surface.Clear();
      mortar::EvaluateShapeFunction(
          N_mortar, projected_gauss_point.GetEta(), std::integral_constant<unsigned int, 1>{});
      beam::EvaluateShapeFunction(N_beam, projected_gauss_point.GetEta(),
          std::integral_constant<unsigned int, 1>{}, this->Element1());
      surface::EvaluateShapeFunction(N_surface, projected_gauss_point.GetXi(),
          std::integral_constant<unsigned int, 2>{}, this->GeometryPair()->Element2());

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
  LINALG::Matrix<beam::n_dof_, 1, double> beam_coupling_dof(true);
  LINALG::Matrix<surface::n_dof_, 1, double> surface_coupling_dof(true);
  switch (this->Params()->BeamToSolidSurfaceMeshtyingParams()->GetCouplingType())
  {
    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::reference_configuration_forced_to_zero:
    {
      for (unsigned int i_beam = 0; i_beam < beam::n_dof_; i_beam++)
        beam_coupling_dof(i_beam) = FADUTILS::CastToDouble(this->ele1pos_(i_beam));
      for (unsigned int i_surface = 0; i_surface < surface::n_dof_; i_surface++)
        surface_coupling_dof(i_surface) =
            FADUTILS::CastToDouble(this->face_element_->GetFacePosition()(i_surface));
      break;
    }
    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::displacement:
    {
      for (unsigned int i_beam = 0; i_beam < beam::n_dof_; i_beam++)
        beam_coupling_dof(i_beam) =
            FADUTILS::CastToDouble(this->ele1pos_(i_beam)) - this->ele1posref_(i_beam);
      for (unsigned int i_surface = 0; i_surface < surface::n_dof_; i_surface++)
        surface_coupling_dof(i_surface) =
            FADUTILS::CastToDouble(this->face_element_->GetFacePosition()(i_surface)) -
            this->face_element_->GetFaceReferencePosition()(i_surface);
      break;
    }
    default:
      dserror("Wrong coupling type.");
  }
  for (unsigned int i_lambda = 0; i_lambda < mortar::n_dof_; i_lambda++)
  {
    for (unsigned int i_beam = 0; i_beam < beam::n_dof_; i_beam++)
      local_constraint(i_lambda) +=
          local_D(i_lambda, i_beam) * FADUTILS::CastToDouble(beam_coupling_dof(i_beam));
    for (unsigned int i_surface = 0; i_surface < surface::n_dof_; i_surface++)
      local_constraint(i_lambda) -=
          local_M(i_lambda, i_surface) * FADUTILS::CastToDouble(surface_coupling_dof(i_surface));
  }
}


/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarFactory(
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
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri3, t_line2>());
        case DRT::Element::tri6:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri6, t_line2>());
        case DRT::Element::quad4:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad4, t_line2>());
        case DRT::Element::quad8:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad8, t_line2>());
        case DRT::Element::quad9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad9, t_line2>());
        case DRT::Element::nurbs9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_nurbs9, t_line2>());
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
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri3, t_line3>());
        case DRT::Element::tri6:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri6, t_line3>());
        case DRT::Element::quad4:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad4, t_line3>());
        case DRT::Element::quad8:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad8, t_line3>());
        case DRT::Element::quad9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad9, t_line3>());
        case DRT::Element::nurbs9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_nurbs9, t_line3>());
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
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri3, t_line4>());
        case DRT::Element::tri6:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri6, t_line4>());
        case DRT::Element::quad4:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad4, t_line4>());
        case DRT::Element::quad8:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad8, t_line4>());
        case DRT::Element::quad9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad9, t_line4>());
        case DRT::Element::nurbs9:
          return Teuchos::rcp(
              new BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_nurbs9, t_line4>());
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
