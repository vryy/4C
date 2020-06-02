/*----------------------------------------------------------------------*/
/*! \file

\brief Mortar mesh tying element for between a 3D beam and a surface element.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_surface_meshtying_pair_mortar.H"

#include "beam_contact_params.H"
#include "beam_to_solid_surface_meshtying_params.H"
#include "../drt_geometry_pair/geometry_pair_line_to_surface.H"
#include "../drt_geometry_pair/geometry_pair_element_functions.H"
#include "../drt_geometry_pair/geometry_pair_element_faces.H"


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
bool BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortar<beam, surface, mortar>::EvaluateDM(
    LINALG::SerialDenseMatrix& local_D, LINALG::SerialDenseMatrix& local_M,
    LINALG::SerialDenseVector& local_kappa, LINALG::SerialDenseVector& local_constraint_offset)
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
  if (this->line_to_3D_segments_.size() == 0) return false;

  // Initialize variables for local mortar matrices.
  LINALG::Matrix<mortar::n_dof_, beam::n_dof_, double> D(true);
  LINALG::Matrix<mortar::n_dof_, surface::n_dof_, double> M(true);
  LINALG::Matrix<mortar::n_dof_, 1, double> kappa(true);

  // Initialize variables for shape function values.
  LINALG::Matrix<1, mortar::n_nodes_ * mortar::n_val_, double> N_mortar(true);
  LINALG::Matrix<1, beam::n_nodes_ * beam::n_val_, double> N_beam(true);
  LINALG::Matrix<1, surface::n_nodes_ * surface::n_val_, double> N_surface(true);

  // Initialize variable for beam position derivative.
  LINALG::Matrix<3, 1, double> dr_beam_ref(true);

  // Initialize scalar variables.Clear
  double segment_jacobian, beam_segmentation_factor;

  // Calculate the mortar matrices.
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
                D(i_mortar_node * mortar::n_val_ * 3 + i_mortar_val * 3 + i_dim,
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
                M(i_mortar_node * mortar::n_val_ * 3 + i_mortar_val * 3 + i_dim,
                    i_surface_node * surface::n_val_ * 3 + i_surface_val * 3 + i_dim) +=
                    N_mortar(i_mortar_node * mortar::n_val_ + i_mortar_val) *
                    N_surface(i_surface_node * surface::n_val_ + i_surface_val) *
                    projected_gauss_point.GetGaussWeight() * segment_jacobian;

      // Fill in the local templated mortar scaling vector kappa.
      for (unsigned int i_mortar_node = 0; i_mortar_node < mortar::n_nodes_; i_mortar_node++)
        for (unsigned int i_mortar_val = 0; i_mortar_val < mortar::n_val_; i_mortar_val++)
          for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
            kappa(i_mortar_node * mortar::n_val_ * 3 + i_mortar_val * 3 + i_dim) +=
                N_mortar(i_mortar_node * mortar::n_val_ + i_mortar_val) *
                projected_gauss_point.GetGaussWeight() * segment_jacobian;
    }
  }

  // Put the values of the templated matrices into the local matrices that are returned.
  local_D.Shape(mortar::n_dof_, beam::n_dof_);
  local_M.Shape(mortar::n_dof_, surface::n_dof_);
  local_kappa.Shape(mortar::n_dof_, 1);
  for (unsigned int i_row = 0; i_row < mortar::n_dof_; i_row++)
    for (unsigned int i_col = 0; i_col < beam::n_dof_; i_col++)
      local_D(i_row, i_col) = D(i_row, i_col);
  for (unsigned int i_row = 0; i_row < mortar::n_dof_; i_row++)
    for (unsigned int i_col = 0; i_col < surface::n_dof_; i_col++)
      local_M(i_row, i_col) = M(i_row, i_col);
  for (unsigned int i_row = 0; i_row < mortar::n_dof_; i_row++) local_kappa(i_row) = kappa(i_row);

  // Calculate the constraint offset vector.
  local_constraint_offset.Size(mortar::n_dof_);
  switch (this->Params()->BeamToSolidSurfaceMeshtyingParams()->GetCouplingType())
  {
    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::reference_configuration_forced_to_zero:
    {
      // Add the reference offset values.
      for (unsigned int i_dof_lambda = 0; i_dof_lambda < mortar::n_dof_; i_dof_lambda++)
      {
        double row_value = 0.0;
        for (unsigned int i_dof_beam = 0; i_dof_beam < beam::n_dof_; i_dof_beam++)
          row_value += local_D(i_dof_lambda, i_dof_beam) * this->ele1posref_(i_dof_beam);
        for (unsigned int i_dof_surface = 0; i_dof_surface < surface::n_dof_; i_dof_surface++)
          row_value -= local_M(i_dof_lambda, i_dof_surface) *
                       this->face_element_->GetFaceReferencePosition()(i_dof_surface);
        local_constraint_offset(i_dof_lambda) = row_value;
      }
      break;
    }
    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::displacement:
    {
      // In this case we do not need to add a constraint offset.
      break;
    }
    default:
      dserror("Wrong coupling type.");
  }

  // If we get to this point, the pair has a mortar contribution.
  return true;
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
