/*!

\brief Meshtying element for meshtying between a 3D beam and a 3D solid element using mortar shape
functions for the traction.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_volume_meshtying_pair_mortar.H"

#include "beam_to_solid_vtu_output_writer_base.H"
#include "beam_to_solid_vtu_output_writer_visualization.H"
#include "beam_to_solid_volume_meshtying_vtk_output_params.H"
#include "beam_to_solid_mortar_manager.H"

#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_geometry_pair/geometry_pair_element_types.H"
#include "../drt_geometry_pair/geometry_pair_line_to_volume.H"


/**
 *
 */
template <typename beam, typename solid, typename mortar>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam, solid,
    mortar>::BeamToSolidVolumeMeshtyingPairMortar()
    : BeamToSolidVolumeMeshtyingPairBase<beam, solid>()
{
  // Empty constructor.
}


/**
 *
 */
template <typename beam, typename solid, typename mortar>
bool BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam, solid, mortar>::EvaluateDM(
    LINALG::SerialDenseMatrix& local_D, LINALG::SerialDenseMatrix& local_M,
    LINALG::SerialDenseVector& local_kappa)
{
  // Call Evaluate on the geometry Pair. Only do this once for meshtying.
  if (!this->meshtying_is_evaluated_)
  {
    this->CastGeometryPair()->Evaluate(
        this->ele1posref_, this->ele2posref_, this->line_to_volume_segments_);
    this->meshtying_is_evaluated_ = true;
  }

  // If there are no intersection segments, return no contact status.
  if (this->line_to_volume_segments_.size() == 0) return false;

  // Initialize variables for local mortar matrices.
  LINALG::TMatrix<double, mortar::n_dof_, beam::n_dof_> D(true);
  LINALG::TMatrix<double, mortar::n_dof_, solid::n_dof_> M(true);
  LINALG::TMatrix<double, mortar::n_dof_, 1> kappa(true);

  // Initialize variables for shape function values.
  LINALG::TMatrix<double, 1, mortar::n_nodes_ * mortar::n_val_> N_mortar(true);
  LINALG::TMatrix<double, 1, beam::n_nodes_ * beam::n_val_> N_beam(true);
  LINALG::TMatrix<double, 1, solid::n_nodes_ * solid::n_val_> N_solid(true);

  // Initialize variable for beam position derivative.
  LINALG::TMatrix<double, 3, 1> dr_beam_ref(true);

  // Initialize scalar variables.Clear
  double segment_jacobian, beam_segmentation_factor;

  // Calculate the mortar matrices.
  // Loop over segments.
  for (unsigned int i_segment = 0; i_segment < this->line_to_volume_segments_.size(); i_segment++)
  {
    // Factor to account for a segment length not from -1 to 1.
    beam_segmentation_factor = 0.5 * this->line_to_volume_segments_[i_segment].GetSegmentLength();

    // Gauss point loop.
    for (unsigned int i_gp = 0;
         i_gp < this->line_to_volume_segments_[i_segment].GetProjectionPoints().size(); i_gp++)
    {
      // Get the current Gauss point.
      const GEOMETRYPAIR::ProjectionPointLineToVolume<double>& projected_gauss_point =
          this->line_to_volume_segments_[i_segment].GetProjectionPoints()[i_gp];

      // Get the jacobian in the reference configuration.
      GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(
          projected_gauss_point.GetEta(), this->ele1posref_, dr_beam_ref, this->Element1());

      // Jacobian including the segment length.
      segment_jacobian = dr_beam_ref.Norm2() * beam_segmentation_factor;

      // Get the shape function matrices.
      N_mortar.Clear();
      N_beam.Clear();
      N_solid.Clear();
      mortar::EvaluateShapeFunction(
          N_mortar, projected_gauss_point.GetEta(), std::integral_constant<unsigned int, 1>{});
      beam::EvaluateShapeFunction(N_beam, projected_gauss_point.GetEta(),
          std::integral_constant<unsigned int, 1>{}, this->Element1());
      solid::EvaluateShapeFunction(N_solid, projected_gauss_point.GetXi(),
          std::integral_constant<unsigned int, 3>{}, this->Element2());

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
          for (unsigned int i_solid_node = 0; i_solid_node < solid::n_nodes_; i_solid_node++)
            for (unsigned int i_solid_val = 0; i_solid_val < solid::n_val_; i_solid_val++)
              for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
                M(i_mortar_node * mortar::n_val_ * 3 + i_mortar_val * 3 + i_dim,
                    i_solid_node * solid::n_val_ * 3 + i_solid_val * 3 + i_dim) +=
                    N_mortar(i_mortar_node * mortar::n_val_ + i_mortar_val) *
                    N_solid(i_solid_node * solid::n_val_ + i_solid_val) *
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
  local_M.Shape(mortar::n_dof_, solid::n_dof_);
  local_kappa.Shape(mortar::n_dof_, 1);
  for (unsigned int i_row = 0; i_row < mortar::n_dof_; i_row++)
    for (unsigned int i_col = 0; i_col < beam::n_dof_; i_col++)
      local_D(i_row, i_col) = D(i_row, i_col);
  for (unsigned int i_row = 0; i_row < mortar::n_dof_; i_row++)
    for (unsigned int i_col = 0; i_col < solid::n_dof_; i_col++)
      local_M(i_row, i_col) = M(i_row, i_col);
  for (unsigned int i_row = 0; i_row < mortar::n_dof_; i_row++) local_kappa(i_row) = kappa(i_row);

  // If we get to this point, the pair has a mortar contribution.
  return true;
}

/**
 *
 */
template <typename beam, typename solid, typename mortar>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam, solid,
    mortar>::GetPairVisualization(Teuchos::RCP<BeamToSolidVtuOutputWriterBase> visualization_writer,
    const Teuchos::ParameterList& visualization_params) const
{
  // Get visualization of base method.
  BeamToSolidVolumeMeshtyingPairBase<beam, solid>::GetPairVisualization(
      visualization_writer, visualization_params);


  Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_discret =
      visualization_writer->GetVisualizationWriter("mortar");
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_continuous =
      visualization_writer->GetVisualizationWriter("mortar-continuous");
  if (visualization_discret != Teuchos::null || visualization_continuous != Teuchos::null)
  {
    // Setup variables.
    LINALG::TMatrix<double, mortar::n_dof_, 1> q_lambda;
    LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> X;
    LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> r;
    LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> u;
    LINALG::TMatrix<double, 3, 1> lambda_discret;
    LINALG::TMatrix<double, 3, 1> xi_mortar_node;

    // Get the mortar manager and the global lambda vector, those objects will be used to get the
    // discrete Lagrange multiplier values for this pair.
    Teuchos::RCP<const BEAMINTERACTION::BeamToSolidMortarManager> mortar_manager =
        visualization_params.get<Teuchos::RCP<const BEAMINTERACTION::BeamToSolidMortarManager>>(
            "mortar_manager");
    Teuchos::RCP<Epetra_Vector> lambda =
        visualization_params.get<Teuchos::RCP<Epetra_Vector>>("lambda");

    // Get the lambda GIDs of this pair.
    Teuchos::RCP<const BeamContactPair> this_rcp = Teuchos::rcp(this, false);
    std::vector<int> lambda_row;
    std::vector<double> lambda_pair;
    mortar_manager->LocationVector(this_rcp, lambda_row);
    DRT::UTILS::ExtractMyValues(*lambda, lambda_pair, lambda_row);
    for (unsigned int i_dof; i_dof < mortar::n_dof_; i_dof++) q_lambda(i_dof) = lambda_pair[i_dof];

    // Add the discrete values of the Lagrange multipliers.
    if (visualization_discret != Teuchos::null)
    {
      // Get the visualization vectors.
      std::vector<double>& point_coordinates =
          visualization_discret->GetMutablePointCoordinateVector();
      std::vector<double>& displacement =
          visualization_discret->GetMutablePointDataVector("displacement");
      std::vector<double>& lambda_vis = visualization_discret->GetMutablePointDataVector("lambda");

      for (unsigned int i_node = 0; i_node < mortar::n_nodes_; i_node++)
      {
        // Get the local coordinate of this node.
        xi_mortar_node = DRT::UTILS::getNodeCoordinates(i_node, mortar::discretization_);

        // Get position and displacement of the mortar node.
        GEOMETRYPAIR::EvaluatePosition<beam>(
            xi_mortar_node(0), this->ele1pos_, r, this->Element1());
        GEOMETRYPAIR::EvaluatePosition<beam>(
            xi_mortar_node(0), this->ele1posref_, X, this->Element1());
        u = r;
        u -= X;

        // Get the discrete Lagrangian multiplier.
        GEOMETRYPAIR::EvaluatePosition<mortar>(xi_mortar_node(0), q_lambda, lambda_discret);

        // Add to output data.
        for (unsigned int dim = 0; dim < 3; dim++)
        {
          point_coordinates.push_back(FADUTILS::CastToDouble(X(dim)));
          displacement.push_back(FADUTILS::CastToDouble(u(dim)));
          lambda_vis.push_back(FADUTILS::CastToDouble(lambda_discret(dim)));
        }
      }
    }


    // Add the continuous values for the Lagrange multipliers.
    if (visualization_continuous != Teuchos::null)
    {
      const unsigned int mortar_segments =
          visualization_params
              .get<Teuchos::RCP<const BeamToSolidVolumeMeshtyingVtkOutputParams>>(
                  "output_params_ptr")
              ->GetMortarLambdaContinuousSegments();
      double xi;
      std::vector<double>& point_coordinates =
          visualization_continuous->GetMutablePointCoordinateVector(
              (mortar_segments + 1) * 3 * this->line_to_volume_segments_.size());
      std::vector<double>& displacement = visualization_continuous->GetMutablePointDataVector(
          "displacement", (mortar_segments + 1) * 3 * this->line_to_volume_segments_.size());
      std::vector<double>& lambda_vis = visualization_continuous->GetMutablePointDataVector(
          "lambda", (mortar_segments + 1) * 3 * this->line_to_volume_segments_.size());
      std::vector<uint8_t>& cell_type = visualization_continuous->GetMutableCellTypeVector();
      std::vector<int32_t>& cell_offset = visualization_continuous->GetMutableCellOffsetVector();

      for (const auto& segment : this->line_to_volume_segments_)
      {
        for (unsigned int i_curve_segment = 0; i_curve_segment <= mortar_segments;
             i_curve_segment++)
        {
          // Get the position, displacement and lambda value at the current point.
          xi = segment.GetEtaA() +
               i_curve_segment * (segment.GetEtaB() - segment.GetEtaA()) / (double)mortar_segments;
          GEOMETRYPAIR::EvaluatePosition<beam>(xi, this->ele1pos_, r, this->Element1());
          GEOMETRYPAIR::EvaluatePosition<beam>(xi, this->ele1posref_, X, this->Element1());
          u = r;
          u -= X;
          GEOMETRYPAIR::EvaluatePosition<mortar>(xi, q_lambda, lambda_discret);

          // Add to output data.
          for (unsigned int dim = 0; dim < 3; dim++)
          {
            point_coordinates.push_back(FADUTILS::CastToDouble(X(dim)));
            displacement.push_back(FADUTILS::CastToDouble(u(dim)));
            lambda_vis.push_back(FADUTILS::CastToDouble(lambda_discret(dim)));
          }
        }

        // Add the cell for this segment (poly line).
        cell_type.push_back(4);
        cell_offset.push_back(point_coordinates.size() / 3);
      }
    }
  }
}


/**
 *
 */
template <typename beam, typename solid, typename mortar>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam, solid,
    mortar>::EvaluatePenaltyForce(const LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1>& r_beam,
    const LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1>& r_solid,
    LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1>& force) const
{
  force.PutScalar(0.);
}


/**
 * Explicit template initialization of template class.
 */
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line2>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex20, GEOMETRYPAIR::t_line2>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex27, GEOMETRYPAIR::t_line2>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet4, GEOMETRYPAIR::t_line2>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet10, GEOMETRYPAIR::t_line2>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_nurbs27, GEOMETRYPAIR::t_line2>;

template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line3>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex20, GEOMETRYPAIR::t_line3>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex27, GEOMETRYPAIR::t_line3>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet4, GEOMETRYPAIR::t_line3>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet10, GEOMETRYPAIR::t_line3>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_nurbs27, GEOMETRYPAIR::t_line3>;

template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line4>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex20, GEOMETRYPAIR::t_line4>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex27, GEOMETRYPAIR::t_line4>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet4, GEOMETRYPAIR::t_line4>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet10, GEOMETRYPAIR::t_line4>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_nurbs27, GEOMETRYPAIR::t_line4>;
