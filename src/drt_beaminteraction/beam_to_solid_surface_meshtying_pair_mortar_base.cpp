/*----------------------------------------------------------------------*/
/*! \file

\brief Base mortar mesh tying element for between a 3D beam and a surface element.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_surface_meshtying_pair_mortar_base.H"

#include "beam_to_solid_vtu_output_writer_base.H"
#include "beam_to_solid_vtu_output_writer_visualization.H"
#include "beam_to_solid_surface_vtk_output_params.H"
#include "beam_to_solid_mortar_manager.H"
#include "../drt_geometry_pair/geometry_pair_line_to_surface.H"
#include "../drt_geometry_pair/geometry_pair_element_functions.H"

#include <unordered_set>


/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarBase<scalar_type, beam, surface,
    mortar>::BeamToSolidSurfaceMeshtyingPairMortarBase()
    : base_class()
{
  // Empty constructor.
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarBase<scalar_type, beam, surface,
    mortar>::GetPairVisualization(Teuchos::RCP<BeamToSolidVtuOutputWriterBase> visualization_writer,
    Teuchos::ParameterList& visualization_params) const
{
  // Get visualization of base method.
  base_class::GetPairVisualization(visualization_writer, visualization_params);


  Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_discret =
      visualization_writer->GetVisualizationWriter("btssc-mortar");
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> visualization_continuous =
      visualization_writer->GetVisualizationWriter("btssc-mortar-continuous");
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization>
      visualization_nodal_forces =
          visualization_writer->GetVisualizationWriter("btssc-nodal-forces");
  if (visualization_discret != Teuchos::null or visualization_continuous != Teuchos::null or
      visualization_nodal_forces != Teuchos::null)
  {
    // Setup variables.
    LINALG::Matrix<mortar::n_dof_, 1, double> q_lambda;
    LINALG::Matrix<3, 1, scalar_type> X;
    LINALG::Matrix<3, 1, scalar_type> r;
    LINALG::Matrix<3, 1, scalar_type> u;
    LINALG::Matrix<3, 1, double> lambda_discret;
    LINALG::Matrix<3, 1, double> xi_mortar_node;

    // Get the mortar manager and the global lambda vector, those objects will be used to get the
    // discrete Lagrange multiplier values for this pair.
    Teuchos::RCP<const BEAMINTERACTION::BeamToSolidMortarManager> mortar_manager =
        visualization_params.get<Teuchos::RCP<const BEAMINTERACTION::BeamToSolidMortarManager>>(
            "mortar_manager");
    Teuchos::RCP<Epetra_Vector> lambda =
        visualization_params.get<Teuchos::RCP<Epetra_Vector>>("lambda");

    // Get the lambda GIDs of this pair.
    std::vector<int> lambda_row;
    std::vector<double> lambda_pair;
    mortar_manager->LocationVector(this, lambda_row);
    DRT::UTILS::ExtractMyValues(*lambda, lambda_pair, lambda_row);
    for (unsigned int i_dof = 0; i_dof < mortar::n_dof_; i_dof++)
      q_lambda(i_dof) = lambda_pair[i_dof];

    // Add the discrete values of the Lagrange multipliers.
    if (visualization_discret != Teuchos::null)
    {
      // Check if data for this beam was already written.
      Teuchos::RCP<std::unordered_set<int>> beam_tracker =
          visualization_params.get<Teuchos::RCP<std::unordered_set<int>>>("beam_tracker");

      auto it = beam_tracker->find(this->Element1()->Id());
      if (it == beam_tracker->end())
      {
        // Only do something if this beam element did not write any output yet.

        // Add this element Id to the tracker.
        beam_tracker->insert(this->Element1()->Id());

        // Get the visualization vectors.
        std::vector<double>& point_coordinates =
            visualization_discret->GetMutablePointCoordinateVector();
        std::vector<double>& displacement =
            visualization_discret->GetMutablePointDataVector("displacement");
        std::vector<double>& lambda_vis =
            visualization_discret->GetMutablePointDataVector("lambda");

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
    }


    // Add the continuous values for the Lagrange multipliers.
    if (visualization_continuous != Teuchos::null and this->line_to_3D_segments_.size() > 0)
    {
      const unsigned int mortar_segments =
          visualization_params
              .get<Teuchos::RCP<const BeamToSolidSurfaceVtkOutputParams>>("btssc-output_params_ptr")
              ->GetMortarLambdaContinuousSegments();
      double xi;
      std::vector<double>& point_coordinates =
          visualization_continuous->GetMutablePointCoordinateVector(
              (mortar_segments + 1) * 3 * this->line_to_3D_segments_.size());
      std::vector<double>& displacement = visualization_continuous->GetMutablePointDataVector(
          "displacement", (mortar_segments + 1) * 3 * this->line_to_3D_segments_.size());
      std::vector<double>& lambda_vis = visualization_continuous->GetMutablePointDataVector(
          "lambda", (mortar_segments + 1) * 3 * this->line_to_3D_segments_.size());
      std::vector<uint8_t>& cell_type = visualization_continuous->GetMutableCellTypeVector();
      std::vector<int32_t>& cell_offset = visualization_continuous->GetMutableCellOffsetVector();

      for (const auto& segment : this->line_to_3D_segments_)
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


    // Calculate the global moment of the coupling load.
    if (visualization_nodal_forces != Teuchos::null)
    {
      // Get the global moment vector.
      auto line_load_moment_origin =
          visualization_params.get<Teuchos::RCP<LINALG::Matrix<3, 1, double>>>(
              "global_coupling_moment_origin");

      // Initialize variables for local values.
      LINALG::Matrix<3, 1, double> dr_beam_ref(true);
      LINALG::Matrix<3, 1, double> lambda_gauss_point(true);
      LINALG::Matrix<3, 1, double> r_gauss_point(true);
      LINALG::Matrix<3, 1, double> temp_moment(true);

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

          // Evaluate the coupling load at this point.
          GEOMETRYPAIR::EvaluatePosition<mortar>(
              projected_gauss_point.GetEta(), q_lambda, lambda_gauss_point);

          // Get the position at this Gauss point.
          GEOMETRYPAIR::EvaluatePosition<beam>(projected_gauss_point.GetEta(),
              FADUTILS::CastToDouble(this->ele1pos_), r_gauss_point, this->Element1());

          // Calculate moment around origin.
          temp_moment.CrossProduct(r_gauss_point, lambda_gauss_point);
          temp_moment.Scale(projected_gauss_point.GetGaussWeight() * segment_jacobian);
          (*line_load_moment_origin) += temp_moment;
        }
      }
    }
  }
}


/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_tri3>, t_hermite, t_tri3, t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_tri6>, t_hermite, t_tri6, t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_quad4>, t_hermite, t_quad4, t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_quad8>, t_hermite, t_quad8, t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_quad9>, t_hermite, t_quad9, t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_nurbs9>, t_hermite, t_nurbs9, t_line2>;

  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_tri3>, t_hermite, t_tri3, t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_tri6>, t_hermite, t_tri6, t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_quad4>, t_hermite, t_quad4, t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_quad8>, t_hermite, t_quad8, t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_quad9>, t_hermite, t_quad9, t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_nurbs9>, t_hermite, t_nurbs9, t_line3>;

  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_tri3>, t_hermite, t_tri3, t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_tri6>, t_hermite, t_tri6, t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_quad4>, t_hermite, t_quad4, t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_quad8>, t_hermite, t_quad8, t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_quad9>, t_hermite, t_quad9, t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_scalar_type<t_hermite, t_nurbs9>, t_hermite, t_nurbs9, t_line4>;


  template class BeamToSolidSurfaceMeshtyingPairMortarBase<line_to_surface_patch_scalar_type,
      t_hermite, t_tri3, t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<line_to_surface_patch_scalar_type,
      t_hermite, t_tri6, t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<line_to_surface_patch_scalar_type,
      t_hermite, t_quad4, t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<line_to_surface_patch_scalar_type,
      t_hermite, t_quad8, t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<line_to_surface_patch_scalar_type,
      t_hermite, t_quad9, t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_patch_nurbs_scalar_type<t_hermite, t_nurbs9>, t_hermite, t_nurbs9, t_line2>;

  template class BeamToSolidSurfaceMeshtyingPairMortarBase<line_to_surface_patch_scalar_type,
      t_hermite, t_tri3, t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<line_to_surface_patch_scalar_type,
      t_hermite, t_tri6, t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<line_to_surface_patch_scalar_type,
      t_hermite, t_quad4, t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<line_to_surface_patch_scalar_type,
      t_hermite, t_quad8, t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<line_to_surface_patch_scalar_type,
      t_hermite, t_quad9, t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_patch_nurbs_scalar_type<t_hermite, t_nurbs9>, t_hermite, t_nurbs9, t_line3>;

  template class BeamToSolidSurfaceMeshtyingPairMortarBase<line_to_surface_patch_scalar_type,
      t_hermite, t_tri3, t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<line_to_surface_patch_scalar_type,
      t_hermite, t_tri6, t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<line_to_surface_patch_scalar_type,
      t_hermite, t_quad4, t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<line_to_surface_patch_scalar_type,
      t_hermite, t_quad8, t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<line_to_surface_patch_scalar_type,
      t_hermite, t_quad9, t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_patch_nurbs_scalar_type<t_hermite, t_nurbs9>, t_hermite, t_nurbs9, t_line4>;
}  // namespace BEAMINTERACTION
