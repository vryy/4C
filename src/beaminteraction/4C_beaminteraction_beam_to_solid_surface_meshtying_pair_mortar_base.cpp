/*----------------------------------------------------------------------*/
/*! \file

\brief Base mortar mesh tying element for between a 3D beam and a surface element.

\level 3
*/


#include "4C_beaminteraction_beam_to_solid_surface_meshtying_pair_mortar_base.hpp"

#include "4C_beaminteraction_beam_to_solid_mortar_manager.hpp"
#include "4C_beaminteraction_beam_to_solid_surface_visualization_output_params.hpp"
#include "4C_beaminteraction_beam_to_solid_utils.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_base.hpp"
#include "4C_beaminteraction_beam_to_solid_visualization_output_writer_visualization.hpp"
#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_line_to_surface.hpp"

#include <unordered_set>

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarBase<scalar_type, beam, surface,
    mortar>::BeamToSolidSurfaceMeshtyingPairMortarBase()
    : base_class(), n_mortar_rot_(0)
{
  // Empty constructor.
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface, typename mortar>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortarBase<scalar_type, beam, surface,
    mortar>::GetPairVisualization(Teuchos::RCP<BeamToSolidVisualizationOutputWriterBase>
                                      visualization_writer,
    Teuchos::ParameterList& visualization_params) const
{
  // Get visualization of base method.
  base_class::GetPairVisualization(visualization_writer, visualization_params);

  Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_discret =
      visualization_writer->GetVisualizationWriter("btssc-mortar");
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_continuous =
      visualization_writer->GetVisualizationWriter("btssc-mortar-continuous");
  Teuchos::RCP<BEAMINTERACTION::BeamToSolidOutputWriterVisualization> visualization_nodal_forces =
      visualization_writer->GetVisualizationWriter("btssc-nodal-forces");
  if (visualization_discret.is_null() and visualization_continuous.is_null() and
      visualization_nodal_forces.is_null())
    return;

  const Teuchos::RCP<const BeamToSolidSurfaceVisualizationOutputParams>& output_params_ptr =
      visualization_params.get<Teuchos::RCP<const BeamToSolidSurfaceVisualizationOutputParams>>(
          "btssc-output_params_ptr");
  const bool write_unique_ids = output_params_ptr->GetWriteUniqueIDsFlag();

  if (visualization_discret != Teuchos::null or visualization_continuous != Teuchos::null or
      visualization_nodal_forces != Teuchos::null)
  {
    // Setup variables.
    CORE::LINALG::Matrix<3, 1, scalar_type> X;
    CORE::LINALG::Matrix<3, 1, scalar_type> r;
    CORE::LINALG::Matrix<3, 1, scalar_type> u;
    CORE::LINALG::Matrix<3, 1, double> lambda_discret;
    CORE::LINALG::Matrix<3, 1, double> xi_mortar_node;

    // Get the mortar manager and the global lambda vector, those objects will be used to get the
    // discrete Lagrange multiplier values for this pair.
    Teuchos::RCP<const BEAMINTERACTION::BeamToSolidMortarManager> mortar_manager =
        visualization_params.get<Teuchos::RCP<const BEAMINTERACTION::BeamToSolidMortarManager>>(
            "mortar_manager");
    Teuchos::RCP<Epetra_Vector> lambda =
        visualization_params.get<Teuchos::RCP<Epetra_Vector>>("lambda");

    // Get the lambda GIDs of this pair.
    auto q_lambda = GEOMETRYPAIR::InitializeElementData<mortar, double>::Initialize(nullptr);
    const auto& [lambda_row_pos, dummy] = mortar_manager->LocationVector(*this);
    std::vector<double> lambda_pair;
    CORE::FE::ExtractMyValues(*lambda, lambda_pair, lambda_row_pos);
    for (unsigned int i_dof = 0; i_dof < mortar::n_dof_; i_dof++)
      q_lambda.element_position_(i_dof) = lambda_pair[i_dof];

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
        auto& visualization_data = visualization_discret->GetVisualizationData();
        std::vector<double>& point_coordinates = visualization_data.GetPointCoordinates();
        std::vector<double>& displacement = visualization_data.GetPointData<double>("displacement");
        std::vector<double>& lambda_vis = visualization_data.GetPointData<double>("lambda");

        std::vector<double>* pair_beam_id = nullptr;
        std::vector<double>* pair_solid_id = nullptr;
        if (write_unique_ids)
        {
          pair_beam_id = &(visualization_data.GetPointData<double>("uid_0_pair_beam_id"));
          pair_solid_id = &(visualization_data.GetPointData<double>("uid_1_pair_solid_id"));
        }

        for (unsigned int i_node = 0; i_node < mortar::n_nodes_; i_node++)
        {
          // Get the local coordinate of this node.
          xi_mortar_node = CORE::FE::GetNodeCoordinates(i_node, mortar::discretization_);

          // Get position and displacement of the mortar node.
          GEOMETRYPAIR::EvaluatePosition<beam>(xi_mortar_node(0), this->ele1pos_, r);
          GEOMETRYPAIR::EvaluatePosition<beam>(xi_mortar_node(0), this->ele1posref_, X);
          u = r;
          u -= X;

          // Get the discrete Lagrangian multiplier.
          GEOMETRYPAIR::EvaluatePosition<mortar>(xi_mortar_node(0), q_lambda, lambda_discret);

          // Add to output data.
          for (unsigned int dim = 0; dim < 3; dim++)
          {
            point_coordinates.push_back(CORE::FADUTILS::CastToDouble(X(dim)));
            displacement.push_back(CORE::FADUTILS::CastToDouble(u(dim)));
            lambda_vis.push_back(CORE::FADUTILS::CastToDouble(lambda_discret(dim)));
          }

          if (write_unique_ids)
          {
            pair_beam_id->push_back(this->Element1()->Id());
            pair_solid_id->push_back(this->Element2()->Id());
          }
        }
      }
    }


    // Add the continuous values for the Lagrange multipliers.
    if (visualization_continuous != Teuchos::null and this->line_to_3D_segments_.size() > 0)
    {
      const unsigned int mortar_segments =
          visualization_params
              .get<Teuchos::RCP<const BeamToSolidSurfaceVisualizationOutputParams>>(
                  "btssc-output_params_ptr")
              ->GetMortarLambdaContinuousSegments();
      double xi;
      auto& visualization_data = visualization_continuous->GetVisualizationData();
      std::vector<double>& point_coordinates = visualization_data.GetPointCoordinates(
          (mortar_segments + 1) * 3 * this->line_to_3D_segments_.size());
      std::vector<double>& displacement = visualization_data.GetPointData<double>(
          "displacement", (mortar_segments + 1) * 3 * this->line_to_3D_segments_.size());
      std::vector<double>& lambda_vis = visualization_data.GetPointData<double>(
          "lambda", (mortar_segments + 1) * 3 * this->line_to_3D_segments_.size());
      std::vector<uint8_t>& cell_types = visualization_data.GetCellTypes();
      std::vector<int32_t>& cell_offsets = visualization_data.GetCellOffsets();

      std::vector<double>* pair_point_beam_id = nullptr;
      std::vector<double>* pair_point_solid_id = nullptr;
      std::vector<double>* pair_cell_beam_id = nullptr;
      std::vector<double>* pair_cell_solid_id = nullptr;
      if (write_unique_ids)
      {
        pair_point_beam_id = &(visualization_data.GetPointData<double>("uid_0_pair_beam_id"));
        pair_point_solid_id = &(visualization_data.GetPointData<double>("uid_1_pair_solid_id"));
        pair_cell_beam_id = &(visualization_data.GetCellData<double>("uid_0_pair_beam_id"));
        pair_cell_solid_id = &(visualization_data.GetCellData<double>("uid_1_pair_solid_id"));
      }

      for (const auto& segment : this->line_to_3D_segments_)
      {
        for (unsigned int i_curve_segment = 0; i_curve_segment <= mortar_segments;
             i_curve_segment++)
        {
          // Get the position, displacement and lambda value at the current point.
          xi = segment.GetEtaA() +
               i_curve_segment * (segment.GetEtaB() - segment.GetEtaA()) / (double)mortar_segments;
          GEOMETRYPAIR::EvaluatePosition<beam>(xi, this->ele1pos_, r);
          GEOMETRYPAIR::EvaluatePosition<beam>(xi, this->ele1posref_, X);
          u = r;
          u -= X;
          GEOMETRYPAIR::EvaluatePosition<mortar>(xi, q_lambda, lambda_discret);

          // Add to output data.
          for (unsigned int dim = 0; dim < 3; dim++)
          {
            point_coordinates.push_back(CORE::FADUTILS::CastToDouble(X(dim)));
            displacement.push_back(CORE::FADUTILS::CastToDouble(u(dim)));
            lambda_vis.push_back(CORE::FADUTILS::CastToDouble(lambda_discret(dim)));
          }
        }

        // Add the cell for this segment (poly line).
        cell_types.push_back(4);
        cell_offsets.push_back(point_coordinates.size() / 3);

        if (write_unique_ids)
        {
          pair_cell_beam_id->push_back(this->Element1()->Id());
          pair_cell_solid_id->push_back(this->Element2()->Id());
          for (unsigned int i_curve_segment = 0; i_curve_segment <= mortar_segments;
               i_curve_segment++)
          {
            pair_point_beam_id->push_back(this->Element1()->Id());
            pair_point_solid_id->push_back(this->Element2()->Id());
          }
        }
      }
    }


    // Calculate the global moment of the coupling load.
    if (visualization_nodal_forces != Teuchos::null)
    {
      // Get the global moment vector.
      auto line_load_moment_origin =
          visualization_params.get<Teuchos::RCP<CORE::LINALG::Matrix<3, 1, double>>>(
              "global_coupling_moment_origin");

      // Initialize variables for local values.
      CORE::LINALG::Matrix<3, 1, double> dr_beam_ref(true);
      CORE::LINALG::Matrix<3, 1, double> lambda_gauss_point(true);
      CORE::LINALG::Matrix<3, 1, double> r_gauss_point(true);
      CORE::LINALG::Matrix<3, 1, double> temp_moment(true);

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
        const unsigned int n_gp =
            this->line_to_3D_segments_[i_segment].GetProjectionPoints().size();
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

          // Evaluate the coupling load at this point.
          GEOMETRYPAIR::EvaluatePosition<mortar>(
              projected_gauss_point.GetEta(), q_lambda, lambda_gauss_point);

          // Get the position at this Gauss point.
          GEOMETRYPAIR::EvaluatePosition<beam>(projected_gauss_point.GetEta(),
              GEOMETRYPAIR::ElementDataToDouble<beam>::ToDouble(this->ele1pos_), r_gauss_point);

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
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite, t_nurbs9,
      t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex8>, t_hermite, t_quad4, t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex20>, t_hermite, t_quad8,
      t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex27>, t_hermite, t_quad9,
      t_line2>;

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
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite, t_nurbs9,
      t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex8>, t_hermite, t_quad4, t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex20>, t_hermite, t_quad8,
      t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex27>, t_hermite, t_quad9,
      t_line3>;

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
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite, t_nurbs9,
      t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex8>, t_hermite, t_quad4, t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex20>, t_hermite, t_quad8,
      t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortarBase<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_hex27>, t_hermite, t_quad9,
      t_line4>;
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE
