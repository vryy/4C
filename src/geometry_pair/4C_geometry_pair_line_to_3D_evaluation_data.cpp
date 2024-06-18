/*----------------------------------------------------------------------*/
/*! \file

\brief Container for parameters for line to 3D pairs, as well as global evaluation data.

\level 1
*/


#include "4C_geometry_pair_line_to_3D_evaluation_data.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
GEOMETRYPAIR::LineTo3DEvaluationData::LineTo3DEvaluationData(
    const Teuchos::ParameterList& input_parameter_list)
    : GeometryEvaluationDataBase(input_parameter_list),
      strategy_(Inpar::GEOMETRYPAIR::LineTo3DStrategy::none),
      gauss_rule_(Core::FE::GaussRule1D::undefined),
      integration_points_circumference_(-1),
      gauss_point_projection_tracker_(),
      n_search_points_(0),
      not_all_gauss_points_project_valid_action_(
          Inpar::GEOMETRYPAIR::NotAllGaussPointsProjectValidAction::fail),
      segment_tracker_()
{
  // Get parameters from the input file.
  {
    strategy_ = Teuchos::getIntegralValue<Inpar::GEOMETRYPAIR::LineTo3DStrategy>(
        input_parameter_list, "GEOMETRY_PAIR_STRATEGY");

    n_search_points_ = input_parameter_list.get<int>("GEOMETRY_PAIR_SEGMENTATION_SEARCH_POINTS");
    not_all_gauss_points_project_valid_action_ =
        Teuchos::getIntegralValue<Inpar::GEOMETRYPAIR::NotAllGaussPointsProjectValidAction>(
            input_parameter_list,
            "GEOMETRY_PAIR_SEGMENTATION_NOT_ALL_GAUSS_POINTS_PROJECT_VALID_ACTION");

    gauss_rule_ =
        Inpar::GEOMETRYPAIR::IntToGaussRule1D(input_parameter_list.get<int>("GAUSS_POINTS"));

    integration_points_circumference_ =
        input_parameter_list.get<int>("INTEGRATION_POINTS_CIRCUMFERENCE");
  }

  // Initialize evaluation data structures.
  clear();
}

/**
 *
 */
void GEOMETRYPAIR::LineTo3DEvaluationData::clear()
{
  // Call reset on the base method.
  GeometryEvaluationDataBase::clear();

  // Initialize evaluation data structures.
  {
    // Tracker for gauss point projection method.
    gauss_point_projection_tracker_.clear();

    // Segment tracker for segmentation.
    segment_tracker_.clear();
  }
}

/**
 *
 */
void GEOMETRYPAIR::LineTo3DEvaluationData::ResetTracker()
{
  for (auto& data : gauss_point_projection_tracker_)
    std::fill(data.second.begin(), data.second.end(), false);

  for (auto& data : segment_tracker_) data.second.clear();
}

FOUR_C_NAMESPACE_CLOSE
