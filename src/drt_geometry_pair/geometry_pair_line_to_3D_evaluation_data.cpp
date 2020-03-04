/*----------------------------------------------------------------------*/
/*! \file

\brief Container for parameters for line to 3D pairs, as well as global evaluation data.

\level 1
\maintainer Ivo Steinbrecher
*/


#include "geometry_pair_line_to_3D_evaluation_data.H"

#include "../drt_inpar/inpar_beam_to_solid.H"


/**
 *
 */
GEOMETRYPAIR::LineTo3DEvaluationData::LineTo3DEvaluationData(
    const Teuchos::ParameterList& input_parameter_list)
    : GeometryEvaluationDataBase(input_parameter_list),
      strategy_(INPAR::GEOMETRYPAIR::LineTo3DStrategy::none),
      gauss_rule_(DRT::UTILS::GaussRule1D::intrule1D_undefined),
      integration_points_circumference_(-1),
      gauss_point_projection_tracker_(),
      n_search_points_(0),
      segment_tracker_()
{
  // Get parameters from the input file.
  {
    strategy_ = Teuchos::getIntegralValue<INPAR::GEOMETRYPAIR::LineTo3DStrategy>(
        input_parameter_list, "GEOMETRY_PAIR_STRATEGY");

    n_search_points_ = input_parameter_list.get<int>("GEOMETRY_PAIR_SEARCH_POINTS");

    gauss_rule_ =
        INPAR::GEOMETRYPAIR::IntToGaussRule1D(input_parameter_list.get<int>("GAUSS_POINTS"));

    integration_points_circumference_ =
        input_parameter_list.get<int>("INTEGRATION_POINTS_CIRCUMFERENCE");
  }

  // Initialize evaluation data structures.
  Clear();
}

/**
 *
 */
void GEOMETRYPAIR::LineTo3DEvaluationData::Clear()
{
  // Call reset on the base method.
  GeometryEvaluationDataBase::Clear();

  // Initialize evaluation data structures.
  {
    // Tracker for gauss point projection method.
    gauss_point_projection_tracker_.clear();

    // Segment tracker for segmentation.
    segment_tracker_.clear();
  }
}
