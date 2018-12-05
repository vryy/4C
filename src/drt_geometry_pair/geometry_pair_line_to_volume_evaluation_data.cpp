/*!
\file geometry_pair_line_to_volume_evaluation_data.cpp

\brief container for parameters for line to volume pairs, as well as global evaluation data.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "geometry_pair_line_to_volume_evaluation_data.H"

#include "../drt_lib/drt_globalproblem.H"


/**
 *
 */
GEOMETRYPAIR::LineToVolumeEvaluationData::LineToVolumeEvaluationData()
    : strategy_(INPAR::GEOMETRYPAIR::LineToVolumeStrategy::none),
      gauss_rule_(DRT::UTILS::GaussRule1D::intrule1D_undefined),
      gauss_point_projection_tracker_(),
      n_search_points_(0),
      segment_tracker_()
{
  // Empty Constructor
}


/**
 *
 */
void GEOMETRYPAIR::LineToVolumeEvaluationData::Init()
{
  // Call init in base class first.
  GEOMETRYPAIR::GeometryEvaluationDataBase::Init();

  // Get parameters from the input file.
  {
    // Teuchos parameter list for line to volume pairs.
    const Teuchos::ParameterList& line_to_volume_params_list =
        DRT::Problem::Instance()->GeometryPairParams().sublist("LINE TO VOLUME");

    // Get the segmentation strategy.
    strategy_ = Teuchos::getIntegralValue<INPAR::GEOMETRYPAIR::LineToVolumeStrategy>(
        line_to_volume_params_list, "STRATEGY");

    // Get number of search points.
    n_search_points_ = line_to_volume_params_list.get<int>("SEARCH_POINTS");
  }

  // Initialize evaluation data structures.
  {
    // Tracker for gauss point projection method.
    gauss_point_projection_tracker_.clear();

    // Segment tracker for segmentation.
    segment_tracker_.clear();
  }

  // Set init flag.
  isinit_ = true;
}


/**
 *
 */
void GEOMETRYPAIR::LineToVolumeEvaluationData::Setup()
{
  // Call  setup in base class.
  GEOMETRYPAIR::GeometryEvaluationDataBase::Setup();

  issetup_ = true;
}
