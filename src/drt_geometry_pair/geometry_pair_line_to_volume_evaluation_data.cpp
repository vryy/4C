/*----------------------------------------------------------------------*/
/*! \file

\brief container for parameters for line to volume pairs, as well as global evaluation data.

\level 1
\maintainer Ivo Steinbrecher
*/


#include "geometry_pair_line_to_volume_evaluation_data.H"

#include "../drt_inpar/inpar_beam_to_solid.H"
#include "../drt_lib/drt_globalproblem.H"


/**
 *
 */
GEOMETRYPAIR::LineToVolumeEvaluationData::LineToVolumeEvaluationData()
    : strategy_(INPAR::GEOMETRYPAIR::LineTo3DStrategy::none),
      gauss_rule_(DRT::UTILS::GaussRule1D::intrule1D_undefined),
      integration_points_circumfence_(-1),
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
    const Teuchos::ParameterList& line_to_volume_params_list =
        DRT::Problem::Instance()->BeamInteractionParams().sublist("BEAM TO SOLID VOLUME MESHTYING");

    strategy_ = Teuchos::getIntegralValue<INPAR::GEOMETRYPAIR::LineTo3DStrategy>(
        line_to_volume_params_list, "GEOMETRY_PAIR_STRATEGY");

    n_search_points_ = line_to_volume_params_list.get<int>("GEOMETRY_PAIR_SEARCH_POINTS");

    gauss_rule_ =
        INPAR::BEAMTOSOLID::IntToGaussRule1D(line_to_volume_params_list.get<int>("GAUSS_POINTS"));

    integration_points_circumfence_ =
        line_to_volume_params_list.get<int>("INTEGRATION_POINTS_CIRCUMFENCE");
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
