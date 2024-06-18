/*----------------------------------------------------------------------*/
/*! \file

\brief Object to store the beam to solid surface output (visualization) parameters.

\level 3

*/


#include "4C_beaminteraction_beam_to_solid_surface_visualization_output_params.hpp"

#include "4C_global_data.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputParams::
    BeamToSolidSurfaceVisualizationOutputParams()
    : isinit_(false),
      issetup_(false),
      output_interval_steps_(-1),
      output_every_iteration_(false),
      output_flag_(false),
      nodal_forces_(false),
      averaged_normals_(false),
      mortar_lambda_discret_(false),
      mortar_lambda_continuous_(false),
      mortar_lambda_continuous_segments_(0),
      segmentation_(false),
      integration_points_(false),
      write_unique_ids_(false)
{
  // empty constructor
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputParams::init()
{
  issetup_ = false;
  isinit_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputParams::setup()
{
  check_init();

  // Teuchos parameter lists from input file.
  const Teuchos::ParameterList& beam_to_solid_volume_meshtying_visualization_output_paramslist =
      Global::Problem::Instance()
          ->beam_interaction_params()
          .sublist("BEAM TO SOLID SURFACE")
          .sublist("RUNTIME VTK OUTPUT");
  const Teuchos::ParameterList& global_visualization_output_paramslist =
      Global::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT");

  // Get global parameters.
  output_interval_steps_ = global_visualization_output_paramslist.get<int>("INTERVAL_STEPS");
  output_every_iteration_ = (bool)Core::UTILS::IntegralValue<int>(
      global_visualization_output_paramslist, "EVERY_ITERATION");

  // Get beam to solid surface specific parameters.
  output_flag_ = (bool)Core::UTILS::IntegralValue<int>(
      beam_to_solid_volume_meshtying_visualization_output_paramslist, "WRITE_OUTPUT");

  nodal_forces_ = (bool)Core::UTILS::IntegralValue<int>(
      beam_to_solid_volume_meshtying_visualization_output_paramslist, "NODAL_FORCES");

  averaged_normals_ = (bool)Core::UTILS::IntegralValue<int>(
      beam_to_solid_volume_meshtying_visualization_output_paramslist, "AVERAGED_NORMALS");

  mortar_lambda_discret_ = (bool)Core::UTILS::IntegralValue<int>(
      beam_to_solid_volume_meshtying_visualization_output_paramslist, "MORTAR_LAMBDA_DISCRET");

  mortar_lambda_continuous_ = (bool)Core::UTILS::IntegralValue<int>(
      beam_to_solid_volume_meshtying_visualization_output_paramslist, "MORTAR_LAMBDA_CONTINUOUS");

  mortar_lambda_continuous_segments_ =
      beam_to_solid_volume_meshtying_visualization_output_paramslist.get<int>(
          "MORTAR_LAMBDA_CONTINUOUS_SEGMENTS");

  segmentation_ = (bool)Core::UTILS::IntegralValue<int>(
      beam_to_solid_volume_meshtying_visualization_output_paramslist, "SEGMENTATION");

  integration_points_ = (bool)Core::UTILS::IntegralValue<int>(
      beam_to_solid_volume_meshtying_visualization_output_paramslist, "INTEGRATION_POINTS");

  write_unique_ids_ = (bool)Core::UTILS::IntegralValue<int>(
      beam_to_solid_volume_meshtying_visualization_output_paramslist, "UNIQUE_IDS");

  // Set the setup flag.
  issetup_ = true;
}

FOUR_C_NAMESPACE_CLOSE
