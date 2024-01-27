/*----------------------------------------------------------------------*/
/*! \file

\brief Object to store the beam to solid surface output (visualization) parameters.

\level 3

*/


#include "baci_beaminteraction_beam_to_solid_surface_visualization_output_params.H"

#include "baci_global_data.H"
#include "baci_inpar_parameterlist_utils.H"

BACI_NAMESPACE_OPEN


/**
 *
 */
BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputParams::
    BeamToSolidSurfaceVisualizationOutputParams()
    : isinit_(false),
      issetup_(false),
      visualization_parameters_(IO::VisualizationParametersFactory(
          GLOBAL::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"))),
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
void BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputParams::Init()
{
  issetup_ = false;
  isinit_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputParams::Setup()
{
  CheckInit();

  // Teuchos parameter lists from input file.
  const Teuchos::ParameterList& beam_to_solid_volume_meshtying_visualization_output_paramslist =
      GLOBAL::Problem::Instance()
          ->BeamInteractionParams()
          .sublist("BEAM TO SOLID SURFACE")
          .sublist("RUNTIME VTK OUTPUT");
  const Teuchos::ParameterList& global_visualization_output_paramslist =
      GLOBAL::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT");

  // Get global parameters.
  output_interval_steps_ = global_visualization_output_paramslist.get<int>("INTERVAL_STEPS");
  output_every_iteration_ =
      (bool)INPUT::IntegralValue<int>(global_visualization_output_paramslist, "EVERY_ITERATION");
  visualization_parameters_.every_iteration_ = output_every_iteration_;

  // Get beam to solid surface specific parameters.
  output_flag_ = (bool)INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_visualization_output_paramslist, "WRITE_OUTPUT");

  nodal_forces_ = (bool)INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_visualization_output_paramslist, "NODAL_FORCES");

  averaged_normals_ = (bool)INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_visualization_output_paramslist, "AVERAGED_NORMALS");

  mortar_lambda_discret_ = (bool)INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_visualization_output_paramslist, "MORTAR_LAMBDA_DISCRET");

  mortar_lambda_continuous_ = (bool)INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_visualization_output_paramslist, "MORTAR_LAMBDA_CONTINUOUS");

  mortar_lambda_continuous_segments_ =
      beam_to_solid_volume_meshtying_visualization_output_paramslist.get<int>(
          "MORTAR_LAMBDA_CONTINUOUS_SEGMENTS");

  segmentation_ = (bool)INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_visualization_output_paramslist, "SEGMENTATION");

  integration_points_ = (bool)INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_visualization_output_paramslist, "INTEGRATION_POINTS");

  write_unique_ids_ = (bool)INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_visualization_output_paramslist, "UNIQUE_IDS");

  // Set the setup flag.
  issetup_ = true;
}

BACI_NAMESPACE_CLOSE
