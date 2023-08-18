/*----------------------------------------------------------------------*/
/*! \file

\brief Object to store the beam to solid volume meshtying output (visualization) parameters.

\level 3

*/


#include "baci_beaminteraction_beam_to_solid_volume_meshtying_vtk_output_params.H"

#include "baci_lib_globalproblem.H"


/**
 *
 */
BEAMINTERACTION::BeamToSolidVolumeMeshtyingVtkOutputParams::
    BeamToSolidVolumeMeshtyingVtkOutputParams()
    : isinit_(false),
      issetup_(false),
      visualization_parameters_(IO::VisualizationParametersFactory(
          DRT::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"))),
      output_interval_steps_(-1),
      output_every_iteration_(false),
      output_flag_(false),
      nodal_forces_(false),
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
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingVtkOutputParams::Init()
{
  issetup_ = false;
  isinit_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingVtkOutputParams::Setup()
{
  CheckInit();

  // Teuchos parameter lists from input file.
  const Teuchos::ParameterList& beam_to_solid_volume_meshtying_vtk_paramslist =
      DRT::Problem::Instance()
          ->BeamInteractionParams()
          .sublist("BEAM TO SOLID VOLUME MESHTYING")
          .sublist("RUNTIME VTK OUTPUT");
  const Teuchos::ParameterList& global_vtk_paramslist =
      DRT::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT");

  // Get global parameters.
  output_interval_steps_ = global_vtk_paramslist.get<int>("INTERVAL_STEPS");
  output_every_iteration_ =
      (bool)DRT::INPUT::IntegralValue<int>(global_vtk_paramslist, "EVERY_ITERATION");
  visualization_parameters_.every_iteration_ = output_every_iteration_;

  // Get beam to solid volume mesh tying specific parameters.
  output_flag_ = (bool)DRT::INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_vtk_paramslist, "WRITE_OUTPUT");

  nodal_forces_ = (bool)DRT::INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_vtk_paramslist, "NODAL_FORCES");

  mortar_lambda_discret_ = (bool)DRT::INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_vtk_paramslist, "MORTAR_LAMBDA_DISCRET");

  mortar_lambda_continuous_ = (bool)DRT::INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_vtk_paramslist, "MORTAR_LAMBDA_CONTINUOUS");

  mortar_lambda_continuous_segments_ =
      beam_to_solid_volume_meshtying_vtk_paramslist.get<int>("MORTAR_LAMBDA_CONTINUOUS_SEGMENTS");

  segmentation_ = (bool)DRT::INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_vtk_paramslist, "SEGMENTATION");

  integration_points_ = (bool)DRT::INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_vtk_paramslist, "INTEGRATION_POINTS");

  write_unique_ids_ = (bool)DRT::INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_vtk_paramslist, "UNIQUE_IDS");

  // Set the setup flag.
  issetup_ = true;
}
