/*----------------------------------------------------------------------*/
/*! \file

\brief Object to store the beam to solid surface output (visualization) parameters.

\level 3

*/


#include "beam_to_solid_surface_vtk_output_params.H"

#include "../drt_lib/drt_globalproblem.H"


/**
 *
 */
BEAMINTERACTION::BeamToSolidSurfaceVtkOutputParams::BeamToSolidSurfaceVtkOutputParams()
    : isinit_(false),
      issetup_(false),
      output_data_format_(INPAR::IO_RUNTIME_VTK::vague),
      output_interval_steps_(-1),
      output_every_iteration_(false),
      output_flag_(false),
      nodal_forces_(false),
      averaged_normals_(false),
      mortar_lambda_discret_(false),
      mortar_lambda_continuous_(false),
      mortar_lambda_continuous_segments_(0),
      segmentation_(false),
      integration_points_(false)
{
  // empty constructor
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceVtkOutputParams::Init()
{
  issetup_ = false;
  isinit_ = true;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceVtkOutputParams::Setup()
{
  CheckInit();

  // Teuchos parameter lists from input file.
  const Teuchos::ParameterList& beam_to_solid_volume_meshtying_vtk_paramslist =
      DRT::Problem::Instance()
          ->BeamInteractionParams()
          .sublist("BEAM TO SOLID SURFACE")
          .sublist("RUNTIME VTK OUTPUT");
  const Teuchos::ParameterList& global_vtk_paramslist =
      DRT::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT");

  // Get global parameters.
  output_data_format_ = DRT::INPUT::IntegralValue<INPAR::IO_RUNTIME_VTK::OutputDataFormat>(
      global_vtk_paramslist, "OUTPUT_DATA_FORMAT");
  output_interval_steps_ = global_vtk_paramslist.get<int>("INTERVAL_STEPS");
  output_every_iteration_ =
      (bool)DRT::INPUT::IntegralValue<int>(global_vtk_paramslist, "EVERY_ITERATION");

  // Get beam to solid surface specific parameters.
  output_flag_ = (bool)DRT::INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_vtk_paramslist, "WRITE_OUTPUT");

  nodal_forces_ = (bool)DRT::INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_vtk_paramslist, "NODAL_FORCES");

  averaged_normals_ = (bool)DRT::INPUT::IntegralValue<int>(
      beam_to_solid_volume_meshtying_vtk_paramslist, "AVERAGED_NORMALS");

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

  // Set the setup flag.
  issetup_ = true;
}
