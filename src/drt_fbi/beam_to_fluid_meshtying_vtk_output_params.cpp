/*----------------------------------------------------------------------*/
/*! \file

\brief Object to store the beam to fluid meshtying output (visualization) parameters.

\level 2

\maintainer Nora Hagmeyer
*/


#include "beam_to_fluid_meshtying_vtk_output_params.H"

#include "../drt_inpar/inpar_IO_runtime_vtk_output.H"
#include "../drt_lib/drt_globalproblem.H"

FBI::BeamToFluidMeshtyingVtkOutputParams::BeamToFluidMeshtyingVtkOutputParams()
    : BEAMINTERACTION::BeamToSolidVolumeMeshtyingVtkOutputParams(), constraint_violation_(false)
{
  // empty constructor
}
/*----------------------------------------------------------------------------------------------------*/
void FBI::BeamToFluidMeshtyingVtkOutputParams::Setup()
{
  CheckInit();

  // Teuchos parameter lists from input file.
  const Teuchos::ParameterList& beam_to_fluid_meshtying_vtk_paramslist =
      DRT::Problem::Instance()
          ->FBIParams()
          .sublist("BEAM TO FLUID MESHTYING")
          .sublist("RUNTIME VTK OUTPUT");
  const Teuchos::ParameterList& global_vtk_paramslist =
      DRT::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT");

  // Get global parameters.
  output_data_format_ = DRT::INPUT::IntegralValue<INPAR::IO_RUNTIME_VTK::OutputDataFormat>(
      global_vtk_paramslist, "OUTPUT_DATA_FORMAT");
  output_interval_steps_ = global_vtk_paramslist.get<int>("INTERVAL_STEPS");
  output_every_iteration_ =
      (bool)DRT::INPUT::IntegralValue<int>(global_vtk_paramslist, "EVERY_ITERATION");

  // Get beam to fluid mesh tying specific parameters.
  output_flag_ =
      (bool)DRT::INPUT::IntegralValue<int>(beam_to_fluid_meshtying_vtk_paramslist, "WRITE_OUTPUT");

  nodal_forces_ =
      (bool)DRT::INPUT::IntegralValue<int>(beam_to_fluid_meshtying_vtk_paramslist, "NODAL_FORCES");

  segmentation_ =
      (bool)DRT::INPUT::IntegralValue<int>(beam_to_fluid_meshtying_vtk_paramslist, "SEGMENTATION");

  integration_points_ = (bool)DRT::INPUT::IntegralValue<int>(
      beam_to_fluid_meshtying_vtk_paramslist, "INTEGRATION_POINTS");

  constraint_violation_ = (bool)DRT::INPUT::IntegralValue<int>(
      beam_to_fluid_meshtying_vtk_paramslist, "CONSTRAINT_VIOLATION");

  // Set the setup flag.
  issetup_ = true;
}
