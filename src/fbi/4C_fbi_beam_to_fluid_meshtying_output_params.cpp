/*----------------------------------------------------------------------*/
/*! \file

\brief Object to store the beam to fluid meshtying output (visualization) parameters.

\level 2

*/


#include "4C_fbi_beam_to_fluid_meshtying_output_params.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_IO_runtime_output.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

FBI::BeamToFluidMeshtyingVtkOutputParams::BeamToFluidMeshtyingVtkOutputParams()
    : BEAMINTERACTION::BeamToSolidVolumeMeshtyingVisualizationOutputParams(),
      constraint_violation_(false)
{
  // empty constructor
}
/*----------------------------------------------------------------------------------------------------*/
void FBI::BeamToFluidMeshtyingVtkOutputParams::Setup()
{
  CheckInit();

  // Teuchos parameter lists from input file.
  const Teuchos::ParameterList& beam_to_fluid_meshtying_visualization_output_paramslist =
      GLOBAL::Problem::Instance()
          ->FBIParams()
          .sublist("BEAM TO FLUID MESHTYING")
          .sublist("RUNTIME VTK OUTPUT");
  const Teuchos::ParameterList& global_visualization_output_paramslist =
      GLOBAL::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT");

  // Get global parameters.
  output_interval_steps_ = global_visualization_output_paramslist.get<int>("INTERVAL_STEPS");
  output_every_iteration_ = (bool)CORE::UTILS::IntegralValue<int>(
      global_visualization_output_paramslist, "EVERY_ITERATION");

  // Get beam to fluid mesh tying specific parameters.
  output_flag_ = (bool)CORE::UTILS::IntegralValue<int>(
      beam_to_fluid_meshtying_visualization_output_paramslist, "WRITE_OUTPUT");

  nodal_forces_ = (bool)CORE::UTILS::IntegralValue<int>(
      beam_to_fluid_meshtying_visualization_output_paramslist, "NODAL_FORCES");

  segmentation_ = (bool)CORE::UTILS::IntegralValue<int>(
      beam_to_fluid_meshtying_visualization_output_paramslist, "SEGMENTATION");

  integration_points_ = (bool)CORE::UTILS::IntegralValue<int>(
      beam_to_fluid_meshtying_visualization_output_paramslist, "INTEGRATION_POINTS");

  constraint_violation_ = (bool)CORE::UTILS::IntegralValue<int>(
      beam_to_fluid_meshtying_visualization_output_paramslist, "CONSTRAINT_VIOLATION");

  // Set the setup flag.
  issetup_ = true;
}

FOUR_C_NAMESPACE_CLOSE
