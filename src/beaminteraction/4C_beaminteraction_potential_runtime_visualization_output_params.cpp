/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container for input parameters for vtk-based visualization of potential-based beam
       interactions

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#include "4C_beaminteraction_potential_runtime_visualization_output_params.hpp"

#include "4C_global_data.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BEAMINTERACTION::BeamToBeamPotentialRuntimeOutputParams::BeamToBeamPotentialRuntimeOutputParams(
    const double restart_time)
    : isinit_(false),
      issetup_(false),
      visualization_parameters_(Core::IO::VisualizationParametersFactory(
          Global::Problem::instance()->io_params().sublist("RUNTIME VTK OUTPUT"),
          *Global::Problem::instance()->output_control_file(), restart_time)),
      output_interval_steps_(-1),
      output_every_iteration_(false),
      output_forces_(false),
      output_moments_(false),
      write_force_moment_per_elepair_(false)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamPotentialRuntimeOutputParams::init(
    const Teuchos::ParameterList& beam_contact_visualization_output_paramslist)
{
  issetup_ = false;


  /****************************************************************************/
  // get and check required parameters
  /****************************************************************************/
  output_interval_steps_ = beam_contact_visualization_output_paramslist.get<int>("INTERVAL_STEPS");

  output_every_iteration_ = (bool)Core::UTILS::IntegralValue<int>(
      beam_contact_visualization_output_paramslist, "EVERY_ITERATION");
  visualization_parameters_.every_iteration_ = output_every_iteration_;

  /****************************************************************************/
  output_forces_ =
      (bool)Core::UTILS::IntegralValue<int>(beam_contact_visualization_output_paramslist, "FORCES");

  /****************************************************************************/
  output_moments_ = (bool)Core::UTILS::IntegralValue<int>(
      beam_contact_visualization_output_paramslist, "MOMENTS");

  /****************************************************************************/
  write_force_moment_per_elepair_ = (bool)Core::UTILS::IntegralValue<int>(
      beam_contact_visualization_output_paramslist, "WRITE_FORCE_MOMENT_PER_ELEMENTPAIR");

  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamPotentialRuntimeOutputParams::setup()
{
  throw_error_if_not_init();

  // empty for now

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamPotentialRuntimeOutputParams::throw_error_if_not_init_and_setup()
    const
{
  if (!is_init() or !is_setup()) FOUR_C_THROW("Call init() and setup() first!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamPotentialRuntimeOutputParams::throw_error_if_not_init() const
{
  if (!is_init()) FOUR_C_THROW("init() has not been called, yet!");
}

FOUR_C_NAMESPACE_CLOSE
