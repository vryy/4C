#include "4C_beaminteraction_contact_runtime_visualization_output_params.hpp"

#include "4C_global_data.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BEAMINTERACTION::BeamContactRuntimeVisualizationOutputParams::
    BeamContactRuntimeVisualizationOutputParams(const double restart_time)
    : isinit_(false),
      issetup_(false),
      visualization_parameters_(Core::IO::visualization_parameters_factory(
          Global::Problem::instance()->io_params().sublist("RUNTIME VTK OUTPUT"),
          *Global::Problem::instance()->output_control_file(), restart_time)),
      output_interval_steps_(-1),
      output_every_iteration_(false),
      output_forces_(false),
      output_gaps_(false)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactRuntimeVisualizationOutputParams::init()
{
  issetup_ = false;
  // empty for now

  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactRuntimeVisualizationOutputParams::setup()
{
  throw_error_if_not_init();

  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_contact_visualization_output_paramslist =
      Global::Problem::instance()->beam_contact_params().sublist("RUNTIME VTK OUTPUT");

  /****************************************************************************/
  // get and check required parameters
  /****************************************************************************/
  output_interval_steps_ = beam_contact_visualization_output_paramslist.get<int>("INTERVAL_STEPS");

  output_every_iteration_ =
      beam_contact_visualization_output_paramslist.get<bool>("EVERY_ITERATION");
  visualization_parameters_.every_iteration_ = output_every_iteration_;

  /****************************************************************************/
  output_forces_ = beam_contact_visualization_output_paramslist.get<bool>("CONTACT_FORCES");

  /****************************************************************************/
  output_gaps_ = beam_contact_visualization_output_paramslist.get<bool>("GAPS");


  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactRuntimeVisualizationOutputParams::
    throw_error_if_not_init_and_setup() const
{
  if (!is_init() or !is_setup()) FOUR_C_THROW("Call init() and setup() first!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactRuntimeVisualizationOutputParams::throw_error_if_not_init() const
{
  if (!is_init()) FOUR_C_THROW("init() has not been called, yet!");
}

FOUR_C_NAMESPACE_CLOSE
