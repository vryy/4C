/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container holding all input parameters for vtk-based visualization of beam contact

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#include "baci_beaminteraction_contact_runtime_visualization_output_params.hpp"

#include "baci_global_data.hpp"
#include "baci_utils_exceptions.hpp"

BACI_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BEAMINTERACTION::BeamContactRuntimeVisualizationOutputParams::
    BeamContactRuntimeVisualizationOutputParams(const double restart_time)
    : isinit_(false),
      issetup_(false),
      visualization_parameters_(IO::VisualizationParametersFactory(
          GLOBAL::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"),
          *GLOBAL::Problem::Instance()->OutputControlFile(), restart_time)),
      output_interval_steps_(-1),
      output_every_iteration_(false),
      output_forces_(false),
      output_gaps_(false)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactRuntimeVisualizationOutputParams::Init()
{
  issetup_ = false;
  // empty for now

  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactRuntimeVisualizationOutputParams::Setup()
{
  ThrowErrorIfNotInit();

  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_contact_visualization_output_paramslist =
      GLOBAL::Problem::Instance()->BeamContactParams().sublist("RUNTIME VTK OUTPUT");

  /****************************************************************************/
  // get and check required parameters
  /****************************************************************************/
  output_interval_steps_ = beam_contact_visualization_output_paramslist.get<int>("INTERVAL_STEPS");

  output_every_iteration_ = (bool)CORE::UTILS::IntegralValue<int>(
      beam_contact_visualization_output_paramslist, "EVERY_ITERATION");
  visualization_parameters_.every_iteration_ = output_every_iteration_;

  /****************************************************************************/
  output_forces_ = (bool)CORE::UTILS::IntegralValue<int>(
      beam_contact_visualization_output_paramslist, "CONTACT_FORCES");

  /****************************************************************************/
  output_gaps_ =
      (bool)CORE::UTILS::IntegralValue<int>(beam_contact_visualization_output_paramslist, "GAPS");


  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactRuntimeVisualizationOutputParams::ThrowErrorIfNotInitAndSetup()
    const
{
  if (!IsInit() or !IsSetup()) dserror("Call Init() and Setup() first!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactRuntimeVisualizationOutputParams::ThrowErrorIfNotInit() const
{
  if (!IsInit()) dserror("Init() has not been called, yet!");
}

BACI_NAMESPACE_CLOSE
