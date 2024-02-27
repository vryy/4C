/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container for input parameters for vtk-based visualization of potential-based beam
       interactions

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#include "baci_beaminteraction_potential_runtime_visualization_output_params.hpp"

#include "baci_global_data.hpp"
#include "baci_utils_exceptions.hpp"

BACI_NAMESPACE_OPEN

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BEAMINTERACTION::BeamToBeamPotentialRuntimeOutputParams::BeamToBeamPotentialRuntimeOutputParams(
    const double restart_time)
    : isinit_(false),
      issetup_(false),
      visualization_parameters_(IO::VisualizationParametersFactory(
          GLOBAL::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"),
          *GLOBAL::Problem::Instance()->OutputControlFile(), restart_time)),
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
void BEAMINTERACTION::BeamToBeamPotentialRuntimeOutputParams::Init(
    const Teuchos::ParameterList& beam_contact_visualization_output_paramslist)
{
  issetup_ = false;


  /****************************************************************************/
  // get and check required parameters
  /****************************************************************************/
  output_interval_steps_ = beam_contact_visualization_output_paramslist.get<int>("INTERVAL_STEPS");

  output_every_iteration_ = (bool)CORE::UTILS::IntegralValue<int>(
      beam_contact_visualization_output_paramslist, "EVERY_ITERATION");
  visualization_parameters_.every_iteration_ = output_every_iteration_;

  /****************************************************************************/
  output_forces_ =
      (bool)CORE::UTILS::IntegralValue<int>(beam_contact_visualization_output_paramslist, "FORCES");

  /****************************************************************************/
  output_moments_ = (bool)CORE::UTILS::IntegralValue<int>(
      beam_contact_visualization_output_paramslist, "MOMENTS");

  /****************************************************************************/
  write_force_moment_per_elepair_ = (bool)CORE::UTILS::IntegralValue<int>(
      beam_contact_visualization_output_paramslist, "WRITE_FORCE_MOMENT_PER_ELEMENTPAIR");

  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamPotentialRuntimeOutputParams::Setup()
{
  ThrowErrorIfNotInit();

  // empty for now

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamPotentialRuntimeOutputParams::ThrowErrorIfNotInitAndSetup() const
{
  if (!IsInit() or !IsSetup()) dserror("Call Init() and Setup() first!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamPotentialRuntimeOutputParams::ThrowErrorIfNotInit() const
{
  if (!IsInit()) dserror("Init() has not been called, yet!");
}

BACI_NAMESPACE_CLOSE
