/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_contact_runtime_vtk_output_params.cpp

\brief data container holding all input parameters for vtk-based visualization of beam contact

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam_contact_runtime_vtk_output_params.H"

#include "../drt_lib/drt_dserror.H"

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BEAMINTERACTION::BeamToBeamContactRuntimeVtkParams::BeamToBeamContactRuntimeVtkParams()
: isinit_(false),
  issetup_(false),
  output_data_format_( INPAR::BEAMCONTACT::vague ),
  output_interval_steps_(-1),
  output_every_iteration_(false),
  output_forces_(false),
  output_gaps_(false)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamContactRuntimeVtkParams::Init(
    const Teuchos::ParameterList& beam_contact_vtk_paramslist )
{
  issetup_ = false;


  /****************************************************************************/
  // get and check required parameters
  /****************************************************************************/
  output_data_format_ =
      DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::OutputDataFormat>(
          beam_contact_vtk_paramslist, "OUTPUT_DATA_FORMAT" );

  output_interval_steps_ = beam_contact_vtk_paramslist.get<int>("INTERVAL_STEPS");

  output_every_iteration_ =
      (bool) DRT::INPUT::IntegralValue<int>(beam_contact_vtk_paramslist, "EVERY_ITERATION");

  if ( output_every_iteration_ )
    dserror("not implemented yet!");


  /****************************************************************************/
  output_forces_ =
      (bool) DRT::INPUT::IntegralValue<int>(beam_contact_vtk_paramslist, "CONTACT_FORCES");

  /****************************************************************************/
  output_gaps_ =
      (bool) DRT::INPUT::IntegralValue<int>(beam_contact_vtk_paramslist, "GAPS");


  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamContactRuntimeVtkParams::Setup()
{
  ThrowErrorIfNotInit();

  // empty for now

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamContactRuntimeVtkParams::ThrowErrorIfNotInitAndSetup() const
{
  if (!IsInit() or !IsSetup())
    dserror("Call Init() and Setup() first!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamContactRuntimeVtkParams::ThrowErrorIfNotInit() const
{
  if (!IsInit())
    dserror("Init() has not been called, yet!");
}
