/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container holding all input parameters for vtk-based visualization of beam contact

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam_contact_runtime_vtk_output_params.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BEAMINTERACTION::BeamContactRuntimeVtkOutputParams::BeamContactRuntimeVtkOutputParams()
    : isinit_(false),
      issetup_(false),
      output_data_format_(INPAR::BEAMCONTACT::vague),
      output_interval_steps_(-1),
      output_every_iteration_(false),
      output_forces_(false),
      output_gaps_(false)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactRuntimeVtkOutputParams::Init()
{
  issetup_ = false;
  // empty for now

  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactRuntimeVtkOutputParams::Setup()
{
  ThrowErrorIfNotInit();

  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_contact_vtk_paramslist =
      DRT::Problem::Instance()->BeamContactParams().sublist("RUNTIME VTK OUTPUT");

  /****************************************************************************/
  // get and check required parameters
  /****************************************************************************/
  output_data_format_ = DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::OutputDataFormat>(
      beam_contact_vtk_paramslist, "OUTPUT_DATA_FORMAT");

  output_interval_steps_ = beam_contact_vtk_paramslist.get<int>("INTERVAL_STEPS");

  output_every_iteration_ =
      (bool)DRT::INPUT::IntegralValue<int>(beam_contact_vtk_paramslist, "EVERY_ITERATION");

  /****************************************************************************/
  output_forces_ =
      (bool)DRT::INPUT::IntegralValue<int>(beam_contact_vtk_paramslist, "CONTACT_FORCES");

  /****************************************************************************/
  output_gaps_ = (bool)DRT::INPUT::IntegralValue<int>(beam_contact_vtk_paramslist, "GAPS");


  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactRuntimeVtkOutputParams::ThrowErrorIfNotInitAndSetup() const
{
  if (!IsInit() or !IsSetup()) dserror("Call Init() and Setup() first!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactRuntimeVtkOutputParams::ThrowErrorIfNotInit() const
{
  if (!IsInit()) dserror("Init() has not been called, yet!");
}
