/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container for input parameters for vtk-based visualization of potential-based beam
       interactions

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam_potential_runtime_vtk_output_params.H"

#include "../drt_lib/drt_dserror.H"

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BEAMINTERACTION::BeamToBeamPotentialRuntimeVtkParams::BeamToBeamPotentialRuntimeVtkParams()
    : isinit_(false),
      issetup_(false),
      output_data_format_(INPAR::BEAMPOTENTIAL::vague),
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
void BEAMINTERACTION::BeamToBeamPotentialRuntimeVtkParams::Init(
    const Teuchos::ParameterList& beam_contact_vtk_paramslist)
{
  issetup_ = false;


  /****************************************************************************/
  // get and check required parameters
  /****************************************************************************/
  output_data_format_ = DRT::INPUT::IntegralValue<INPAR::BEAMPOTENTIAL::OutputDataFormat>(
      beam_contact_vtk_paramslist, "OUTPUT_DATA_FORMAT");

  output_interval_steps_ = beam_contact_vtk_paramslist.get<int>("INTERVAL_STEPS");

  output_every_iteration_ =
      (bool)DRT::INPUT::IntegralValue<int>(beam_contact_vtk_paramslist, "EVERY_ITERATION");

  /****************************************************************************/
  output_forces_ = (bool)DRT::INPUT::IntegralValue<int>(beam_contact_vtk_paramslist, "FORCES");

  /****************************************************************************/
  output_moments_ = (bool)DRT::INPUT::IntegralValue<int>(beam_contact_vtk_paramslist, "MOMENTS");

  /****************************************************************************/
  write_force_moment_per_elepair_ = (bool)DRT::INPUT::IntegralValue<int>(
      beam_contact_vtk_paramslist, "WRITE_FORCE_MOMENT_PER_ELEMENTPAIR");

  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamPotentialRuntimeVtkParams::Setup()
{
  ThrowErrorIfNotInit();

  // empty for now

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamPotentialRuntimeVtkParams::ThrowErrorIfNotInitAndSetup() const
{
  if (!IsInit() or !IsSetup()) dserror("Call Init() and Setup() first!");
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToBeamPotentialRuntimeVtkParams::ThrowErrorIfNotInit() const
{
  if (!IsInit()) dserror("Init() has not been called, yet!");
}
