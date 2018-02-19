/*----------------------------------------------------------------------------*/
/*!
\file beam_to_solid_volume_meshtying_params.cpp

\brief data container holding all beam to solid volume meshtying input parameters

\level 3

\maintainer Alexander Popp
*/
/*----------------------------------------------------------------------------*/


#include "beam_to_solid_volume_meshtying_params.H"

#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams::BeamToSolidVolumeMeshtyingParams():
isinit_(false),
issetup_(false),
BTSVOLMT_penalty_param_(-1.0)
{
  // Empty Constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams::Init()
{
  issetup_ = false;

  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_to_solid_contact_params_list =
      DRT::Problem::Instance()->BeamInteractionParams().sublist("BEAM TO SOLID CONTACT");

  /****************************************************************************/
  // get and check required parameters
  /****************************************************************************/

  /****************************************************************************/
  // get penalty parameter
  BTSVOLMT_penalty_param_ =
      beam_to_solid_contact_params_list.get<double>("BEAMS_BTSVOLMTPENALTYPARAM");

  if (BTSVOLMT_penalty_param_ < 0.0)
    dserror("beam-to-volume-meshtying penalty parameter must not be negative!");


  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams::Setup()
{
  CheckInit();

  // empty for now

  issetup_ = true;
}
