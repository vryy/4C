/*-----------------------------------------------------------------------------------------------*/
/*!

\brief data container holding all beam to sphere contact input parameters

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam_to_sphere_contact_params.H"

#include "../drt_lib/drt_globalproblem.H"


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BEAMINTERACTION::BeamToSphereContactParams::BeamToSphereContactParams()
    : isinit_(false), issetup_(false), BTSPH_penalty_param_(-1.0)
{
  // empty constructor
}
/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/

void BEAMINTERACTION::BeamToSphereContactParams::Init()
{
  issetup_ = false;

  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_to_beam_contact_params_list =
      DRT::Problem::Instance()->BeamContactParams();

  BTSPH_penalty_param_ = beam_to_beam_contact_params_list.get<double>("BEAMS_BTSPH_PENALTYPARAM");

  if (BTSPH_penalty_param_ < 0.0) dserror("beam-to-sphere penalty parameter must not be negative!");


  isinit_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamToSphereContactParams::Setup()
{
  CheckInit();

  // empty for now

  issetup_ = true;
}
