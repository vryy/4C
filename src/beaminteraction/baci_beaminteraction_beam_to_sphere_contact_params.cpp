/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief data container holding all beam to sphere contact input parameters

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#include "baci_beaminteraction_beam_to_sphere_contact_params.H"

#include "baci_lib_globalproblem.H"

BACI_NAMESPACE_OPEN


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
BEAMINTERACTION::BeamToSphereContactParams::BeamToSphereContactParams()
    : isinit_(false), issetup_(false), penalty_parameter_(-1.0)
{
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/

void BEAMINTERACTION::BeamToSphereContactParams::Init()
{
  issetup_ = false;

  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_to_sphere_contact_params_list =
      DRT::Problem::Instance()->BeamInteractionParams().sublist("BEAM TO SPHERE CONTACT");

  penalty_parameter_ = beam_to_sphere_contact_params_list.get<double>("PENALTY_PARAMETER");

  if (penalty_parameter_ < 0.0) dserror("beam-to-sphere penalty parameter must not be negative!");


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

BACI_NAMESPACE_CLOSE
