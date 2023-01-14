/*----------------------------------------------------------------------*/
/*! \file

\brief Data container holding all beam to solid surface contact input parameters.

\level 3
*/


#include "beaminteraction_beam_to_solid_surface_contact_params.H"

#include "beaminteraction_beam_to_solid_surface_vtk_output_params.H"
#include "lib_globalproblem.H"
#include "inpar_geometry_pair.H"


/**
 *
 */
BEAMINTERACTION::BeamToSolidSurfaceContactParams::BeamToSolidSurfaceContactParams()
    : BeamToSolidParamsBase(),
      contact_type_(INPAR::BEAMTOSOLID::BeamToSolidSurfaceContact::gap_variation),
      penalty_law_(INPAR::BEAMTOSOLID::BeamToSolidSurfaceContactPenaltyLaw::none),
      penalty_parameter_g0_(0.0),
      output_params_ptr_(Teuchos::null)
{
  // Empty Constructor.
}


/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceContactParams::Init()
{
  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_to_solid_contact_params_list =
      DRT::Problem::Instance()->BeamInteractionParams().sublist("BEAM TO SOLID SURFACE CONTACT");

  // Set the common beam-to-solid parameters.
  SetBaseParams(beam_to_solid_contact_params_list);

  // Get parameters form input file.
  {
    contact_type_ = Teuchos::getIntegralValue<INPAR::BEAMTOSOLID::BeamToSolidSurfaceContact>(
        beam_to_solid_contact_params_list, "CONTACT_TYPE");

    penalty_law_ =
        Teuchos::getIntegralValue<INPAR::BEAMTOSOLID::BeamToSolidSurfaceContactPenaltyLaw>(
            beam_to_solid_contact_params_list, "PENALTY_LAW");

    penalty_parameter_g0_ = beam_to_solid_contact_params_list.get<double>("PENALTY_PARAMETER_G0");
  }

  // Setup the output parameter object.
  {
    output_params_ptr_ =
        Teuchos::rcp<BeamToSolidSurfaceVtkOutputParams>(new BeamToSolidSurfaceVtkOutputParams());
    output_params_ptr_->Init();
    output_params_ptr_->Setup();
  }

  isinit_ = true;
}


/**
 *
 */
int BEAMINTERACTION::BeamToSolidSurfaceContactParams::GetFADOrder() const

{
  switch (GetContactType())
  {
    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceContact::gap_variation:
      return 1;
      break;
    case INPAR::BEAMTOSOLID::BeamToSolidSurfaceContact::potential:
      return 2;
      break;
    default:
      dserror("Got unexpected contact type.");
      return 0;
  }
}
