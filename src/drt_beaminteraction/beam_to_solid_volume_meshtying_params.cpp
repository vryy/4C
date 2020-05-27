/*----------------------------------------------------------------------*/
/*! \file

\brief Data container holding all beam to solid volume meshtying input parameters.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_volume_meshtying_params.H"

#include "beam_to_solid_volume_meshtying_vtk_output_params.H"

#include "../drt_inpar/inpar_geometry_pair.H"
#include "../drt_lib/drt_globalproblem.H"


/**
 *
 */
BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams::BeamToSolidVolumeMeshtyingParams()
    : BeamToSolidParamsBase(),
      integration_points_circumference_(0),
      output_params_ptr_(Teuchos::null),
      couple_restart_state_(false)
{
  // Empty Constructor.
}


/**
 *
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams::Init()
{
  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_to_solid_contact_params_list =
      DRT::Problem::Instance()->BeamInteractionParams().sublist("BEAM TO SOLID VOLUME MESHTYING");

  // Set the common beam-to-solid parameters.
  SetBaseParams(beam_to_solid_contact_params_list);

  // Get parameters form input file.
  {
    // Number of integrations points along the circumference of the cross section.
    integration_points_circumference_ =
        beam_to_solid_contact_params_list.get<int>("INTEGRATION_POINTS_CIRCUMFERENCE");

    // If the restart configuration should be coupled.
    couple_restart_state_ = (bool)DRT::INPUT::IntegralValue<int>(
        beam_to_solid_contact_params_list, "COUPLE_RESTART_STATE");
  }

  // Setup the output parameter object.
  {
    output_params_ptr_ = Teuchos::rcp<BeamToSolidVolumeMeshtyingVtkOutputParams>(
        new BeamToSolidVolumeMeshtyingVtkOutputParams());
    output_params_ptr_->Init();
    output_params_ptr_->Setup();
  }

  isinit_ = true;
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamToSolidVolumeMeshtyingVtkOutputParams>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams::GetVtkOuputParamsPtr()
{
  return output_params_ptr_;
};
