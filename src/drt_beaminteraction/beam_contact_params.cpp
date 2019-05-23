/*----------------------------------------------------------------------------*/
/*!

\brief data container holding pointers to all subcontainers that in turn hold
       all input parameters specific to their problem type

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/

#include "beam_contact_params.H"
#include "beam_to_beam_contact_params.H"
#include "beam_to_sphere_contact_params.H"
#include "beam_to_solid_volume_meshtying_params.H"
#include "beam_contact_runtime_vtk_output_params.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
BEAMINTERACTION::BeamContactParams::BeamContactParams()
    : beam_to_beam_contact_params_(Teuchos::null),
      beam_to_sphere_contact_params_(Teuchos::null),
      beam_to_solid_volume_meshtying_params_(Teuchos::null),
      beam_contact_runtime_vtk_output_params_(Teuchos::null)
{
  // empty constructor
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactParams::BuildBeamToBeamContactParams()
{
  beam_to_beam_contact_params_ = Teuchos::rcp(new BEAMINTERACTION::BeamToBeamContactParams());
  beam_to_beam_contact_params_->Init();
  beam_to_beam_contact_params_->Setup();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactParams::BuildBeamContactRuntimeVtkOutputParams()
{
  beam_contact_runtime_vtk_output_params_ =
      Teuchos::rcp(new BEAMINTERACTION::BeamContactRuntimeVtkOutputParams());
  beam_contact_runtime_vtk_output_params_->Init();
  beam_contact_runtime_vtk_output_params_->Setup();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactParams::BuildBeamToSphereContactParams()
{
  beam_to_sphere_contact_params_ = Teuchos::rcp(new BEAMINTERACTION::BeamToSphereContactParams());
  beam_to_sphere_contact_params_->Init();
  beam_to_sphere_contact_params_->Setup();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void BEAMINTERACTION::BeamContactParams::BuildBeamToSolidVolumeMeshtyingParams()
{
  beam_to_solid_volume_meshtying_params_ =
      Teuchos::rcp(new BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams());
  beam_to_solid_volume_meshtying_params_->Init();
  beam_to_solid_volume_meshtying_params_->Setup();
}
