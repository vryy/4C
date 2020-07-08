/*----------------------------------------------------------------------*/
/*! \file

\brief Data container holding all beam to solid volume meshtying input parameters.

\level 3
*/


#include "beam_to_solid_surface_meshtying_params.H"

#include "beam_to_solid_surface_vtk_output_params.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_geometry_pair.H"


/**
 *
 */
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams::BeamToSolidSurfaceMeshtyingParams()
    : BeamToSolidParamsBase(),
      coupling_type_(INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::none),
      output_params_ptr_(Teuchos::null)
{
  // Empty Constructor.
}


/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams::Init()
{
  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_to_solid_contact_params_list =
      DRT::Problem::Instance()->BeamInteractionParams().sublist("BEAM TO SOLID SURFACE MESHTYING");

  // Set the common beam-to-solid parameters.
  SetBaseParams(beam_to_solid_contact_params_list);

  // Get parameters form input file.
  {
    // Type of coupling evaluation to be used.
    coupling_type_ = Teuchos::getIntegralValue<INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling>(
        beam_to_solid_contact_params_list, "COUPLING_TYPE");
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
Teuchos::RCP<BEAMINTERACTION::BeamToSolidSurfaceVtkOutputParams>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams::GetVtkOuputParamsPtr()
{
  return output_params_ptr_;
};
