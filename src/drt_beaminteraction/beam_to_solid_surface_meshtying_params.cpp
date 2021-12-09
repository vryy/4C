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
      output_params_ptr_(Teuchos::null),
      rotational_coupling_(false),
      rotational_coupling_penalty_parameter_(-1.0),
      rotational_coupling_triad_construction_(
          INPAR::BEAMTOSOLID::BeamToSolidSurfaceRotationCoupling::none)
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

    // Parameters for rotational coupling.
    rotational_coupling_ = (bool)DRT::INPUT::IntegralValue<int>(
        beam_to_solid_contact_params_list, "ROTATIONAL_COUPLING");
    rotational_coupling_penalty_parameter_ =
        beam_to_solid_contact_params_list.get<double>("ROTATIONAL_COUPLING_PENALTY_PARAMETER");
    rotational_coupling_triad_construction_ =
        Teuchos::getIntegralValue<INPAR::BEAMTOSOLID::BeamToSolidSurfaceRotationCoupling>(
            beam_to_solid_contact_params_list, "ROTATIONAL_COUPLING_SURFACE_TRIAD");

    if (rotational_coupling_)
    {
      switch (coupling_type_)
      {
        case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::consistent_fad:
        case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::displacement_fad:
        case INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::
            reference_configuration_forced_to_zero_fad:
          break;
        default:
          dserror(
              "Beam-to-solid surface coupling with rotational coupling is only implemented in "
              "combination with the \"fad\" variants of surface coupling.");
          break;
      }
    }
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
