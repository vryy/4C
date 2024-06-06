/*----------------------------------------------------------------------*/
/*! \file

\brief Data container holding all beam to solid volume meshtying input parameters.

\level 3
*/


#include "4C_beaminteraction_beam_to_solid_surface_meshtying_params.hpp"

#include "4C_beaminteraction_beam_to_solid_surface_visualization_output_params.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_geometry_pair.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams::BeamToSolidSurfaceMeshtyingParams()
    : BeamToSolidParamsBase(),
      coupling_type_(Inpar::BeamToSolid::BeamToSolidSurfaceCoupling::none),
      output_params_ptr_(Teuchos::null),
      rotational_coupling_penalty_parameter_(-1.0),
      rotational_coupling_triad_construction_(
          Inpar::BeamToSolid::BeamToSolidSurfaceRotationCoupling::none)
{
  // Empty Constructor.
}


/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams::Init()
{
  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_to_solid_surface_meshtying_params_list =
      Global::Problem::Instance()->beam_interaction_params().sublist(
          "BEAM TO SOLID SURFACE MESHTYING");

  // Set the common beam-to-solid parameters.
  SetBaseParams(beam_to_solid_surface_meshtying_params_list);

  // Get parameters form input file.
  {
    // Type of coupling evaluation to be used.
    coupling_type_ = Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidSurfaceCoupling>(
        beam_to_solid_surface_meshtying_params_list, "COUPLING_TYPE");

    // Parameters for rotational coupling.
    rotational_coupling_ = (bool)Core::UTILS::IntegralValue<int>(
        beam_to_solid_surface_meshtying_params_list, "ROTATIONAL_COUPLING");
    rotational_coupling_penalty_parameter_ =
        beam_to_solid_surface_meshtying_params_list.get<double>(
            "ROTATIONAL_COUPLING_PENALTY_PARAMETER");
    rotational_coupling_triad_construction_ =
        Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidSurfaceRotationCoupling>(
            beam_to_solid_surface_meshtying_params_list, "ROTATIONAL_COUPLING_SURFACE_TRIAD");

    if (rotational_coupling_)
    {
      switch (coupling_type_)
      {
        case Inpar::BeamToSolid::BeamToSolidSurfaceCoupling::consistent_fad:
        case Inpar::BeamToSolid::BeamToSolidSurfaceCoupling::displacement_fad:
        case Inpar::BeamToSolid::BeamToSolidSurfaceCoupling::
            reference_configuration_forced_to_zero_fad:
          break;
        default:
          FOUR_C_THROW(
              "Beam-to-solid surface coupling with rotational coupling is only implemented in "
              "combination with the \"fad\" variants of surface coupling.");
      }
    }
  }

  // Setup the output parameter object.
  {
    output_params_ptr_ = Teuchos::rcp<BeamToSolidSurfaceVisualizationOutputParams>(
        new BeamToSolidSurfaceVisualizationOutputParams());
    output_params_ptr_->Init();
    output_params_ptr_->Setup();
  }

  isinit_ = true;
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamToSolidSurfaceVisualizationOutputParams>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams::get_visualization_output_params_ptr()
{
  return output_params_ptr_;
};

FOUR_C_NAMESPACE_CLOSE
