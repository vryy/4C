/*----------------------------------------------------------------------*/
/*! \file

\brief Data container holding all beam to solid volume meshtying input parameters.

\level 3
*/


#include "4C_beaminteraction_beam_to_solid_volume_meshtying_params.hpp"

#include "4C_beaminteraction_beam_to_solid_volume_meshtying_visualization_output_params.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_geometry_pair.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams::BeamToSolidVolumeMeshtyingParams()
    : BeamToSolidParamsBase(),
      integration_points_circumference_(0),
      rotational_coupling_triad_construction_(
          Inpar::BeamToSolid::BeamToSolidRotationCoupling::none),
      rotational_coupling_penalty_parameter_(0.0),
      output_params_ptr_(Teuchos::null),
      couple_restart_state_(false)
{
  // Empty Constructor.
}


/**
 *
 */
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams::init()
{
  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_to_solid_contact_params_list =
      Global::Problem::Instance()->beam_interaction_params().sublist(
          "BEAM TO SOLID VOLUME MESHTYING");

  // Set the common beam-to-solid parameters.
  SetBaseParams(beam_to_solid_contact_params_list);

  // Get parameters form input file.
  {
    // Number of integrations points along the circumference of the cross section.
    integration_points_circumference_ =
        beam_to_solid_contact_params_list.get<int>("INTEGRATION_POINTS_CIRCUMFERENCE");

    // Type of rotational coupling.
    rotational_coupling_triad_construction_ =
        Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidRotationCoupling>(
            beam_to_solid_contact_params_list, "ROTATION_COUPLING");
    rotational_coupling_ = rotational_coupling_triad_construction_ !=
                           Inpar::BeamToSolid::BeamToSolidRotationCoupling::none;

    // Mortar contact discretization to be used.
    mortar_shape_function_rotation_ =
        Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidMortarShapefunctions>(
            beam_to_solid_contact_params_list, "ROTATION_COUPLING_MORTAR_SHAPE_FUNCTION");
    if (get_contact_discretization() ==
            Inpar::BeamToSolid::BeamToSolidContactDiscretization::mortar and
        rotational_coupling_ and
        mortar_shape_function_rotation_ ==
            Inpar::BeamToSolid::BeamToSolidMortarShapefunctions::none)
      FOUR_C_THROW(
          "If mortar coupling and rotational coupling are activated, the shape function type for "
          "rotational coupling has to be explicitly given.");

    // Penalty parameter for rotational coupling.
    rotational_coupling_penalty_parameter_ =
        beam_to_solid_contact_params_list.get<double>("ROTATION_COUPLING_PENALTY_PARAMETER");

    // If the restart configuration should be coupled.
    couple_restart_state_ = (bool)Core::UTILS::IntegralValue<int>(
        beam_to_solid_contact_params_list, "COUPLE_RESTART_STATE");
  }

  // Setup the output parameter object.
  {
    output_params_ptr_ = Teuchos::rcp<BeamToSolidVolumeMeshtyingVisualizationOutputParams>(
        new BeamToSolidVolumeMeshtyingVisualizationOutputParams());
    output_params_ptr_->init();
    output_params_ptr_->setup();
  }

  // Sanity checks.
  if (rotational_coupling_ and couple_restart_state_)
    FOUR_C_THROW(
        "Coupling restart state combined with rotational coupling is not yet implemented!");

  isinit_ = true;
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamToSolidVolumeMeshtyingVisualizationOutputParams>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams::get_visualization_output_params_ptr()
{
  return output_params_ptr_;
};

FOUR_C_NAMESPACE_CLOSE
