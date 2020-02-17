/*----------------------------------------------------------------------*/
/*! \file

\brief Data container holding all beam to solid volume meshtying input parameters.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_volume_meshtying_params.H"

#include "beam_to_solid_volume_meshtying_vtk_output_params.H"

#include "../drt_lib/drt_globalproblem.H"


/**
 *
 */
BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams::BeamToSolidVolumeMeshtyingParams()
    : isinit_(false),
      issetup_(false),
      constraint_enforcement_(INPAR::BEAMTOSOLID::BeamToSolidConstraintEnforcement::none),
      contact_discretization_(INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization::none),
      mortar_shape_function_(INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::none),
      penalty_parameter_(-1.0),
      gauss_rule_(DRT::UTILS::GaussRule1D::intrule1D_undefined),
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

  // Get parameters form input file.
  {
    // Constraint enforcement.
    constraint_enforcement_ =
        Teuchos::getIntegralValue<INPAR::BEAMTOSOLID::BeamToSolidConstraintEnforcement>(
            beam_to_solid_contact_params_list, "CONSTRAINT_STRATEGY");

    // Contact discretization to be used.
    contact_discretization_ =
        Teuchos::getIntegralValue<INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization>(
            beam_to_solid_contact_params_list, "CONTACT_DISCRETIZATION");

    // Contact discretization to be used.
    mortar_shape_function_ =
        Teuchos::getIntegralValue<INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions>(
            beam_to_solid_contact_params_list, "MORTAR_SHAPE_FUNCTION");

    // Penalty parameter.
    penalty_parameter_ = beam_to_solid_contact_params_list.get<double>("PENALTY_PARAMETER");
    if (penalty_parameter_ < 0.0)
      dserror("beam-to-volume-meshtying penalty parameter must not be negative!");

    // Gauss rule for integration along the beam (segments).
    gauss_rule_ = INPAR::BEAMTOSOLID::IntToGaussRule1D(
        beam_to_solid_contact_params_list.get<int>("GAUSS_POINTS"));

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
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams::Setup()
{
  CheckInit();

  // Empty for now.

  issetup_ = true;
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamToSolidVolumeMeshtyingVtkOutputParams>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingParams::GetVtkOuputParamsPtr()
{
  return output_params_ptr_;
};
