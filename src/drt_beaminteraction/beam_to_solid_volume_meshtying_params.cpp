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
      constraint_enforcement_(INPAR::BEAMINTERACTION::BeamToSolidVolumeConstraintEnforcement::none),
      contact_discretization_(INPAR::BEAMINTERACTION::BeamToSolidVolumeContactDiscretization::none),
      mortar_shape_function_(INPAR::BEAMINTERACTION::BeamToSolidVolumeMortarShapefunctions::none),
      penalty_parameter_(-1.0),
      gauss_rule_(DRT::UTILS::GaussRule1D::intrule1D_undefined)
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
        Teuchos::getIntegralValue<INPAR::BEAMINTERACTION::BeamToSolidVolumeConstraintEnforcement>(
            beam_to_solid_contact_params_list, "CONSTRAINT_STRATEGY");

    // Contact discretization to be used.
    contact_discretization_ =
        Teuchos::getIntegralValue<INPAR::BEAMINTERACTION::BeamToSolidVolumeContactDiscretization>(
            beam_to_solid_contact_params_list, "CONTACT_DISCRETIZATION");

    // Contact discretization to be used.
    mortar_shape_function_ =
        Teuchos::getIntegralValue<INPAR::BEAMINTERACTION::BeamToSolidVolumeMortarShapefunctions>(
            beam_to_solid_contact_params_list, "MORTAR_SHAPE_FUNCTION");

    // Penalty parameter.
    penalty_parameter_ = beam_to_solid_contact_params_list.get<double>("PENALTY_PARAMETER");
    if (penalty_parameter_ < 0.0)
      dserror("beam-to-volume-meshtying penalty parameter must not be negative!");

    // Gauss rule for integration along the beam (segments).
    gauss_rule_ = INPAR::BEAMINTERACTION::IntToGaussRule1D(
        beam_to_solid_contact_params_list.get<int>("GAUSS_POINTS"));

    // Number of integrations points along the circumfence of the cross section.
    integration_points_circumfence_ =
        beam_to_solid_contact_params_list.get<int>("INTEGRATION_POINTS_CIRCUMFENCE");
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
