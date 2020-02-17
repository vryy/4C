/*----------------------------------------------------------------------*/
/*! \file

\brief Data container holding all beam to solid volume meshtying input parameters.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_surface_meshtying_params.H"

#include "../drt_lib/drt_globalproblem.H"


/**
 *
 */
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams::BeamToSolidSurfaceMeshtyingParams()
    : isinit_(false),
      issetup_(false),
      constraint_enforcement_(INPAR::BEAMTOSOLID::BeamToSolidConstraintEnforcement::none),
      contact_discretization_(INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization::none),
      coupling_type_(INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling::none),
      mortar_shape_function_(INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::none),
      penalty_parameter_(-1.0),
      gauss_rule_(DRT::UTILS::GaussRule1D::intrule1D_undefined)
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

    // Type of coupling evaluation to be used.
    coupling_type_ = Teuchos::getIntegralValue<INPAR::BEAMTOSOLID::BeamToSolidSurfaceCoupling>(
        beam_to_solid_contact_params_list, "COUPLING_TYPE");

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
  }

  isinit_ = true;
}


/**
 *
 */
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingParams::Setup()
{
  CheckInit();

  // Empty for now.

  issetup_ = true;
}
