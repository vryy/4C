/*----------------------------------------------------------------------*/
/*! \file

\brief Base data container holding data for beam-to-solid interactions.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_params_base.H"

#include "../drt_inpar/inpar_geometry_pair.H"
#include "../drt_lib/drt_globalproblem.H"


/**
 *
 */
BEAMINTERACTION::BeamToSolidParamsBase::BeamToSolidParamsBase()
    : isinit_(false),
      issetup_(false),
      constraint_enforcement_(INPAR::BEAMTOSOLID::BeamToSolidConstraintEnforcement::none),
      contact_discretization_(INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization::none),
      mortar_shape_function_(INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions::none),
      penalty_parameter_(-1.0),
      gauss_rule_(DRT::UTILS::GaussRule1D::intrule1D_undefined)
{
  // Empty Constructor.
}


/**
 *
 */
void BEAMINTERACTION::BeamToSolidParamsBase::SetBaseParams(
    const Teuchos::ParameterList& beam_to_solid_params_list)
{
  // Get parameters form input file.
  {
    // Constraint enforcement.
    constraint_enforcement_ =
        Teuchos::getIntegralValue<INPAR::BEAMTOSOLID::BeamToSolidConstraintEnforcement>(
            beam_to_solid_params_list, "CONSTRAINT_STRATEGY");

    // Contact discretization to be used.
    contact_discretization_ =
        Teuchos::getIntegralValue<INPAR::BEAMTOSOLID::BeamToSolidContactDiscretization>(
            beam_to_solid_params_list, "CONTACT_DISCRETIZATION");

    // Contact discretization to be used.
    mortar_shape_function_ =
        Teuchos::getIntegralValue<INPAR::BEAMTOSOLID::BeamToSolidMortarShapefunctions>(
            beam_to_solid_params_list, "MORTAR_SHAPE_FUNCTION");

    // Penalty parameter.
    penalty_parameter_ = beam_to_solid_params_list.get<double>("PENALTY_PARAMETER");
    if (penalty_parameter_ < 0.0)
      dserror("beam-to-volume-meshtying penalty parameter must not be negative!");

    // Gauss rule for integration along the beam (segments).
    gauss_rule_ =
        INPAR::GEOMETRYPAIR::IntToGaussRule1D(beam_to_solid_params_list.get<int>("GAUSS_POINTS"));
  }

  isinit_ = true;
}


/**
 *
 */
void BEAMINTERACTION::BeamToSolidParamsBase::Setup()
{
  CheckInit();

  // Empty for now.

  issetup_ = true;
}
