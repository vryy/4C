/*----------------------------------------------------------------------*/
/*! \file

\brief Base data container holding data for beam-to-solid interactions.

\level 3
*/


#include "4C_beaminteraction_beam_to_solid_params_base.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_geometry_pair.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BEAMINTERACTION::BeamToSolidParamsBase::BeamToSolidParamsBase()
    : isinit_(false),
      issetup_(false),
      constraint_enforcement_(Inpar::BeamToSolid::BeamToSolidConstraintEnforcement::none),
      contact_discretization_(Inpar::BeamToSolid::BeamToSolidContactDiscretization::none),
      mortar_shape_function_(Inpar::BeamToSolid::BeamToSolidMortarShapefunctions::none),
      penalty_parameter_(-1.0),
      gauss_rule_(Core::FE::GaussRule1D::undefined),
      rotational_coupling_(false)
{
  // Empty Constructor.
}


/**
 *
 */
void BEAMINTERACTION::BeamToSolidParamsBase::set_base_params(
    const Teuchos::ParameterList& beam_to_solid_params_list)
{
  // Get parameters form input file.
  {
    // Constraint enforcement.
    constraint_enforcement_ =
        Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidConstraintEnforcement>(
            beam_to_solid_params_list, "CONSTRAINT_STRATEGY");

    // Contact discretization to be used.
    contact_discretization_ =
        Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidContactDiscretization>(
            beam_to_solid_params_list, "CONTACT_DISCRETIZATION");

    // Contact discretization to be used.
    mortar_shape_function_ =
        Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidMortarShapefunctions>(
            beam_to_solid_params_list, "MORTAR_SHAPE_FUNCTION");

    // Penalty parameter.
    penalty_parameter_ = beam_to_solid_params_list.get<double>("PENALTY_PARAMETER");
    if (penalty_parameter_ < 0.0)
      FOUR_C_THROW("beam-to-volume-meshtying penalty parameter must not be negative!");

    // Gauss rule for integration along the beam (segments).
    gauss_rule_ =
        Inpar::GEOMETRYPAIR::IntToGaussRule1D(beam_to_solid_params_list.get<int>("GAUSS_POINTS"));
  }

  isinit_ = true;
}


/**
 *
 */
void BEAMINTERACTION::BeamToSolidParamsBase::setup()
{
  check_init();

  // Empty for now.

  issetup_ = true;
}

FOUR_C_NAMESPACE_CLOSE
