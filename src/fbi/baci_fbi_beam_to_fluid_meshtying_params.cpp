/*----------------------------------------------------------------------*/
/*! \file

\brief Data container holding all beam to fluid volume meshtying input parameters.

\level 1
*/


#include "baci_fbi_beam_to_fluid_meshtying_params.H"

#include "baci_fbi_beam_to_fluid_meshtying_output_params.H"
#include "baci_inpar_fbi.H"
#include "baci_inpar_geometry_pair.H"
#include "baci_lib_globalproblem.H"
#include "baci_utils_exceptions.H"

BACI_NAMESPACE_OPEN

/**
 *
 */
FBI::BeamToFluidMeshtyingParams::BeamToFluidMeshtyingParams()
    : isinit_(false),
      issetup_(false),
      constraint_enforcement_(INPAR::FBI::BeamToFluidConstraintEnforcement::none),
      meshtying_discretization_(INPAR::FBI::BeamToFluidDiscretization::none),
      penalty_parameter_(-1.0),
      gauss_rule_(CORE::FE::GaussRule1D::undefined),
      calcfluidweakdirichletforce_(false),
      mortar_shape_function_(INPAR::FBI::BeamToFluidMeshtingMortarShapefunctions::none)
{
  // Empty Constructor.
}

/**
 *
 */
void FBI::BeamToFluidMeshtyingParams::Init()
{
  // Teuchos parameter list for beam contact
  const Teuchos::ParameterList& beam_to_fluid_meshtying_params_list =
      DRT::Problem::Instance()->getParameterList()->sublist("FLUID BEAM INTERACTION");

  // Get parameters form input file.
  // Constraint enforcement.
  constraint_enforcement_ = Teuchos::getIntegralValue<INPAR::FBI::BeamToFluidConstraintEnforcement>(
      beam_to_fluid_meshtying_params_list.sublist("BEAM TO FLUID MESHTYING"),
      "CONSTRAINT_STRATEGY");

  // Constraint enforcement.
  mortar_shape_function_ =
      Teuchos::getIntegralValue<INPAR::FBI::BeamToFluidMeshtingMortarShapefunctions>(
          beam_to_fluid_meshtying_params_list.sublist("BEAM TO FLUID MESHTYING"),
          "MORTAR_SHAPE_FUNCTION");

  // Contact discretization to be used.
  meshtying_discretization_ = Teuchos::getIntegralValue<INPAR::FBI::BeamToFluidDiscretization>(
      beam_to_fluid_meshtying_params_list.sublist("BEAM TO FLUID MESHTYING"),
      "MESHTYING_DISCRETIZATION");

  // Penalty parameter.
  penalty_parameter_ = (beam_to_fluid_meshtying_params_list.sublist("BEAM TO FLUID MESHTYING"))
                           .get<double>("PENALTY_PARAMETER");
  if (penalty_parameter_ < 0.0)
    dserror("beam-to-volume-meshtying penalty parameter must not be negative!");

  // Gauss rule for integration along the beam (segments).
  gauss_rule_ = INPAR::GEOMETRYPAIR::IntToGaussRule1D(
      beam_to_fluid_meshtying_params_list.sublist("BEAM TO FLUID MESHTYING")
          .get<int>("GAUSS_POINTS"));
  isinit_ = true;

  // Create and get visualization output parameter
  output_params_ = Teuchos::rcp<FBI::BeamToFluidMeshtyingVtkOutputParams>(
      new FBI::BeamToFluidMeshtyingVtkOutputParams());
  output_params_->Init();
  output_params_->Setup();
}


/**
 *
 */
void FBI::BeamToFluidMeshtyingParams::Setup()
{
  CheckInit();

  // Empty for now.

  issetup_ = true;
}

BACI_NAMESPACE_CLOSE
