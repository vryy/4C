/*----------------------------------------------------------------------*/
/*! \file
\file beam_to_fluid_meshtying_params.cpp

\brief Data container holding all beam to fluid volume meshtying input parameters.

\level 3
\maintainer Nora Hagmeyer
*/


#include "../drt_lib/drt_globalproblem.H"
#include "beam_to_fluid_meshtying_params.H"
#include "../drt_inpar/inpar_fbi.H"


/**
 *
 */
FBI::BeamToFluidMeshtyingParams::BeamToFluidMeshtyingParams()
    : isinit_(false),
      issetup_(false),
      constraint_enforcement_(INPAR::FBI::BeamToFluidConstraintEnforcement::none),
      meshtying_discretization_(INPAR::FBI::BeamToFluidDiscretization::none),
      penalty_parameter_(-1.0),
      gauss_rule_(DRT::UTILS::GaussRule1D::intrule1D_undefined)
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
  gauss_rule_ = INPAR::FBI::IntToGaussRule1D(
      beam_to_fluid_meshtying_params_list.sublist("BEAM TO FLUID MESHTYING")
          .get<int>("GAUSS_POINTS"));
  isinit_ = true;
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
