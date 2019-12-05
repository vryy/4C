/*-----------------------------------------------------------*/
/*! \file
\file inpar_fbi.cpp

\brief input parameter for Fluid-Beam Interaction

\maintainer Nora Hagmeyer

\level 3

*/
/*-----------------------------------------------------------*/


#include "inpar_fbi.H"
#include "drt_validparameters.H"
#include "../drt_lib/drt_conditiondefinition.H"
#include "inpar_geometry_pair.H"

void INPAR::FBI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::Array<std::string> yesnotuple =
      tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

  Teuchos::ParameterList& fbi = list->sublist("FLUID BEAM INTERACTION", false, "");

  /*----------------------------------------------------------------------*/
  /* parameters for beam to fluid meshtying */

  Teuchos::ParameterList& beam_to_fluid_meshtying =
      fbi.sublist("BEAM TO FLUID MESHTYING", false, "");

  setStringToIntegralParameter<int>("COUPLING", "two-way", "Type of coupling for FBI",
      tuple<std::string>("two-way", "fluid", "solid"),
      tuple<int>(
          BeamToFluidCoupling::twoway, BeamToFluidCoupling::fluid, BeamToFluidCoupling::solid),
      &fbi);

  setStringToIntegralParameter<BeamToFluidDiscretization>("MESHTYING_DISCRETIZATION", "none",
      "Type of employed meshtying discretization",
      tuple<std::string>("none", "gauss_point_to_segment"),
      tuple<FBI::BeamToFluidDiscretization>(
          BeamToFluidDiscretization::none, BeamToFluidDiscretization::gauss_point_to_segment),
      &beam_to_fluid_meshtying);

  setStringToIntegralParameter<BeamToFluidConstraintEnforcement>("CONSTRAINT_STRATEGY", "none",
      "Type of employed constraint enforcement strategy", tuple<std::string>("none", "penalty"),
      tuple<BeamToFluidConstraintEnforcement>(
          BeamToFluidConstraintEnforcement::none, BeamToFluidConstraintEnforcement::penalty),
      &beam_to_fluid_meshtying);

  DoubleParameter("PENALTY_PARAMETER", 0.0, "Penalty parameter for beam-to-Fluid volume meshtying",
      &beam_to_fluid_meshtying);

  IntParameter("GAUSS_POINTS", 6, "Number of Gauss Points for the integral evaluations",
      &beam_to_fluid_meshtying);

  DoubleParameter("SEARCH_RADIUS", 1000, "Search radius for beam-to-fluid volume meshtying",
      &beam_to_fluid_meshtying);

  // Add the geometry pair input parameters.
  INPAR::GEOMETRYPAIR::SetValidParametersLineTo3D(beam_to_fluid_meshtying);

  /*----------------------------------------------------------------------*/

  // Create subsection for runtime vtk output.
  Teuchos::ParameterList& beam_to_fluid_meshtying_vtk =
      beam_to_fluid_meshtying.sublist("RUNTIME VTK OUTPUT", false, "");

  // Whether to write vtp output at all for beam to fluid meshtying.
  setStringToIntegralParameter<int>("WRITE_OUTPUT", "No",
      "Enable / disable beam-to-fluid mesh tying output.", yesnotuple, yesnovalue,
      &beam_to_fluid_meshtying_vtk);

  setStringToIntegralParameter<int>("NODAL_FORCES", "No",
      "Enable / disable output of the resulting nodal forces due to beam to Fluid interaction.",
      yesnotuple, yesnovalue, &beam_to_fluid_meshtying_vtk);

  setStringToIntegralParameter<int>("SEGMENTATION", "No",
      "Enable / disable output of segmentation points.", yesnotuple, yesnovalue,
      &beam_to_fluid_meshtying_vtk);

  setStringToIntegralParameter<int>("INTEGRATION_POINTS", "No",
      "Enable / disable output of used integration points. If the meshtying method has 'forces' at "
      "the integration point, they will also be output.",
      yesnotuple, yesnovalue, &beam_to_fluid_meshtying_vtk);

  // ...
}

void INPAR::FBI::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  /*-------------------------------------------------------------------*/
}
