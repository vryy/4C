/*-----------------------------------------------------------*/
/*! \file

\brief input parameter for beaminteraction


\level 2

*/
/*-----------------------------------------------------------*/


#include "4C_inpar_beaminteraction.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_inpar_beam_to_solid.hpp"
#include "4C_inpar_beamcontact.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


void Inpar::BEAMINTERACTION::BeamInteractionConditionsGetAll(
    std::vector<Inpar::BEAMINTERACTION::BeamInteractionConditions>& interactions)
{
  interactions = {Inpar::BEAMINTERACTION::BeamInteractionConditions::beam_to_beam_contact,
      Inpar::BEAMINTERACTION::BeamInteractionConditions::beam_to_beam_point_coupling,
      Inpar::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_volume_meshtying,
      Inpar::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_surface_meshtying,
      Inpar::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_surface_contact};
}

void Inpar::BEAMINTERACTION::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& beaminteraction = list->sublist("BEAM INTERACTION", false, "");

  setStringToIntegralParameter<int>("REPARTITIONSTRATEGY", "Adaptive",
      "Type of employed repartitioning strategy",
      tuple<std::string>("Adaptive", "adaptive", "Everydt", "everydt"),
      tuple<int>(repstr_adaptive, repstr_adaptive, repstr_everydt, repstr_everydt),
      &beaminteraction);

  setStringToIntegralParameter<SearchStrategy>("SEARCH_STRATEGY", "bruteforce_with_binning",
      "Type of search strategy used for finding coupling pairs",
      tuple<std::string>("bruteforce_with_binning", "bounding_volume_hierarchy"),
      tuple<SearchStrategy>(
          SearchStrategy::bruteforce_with_binning, SearchStrategy::bounding_volume_hierarchy),
      &beaminteraction);

  /*----------------------------------------------------------------------*/
  /* parameters for crosslinking submodel */

  Teuchos::ParameterList& crosslinking = beaminteraction.sublist("CROSSLINKING", false, "");

  // remove this some day
  Core::UTILS::BoolParameter("CROSSLINKER", "No", "Crosslinker in problem", &crosslinking);

  // bounding box for initial random crosslinker position
  Core::UTILS::StringParameter("INIT_LINKER_BOUNDINGBOX", "1e12 1e12 1e12 1e12 1e12 1e12",
      "Linker are initially set randomly within this bounding box", &crosslinking);

  // time step for stochastic events concerning crosslinking
  Core::UTILS::DoubleParameter("TIMESTEP", -1.0,
      "time step for stochastic events concerning crosslinking (e.g. diffusion, p_link, p_unlink) ",
      &crosslinking);
  // Reading double parameter for viscosity of background fluid
  Core::UTILS::DoubleParameter("VISCOSITY", 0.0, "viscosity", &crosslinking);
  // Reading double parameter for thermal energy in background fluid (temperature * Boltzmann
  // constant)
  Core::UTILS::DoubleParameter("KT", 0.0, "thermal energy", &crosslinking);
  // number of initial (are set right in the beginning) crosslinker of certain type
  setNumericStringParameter("MAXNUMINITCROSSLINKERPERTYPE", "0",
      "number of initial crosslinker of certain type (additional to NUMCROSSLINKERPERTYPE) ",
      &crosslinking);
  // number of crosslinker of certain type
  setNumericStringParameter(
      "NUMCROSSLINKERPERTYPE", "0", "number of crosslinker of certain type ", &crosslinking);
  // material number characterizing crosslinker type
  setNumericStringParameter("MATCROSSLINKERPERTYPE", "-1",
      "material number characterizing crosslinker type ", &crosslinking);
  // maximal number of binding partner per filament binding spot for each binding spot type
  setNumericStringParameter("MAXNUMBONDSPERFILAMENTBSPOT", "1",
      "maximal number of bonds per filament binding spot", &crosslinking);
  // distance between two binding spots on a filament (same on all filaments)
  setNumericStringParameter("FILAMENTBSPOTINTERVALGLOBAL", "-1.0",
      "distance between two binding spots on all filaments", &crosslinking);
  // distance between two binding spots on a filament (as percentage of current filament length)
  setNumericStringParameter("FILAMENTBSPOTINTERVALLOCAL", "-1.0",
      "distance between two binding spots on current filament", &crosslinking);
  // start and end for bspots on a filament in arc parameter (same on each filament independent of
  // their length)
  Core::UTILS::StringParameter("FILAMENTBSPOTRANGEGLOBAL", "-1.0 -1.0",
      "Lower and upper arc parameter bound for binding spots on a filament", &crosslinking);
  // start and end for bspots on a filament in percent of reference filament length
  Core::UTILS::StringParameter("FILAMENTBSPOTRANGELOCAL", "0.0 1.0",
      "Lower and upper arc parameter bound for binding spots on a filament", &crosslinking);


  /*----------------------------------------------------------------------*/
  /* parameters for sphere beam link submodel */

  Teuchos::ParameterList& spherebeamlink = beaminteraction.sublist("SPHERE BEAM LINK", false, "");

  Core::UTILS::BoolParameter("SPHEREBEAMLINKING", "No", "Integrins in problem", &spherebeamlink);

  // Reading double parameter for contraction rate for active linker
  Core::UTILS::DoubleParameter("CONTRACTIONRATE", 0.0,
      "contraction rate of cell (integrin linker) in [microm/s]", &spherebeamlink);
  // time step for stochastic events concerning sphere beam linking
  Core::UTILS::DoubleParameter("TIMESTEP", -1.0,
      "time step for stochastic events concerning sphere beam linking (e.g. catch-slip-bond "
      "behavior) ",
      &spherebeamlink);
  setNumericStringParameter(
      "MAXNUMLINKERPERTYPE", "0", "number of crosslinker of certain type ", &spherebeamlink);
  // material number characterizing crosslinker type
  setNumericStringParameter("MATLINKERPERTYPE", "-1",
      "material number characterizing crosslinker type ", &spherebeamlink);
  // distance between two binding spots on a filament (same on all filaments)
  setNumericStringParameter("FILAMENTBSPOTINTERVALGLOBAL", "-1.0",
      "distance between two binding spots on all filaments", &spherebeamlink);
  // distance between two binding spots on a filament (as percentage of current filament length)
  setNumericStringParameter("FILAMENTBSPOTINTERVALLOCAL", "-1.0",
      "distance between two binding spots on current filament", &spherebeamlink);
  // start and end for bspots on a filament in arc parameter (same on each filament independent of
  // their length)
  Core::UTILS::StringParameter("FILAMENTBSPOTRANGEGLOBAL", "-1.0 -1.0",
      "Lower and upper arc parameter bound for binding spots on a filament", &spherebeamlink);
  // start and end for bspots on a filament in percent of reference filament length
  Core::UTILS::StringParameter("FILAMENTBSPOTRANGELOCAL", "0.0 1.0",
      "Lower and upper arc parameter bound for binding spots on a filament", &spherebeamlink);

  /*----------------------------------------------------------------------*/
  /* parameters for beam to ? contact submodel*/
  /*----------------------------------------------------------------------*/

  /*----------------------------------------------------------------------*/
  /* parameters for beam to beam contact */

  Teuchos::ParameterList& beamtobeamcontact =
      beaminteraction.sublist("BEAM TO BEAM CONTACT", false, "");

  setStringToIntegralParameter<int>("STRATEGY", "None", "Type of employed solving strategy",
      tuple<std::string>("None", "none", "Penalty", "penalty"),
      tuple<int>(bstr_none, bstr_none, bstr_penalty, bstr_penalty), &beamtobeamcontact);

  // ...

  /*----------------------------------------------------------------------*/
  /* parameters for beam to sphere contact */

  Teuchos::ParameterList& beamtospherecontact =
      beaminteraction.sublist("BEAM TO SPHERE CONTACT", false, "");

  setStringToIntegralParameter<int>("STRATEGY", "None", "Type of employed solving strategy",
      tuple<std::string>("None", "none", "Penalty", "penalty"),
      tuple<int>(bstr_none, bstr_none, bstr_penalty, bstr_penalty), &beamtospherecontact);

  Core::UTILS::DoubleParameter("PENALTY_PARAMETER", 0.0,
      "Penalty parameter for beam-to-rigidsphere contact", &beamtospherecontact);

  // ...

  /*----------------------------------------------------------------------*/
  /* parameters for beam to solid contact */
  BeamToSolid::SetValidParameters(list);
}

void Inpar::BEAMINTERACTION::SetValidConditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  /*-------------------------------------------------------------------*/
  // beam potential interaction: atom/charge density per unit length on LINE
  Teuchos::RCP<Core::Conditions::ConditionDefinition> beam_filament_condition =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN LINE BEAM FILAMENT CONDITIONS",
          "BeamLineFilamentCondition", "Beam_Line_Filament_Condition",
          Core::Conditions::FilamentBeamLineCondition, false,
          Core::Conditions::geometry_type_line));

  beam_filament_condition->AddComponent(Teuchos::rcp(new Input::SeparatorComponent("ID")));
  beam_filament_condition->AddComponent(Teuchos::rcp(new Input::IntComponent("FilamentId")));
  beam_filament_condition->AddComponent(
      Teuchos::rcp(new Input::SeparatorComponent("TYPE", "", true)));
  beam_filament_condition->AddComponent(
      Teuchos::rcp(new Input::SelectionComponent("Type", "Arbitrary",
          Teuchos::tuple<std::string>(
              "Arbitrary", "arbitrary", "Actin", "actin", "Collagen", "collagen"),
          Teuchos::tuple<std::string>(
              "Arbitrary", "arbitrary", "Actin", "actin", "Collagen", "collagen"),
          true)));

  condlist.push_back(beam_filament_condition);

  /*-------------------------------------------------------------------*/
  Teuchos::RCP<Core::Conditions::ConditionDefinition> penalty_coupling_condition = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN POINT PENALTY COUPLING CONDITIONS",
          "PenaltyPointCouplingCondition", "Couples beam nodes that lie on the same position",
          Core::Conditions::PenaltyPointCouplingCondition, false,
          Core::Conditions::geometry_type_point));

  Input::AddNamedReal(penalty_coupling_condition, "POSITIONAL_PENALTY_PARAMETER");
  Input::AddNamedReal(penalty_coupling_condition, "ROTATIONAL_PENALTY_PARAMETER");

  condlist.push_back(penalty_coupling_condition);

  // beam-to-beam interactions
  Inpar::BEAMCONTACT::SetValidConditions(condlist);

  // beam-to-solid interactions
  Inpar::BeamToSolid::SetValidConditions(condlist);
}

FOUR_C_NAMESPACE_CLOSE
