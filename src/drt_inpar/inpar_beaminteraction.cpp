/*-----------------------------------------------------------*/
/*! \file

\brief input parameter for beaminteraction


\level 2

*/
/*-----------------------------------------------------------*/


#include "inpar_beaminteraction.H"
#include "drt_validparameters.H"
#include "../drt_lib/drt_conditiondefinition.H"
#include "inpar_beam_to_solid.H"


void INPAR::BEAMINTERACTION::BeamInteractionConditionsGetAll(
    std::vector<INPAR::BEAMINTERACTION::BeamInteractionConditions>& interactions)
{
  interactions = {INPAR::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_volume_meshtying,
      INPAR::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_surface_meshtying};
}

void INPAR::BEAMINTERACTION::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::Array<std::string> yesnotuple =
      tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

  Teuchos::ParameterList& beaminteraction = list->sublist("BEAM INTERACTION", false, "");

  setStringToIntegralParameter<int>("REPARTITIONSTRATEGY", "Adaptive",
      "Type of employed repartitioning strategy",
      tuple<std::string>("Adaptive", "adaptive", "Everydt", "everydt"),
      tuple<int>(repstr_adaptive, repstr_adaptive, repstr_everydt, repstr_everydt),
      &beaminteraction);


  /*----------------------------------------------------------------------*/
  /* parameters for crosslinking submodel */

  Teuchos::ParameterList& crosslinking = beaminteraction.sublist("CROSSLINKING", false, "");

  // remove this some day
  setStringToIntegralParameter<int>(
      "CROSSLINKER", "No", "Crosslinker in problem", yesnotuple, yesnovalue, &crosslinking);

  // bounding box for initial random crosslinker position
  StringParameter("INIT_LINKER_BOUNDINGBOX", "1e12 1e12 1e12 1e12 1e12 1e12",
      "Linker are initially set randomly within this bounding box", &crosslinking);

  // time step for stochastic events concerning crosslinking
  DoubleParameter("TIMESTEP", -1.0,
      "time step for stochastic events concerning crosslinking (e.g. diffusion, p_link, p_unlink) ",
      &crosslinking);
  // Reading double parameter for viscosity of background fluid
  DoubleParameter("VISCOSITY", 0.0, "viscosity", &crosslinking);
  // Reading double parameter for thermal energy in background fluid (temperature * Boltzmann
  // constant)
  DoubleParameter("KT", 0.0, "thermal energy", &crosslinking);
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
  StringParameter("FILAMENTBSPOTRANGEGLOBAL", "-1.0 -1.0",
      "Lower and upper arc parameter bound for binding spots on a filament", &crosslinking);
  // start and end for bspots on a filament in percent of reference filament length
  StringParameter("FILAMENTBSPOTRANGELOCAL", "0.0 1.0",
      "Lower and upper arc parameter bound for binding spots on a filament", &crosslinking);


  /*----------------------------------------------------------------------*/
  /* parameters for sphere beam link submodel */

  Teuchos::ParameterList& spherebeamlink = beaminteraction.sublist("SPHERE BEAM LINK", false, "");

  setStringToIntegralParameter<int>(
      "SPHEREBEAMLINKING", "No", "Integrins in problem", yesnotuple, yesnovalue, &spherebeamlink);

  // Reading double parameter for contraction rate for active linker
  DoubleParameter("CONTRACTIONRATE", 0.0, "contraction rate of cell (integrin linker) in [μm/s]",
      &spherebeamlink);
  // time step for stochastic events concerning sphere beam linking
  DoubleParameter("TIMESTEP", -1.0,
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
  StringParameter("FILAMENTBSPOTRANGEGLOBAL", "-1.0 -1.0",
      "Lower and upper arc parameter bound for binding spots on a filament", &spherebeamlink);
  // start and end for bspots on a filament in percent of reference filament length
  StringParameter("FILAMENTBSPOTRANGELOCAL", "0.0 1.0",
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

  // ...

  /*----------------------------------------------------------------------*/
  /* parameters for beam to solid contact */
  BEAMTOSOLID::SetValidParameters(list);
}

void INPAR::BEAMINTERACTION::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  /*-------------------------------------------------------------------*/
  // beam potential interaction: atom/charge density per unit length on LINE
  Teuchos::RCP<ConditionDefinition> beam_filament_condition =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE BEAM FILAMENT CONDITIONS",
          "BeamLineFilamentCondition", "Beam_Line_Filament_Condition",
          DRT::Condition::FilamentBeamLineCondition, false, DRT::Condition::Line));

  beam_filament_condition->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("ID")));
  beam_filament_condition->AddComponent(
      Teuchos::rcp(new IntConditionComponent("FilamentId", false, false)));
  beam_filament_condition->AddComponent(
      Teuchos::rcp(new SeparatorConditionComponent("TYPE", true)));
  beam_filament_condition->AddComponent(
      Teuchos::rcp(new StringConditionComponent("Type", "Arbitrary",
          Teuchos::tuple<std::string>(
              "Arbitrary", "arbitrary", "Actin", "actin", "Collagen", "collagen"),
          Teuchos::tuple<std::string>(
              "Arbitrary", "arbitrary", "Actin", "actin", "Collagen", "collagen"),
          true)));

  condlist.push_back(beam_filament_condition);

  // beam-to-solid interaction
  INPAR::BEAMTOSOLID::SetValidConditions(condlist);
}
