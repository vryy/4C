/*----------------------------------------------------------------------*/
/*! \file
\brief input quantities and globally accessible enumerations for scatra-scatra interface coupling

\level 2


*/
/*----------------------------------------------------------------------*/
#include "4C_inpar_s2i.hpp"

#include "4C_discretization_condition_definition.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*------------------------------------------------------------------------*
 | set valid parameters for scatra-scatra interface coupling   fang 01/16 |
 *------------------------------------------------------------------------*/
void INPAR::S2I::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& s2icoupling =
      list->sublist("SCALAR TRANSPORT DYNAMIC", true)
          .sublist(
              "S2I COUPLING", false, "control parameters for scatra-scatra interface coupling");

  // type of mortar meshtying
  setStringToIntegralParameter<int>("COUPLINGTYPE", "Undefined", "type of mortar meshtying",
      tuple<std::string>("Undefined", "MatchingNodes", "StandardMortar", "SaddlePointMortar_Petrov",
          "SaddlePointMortar_Bubnov", "CondensedMortar_Petrov", "CondensedMortar_Bubnov",
          "StandardNodeToSegment"),
      tuple<int>(coupling_undefined, coupling_matching_nodes, coupling_mortar_standard,
          coupling_mortar_saddlepoint_petrov, coupling_mortar_saddlepoint_bubnov,
          coupling_mortar_condensed_petrov, coupling_mortar_condensed_bubnov,
          coupling_nts_standard),
      &s2icoupling);

  // flag for interface side underlying Lagrange multiplier definition
  setStringToIntegralParameter<int>("LMSIDE", "slave",
      "flag for interface side underlying Lagrange multiplier definition",
      tuple<std::string>("slave", "master"), tuple<int>(side_slave, side_master), &s2icoupling);

  // flag for evaluation of interface linearizations and residuals on slave side only
  CORE::UTILS::BoolParameter("SLAVEONLY", "No",
      "flag for evaluation of interface linearizations and residuals on slave side only",
      &s2icoupling);

  // node-to-segment projection tolerance
  CORE::UTILS::DoubleParameter(
      "NTSPROJTOL", 0.0, "node-to-segment projection tolerance", &s2icoupling);

  // flag for evaluation of scatra-scatra interface coupling involving interface layer growth
  setStringToIntegralParameter<int>("INTLAYERGROWTH_EVALUATION", "none",
      "flag for evaluation of scatra-scatra interface coupling involving interface layer growth",
      tuple<std::string>("none", "monolithic", "semi-implicit"),
      tuple<int>(
          growth_evaluation_none, growth_evaluation_monolithic, growth_evaluation_semi_implicit),
      &s2icoupling);

  // local Newton-Raphson convergence tolerance for scatra-scatra interface coupling involving
  // interface layer growth
  CORE::UTILS::DoubleParameter("INTLAYERGROWTH_CONVTOL", 1.e-12,
      "local Newton-Raphson convergence tolerance for scatra-scatra interface coupling involving "
      "interface layer growth",
      &s2icoupling);

  // maximum number of local Newton-Raphson iterations for scatra-scatra interface coupling
  // involving interface layer growth
  CORE::UTILS::IntParameter("INTLAYERGROWTH_ITEMAX", 5,
      "maximum number of local Newton-Raphson iterations for scatra-scatra interface coupling "
      "involving interface layer growth",
      &s2icoupling);

  // ID of linear solver for monolithic scatra-scatra interface coupling involving interface layer
  // growth
  CORE::UTILS::IntParameter("INTLAYERGROWTH_LINEAR_SOLVER", -1,
      "ID of linear solver for monolithic scatra-scatra interface coupling involving interface "
      "layer growth",
      &s2icoupling);

  // modified time step size for scatra-scatra interface coupling involving interface layer growth
  CORE::UTILS::DoubleParameter("INTLAYERGROWTH_TIMESTEP", -1.,
      "modified time step size for scatra-scatra interface coupling involving interface layer "
      "growth",
      &s2icoupling);

  CORE::UTILS::BoolParameter("MESHTYING_CONDITIONS_INDEPENDENT_SETUP", "No",
      "mesh tying for different conditions should be setup independently", &s2icoupling);

  CORE::UTILS::BoolParameter("OUTPUT_INTERFACE_FLUX", "No",
      "evaluate integral of coupling flux on slave side for each s2i condition and write it to csv "
      "file",
      &s2icoupling);
}


/*------------------------------------------------------------------------*
 | set valid conditions for scatra-scatra interface coupling   fang 01/16 |
 *------------------------------------------------------------------------*/
void INPAR::S2I::SetValidConditions(
    std::vector<Teuchos::RCP<CORE::Conditions::ConditionDefinition>>& condlist)
{
  using namespace INPUT;

  /*--------------------------------------------------------------------*/
  // scatra-scatra interface mesh tying condition
  {
    // definition of scatra-scatra interface mesh tying line condition
    auto s2imeshtyingline = Teuchos::rcp(
        new CORE::Conditions::ConditionDefinition("DESIGN S2I MESHTYING LINE CONDITIONS",
            "S2IMeshtying", "Scatra-scatra line interface mesh tying",
            CORE::Conditions::S2IMeshtying, true, CORE::Conditions::geometry_type_line));

    // definition of scatra-scatra interface mesh tying surface condition
    auto s2imeshtyingsurf = Teuchos::rcp(
        new CORE::Conditions::ConditionDefinition("DESIGN S2I MESHTYING SURF CONDITIONS",
            "S2IMeshtying", "Scatra-scatra surface interface mesh tying",
            CORE::Conditions::S2IMeshtying, true, CORE::Conditions::geometry_type_surface));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<INPUT::LineComponent>> s2imeshtyingcomponents;
    s2imeshtyingcomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent("ConditionID")));
    s2imeshtyingcomponents.emplace_back(Teuchos::rcp(new INPUT::SelectionComponent("interface side",
        "Undefined", Teuchos::tuple<std::string>("Undefined", "Slave", "Master"),
        Teuchos::tuple<int>(
            INPAR::S2I::side_undefined, INPAR::S2I::side_slave, INPAR::S2I::side_master))));
    s2imeshtyingcomponents.emplace_back(
        Teuchos::rcp(new INPUT::SeparatorComponent("S2I_KINETICS_ID")));
    s2imeshtyingcomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent("S2IKineticsID")));

    // insert input file line components into condition definitions
    for (auto& conditioncomponent : s2imeshtyingcomponents)
    {
      s2imeshtyingline->AddComponent(conditioncomponent);
      s2imeshtyingsurf->AddComponent(conditioncomponent);
    }

    condlist.push_back(s2imeshtyingline);
    condlist.push_back(s2imeshtyingsurf);
  }

  /*--------------------------------------------------------------------*/
  // scatra-scatra interface no evaluation condition
  {
    // definition of scatra-scatra interface no evaluation line condition
    auto s2inoevaluationline = Teuchos::rcp(new CORE::Conditions::ConditionDefinition(
        "DESIGN S2I NO EVALUATION LINE CONDITIONS", "S2INoEvaluation",
        "Scatra-scatra interface no evaluation line condition. This condition can be used to "
        "deactivate the evaluation of the corresponding `S2IKinetics` condition. Another usage "
        "is in coupled algorithms, where a specific `S2IKinetics` condition should not be "
        "evaluated within the scalar transport framework, as it already evaluated elsewhere. "
        "One example is the SSI contact.",
        CORE::Conditions::S2INoEvaluation, true, CORE::Conditions::geometry_type_line));

    // definition of scatra-scatra interface no evaluation surface condition
    auto s2inoevaluationsurf = Teuchos::rcp(new CORE::Conditions::ConditionDefinition(
        "DESIGN S2I NO EVALUATION SURF CONDITIONS", "S2INoEvaluation",
        "Scatra-scatra interface no evaluation surface condition. This condition can be used to "
        "deactivate the evaluation of the corresponding `S2IKinetics` condition. Another usage "
        "is in coupled algorithms, where a specific `S2IKinetics` condition should not be "
        "evaluated within the scalar transport framework, as it already evaluated elsewhere. "
        "One example is the SSI contact.",
        CORE::Conditions::S2INoEvaluation, true, CORE::Conditions::geometry_type_surface));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<INPUT::LineComponent>> s2inoevaluationcomponents;
    s2inoevaluationcomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent("ConditionID")));
    s2inoevaluationcomponents.emplace_back(Teuchos::rcp(new INPUT::SelectionComponent(
        "interface side", "Undefined", Teuchos::tuple<std::string>("Undefined", "Slave", "Master"),
        Teuchos::tuple<int>(
            INPAR::S2I::side_undefined, INPAR::S2I::side_slave, INPAR::S2I::side_master))));
    s2inoevaluationcomponents.emplace_back(
        Teuchos::rcp(new INPUT::SeparatorComponent("S2I_KINETICS_ID")));
    s2inoevaluationcomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent("S2IKineticsID")));

    // insert input file line components into condition definitions
    for (auto& conditioncomponent : s2inoevaluationcomponents)
    {
      s2inoevaluationline->AddComponent(conditioncomponent);
      s2inoevaluationsurf->AddComponent(conditioncomponent);
    }

    condlist.push_back(s2inoevaluationline);
    condlist.push_back(s2inoevaluationsurf);
  }

  /*--------------------------------------------------------------------*/
  // scatra-scatra interface kinetics condition
  {
    // definition of scatra-scatra interface kinetics point condition
    auto s2ikineticspoint = Teuchos::rcp(
        new CORE::Conditions::ConditionDefinition("DESIGN S2I KINETICS POINT CONDITIONS",
            "S2IKinetics", "Scatra-scatra line interface kinetics", CORE::Conditions::S2IKinetics,
            true, CORE::Conditions::geometry_type_point));

    // definition of scatra-scatra interface kinetics line condition
    auto s2ikineticsline = Teuchos::rcp(
        new CORE::Conditions::ConditionDefinition("DESIGN S2I KINETICS LINE CONDITIONS",
            "S2IKinetics", "Scatra-scatra line interface kinetics", CORE::Conditions::S2IKinetics,
            true, CORE::Conditions::geometry_type_line));

    // definition of scatra-scatra interface kinetics surface condition
    auto s2ikineticssurf = Teuchos::rcp(
        new CORE::Conditions::ConditionDefinition("DESIGN S2I KINETICS SURF CONDITIONS",
            "S2IKinetics", "Scatra-scatra surface interface kinetics",
            CORE::Conditions::S2IKinetics, true, CORE::Conditions::geometry_type_surface));

    // Macro-micro coupling condition for micro scale in multi-scale scalar transport problems
    auto multiscalecouplingpoint = Teuchos::rcp(new CORE::Conditions::ConditionDefinition(
        "DESIGN SCATRA MULTI-SCALE COUPLING POINT CONDITIONS", "ScatraMultiScaleCoupling",
        "Scalar transport multi-scale coupling condition",
        CORE::Conditions::ScatraMultiScaleCoupling, false, CORE::Conditions::geometry_type_point));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<INPUT::LineComponent>> s2icomponents;
    {
      // interface ID
      s2icomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent("ConditionID")));

      // interface sides for scatra-scatra interface kinetics
      std::map<int, std::pair<std::string, std::vector<Teuchos::RCP<INPUT::LineComponent>>>>
          interface_choices;
      {
        {
          // undefined side
          std::vector<Teuchos::RCP<INPUT::LineComponent>> undefined_side;
          interface_choices.emplace(side_undefined, std::make_pair("Undefined", undefined_side));
        }

        {
          // slave side
          std::vector<Teuchos::RCP<INPUT::LineComponent>> slaveside;

          // Collect the diffent model selection choices in a map.
          std::map<int, std::pair<std::string, std::vector<Teuchos::RCP<INPUT::LineComponent>>>>
              kinetic_model_choices;
          {
            {
              // constant and linear permeability
              std::vector<Teuchos::RCP<INPUT::LineComponent>> constlinperm;

              constlinperm.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("NUMSCAL")));
              constlinperm.emplace_back(Teuchos::rcp(new INPUT::IntComponent("numscal")));
              constlinperm.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("PERMEABILITIES")));
              constlinperm.emplace_back(Teuchos::rcp(new INPUT::RealVectorComponent(
                  "permeabilities", INPUT::LengthFromInt("numscal"))));

              constlinperm.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("IS_PSEUDO_CONTACT")));
              constlinperm.emplace_back(Teuchos::rcp(new INPUT::IntComponent("is_pseudo_contact")));

              kinetic_model_choices.emplace(INPAR::S2I::kinetics_constperm,
                  std::make_pair("ConstantPermeability", constlinperm));
              kinetic_model_choices.emplace(INPAR::S2I::kinetics_linearperm,
                  std::make_pair("LinearPermeability", constlinperm));
            }

            {
              // Butler-Volmer
              std::vector<Teuchos::RCP<INPUT::LineComponent>> butlervolmer;
              // total number of existing scalars
              butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("NUMSCAL")));
              butlervolmer.emplace_back(Teuchos::rcp(new INPUT::IntComponent("numscal")));
              butlervolmer.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("STOICHIOMETRIES")));
              butlervolmer.emplace_back(Teuchos::rcp(new INPUT::IntVectorComponent(
                  "stoichiometries", INPUT::LengthFromInt("numscal"))));
              butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("E-")));
              butlervolmer.emplace_back(Teuchos::rcp(new INPUT::IntComponent("e-")));
              butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("K_R")));
              butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("k_r")));
              butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_A")));
              butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("alpha_a")));
              butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_C")));
              butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("alpha_c")));
              butlervolmer.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("IS_PSEUDO_CONTACT")));
              butlervolmer.emplace_back(Teuchos::rcp(new INPUT::IntComponent("is_pseudo_contact")));

              // same components can be reused for multiple models
              kinetic_model_choices.emplace(
                  kinetics_butlervolmer, std::make_pair("Butler-Volmer", butlervolmer));
              kinetic_model_choices.emplace(kinetics_butlervolmerlinearized,
                  std::make_pair("Butler-Volmer_Linearized", butlervolmer));
              kinetic_model_choices.emplace(kinetics_butlervolmerreduced,
                  std::make_pair("Butler-VolmerReduced", butlervolmer));
              kinetic_model_choices.emplace(kinetics_butlervolmerreducedlinearized,
                  std::make_pair("Butler-VolmerReduced_Linearized", butlervolmer));
            }

            {
              // Butler-Volmer-Peltier
              std::vector<Teuchos::RCP<INPUT::LineComponent>> butlervolmerpeltier;

              butlervolmerpeltier.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("NUMSCAL")));
              butlervolmerpeltier.emplace_back(Teuchos::rcp(new INPUT::IntComponent("numscal")));
              butlervolmerpeltier.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("STOICHIOMETRIES")));
              butlervolmerpeltier.emplace_back(Teuchos::rcp(new INPUT::IntVectorComponent(
                  "stoichiometries", INPUT::LengthFromInt("numscal"))));
              butlervolmerpeltier.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("E-")));
              butlervolmerpeltier.emplace_back(Teuchos::rcp(new INPUT::IntComponent("e-")));
              butlervolmerpeltier.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("K_R")));
              butlervolmerpeltier.emplace_back(Teuchos::rcp(new INPUT::RealComponent("k_r")));
              butlervolmerpeltier.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_A")));
              butlervolmerpeltier.emplace_back(Teuchos::rcp(new INPUT::RealComponent("alpha_a")));
              butlervolmerpeltier.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_C")));
              butlervolmerpeltier.emplace_back(Teuchos::rcp(new INPUT::RealComponent("alpha_c")));
              butlervolmerpeltier.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("IS_PSEUDO_CONTACT")));
              butlervolmerpeltier.emplace_back(
                  Teuchos::rcp(new INPUT::IntComponent("is_pseudo_contact")));
              butlervolmerpeltier.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("PELTIER")));
              butlervolmerpeltier.emplace_back(Teuchos::rcp(new INPUT::RealComponent("peltier")));

              kinetic_model_choices.emplace(kinetics_butlervolmerpeltier,
                  std::make_pair("Butler-Volmer-Peltier", std::move(butlervolmerpeltier)));
            }

            {
              // Butler-Volmer-reduced with interface capacitance
              std::vector<Teuchos::RCP<INPUT::LineComponent>> butlervolmerreducedcapacitance;
              // total number of existing scalars
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("NUMSCAL")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new INPUT::IntComponent("numscal")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("STOICHIOMETRIES")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new INPUT::IntVectorComponent(
                      "stoichiometries", INPUT::LengthFromInt("numscal"))));

              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("E-")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new INPUT::IntComponent("e-")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("K_R")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("k_r")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("CAPACITANCE")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("capacitance")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_A")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("alpha_a")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_C")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("alpha_c")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("IS_PSEUDO_CONTACT")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new INPUT::IntComponent("is_pseudo_contact")));

              kinetic_model_choices.emplace(kinetics_butlervolmerreducedcapacitance,
                  std::make_pair("Butler-VolmerReduced_Capacitance",
                      std::move(butlervolmerreducedcapacitance)));
            }

            {
              // Butler-Volmer-Resistance
              std::vector<Teuchos::RCP<INPUT::LineComponent>> butlervolmerresistance;

              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("NUMSCAL")));
              butlervolmerresistance.emplace_back(Teuchos::rcp(new INPUT::IntComponent("numscal")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("STOICHIOMETRIES")));
              butlervolmerresistance.emplace_back(Teuchos::rcp(new INPUT::IntVectorComponent(
                  "stoichiometries", INPUT::LengthFromInt("numscal"))));

              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("E-")));
              butlervolmerresistance.emplace_back(Teuchos::rcp(new INPUT::IntComponent("e-")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("K_R")));
              butlervolmerresistance.emplace_back(Teuchos::rcp(new INPUT::RealComponent("k_r")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_A")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("alpha_a")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_C")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("alpha_c")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("IS_PSEUDO_CONTACT")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new INPUT::IntComponent("is_pseudo_contact")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("RESISTANCE")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("resistance")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("CONVTOL_IMPLBUTLERVOLMER")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("CONVTOL_IMPLBUTLERVOLMER")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("ITEMAX_IMPLBUTLERVOLMER")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new INPUT::IntComponent("ITEMAX_IMPLBUTLERVOLMER")));

              kinetic_model_choices.emplace(kinetics_butlervolmerresistance,
                  std::make_pair("Butler-Volmer_Resistance", std::move(butlervolmerresistance)));
            }

            {
              // Butler-Volmer-Reduced with resistance
              std::vector<Teuchos::RCP<INPUT::LineComponent>> butlervolmerreducedwithresistance;

              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("NUMSCAL")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::IntComponent("numscal")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("STOICHIOMETRIES")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::IntVectorComponent(
                      "stoichiometries", INPUT::LengthFromInt("numscal"))));

              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("E-")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::IntComponent("e-")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("K_R")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("k_r")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_A")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("alpha_a")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_C")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("alpha_c")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("IS_PSEUDO_CONTACT")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::IntComponent("is_pseudo_contact")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("RESISTANCE")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("resistance")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("CONVTOL_IMPLBUTLERVOLMER")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("CONVTOL_IMPLBUTLERVOLMER")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("ITEMAX_IMPLBUTLERVOLMER")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new INPUT::IntComponent("ITEMAX_IMPLBUTLERVOLMER")));

              kinetic_model_choices.emplace(kinetics_butlervolmerreducedresistance,
                  std::make_pair("Butler-VolmerReduced_Resistance",
                      std::move(butlervolmerreducedwithresistance)));
            }

            {
              // Butler-Volmer-reduced-thermoresistance
              std::vector<Teuchos::RCP<INPUT::LineComponent>> butlervolmerreducedthermo;
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("NUMSCAL")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new INPUT::IntComponent("numscal")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("STOICHIOMETRIES")));
              butlervolmerreducedthermo.emplace_back(Teuchos::rcp(new INPUT::IntVectorComponent(
                  "stoichiometries", INPUT::LengthFromInt("numscal"))));

              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("E-")));
              butlervolmerreducedthermo.emplace_back(Teuchos::rcp(new INPUT::IntComponent("e-")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("K_R")));
              butlervolmerreducedthermo.emplace_back(Teuchos::rcp(new INPUT::RealComponent("k_r")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_A")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("alpha_a")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_C")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("alpha_c")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("IS_PSEUDO_CONTACT")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new INPUT::IntComponent("is_pseudo_contact")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("THERMOPERM")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("thermoperm")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("MOLAR_HEAT_CAPACITY")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("molar_heat_capacity")));

              kinetic_model_choices.emplace(kinetics_butlervolmerreducedthermoresistance,
                  std::make_pair(
                      "Butler-VolmerReduced_ThermoResistance", butlervolmerreducedthermo));
            }

            {
              // constant interface resistance
              std::vector<Teuchos::RCP<INPUT::LineComponent>> constantinterfaceresistance;

              constantinterfaceresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("ONOFF")));
              constantinterfaceresistance.emplace_back(
                  Teuchos::rcp(new INPUT::IntVectorComponent("onoff", 2)));
              constantinterfaceresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("RESISTANCE")));
              constantinterfaceresistance.emplace_back(
                  Teuchos::rcp(new INPUT::RealComponent("resistance")));
              constantinterfaceresistance.emplace_back(new INPUT::SeparatorComponent("E-"));
              constantinterfaceresistance.emplace_back(new INPUT::IntComponent("e-"));
              constantinterfaceresistance.emplace_back(
                  Teuchos::rcp(new INPUT::SeparatorComponent("IS_PSEUDO_CONTACT")));
              constantinterfaceresistance.emplace_back(
                  Teuchos::rcp(new INPUT::IntComponent("is_pseudo_contact")));

              kinetic_model_choices.emplace(kinetics_constantinterfaceresistance,
                  std::make_pair(
                      "ConstantInterfaceResistance", std::move(constantinterfaceresistance)));
            }

            {
              // no interface flux
              std::vector<Teuchos::RCP<INPUT::LineComponent>> nointerfaceflux;
              kinetic_model_choices.emplace(
                  kinetics_nointerfaceflux, std::make_pair("NoInterfaceFlux", nointerfaceflux));
            }
          }  // kinetic models for scatra-scatra interface kinetics

          // insert kinetic models into vector with slave-side condition components
          slaveside.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("KINETIC_MODEL")));
          slaveside.emplace_back(Teuchos::rcp(new INPUT::SwitchComponent(
              "kinetic model", kinetics_butlervolmer, kinetic_model_choices)));

          // add all components from slave side to multi-scale condition
          for (const auto& component : slaveside) multiscalecouplingpoint->AddComponent(component);

          // insert slave-side condition components into vector of interface sides
          interface_choices.emplace(side_slave, std::make_pair("Slave", slaveside));
        }

        {
          // master side
          std::vector<Teuchos::RCP<INPUT::LineComponent>> master_side;
          interface_choices.emplace(side_master, std::make_pair("Master", master_side));
        }
      }  // interface sides for scatra-scatra interface mesh tying

      // insert interface sides into vector with input file line components
      s2icomponents.emplace_back(Teuchos::rcp(
          new INPUT::SwitchComponent("interface side", side_undefined, interface_choices)));
    }

    // insert input file line components into condition definitions
    for (auto& s2icomponent : s2icomponents)
    {
      s2ikineticspoint->AddComponent(s2icomponent);
      s2ikineticsline->AddComponent(s2icomponent);
      s2ikineticssurf->AddComponent(s2icomponent);
    }

    // insert condition definitions into global list of valid condition definitions
    condlist.emplace_back(s2ikineticspoint);
    condlist.emplace_back(s2ikineticsline);
    condlist.emplace_back(s2ikineticssurf);

    condlist.emplace_back(multiscalecouplingpoint);
  }


  /*--------------------------------------------------------------------*/
  // scatra-scatra interface coupling involving interface layer growth
  {
    // definition of scatra-scatra interface coupling line condition involving interface layer
    // growth
    auto s2igrowthline = Teuchos::rcp(
        new CORE::Conditions::ConditionDefinition("DESIGN S2I KINETICS GROWTH LINE CONDITIONS",
            "S2IKineticsGrowth", "Scatra-scatra line interface layer growth kinetics",
            CORE::Conditions::S2IKineticsGrowth, true, CORE::Conditions::geometry_type_line));

    // definition of scatra-scatra interface coupling surface condition involving interface layer
    // growth
    auto s2igrowthsurf = Teuchos::rcp(
        new CORE::Conditions::ConditionDefinition("DESIGN S2I KINETICS GROWTH SURF CONDITIONS",
            "S2IKineticsGrowth", "Scatra-scatra surface interface layer growth kinetics",
            CORE::Conditions::S2IKineticsGrowth, true, CORE::Conditions::geometry_type_surface));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<INPUT::LineComponent>> s2igrowthcomponents;
    {
      // interface ID
      s2igrowthcomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent("ConditionID")));

      // Butler-Volmer
      std::vector<Teuchos::RCP<INPUT::LineComponent>> butlervolmer;
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("NUMSCAL")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::IntComponent("numscal")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("STOICHIOMETRIES")));
      butlervolmer.emplace_back(Teuchos::rcp(
          new INPUT::IntVectorComponent("stoichiometries", INPUT::LengthFromInt("numscal"))));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("E-")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::IntComponent("e-")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("K_R")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("k_r")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_A")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("alpha_a")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_C")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("alpha_c")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("MOLMASS")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("molar mass")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("DENSITY")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("density")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("CONDUCTIVITY")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("conductivity")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("REGTYPE")));
      butlervolmer.emplace_back(
          Teuchos::rcp(new INPUT::SelectionComponent("regularization type", "trigonometrical",
              Teuchos::tuple<std::string>("none", "polynomial", "Hein", "trigonometrical"),
              Teuchos::tuple<int>(INPAR::S2I::regularization_none,
                  INPAR::S2I::regularization_polynomial, INPAR::S2I::regularization_hein,
                  INPAR::S2I::regularization_trigonometrical))));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("REGPAR")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("regularization parameter")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("INITTHICKNESS")));
      butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("initial thickness")));

      // insert kinetic models into vector with input file line components
      s2igrowthcomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("KINETIC_MODEL")));
      s2igrowthcomponents.emplace_back(
          Teuchos::rcp(new INPUT::SwitchComponent("kinetic model", growth_kinetics_butlervolmer,
              {{growth_kinetics_butlervolmer,
                  std::make_pair("Butler-Volmer", std::move(butlervolmer))}})));
    }

    // insert input file line components into condition definitions
    for (auto& s2igrowthcomponent : s2igrowthcomponents)
    {
      s2igrowthline->AddComponent(s2igrowthcomponent);
      s2igrowthsurf->AddComponent(s2igrowthcomponent);
    }

    // insert condition definitions into global list of valid condition definitions
    condlist.emplace_back(s2igrowthline);
    condlist.emplace_back(s2igrowthsurf);
  }

  /*--------------------------------------------------------------------*/
  // scatra-scatra interface with micro-macro coupling for space-charge layers
  {
    auto s2isclcond = Teuchos::rcp(
        new CORE::Conditions::ConditionDefinition("DESIGN S2I SCL COUPLING SURF CONDITIONS",
            "S2ISCLCoupling", "Scatra-scatra surface with SCL micro-macro coupling between",
            CORE::Conditions::S2ISCLCoupling, true, CORE::Conditions::geometry_type_surface));

    s2isclcond->AddComponent(Teuchos::rcp(new INPUT::SelectionComponent("interface side",
        "Undefined", Teuchos::tuple<std::string>("Undefined", "Slave", "Master"),
        Teuchos::tuple<int>(
            INPAR::S2I::side_undefined, INPAR::S2I::side_slave, INPAR::S2I::side_master))));

    condlist.emplace_back(s2isclcond);
  }
}

FOUR_C_NAMESPACE_CLOSE
