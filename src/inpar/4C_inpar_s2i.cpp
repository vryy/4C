#include "4C_inpar_s2i.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_io_linecomponent.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*------------------------------------------------------------------------*
 | set valid parameters for scatra-scatra interface coupling   fang 01/16 |
 *------------------------------------------------------------------------*/
void Inpar::S2I::set_valid_parameters(Teuchos::ParameterList& list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& s2icoupling =
      list.sublist("SCALAR TRANSPORT DYNAMIC", true)
          .sublist(
              "S2I COUPLING", false, "control parameters for scatra-scatra interface coupling");

  // type of mortar meshtying
  setStringToIntegralParameter<CouplingType>("COUPLINGTYPE", "Undefined",
      "type of mortar meshtying",
      tuple<std::string>("Undefined", "MatchingNodes", "StandardMortar", "SaddlePointMortar_Petrov",
          "SaddlePointMortar_Bubnov", "CondensedMortar_Petrov", "CondensedMortar_Bubnov",
          "StandardNodeToSegment"),
      tuple<CouplingType>(coupling_undefined, coupling_matching_nodes, coupling_mortar_standard,
          coupling_mortar_saddlepoint_petrov, coupling_mortar_saddlepoint_bubnov,
          coupling_mortar_condensed_petrov, coupling_mortar_condensed_bubnov,
          coupling_nts_standard),
      &s2icoupling);

  // flag for interface side underlying Lagrange multiplier definition
  setStringToIntegralParameter<InterfaceSides>("LMSIDE", "slave",
      "flag for interface side underlying Lagrange multiplier definition",
      tuple<std::string>("slave", "master"), tuple<InterfaceSides>(side_slave, side_master),
      &s2icoupling);

  // flag for evaluation of interface linearizations and residuals on slave side only
  Core::Utils::bool_parameter("SLAVEONLY", "No",
      "flag for evaluation of interface linearizations and residuals on slave side only",
      &s2icoupling);

  // node-to-segment projection tolerance
  Core::Utils::double_parameter(
      "NTSPROJTOL", 0.0, "node-to-segment projection tolerance", &s2icoupling);

  // flag for evaluation of scatra-scatra interface coupling involving interface layer growth
  setStringToIntegralParameter<GrowthEvaluation>("INTLAYERGROWTH_EVALUATION", "none",
      "flag for evaluation of scatra-scatra interface coupling involving interface layer growth",
      tuple<std::string>("none", "monolithic", "semi-implicit"),
      tuple<GrowthEvaluation>(
          growth_evaluation_none, growth_evaluation_monolithic, growth_evaluation_semi_implicit),
      &s2icoupling);

  // local Newton-Raphson convergence tolerance for scatra-scatra interface coupling involving
  // interface layer growth
  Core::Utils::double_parameter("INTLAYERGROWTH_CONVTOL", 1.e-12,
      "local Newton-Raphson convergence tolerance for scatra-scatra interface coupling involving "
      "interface layer growth",
      &s2icoupling);

  // maximum number of local Newton-Raphson iterations for scatra-scatra interface coupling
  // involving interface layer growth
  Core::Utils::int_parameter("INTLAYERGROWTH_ITEMAX", 5,
      "maximum number of local Newton-Raphson iterations for scatra-scatra interface coupling "
      "involving interface layer growth",
      &s2icoupling);

  // ID of linear solver for monolithic scatra-scatra interface coupling involving interface layer
  // growth
  Core::Utils::int_parameter("INTLAYERGROWTH_LINEAR_SOLVER", -1,
      "ID of linear solver for monolithic scatra-scatra interface coupling involving interface "
      "layer growth",
      &s2icoupling);

  // modified time step size for scatra-scatra interface coupling involving interface layer growth
  Core::Utils::double_parameter("INTLAYERGROWTH_TIMESTEP", -1.,
      "modified time step size for scatra-scatra interface coupling involving interface layer "
      "growth",
      &s2icoupling);

  Core::Utils::bool_parameter("MESHTYING_CONDITIONS_INDEPENDENT_SETUP", "No",
      "mesh tying for different conditions should be setup independently", &s2icoupling);

  Core::Utils::bool_parameter("OUTPUT_INTERFACE_FLUX", "No",
      "evaluate integral of coupling flux on slave side for each s2i condition and write it to csv "
      "file",
      &s2icoupling);
}


/*------------------------------------------------------------------------*
 | set valid conditions for scatra-scatra interface coupling   fang 01/16 |
 *------------------------------------------------------------------------*/
void Inpar::S2I::set_valid_conditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  /*--------------------------------------------------------------------*/
  // scatra-scatra interface mesh tying condition
  {
    // definition of scatra-scatra interface mesh tying line condition
    auto s2imeshtyingline = Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
        "DESIGN S2I MESHTYING LINE CONDITIONS", "S2IMeshtying",
        "Scatra-scatra line interface mesh tying", Core::Conditions::S2IMeshtying, true,
        Core::Conditions::geometry_type_line);

    // definition of scatra-scatra interface mesh tying surface condition
    auto s2imeshtyingsurf = Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
        "DESIGN S2I MESHTYING SURF CONDITIONS", "S2IMeshtying",
        "Scatra-scatra surface interface mesh tying", Core::Conditions::S2IMeshtying, true,
        Core::Conditions::geometry_type_surface);

    // insert input file line components into condition definitions
    for (const auto& cond : {s2imeshtyingline, s2imeshtyingsurf})
    {
      cond->add_component(Teuchos::make_rcp<Input::IntComponent>("ConditionID"));
      cond->add_component(Teuchos::make_rcp<Input::SelectionComponent>("interface side",
          "Undefined", Teuchos::tuple<std::string>("Undefined", "Slave", "Master"),
          Teuchos::tuple<int>(
              Inpar::S2I::side_undefined, Inpar::S2I::side_slave, Inpar::S2I::side_master)));
      add_named_int(cond, "S2I_KINETICS_ID");
      condlist.push_back(cond);
    }
  }

  /*--------------------------------------------------------------------*/
  // scatra-scatra interface no evaluation condition
  {
    // definition of scatra-scatra interface no evaluation line condition
    auto s2inoevaluationline = Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
        "DESIGN S2I NO EVALUATION LINE CONDITIONS", "S2INoEvaluation",
        "Scatra-scatra interface no evaluation line condition. This condition can be used to "
        "deactivate the evaluation of the corresponding `S2IKinetics` condition. Another usage "
        "is in coupled algorithms, where a specific `S2IKinetics` condition should not be "
        "evaluated within the scalar transport framework, as it already evaluated elsewhere. "
        "One example is the SSI contact.",
        Core::Conditions::S2INoEvaluation, true, Core::Conditions::geometry_type_line);

    // definition of scatra-scatra interface no evaluation surface condition
    auto s2inoevaluationsurf = Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
        "DESIGN S2I NO EVALUATION SURF CONDITIONS", "S2INoEvaluation",
        "Scatra-scatra interface no evaluation surface condition. This condition can be used to "
        "deactivate the evaluation of the corresponding `S2IKinetics` condition. Another usage "
        "is in coupled algorithms, where a specific `S2IKinetics` condition should not be "
        "evaluated within the scalar transport framework, as it already evaluated elsewhere. "
        "One example is the SSI contact.",
        Core::Conditions::S2INoEvaluation, true, Core::Conditions::geometry_type_surface);

    // insert input file line components into condition definitions
    for (const auto& cond : {s2inoevaluationline, s2inoevaluationsurf})
    {
      cond->add_component(Teuchos::make_rcp<Input::IntComponent>("ConditionID"));
      cond->add_component(Teuchos::make_rcp<Input::SelectionComponent>("interface side",
          "Undefined", Teuchos::tuple<std::string>("Undefined", "Slave", "Master"),
          Teuchos::tuple<int>(
              Inpar::S2I::side_undefined, Inpar::S2I::side_slave, Inpar::S2I::side_master)));
      add_named_int(cond, "S2I_KINETICS_ID");
      condlist.push_back(cond);
    }
  }

  /*--------------------------------------------------------------------*/
  // scatra-scatra interface kinetics condition
  {
    // definition of scatra-scatra interface kinetics point condition
    auto s2ikineticspoint = Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
        "DESIGN S2I KINETICS POINT CONDITIONS", "S2IKinetics",
        "Scatra-scatra line interface kinetics", Core::Conditions::S2IKinetics, true,
        Core::Conditions::geometry_type_point);

    // definition of scatra-scatra interface kinetics line condition
    auto s2ikineticsline = Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
        "DESIGN S2I KINETICS LINE CONDITIONS", "S2IKinetics",
        "Scatra-scatra line interface kinetics", Core::Conditions::S2IKinetics, true,
        Core::Conditions::geometry_type_line);

    // definition of scatra-scatra interface kinetics surface condition
    auto s2ikineticssurf = Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
        "DESIGN S2I KINETICS SURF CONDITIONS", "S2IKinetics",
        "Scatra-scatra surface interface kinetics", Core::Conditions::S2IKinetics, true,
        Core::Conditions::geometry_type_surface);

    // Macro-micro coupling condition for micro scale in multi-scale scalar transport problems
    auto multiscalecouplingpoint = Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
        "DESIGN SCATRA MULTI-SCALE COUPLING POINT CONDITIONS", "ScatraMultiScaleCoupling",
        "Scalar transport multi-scale coupling condition",
        Core::Conditions::ScatraMultiScaleCoupling, false, Core::Conditions::geometry_type_point);

    // prepare interface sides for scatra-scatra interface kinetics
    std::map<int, std::pair<std::string, std::vector<Teuchos::RCP<Input::LineComponent>>>>
        interface_choices;
    {
      {
        // undefined side
        std::vector<Teuchos::RCP<Input::LineComponent>> undefined_side;
        interface_choices.emplace(side_undefined, std::make_pair("Undefined", undefined_side));
      }

      {
        // slave side
        std::vector<Teuchos::RCP<Input::LineComponent>> slaveside;

        // Collect the diffent model selection choices in a map.
        std::map<int, std::pair<std::string, std::vector<Teuchos::RCP<Input::LineComponent>>>>
            kinetic_model_choices;
        {
          {
            // constant and linear permeability
            std::vector<Teuchos::RCP<Input::LineComponent>> constlinperm;

            constlinperm.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("NUMSCAL"));
            constlinperm.emplace_back(Teuchos::make_rcp<Input::IntComponent>("NUMSCAL"));
            constlinperm.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("PERMEABILITIES"));
            constlinperm.emplace_back(Teuchos::make_rcp<Input::RealVectorComponent>(
                "PERMEABILITIES", Input::LengthFromInt("NUMSCAL")));

            constlinperm.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("IS_PSEUDO_CONTACT"));
            constlinperm.emplace_back(Teuchos::make_rcp<Input::IntComponent>("IS_PSEUDO_CONTACT"));

            kinetic_model_choices.emplace(Inpar::S2I::kinetics_constperm,
                std::make_pair("ConstantPermeability", constlinperm));
            kinetic_model_choices.emplace(Inpar::S2I::kinetics_linearperm,
                std::make_pair("LinearPermeability", constlinperm));
          }

          {
            // Butler-Volmer
            std::vector<Teuchos::RCP<Input::LineComponent>> butlervolmer;
            // total number of existing scalars
            butlervolmer.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("NUMSCAL"));
            butlervolmer.emplace_back(Teuchos::make_rcp<Input::IntComponent>("NUMSCAL"));
            butlervolmer.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("STOICHIOMETRIES"));
            butlervolmer.emplace_back(Teuchos::make_rcp<Input::IntVectorComponent>(
                "STOICHIOMETRIES", Input::LengthFromInt("NUMSCAL")));
            butlervolmer.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("E-"));
            butlervolmer.emplace_back(Teuchos::make_rcp<Input::IntComponent>("E-"));
            butlervolmer.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("K_R"));
            butlervolmer.emplace_back(Teuchos::make_rcp<Input::RealComponent>("K_R"));
            butlervolmer.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("ALPHA_A"));
            butlervolmer.emplace_back(Teuchos::make_rcp<Input::RealComponent>("ALPHA_A"));
            butlervolmer.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("ALPHA_C"));
            butlervolmer.emplace_back(Teuchos::make_rcp<Input::RealComponent>("ALPHA_C"));
            butlervolmer.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("IS_PSEUDO_CONTACT"));
            butlervolmer.emplace_back(Teuchos::make_rcp<Input::IntComponent>("IS_PSEUDO_CONTACT"));

            // same components can be reused for multiple models
            kinetic_model_choices.emplace(
                kinetics_butlervolmer, std::make_pair("Butler-Volmer", butlervolmer));
            kinetic_model_choices.emplace(kinetics_butlervolmerlinearized,
                std::make_pair("Butler-Volmer_Linearized", butlervolmer));
            kinetic_model_choices.emplace(
                kinetics_butlervolmerreduced, std::make_pair("Butler-VolmerReduced", butlervolmer));
            kinetic_model_choices.emplace(kinetics_butlervolmerreducedlinearized,
                std::make_pair("Butler-VolmerReduced_Linearized", butlervolmer));
          }

          {
            // Butler-Volmer-Peltier
            std::vector<Teuchos::RCP<Input::LineComponent>> butlervolmerpeltier;

            butlervolmerpeltier.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("NUMSCAL"));
            butlervolmerpeltier.emplace_back(Teuchos::make_rcp<Input::IntComponent>("NUMSCAL"));
            butlervolmerpeltier.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("STOICHIOMETRIES"));
            butlervolmerpeltier.emplace_back(Teuchos::make_rcp<Input::IntVectorComponent>(
                "STOICHIOMETRIES", Input::LengthFromInt("NUMSCAL")));
            butlervolmerpeltier.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("E-"));
            butlervolmerpeltier.emplace_back(Teuchos::make_rcp<Input::IntComponent>("E-"));
            butlervolmerpeltier.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("K_R"));
            butlervolmerpeltier.emplace_back(Teuchos::make_rcp<Input::RealComponent>("K_R"));
            butlervolmerpeltier.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("ALPHA_A"));
            butlervolmerpeltier.emplace_back(Teuchos::make_rcp<Input::RealComponent>("ALPHA_A"));
            butlervolmerpeltier.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("ALPHA_C"));
            butlervolmerpeltier.emplace_back(Teuchos::make_rcp<Input::RealComponent>("ALPHA_C"));
            butlervolmerpeltier.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("IS_PSEUDO_CONTACT"));
            butlervolmerpeltier.emplace_back(
                Teuchos::make_rcp<Input::IntComponent>("IS_PSEUDO_CONTACT"));
            butlervolmerpeltier.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("PELTIER"));
            butlervolmerpeltier.emplace_back(Teuchos::make_rcp<Input::RealComponent>("PELTIER"));

            kinetic_model_choices.emplace(kinetics_butlervolmerpeltier,
                std::make_pair("Butler-Volmer-Peltier", std::move(butlervolmerpeltier)));
          }

          {
            // Butler-Volmer-reduced with interface capacitance
            std::vector<Teuchos::RCP<Input::LineComponent>> butlervolmerreducedcapacitance;
            // total number of existing scalars
            butlervolmerreducedcapacitance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("NUMSCAL"));
            butlervolmerreducedcapacitance.emplace_back(
                Teuchos::make_rcp<Input::IntComponent>("NUMSCAL"));
            butlervolmerreducedcapacitance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("STOICHIOMETRIES"));
            butlervolmerreducedcapacitance.emplace_back(
                Teuchos::make_rcp<Input::IntVectorComponent>(
                    "STOICHIOMETRIES", Input::LengthFromInt("NUMSCAL")));

            butlervolmerreducedcapacitance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("E-"));
            butlervolmerreducedcapacitance.emplace_back(
                Teuchos::make_rcp<Input::IntComponent>("E-"));
            butlervolmerreducedcapacitance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("K_R"));
            butlervolmerreducedcapacitance.emplace_back(
                Teuchos::make_rcp<Input::RealComponent>("K_R"));
            butlervolmerreducedcapacitance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("CAPACITANCE"));
            butlervolmerreducedcapacitance.emplace_back(
                Teuchos::make_rcp<Input::RealComponent>("CAPACITANCE"));
            butlervolmerreducedcapacitance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("ALPHA_A"));
            butlervolmerreducedcapacitance.emplace_back(
                Teuchos::make_rcp<Input::RealComponent>("ALPHA_A"));
            butlervolmerreducedcapacitance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("ALPHA_C"));
            butlervolmerreducedcapacitance.emplace_back(
                Teuchos::make_rcp<Input::RealComponent>("ALPHA_C"));
            butlervolmerreducedcapacitance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("IS_PSEUDO_CONTACT"));
            butlervolmerreducedcapacitance.emplace_back(
                Teuchos::make_rcp<Input::IntComponent>("IS_PSEUDO_CONTACT"));

            kinetic_model_choices.emplace(kinetics_butlervolmerreducedcapacitance,
                std::make_pair(
                    "Butler-VolmerReduced_Capacitance", std::move(butlervolmerreducedcapacitance)));
          }

          {
            // Butler-Volmer-Resistance
            std::vector<Teuchos::RCP<Input::LineComponent>> butlervolmerresistance;

            butlervolmerresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("NUMSCAL"));
            butlervolmerresistance.emplace_back(Teuchos::make_rcp<Input::IntComponent>("NUMSCAL"));
            butlervolmerresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("STOICHIOMETRIES"));
            butlervolmerresistance.emplace_back(Teuchos::make_rcp<Input::IntVectorComponent>(
                "STOICHIOMETRIES", Input::LengthFromInt("NUMSCAL")));

            butlervolmerresistance.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("E-"));
            butlervolmerresistance.emplace_back(Teuchos::make_rcp<Input::IntComponent>("E-"));
            butlervolmerresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("K_R"));
            butlervolmerresistance.emplace_back(Teuchos::make_rcp<Input::RealComponent>("K_R"));
            butlervolmerresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("ALPHA_A"));
            butlervolmerresistance.emplace_back(Teuchos::make_rcp<Input::RealComponent>("ALPHA_A"));
            butlervolmerresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("ALPHA_C"));
            butlervolmerresistance.emplace_back(Teuchos::make_rcp<Input::RealComponent>("ALPHA_C"));
            butlervolmerresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("IS_PSEUDO_CONTACT"));
            butlervolmerresistance.emplace_back(
                Teuchos::make_rcp<Input::IntComponent>("IS_PSEUDO_CONTACT"));
            butlervolmerresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("RESISTANCE"));
            butlervolmerresistance.emplace_back(
                Teuchos::make_rcp<Input::RealComponent>("RESISTANCE"));
            butlervolmerresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("CONVTOL_IMPLBUTLERVOLMER"));
            butlervolmerresistance.emplace_back(
                Teuchos::make_rcp<Input::RealComponent>("CONVTOL_IMPLBUTLERVOLMER"));
            butlervolmerresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("ITEMAX_IMPLBUTLERVOLMER"));
            butlervolmerresistance.emplace_back(
                Teuchos::make_rcp<Input::IntComponent>("ITEMAX_IMPLBUTLERVOLMER"));

            kinetic_model_choices.emplace(kinetics_butlervolmerresistance,
                std::make_pair("Butler-Volmer_Resistance", std::move(butlervolmerresistance)));
          }

          {
            // Butler-Volmer-Reduced with resistance
            std::vector<Teuchos::RCP<Input::LineComponent>> butlervolmerreducedwithresistance;

            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("NUMSCAL"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::IntComponent>("NUMSCAL"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("STOICHIOMETRIES"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::IntVectorComponent>(
                    "STOICHIOMETRIES", Input::LengthFromInt("NUMSCAL")));

            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("E-"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::IntComponent>("E-"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("K_R"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::RealComponent>("K_R"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("ALPHA_A"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::RealComponent>("ALPHA_A"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("ALPHA_C"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::RealComponent>("ALPHA_C"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("IS_PSEUDO_CONTACT"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::IntComponent>("IS_PSEUDO_CONTACT"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("RESISTANCE"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::RealComponent>("RESISTANCE"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("CONVTOL_IMPLBUTLERVOLMER"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::RealComponent>("CONVTOL_IMPLBUTLERVOLMER"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("ITEMAX_IMPLBUTLERVOLMER"));
            butlervolmerreducedwithresistance.emplace_back(
                Teuchos::make_rcp<Input::IntComponent>("ITEMAX_IMPLBUTLERVOLMER"));

            kinetic_model_choices.emplace(kinetics_butlervolmerreducedresistance,
                std::make_pair("Butler-VolmerReduced_Resistance",
                    std::move(butlervolmerreducedwithresistance)));
          }

          {
            // Butler-Volmer-reduced-thermoresistance
            std::vector<Teuchos::RCP<Input::LineComponent>> butlervolmerreducedthermo;
            butlervolmerreducedthermo.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("NUMSCAL"));
            butlervolmerreducedthermo.emplace_back(
                Teuchos::make_rcp<Input::IntComponent>("NUMSCAL"));
            butlervolmerreducedthermo.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("STOICHIOMETRIES"));
            butlervolmerreducedthermo.emplace_back(Teuchos::make_rcp<Input::IntVectorComponent>(
                "STOICHIOMETRIES", Input::LengthFromInt("NUMSCAL")));

            butlervolmerreducedthermo.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("E-"));
            butlervolmerreducedthermo.emplace_back(Teuchos::make_rcp<Input::IntComponent>("E-"));
            butlervolmerreducedthermo.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("K_R"));
            butlervolmerreducedthermo.emplace_back(Teuchos::make_rcp<Input::RealComponent>("K_R"));
            butlervolmerreducedthermo.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("ALPHA_A"));
            butlervolmerreducedthermo.emplace_back(
                Teuchos::make_rcp<Input::RealComponent>("ALPHA_A"));
            butlervolmerreducedthermo.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("ALPHA_C"));
            butlervolmerreducedthermo.emplace_back(
                Teuchos::make_rcp<Input::RealComponent>("ALPHA_C"));
            butlervolmerreducedthermo.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("IS_PSEUDO_CONTACT"));
            butlervolmerreducedthermo.emplace_back(
                Teuchos::make_rcp<Input::IntComponent>("IS_PSEUDO_CONTACT"));
            butlervolmerreducedthermo.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("THERMOPERM"));
            butlervolmerreducedthermo.emplace_back(
                Teuchos::make_rcp<Input::RealComponent>("THERMOPERM"));
            butlervolmerreducedthermo.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("MOLAR_HEAT_CAPACITY"));
            butlervolmerreducedthermo.emplace_back(
                Teuchos::make_rcp<Input::RealComponent>("MOLAR_HEAT_CAPACITY"));

            kinetic_model_choices.emplace(kinetics_butlervolmerreducedthermoresistance,
                std::make_pair("Butler-VolmerReduced_ThermoResistance", butlervolmerreducedthermo));
          }

          {
            // constant interface resistance
            std::vector<Teuchos::RCP<Input::LineComponent>> constantinterfaceresistance;

            constantinterfaceresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("ONOFF"));
            constantinterfaceresistance.emplace_back(
                Teuchos::make_rcp<Input::IntVectorComponent>("ONOFF", 2));
            constantinterfaceresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("RESISTANCE"));
            constantinterfaceresistance.emplace_back(
                Teuchos::make_rcp<Input::RealComponent>("RESISTANCE"));
            constantinterfaceresistance.emplace_back(new Input::SeparatorComponent("E-"));
            constantinterfaceresistance.emplace_back(new Input::IntComponent("E-"));
            constantinterfaceresistance.emplace_back(
                Teuchos::make_rcp<Input::SeparatorComponent>("IS_PSEUDO_CONTACT"));
            constantinterfaceresistance.emplace_back(
                Teuchos::make_rcp<Input::IntComponent>("IS_PSEUDO_CONTACT"));

            kinetic_model_choices.emplace(
                kinetics_constantinterfaceresistance, std::make_pair("ConstantInterfaceResistance",
                                                          std::move(constantinterfaceresistance)));
          }

          {
            // no interface flux
            std::vector<Teuchos::RCP<Input::LineComponent>> nointerfaceflux;
            kinetic_model_choices.emplace(
                kinetics_nointerfaceflux, std::make_pair("NoInterfaceFlux", nointerfaceflux));
          }
        }  // kinetic models for scatra-scatra interface kinetics

        // insert kinetic models into vector with slave-side condition components
        slaveside.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("KINETIC_MODEL"));
        slaveside.emplace_back(Teuchos::make_rcp<Input::SwitchComponent>(
            "KINETIC_MODEL", kinetics_butlervolmer, kinetic_model_choices));

        // add all components from slave side to multi-scale condition
        for (const auto& component : slaveside) multiscalecouplingpoint->add_component(component);

        // insert slave-side condition components into vector of interface sides
        interface_choices.emplace(side_slave, std::make_pair("Slave", slaveside));
      }

      {
        // master side
        std::vector<Teuchos::RCP<Input::LineComponent>> master_side;
        interface_choices.emplace(side_master, std::make_pair("Master", master_side));
      }
    }  // interface sides for scatra-scatra interface mesh tying

    // add components to conditions
    for (const auto& cond : {s2ikineticspoint, s2ikineticsline, s2ikineticssurf})
    {
      // interface ID
      cond->add_component(Teuchos::make_rcp<Input::IntComponent>("ConditionID"));

      // insert interface sides as line components
      cond->add_component(Teuchos::make_rcp<Input::SwitchComponent>(
          "interface side", side_undefined, interface_choices));

      // insert condition definitions into global list of valid condition definitions
      condlist.emplace_back(cond);
    }

    condlist.emplace_back(multiscalecouplingpoint);
  }



  /*--------------------------------------------------------------------*/
  // scatra-scatra interface coupling involving interface layer growth
  {
    // definition of scatra-scatra interface coupling line condition involving interface layer
    // growth
    auto s2igrowthline = Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
        "DESIGN S2I KINETICS GROWTH LINE CONDITIONS", "S2IKineticsGrowth",
        "Scatra-scatra line interface layer growth kinetics", Core::Conditions::S2IKineticsGrowth,
        true, Core::Conditions::geometry_type_line);

    // definition of scatra-scatra interface coupling surface condition involving interface layer
    // growth
    auto s2igrowthsurf = Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
        "DESIGN S2I KINETICS GROWTH SURF CONDITIONS", "S2IKineticsGrowth",
        "Scatra-scatra surface interface layer growth kinetics",
        Core::Conditions::S2IKineticsGrowth, true, Core::Conditions::geometry_type_surface);

    // Prepare components of Butler-Volmer condition
    std::vector<Teuchos::RCP<Input::LineComponent>> butlervolmer;
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("NUMSCAL"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::IntComponent>("NUMSCAL"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("STOICHIOMETRIES"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::IntVectorComponent>(
        "STOICHIOMETRIES", Input::LengthFromInt("NUMSCAL")));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("E-"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::IntComponent>("E-"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("K_R"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::RealComponent>("K_R"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("ALPHA_A"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::RealComponent>("ALPHA_A"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("ALPHA_C"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::RealComponent>("ALPHA_C"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("MOLMASS"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::RealComponent>("MOLMASS"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("DENSITY"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::RealComponent>("density"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("CONDUCTIVITY"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::RealComponent>("CONDUCTIVITY"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("REGTYPE"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::SelectionComponent>("REGTYPE",
        "trigonometrical",
        Teuchos::tuple<std::string>("none", "polynomial", "Hein", "trigonometrical"),
        Teuchos::tuple<int>(Inpar::S2I::regularization_none, Inpar::S2I::regularization_polynomial,
            Inpar::S2I::regularization_hein, Inpar::S2I::regularization_trigonometrical)));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("REGPAR"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::RealComponent>("REGPAR"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::SeparatorComponent>("INITTHICKNESS"));
    butlervolmer.emplace_back(Teuchos::make_rcp<Input::RealComponent>("INITTHICKNESS"));


    for (const auto& cond : {s2igrowthline, s2igrowthsurf})
    {
      // interface ID
      cond->add_component(Teuchos::make_rcp<Input::IntComponent>("ConditionID"));

      // add kinetic models as input file line components
      cond->add_component(Teuchos::make_rcp<Input::SeparatorComponent>("KINETIC_MODEL"));
      cond->add_component(
          Teuchos::RCP(new Input::SwitchComponent("KINETIC_MODEL", growth_kinetics_butlervolmer,
              {{growth_kinetics_butlervolmer, std::make_pair("Butler-Volmer", butlervolmer)}})));

      // insert condition definitions into global list of valid condition definitions
      condlist.emplace_back(cond);
    }
  }

  /*--------------------------------------------------------------------*/
  // scatra-scatra interface with micro-macro coupling for space-charge layers
  {
    auto s2isclcond = Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
        "DESIGN S2I SCL COUPLING SURF CONDITIONS", "S2ISCLCoupling",
        "Scatra-scatra surface with SCL micro-macro coupling between",
        Core::Conditions::S2ISCLCoupling, true, Core::Conditions::geometry_type_surface);

    s2isclcond->add_component(Teuchos::make_rcp<Input::SelectionComponent>("interface side",
        "Undefined", Teuchos::tuple<std::string>("Undefined", "Slave", "Master"),
        Teuchos::tuple<int>(
            Inpar::S2I::side_undefined, Inpar::S2I::side_slave, Inpar::S2I::side_master)));

    condlist.emplace_back(s2isclcond);
  }
}

FOUR_C_NAMESPACE_CLOSE
