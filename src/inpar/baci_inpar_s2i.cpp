/*----------------------------------------------------------------------*/
/*! \file
\brief input quantities and globally accessible enumerations for scatra-scatra interface coupling

\level 2


*/
/*----------------------------------------------------------------------*/
#include "baci_inpar_s2i.H"

#include "baci_inpar_validparameters.H"
#include "baci_lib_conditiondefinition.H"

/*------------------------------------------------------------------------*
 | set valid parameters for scatra-scatra interface coupling   fang 01/16 |
 *------------------------------------------------------------------------*/
void INPAR::S2I::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
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
  BoolParameter("SLAVEONLY", "No",
      "flag for evaluation of interface linearizations and residuals on slave side only",
      &s2icoupling);

  // node-to-segment projection tolerance
  DoubleParameter("NTSPROJTOL", 0.0, "node-to-segment projection tolerance", &s2icoupling);

  // flag for evaluation of scatra-scatra interface coupling involving interface layer growth
  setStringToIntegralParameter<int>("INTLAYERGROWTH_EVALUATION", "none",
      "flag for evaluation of scatra-scatra interface coupling involving interface layer growth",
      tuple<std::string>("none", "monolithic", "semi-implicit"),
      tuple<int>(
          growth_evaluation_none, growth_evaluation_monolithic, growth_evaluation_semi_implicit),
      &s2icoupling);

  // local Newton-Raphson convergence tolerance for scatra-scatra interface coupling involving
  // interface layer growth
  DoubleParameter("INTLAYERGROWTH_CONVTOL", 1.e-12,
      "local Newton-Raphson convergence tolerance for scatra-scatra interface coupling involving "
      "interface layer growth",
      &s2icoupling);

  // maximum number of local Newton-Raphson iterations for scatra-scatra interface coupling
  // involving interface layer growth
  IntParameter("INTLAYERGROWTH_ITEMAX", 5,
      "maximum number of local Newton-Raphson iterations for scatra-scatra interface coupling "
      "involving interface layer growth",
      &s2icoupling);

  // ID of linear solver for monolithic scatra-scatra interface coupling involving interface layer
  // growth
  IntParameter("INTLAYERGROWTH_LINEAR_SOLVER", -1,
      "ID of linear solver for monolithic scatra-scatra interface coupling involving interface "
      "layer growth",
      &s2icoupling);

  // modified time step size for scatra-scatra interface coupling involving interface layer growth
  DoubleParameter("INTLAYERGROWTH_TIMESTEP", -1.,
      "modified time step size for scatra-scatra interface coupling involving interface layer "
      "growth",
      &s2icoupling);

  BoolParameter("MESHTYING_CONDITIONS_INDEPENDENT_SETUP", "No",
      "mesh tying for different conditions should be setup independently", &s2icoupling);
}


/*------------------------------------------------------------------------*
 | set valid conditions for scatra-scatra interface coupling   fang 01/16 |
 *------------------------------------------------------------------------*/
void INPAR::S2I::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // scatra-scatra interface mesh tying condition
  {
    // definition of scatra-scatra interface mesh tying line condition
    auto s2imeshtyingline =
        Teuchos::rcp(new ConditionDefinition("DESIGN S2I MESHTYING LINE CONDITIONS", "S2IMeshtying",
            "Scatra-scatra line interface mesh tying", DRT::Condition::S2IMeshtying, true,
            DRT::Condition::Line));

    // definition of scatra-scatra interface mesh tying surface condition
    auto s2imeshtyingsurf =
        Teuchos::rcp(new ConditionDefinition("DESIGN S2I MESHTYING SURF CONDITIONS", "S2IMeshtying",
            "Scatra-scatra surface interface mesh tying", DRT::Condition::S2IMeshtying, true,
            DRT::Condition::Surface));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<ConditionComponent>> s2imeshtyingcomponents;
    s2imeshtyingcomponents.emplace_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));
    s2imeshtyingcomponents.emplace_back(Teuchos::rcp(new StringConditionComponent("interface side",
        "Undefined", Teuchos::tuple<std::string>("Undefined", "Slave", "Master"),
        Teuchos::tuple<int>(
            INPAR::S2I::side_undefined, INPAR::S2I::side_slave, INPAR::S2I::side_master))));
    s2imeshtyingcomponents.emplace_back(
        Teuchos::rcp(new SeparatorConditionComponent("S2I_KINETICS_ID")));
    s2imeshtyingcomponents.emplace_back(Teuchos::rcp(new IntConditionComponent("S2IKineticsID")));

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
    auto s2inoevaluationline = Teuchos::rcp(
        new ConditionDefinition("DESIGN S2I NO EVALUATION LINE CONDITIONS", "S2INoEvaluation",
            "Scatra-scatra interface no evaluation line condition. This condition can be used to "
            "deactivate the evaluation of the corresponding `S2IKinetics` condition. Another usage "
            "is in coupled algorithms, where a specific `S2IKinetics` condition should not be "
            "evaluated within the scalar transport framework, as it already evaluated elsewhere. "
            "One example is the SSI contact.",
            DRT::Condition::S2INoEvaluation, true, DRT::Condition::Line));

    // definition of scatra-scatra interface no evaluation surface condition
    auto s2inoevaluationsurf = Teuchos::rcp(new ConditionDefinition(
        "DESIGN S2I NO EVALUATION SURF CONDITIONS", "S2INoEvaluation",
        "Scatra-scatra interface no evaluation surface condition. This condition can be used to "
        "deactivate the evaluation of the corresponding `S2IKinetics` condition. Another usage "
        "is in coupled algorithms, where a specific `S2IKinetics` condition should not be "
        "evaluated within the scalar transport framework, as it already evaluated elsewhere. "
        "One example is the SSI contact.",
        DRT::Condition::S2INoEvaluation, true, DRT::Condition::Surface));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<ConditionComponent>> s2inoevaluationcomponents;
    s2inoevaluationcomponents.emplace_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));
    s2inoevaluationcomponents.emplace_back(Teuchos::rcp(new StringConditionComponent(
        "interface side", "Undefined", Teuchos::tuple<std::string>("Undefined", "Slave", "Master"),
        Teuchos::tuple<int>(
            INPAR::S2I::side_undefined, INPAR::S2I::side_slave, INPAR::S2I::side_master))));
    s2inoevaluationcomponents.emplace_back(
        Teuchos::rcp(new SeparatorConditionComponent("S2I_KINETICS_ID")));
    s2inoevaluationcomponents.emplace_back(
        Teuchos::rcp(new IntConditionComponent("S2IKineticsID")));

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
    auto s2ikineticspoint =
        Teuchos::rcp(new ConditionDefinition("DESIGN S2I KINETICS POINT CONDITIONS", "S2IKinetics",
            "Scatra-scatra line interface kinetics", DRT::Condition::S2IKinetics, true,
            DRT::Condition::Point));

    // definition of scatra-scatra interface kinetics line condition
    auto s2ikineticsline =
        Teuchos::rcp(new ConditionDefinition("DESIGN S2I KINETICS LINE CONDITIONS", "S2IKinetics",
            "Scatra-scatra line interface kinetics", DRT::Condition::S2IKinetics, true,
            DRT::Condition::Line));

    // definition of scatra-scatra interface kinetics surface condition
    auto s2ikineticssurf =
        Teuchos::rcp(new ConditionDefinition("DESIGN S2I KINETICS SURF CONDITIONS", "S2IKinetics",
            "Scatra-scatra surface interface kinetics", DRT::Condition::S2IKinetics, true,
            DRT::Condition::Surface));

    // Macro-micro coupling condition for micro scale in multi-scale scalar transport problems
    auto multiscalecouplingpoint =
        Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA MULTI-SCALE COUPLING POINT CONDITIONS",
            "ScatraMultiScaleCoupling", "Scalar transport multi-scale coupling condition",
            DRT::Condition::ScatraMultiScaleCoupling, false, DRT::Condition::Point));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<ConditionComponent>> s2icomponents;
    {
      // interface ID
      s2icomponents.emplace_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));

      // interface sides for scatra-scatra interface kinetics
      std::vector<Teuchos::RCP<CondCompBundle>> interfacesides;
      {
        {
          // undefined side
          std::vector<Teuchos::RCP<ConditionComponent>> undefinedside;

          // insert undefined-side condition components into vector of interface sides
          interfacesides.emplace_back(Teuchos::rcp(
              new CondCompBundle("Undefined", undefinedside, INPAR::S2I::side_undefined)));
        }

        {
          // slave side
          std::vector<Teuchos::RCP<ConditionComponent>> slaveside;

          // kinetic models for scatra-scatra interface kinetics
          std::vector<Teuchos::RCP<CondCompBundle>> kineticmodels;
          {
            {
              // constant permeability
              std::vector<Teuchos::RCP<ConditionComponent>> constperm;
              // total number of existing scalars
              constperm.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("NUMSCAL")));
              // empty vector --> no separators for integer vectors needed
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp;
              // empty vector --> no integer vectors needed
              std::vector<Teuchos::RCP<IntVectorConditionComponent>> intvectcomp;
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> realsepcomp;
              // string separator in front of real permeability vector in input file line
              realsepcomp.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("PERMEABILITIES")));
              // real vector of constant permeabilities
              std::vector<Teuchos::RCP<RealVectorConditionComponent>> realvectcomp;
              realvectcomp.emplace_back(
                  Teuchos::rcp(new RealVectorConditionComponent("permeabilities", 0)));
              constperm.emplace_back(Teuchos::rcp(new IntRealBundle("permeabilities",
                  Teuchos::rcp(new IntConditionComponent("numscal")), intsepcomp, intvectcomp,
                  realsepcomp, realvectcomp)));
              constperm.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("IS_PSEUDO_CONTACT")));
              constperm.emplace_back(Teuchos::rcp(new IntConditionComponent("is_pseudo_contact")));

              kineticmodels.emplace_back(Teuchos::rcp(new CondCompBundle(
                  "ConstantPermeability", constperm, INPAR::S2I::kinetics_constperm)));
            }

            {
              // Butler-Volmer
              std::vector<Teuchos::RCP<ConditionComponent>> butlervolmer;
              // total number of existing scalars
              butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("NUMSCAL")));
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp;
              // string separator in front of integer stoichiometry vector in input file line
              intsepcomp.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("STOICHIOMETRIES")));
              // integer vector of stoichiometric coefficients
              std::vector<Teuchos::RCP<IntVectorConditionComponent>> intvectcomp;
              intvectcomp.emplace_back(
                  Teuchos::rcp(new IntVectorConditionComponent("stoichiometries", 0)));
              // empty vector --> no separators for real vectors needed
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> realsepcomp;
              // empty vector --> no real vectors needed
              std::vector<Teuchos::RCP<RealVectorConditionComponent>> realvectcomp;
              butlervolmer.emplace_back(Teuchos::rcp(new IntRealBundle("stoichiometries",
                  Teuchos::rcp(new IntConditionComponent("numscal")), intsepcomp, intvectcomp,
                  realsepcomp, realvectcomp)));
              butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("E-")));
              butlervolmer.emplace_back(Teuchos::rcp(new IntConditionComponent("e-")));
              butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("K_R")));
              butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("k_r")));
              butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("ALPHA_A")));
              butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
              butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("ALPHA_C")));
              butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));
              butlervolmer.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("IS_PSEUDO_CONTACT")));
              butlervolmer.emplace_back(
                  Teuchos::rcp(new IntConditionComponent("is_pseudo_contact")));

              kineticmodels.emplace_back(Teuchos::rcp(new CondCompBundle(
                  "Butler-Volmer", butlervolmer, INPAR::S2I::kinetics_butlervolmer)));

              // same components needed for linearized Butler-Volmer
              kineticmodels.emplace_back(Teuchos::rcp(new CondCompBundle("Butler-Volmer_Linearized",
                  butlervolmer, INPAR::S2I::kinetics_butlervolmerlinearized)));

              // same components needed for reduced Butler-Volmer
              kineticmodels.emplace_back(Teuchos::rcp(new CondCompBundle(
                  "Butler-VolmerReduced", butlervolmer, INPAR::S2I::kinetics_butlervolmerreduced)));

              // same components needed for linearized, reduced Butler-Volmer
              kineticmodels.emplace_back(
                  Teuchos::rcp(new CondCompBundle("Butler-VolmerReduced_Linearized", butlervolmer,
                      INPAR::S2I::kinetics_butlervolmerreducedlinearized)));
            }

            {
              // Butler-Volmer-Peltier
              std::vector<Teuchos::RCP<ConditionComponent>> butlervolmerpeltier;
              // total number of existing scalars
              butlervolmerpeltier.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("NUMSCAL")));
              // string separator in front of integer stoichiometry vector in input file line
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp;
              intsepcomp.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("STOICHIOMETRIES")));
              // integer vector of stoichiometric coefficients
              std::vector<Teuchos::RCP<IntVectorConditionComponent>> intvectcomp;
              intvectcomp.emplace_back(
                  Teuchos::rcp(new IntVectorConditionComponent("stoichiometries", 0)));
              // empty vector --> no separators for real vectors needed
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> realsepcomp;
              // empty vector --> no real vectors needed
              std::vector<Teuchos::RCP<RealVectorConditionComponent>> realvectcomp;
              butlervolmerpeltier.emplace_back(Teuchos::rcp(new IntRealBundle("stoichiometries",
                  Teuchos::rcp(new IntConditionComponent("numscal")), intsepcomp, intvectcomp,
                  realsepcomp, realvectcomp)));
              butlervolmerpeltier.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("E-")));
              butlervolmerpeltier.emplace_back(Teuchos::rcp(new IntConditionComponent("e-")));
              butlervolmerpeltier.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("K_R")));
              butlervolmerpeltier.emplace_back(Teuchos::rcp(new RealConditionComponent("k_r")));
              butlervolmerpeltier.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("ALPHA_A")));
              butlervolmerpeltier.emplace_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
              butlervolmerpeltier.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("ALPHA_C")));
              butlervolmerpeltier.emplace_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));
              butlervolmerpeltier.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("IS_PSEUDO_CONTACT")));
              butlervolmerpeltier.emplace_back(
                  Teuchos::rcp(new IntConditionComponent("is_pseudo_contact")));
              butlervolmerpeltier.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("PELTIER")));
              butlervolmerpeltier.emplace_back(Teuchos::rcp(new RealConditionComponent("peltier")));

              kineticmodels.emplace_back(Teuchos::rcp(new CondCompBundle("Butler-Volmer-Peltier",
                  butlervolmerpeltier, INPAR::S2I::kinetics_butlervolmerpeltier)));
            }

            {
              // Butler-Volmer-reduced with interface capacitance
              std::vector<Teuchos::RCP<ConditionComponent>> butlervolmerreducedcapacitance;
              // total number of existing scalars
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("NUMSCAL")));
              // string separator in front of integer stoichiometry vector in input file line
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp;
              intsepcomp.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("STOICHIOMETRIES")));
              // integer vector of stoichiometric coefficients
              std::vector<Teuchos::RCP<IntVectorConditionComponent>> intvectcomp;
              intvectcomp.emplace_back(
                  Teuchos::rcp(new IntVectorConditionComponent("stoichiometries", 0)));
              // empty vector --> no separators for real vectors needed
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> realsepcomp;
              // empty vector --> no real vectors needed
              std::vector<Teuchos::RCP<RealVectorConditionComponent>> realvectcomp;
              butlervolmerreducedcapacitance.emplace_back(Teuchos::rcp(new IntRealBundle(
                  "stoichiometries", Teuchos::rcp(new IntConditionComponent("numscal")), intsepcomp,
                  intvectcomp, realsepcomp, realvectcomp)));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("E-")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new IntConditionComponent("e-")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("K_R")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("k_r")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("CAPACITANCE")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("capacitance")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("ALPHA_A")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("alpha_a")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("ALPHA_C")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("alpha_c")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("IS_PSEUDO_CONTACT")));
              butlervolmerreducedcapacitance.emplace_back(
                  Teuchos::rcp(new IntConditionComponent("is_pseudo_contact")));

              kineticmodels.emplace_back(Teuchos::rcp(new CondCompBundle(
                  "Butler-VolmerReduced_Capacitance", butlervolmerreducedcapacitance,
                  INPAR::S2I::kinetics_butlervolmerreducedcapacitance)));
            }

            {
              // Butler-Volmer-Resistance
              std::vector<Teuchos::RCP<ConditionComponent>> butlervolmerresistance;
              // total number of existing scalars
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("NUMSCAL")));
              // string separator in front of integer stoichiometry vector in input file line
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp;
              intsepcomp.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("STOICHIOMETRIES")));
              // integer vector of stoichiometric coefficients
              std::vector<Teuchos::RCP<IntVectorConditionComponent>> intvectcomp;
              intvectcomp.emplace_back(
                  Teuchos::rcp(new IntVectorConditionComponent("stoichiometries", 0)));
              // empty vector --> no separators for real vectors needed
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> realsepcomp;
              // empty vector --> no real vectors needed
              std::vector<Teuchos::RCP<RealVectorConditionComponent>> realvectcomp;
              butlervolmerresistance.emplace_back(Teuchos::rcp(new IntRealBundle("stoichiometries",
                  Teuchos::rcp(new IntConditionComponent("numscal")), intsepcomp, intvectcomp,
                  realsepcomp, realvectcomp)));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("E-")));
              butlervolmerresistance.emplace_back(Teuchos::rcp(new IntConditionComponent("e-")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("K_R")));
              butlervolmerresistance.emplace_back(Teuchos::rcp(new RealConditionComponent("k_r")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("ALPHA_A")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("alpha_a")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("ALPHA_C")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("alpha_c")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("IS_PSEUDO_CONTACT")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new IntConditionComponent("is_pseudo_contact")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("RESISTANCE")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("resistance")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("CONVTOL_IMPLBUTLERVOLMER")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("CONVTOL_IMPLBUTLERVOLMER")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("ITEMAX_IMPLBUTLERVOLMER")));
              butlervolmerresistance.emplace_back(
                  Teuchos::rcp(new IntConditionComponent("ITEMAX_IMPLBUTLERVOLMER")));

              kineticmodels.emplace_back(Teuchos::rcp(new CondCompBundle("Butler-Volmer_Resistance",
                  butlervolmerresistance, INPAR::S2I::kinetics_butlervolmerresistance)));
            }

            {
              // Butler-Volmer-Reduced with resistance
              std::vector<Teuchos::RCP<ConditionComponent>> butlervolmerreducedwithresistance;
              // total number of existing scalars
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("NUMSCAL")));
              // string separator in front of integer stoichiometry vector in input file line
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp;
              intsepcomp.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("STOICHIOMETRIES")));
              // integer vector of stoichiometric coefficients
              std::vector<Teuchos::RCP<IntVectorConditionComponent>> intvectcomp;
              intvectcomp.emplace_back(
                  Teuchos::rcp(new IntVectorConditionComponent("stoichiometries", 0)));
              // empty vector --> no separators for real vectors needed
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> realsepcomp;
              // empty vector --> no real vectors needed
              std::vector<Teuchos::RCP<RealVectorConditionComponent>> realvectcomp;
              butlervolmerreducedwithresistance.emplace_back(Teuchos::rcp(new IntRealBundle(
                  "stoichiometries", Teuchos::rcp(new IntConditionComponent("numscal")), intsepcomp,
                  intvectcomp, realsepcomp, realvectcomp)));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("E-")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new IntConditionComponent("e-")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("K_R")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("k_r")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("ALPHA_A")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("alpha_a")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("ALPHA_C")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("alpha_c")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("IS_PSEUDO_CONTACT")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new IntConditionComponent("is_pseudo_contact")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("RESISTANCE")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("resistance")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("CONVTOL_IMPLBUTLERVOLMER")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("CONVTOL_IMPLBUTLERVOLMER")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("ITEMAX_IMPLBUTLERVOLMER")));
              butlervolmerreducedwithresistance.emplace_back(
                  Teuchos::rcp(new IntConditionComponent("ITEMAX_IMPLBUTLERVOLMER")));

              kineticmodels.emplace_back(Teuchos::rcp(new CondCompBundle(
                  "Butler-VolmerReduced_Resistance", butlervolmerreducedwithresistance,
                  INPAR::S2I::kinetics_butlervolmerreducedresistance)));
            }

            {
              // Butler-Volmer-reduced-thermoresistance
              std::vector<Teuchos::RCP<ConditionComponent>> butlervolmerreducedthermo;
              // total number of existing scalars
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("NUMSCAL")));
              // string separator in front of integer stoichiometry vector in input file line
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp;
              intsepcomp.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("STOICHIOMETRIES")));
              // integer vector of stoichiometric coefficients
              std::vector<Teuchos::RCP<IntVectorConditionComponent>> intvectcomp;
              intvectcomp.emplace_back(
                  Teuchos::rcp(new IntVectorConditionComponent("stoichiometries", 0)));
              // empty vector --> no separators for real vectors needed
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> realsepcomp;
              // empty vector --> no real vectors needed
              std::vector<Teuchos::RCP<RealVectorConditionComponent>> realvectcomp;
              butlervolmerreducedthermo.emplace_back(Teuchos::rcp(new IntRealBundle(
                  "stoichiometries", Teuchos::rcp(new IntConditionComponent("numscal")), intsepcomp,
                  intvectcomp, realsepcomp, realvectcomp)));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("E-")));
              butlervolmerreducedthermo.emplace_back(Teuchos::rcp(new IntConditionComponent("e-")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("K_R")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("k_r")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("ALPHA_A")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("alpha_a")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("ALPHA_C")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("alpha_c")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("IS_PSEUDO_CONTACT")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new IntConditionComponent("is_pseudo_contact")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("THERMOPERM")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("thermoperm")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("MOLAR_HEAT_CAPACITY")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("molar_heat_capacity")));

              kineticmodels.emplace_back(Teuchos::rcp(new CondCompBundle(
                  "Butler-VolmerReduced_ThermoResistance", butlervolmerreducedthermo,
                  INPAR::S2I::kinetics_butlervolmerreducedthermoresistance)));
            }

            {
              // constant interface resistance
              std::vector<Teuchos::RCP<ConditionComponent>> constantinterfaceresistance;

              constantinterfaceresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
              constantinterfaceresistance.emplace_back(
                  Teuchos::rcp(new IntVectorConditionComponent("onoff", 2)));
              constantinterfaceresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("RESISTANCE")));
              constantinterfaceresistance.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("resistance")));
              constantinterfaceresistance.emplace_back(new SeparatorConditionComponent("E-"));
              constantinterfaceresistance.emplace_back(new IntConditionComponent("e-"));
              constantinterfaceresistance.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("IS_PSEUDO_CONTACT")));
              constantinterfaceresistance.emplace_back(
                  Teuchos::rcp(new IntConditionComponent("is_pseudo_contact")));

              kineticmodels.emplace_back(Teuchos::rcp(
                  new CondCompBundle("ConstantInterfaceResistance", constantinterfaceresistance,
                      INPAR::S2I::kinetics_constantinterfaceresistance)));
            }

            {
              // no interface flux
              std::vector<Teuchos::RCP<ConditionComponent>> nointerfaceflux;

              kineticmodels.emplace_back(Teuchos::rcp(new CondCompBundle(
                  "NoInterfaceFlux", nointerfaceflux, INPAR::S2I::kinetics_nointerfaceflux)));
            }
          }  // kinetic models for scatra-scatra interface kinetics

          // insert kinetic models into vector with slave-side condition components
          slaveside.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("KINETIC_MODEL")));
          slaveside.emplace_back(
              Teuchos::rcp(new CondCompBundleSelector("kinetic model", kineticmodels)));

          // add all components from slave side to multi-scale condition
          for (const auto& component : slaveside) multiscalecouplingpoint->AddComponent(component);

          // insert slave-side condition components into vector of interface sides
          interfacesides.emplace_back(
              Teuchos::rcp(new CondCompBundle("Slave", slaveside, INPAR::S2I::side_slave)));
        }

        {
          // master side
          std::vector<Teuchos::RCP<ConditionComponent>> masterside;

          // insert master-side condition components into vector of interface sides
          interfacesides.emplace_back(
              Teuchos::rcp(new CondCompBundle("Master", masterside, INPAR::S2I::side_master)));
        }
      }  // interface sides for scatra-scatra interface mesh tying

      // insert interface sides into vector with input file line components
      s2icomponents.emplace_back(
          Teuchos::rcp(new CondCompBundleSelector("interface side", interfacesides)));
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
        new ConditionDefinition("DESIGN S2I COUPLING GROWTH LINE CONDITIONS", "S2ICouplingGrowth",
            "Scatra-scatra line interface coupling involving interface layer growth",
            DRT::Condition::S2ICouplingGrowth, true, DRT::Condition::Line));

    // definition of scatra-scatra interface coupling surface condition involving interface layer
    // growth
    auto s2igrowthsurf = Teuchos::rcp(
        new ConditionDefinition("DESIGN S2I COUPLING GROWTH SURF CONDITIONS", "S2ICouplingGrowth",
            "Scatra-scatra surface interface coupling involving interface layer growth",
            DRT::Condition::S2ICouplingGrowth, true, DRT::Condition::Surface));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<ConditionComponent>> s2igrowthcomponents;
    {
      // interface ID
      s2igrowthcomponents.emplace_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));

      // kinetic models for scatra-scatra interface coupling involving interface layer growth
      std::vector<Teuchos::RCP<CondCompBundle>> kineticmodels;
      {
        {
          // Butler-Volmer
          std::vector<Teuchos::RCP<ConditionComponent>> butlervolmer;
          butlervolmer.emplace_back(Teuchos::rcp(
              new SeparatorConditionComponent("NUMSCAL")));  // total number of existing scalars
          std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp;
          intsepcomp.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("STOICHIOMETRIES")));
          std::vector<Teuchos::RCP<IntVectorConditionComponent>>
              intvectcomp;  // string separator in front of integer stoichiometry vector in input
                            // file line
          intvectcomp.emplace_back(Teuchos::rcp(new IntVectorConditionComponent(
              "stoichiometries", 0)));  // integer vector of stoichiometric coefficients
          std::vector<Teuchos::RCP<SeparatorConditionComponent>>
              realsepcomp;  // empty vector --> no separators for real vectors needed
          std::vector<Teuchos::RCP<RealVectorConditionComponent>>
              realvectcomp;  // empty vector --> no real vectors needed
          butlervolmer.emplace_back(Teuchos::rcp(new IntRealBundle("stoichiometries",
              Teuchos::rcp(new IntConditionComponent("numscal")), intsepcomp, intvectcomp,
              realsepcomp, realvectcomp)));
          butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("E-")));
          butlervolmer.emplace_back(Teuchos::rcp(new IntConditionComponent("e-")));
          butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("K_R")));
          butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("k_r")));
          butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("ALPHA_A")));
          butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
          butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("ALPHA_C")));
          butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));
          butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("MOLMASS")));
          butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("molar mass")));
          butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("DENSITY")));
          butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("density")));
          butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("CONDUCTIVITY")));
          butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("conductivity")));
          butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("REGTYPE")));
          butlervolmer.emplace_back(
              Teuchos::rcp(new StringConditionComponent("regularization type", "trigonometrical",
                  Teuchos::tuple<std::string>("none", "polynomial", "Hein", "trigonometrical"),
                  Teuchos::tuple<int>(INPAR::S2I::regularization_none,
                      INPAR::S2I::regularization_polynomial, INPAR::S2I::regularization_hein,
                      INPAR::S2I::regularization_trigonometrical))));
          butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("REGPAR")));
          butlervolmer.emplace_back(
              Teuchos::rcp(new RealConditionComponent("regularization parameter")));
          butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("INITTHICKNESS")));
          butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("initial thickness")));

          kineticmodels.emplace_back(Teuchos::rcp(new CondCompBundle(
              "Butler-Volmer", butlervolmer, INPAR::S2I::growth_kinetics_butlervolmer)));
        }
      }  // kinetic models for scatra-scatra interface coupling involving interface layer growth

      // insert kinetic models into vector with input file line components
      s2igrowthcomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("KINETIC_MODEL")));
      s2igrowthcomponents.emplace_back(
          Teuchos::rcp(new CondCompBundleSelector("kinetic model", kineticmodels)));
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
    auto s2isclcond =
        Teuchos::rcp(new ConditionDefinition("DESIGN S2I SCL COUPLING SURF CONDITIONS",
            "S2ISCLCoupling", "Scatra-scatra surface with SCL micro-macro coupling between",
            DRT::Condition::S2ISCLCoupling, true, DRT::Condition::Surface));

    s2isclcond->AddComponent(Teuchos::rcp(new StringConditionComponent("interface side",
        "Undefined", Teuchos::tuple<std::string>("Undefined", "Slave", "Master"),
        Teuchos::tuple<int>(
            INPAR::S2I::side_undefined, INPAR::S2I::side_slave, INPAR::S2I::side_master))));

    condlist.emplace_back(s2isclcond);
  }
}
