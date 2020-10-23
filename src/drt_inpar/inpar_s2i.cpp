/*----------------------------------------------------------------------*/
/*! \file
\brief input quantities and globally accessible enumerations for scatra-scatra interface coupling

\level 2


*/
/*----------------------------------------------------------------------*/
#include "inpar_s2i.H"

#include "drt_validparameters.H"

#include "../drt_lib/drt_conditiondefinition.H"

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
}


/*------------------------------------------------------------------------*
 | set valid conditions for scatra-scatra interface coupling   fang 01/16 |
 *------------------------------------------------------------------------*/
void INPAR::S2I::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // scatra-scatra interface coupling
  {
    // definition of scatra-scatra interface coupling line condition
    Teuchos::RCP<ConditionDefinition> s2iline =
        Teuchos::rcp(new ConditionDefinition("DESIGN S2I COUPLING LINE CONDITIONS", "S2ICoupling",
            "Scatra-scatra line interface coupling", DRT::Condition::S2ICoupling, true,
            DRT::Condition::Line));

    // definition of scatra-scatra interface coupling surface condition
    Teuchos::RCP<ConditionDefinition> s2isurf =
        Teuchos::rcp(new ConditionDefinition("DESIGN S2I COUPLING SURF CONDITIONS", "S2ICoupling",
            "Scatra-scatra surface interface coupling", DRT::Condition::S2ICoupling, true,
            DRT::Condition::Surface));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<ConditionComponent>> s2icomponents;

    {
      // interface ID
      s2icomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));

      // interface sides for scatra-scatra interface coupling
      std::vector<Teuchos::RCP<CondCompBundle>> interfacesides;
      {
        {
          // undefined side
          std::vector<Teuchos::RCP<ConditionComponent>> undefinedside;

          // insert undefined-side condition components into vector of interface sides
          interfacesides.push_back(Teuchos::rcp(
              new CondCompBundle("Undefined", undefinedside, INPAR::S2I::side_undefined)));
        }

        {
          // slave side
          std::vector<Teuchos::RCP<ConditionComponent>> slaveside;

          // kinetic models for scatra-scatra interface coupling
          std::vector<Teuchos::RCP<CondCompBundle>> kineticmodels;
          {
            {
              // constant permeability
              std::vector<Teuchos::RCP<ConditionComponent>> constperm;
              constperm.push_back(Teuchos::rcp(
                  new SeparatorConditionComponent("numscal")));  // total number of existing scalars
              std::vector<Teuchos::RCP<SeparatorConditionComponent>>
                  intsepcomp;  // empty vector --> no separators for integer vectors needed
              std::vector<Teuchos::RCP<IntVectorConditionComponent>>
                  intvectcomp;  // empty vector --> no integer vectors needed
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> realsepcomp;
              realsepcomp.push_back(Teuchos::rcp(new SeparatorConditionComponent(
                  "permeabilities")));  // string separator in front of real permeability vector in
                                        // input file line
              std::vector<Teuchos::RCP<RealVectorConditionComponent>> realvectcomp;
              realvectcomp.push_back(Teuchos::rcp(new RealVectorConditionComponent(
                  "permeabilities", 0)));  // real vector of constant permeabilities
              constperm.push_back(Teuchos::rcp(new IntRealBundle("permeabilities",
                  Teuchos::rcp(new IntConditionComponent("numscal")), intsepcomp, intvectcomp,
                  realsepcomp, realvectcomp)));

              kineticmodels.push_back(Teuchos::rcp(new CondCompBundle(
                  "ConstantPermeability", constperm, INPAR::S2I::kinetics_constperm)));
            }

            {
              // Butler-Volmer
              std::vector<Teuchos::RCP<ConditionComponent>> butlervolmer;
              butlervolmer.push_back(Teuchos::rcp(
                  new SeparatorConditionComponent("numscal")));  // total number of existing scalars
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp;
              intsepcomp.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("stoichiometries")));
              std::vector<Teuchos::RCP<IntVectorConditionComponent>>
                  intvectcomp;  // string separator in front of integer stoichiometry vector in
                                // input file line
              intvectcomp.push_back(Teuchos::rcp(new IntVectorConditionComponent(
                  "stoichiometries", 0)));  // integer vector of stoichiometric coefficients
              std::vector<Teuchos::RCP<SeparatorConditionComponent>>
                  realsepcomp;  // empty vector --> no separators for real vectors needed
              std::vector<Teuchos::RCP<RealVectorConditionComponent>>
                  realvectcomp;  // empty vector --> no real vectors needed
              butlervolmer.push_back(Teuchos::rcp(new IntRealBundle("stoichiometries",
                  Teuchos::rcp(new IntConditionComponent("numscal")), intsepcomp, intvectcomp,
                  realsepcomp, realvectcomp)));
              butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("e-")));
              butlervolmer.push_back(Teuchos::rcp(new IntConditionComponent("e-")));
              butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("k_r")));
              butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("k_r")));
              butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_a")));
              butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
              butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_c")));
              butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));

              kineticmodels.push_back(Teuchos::rcp(new CondCompBundle(
                  "Butler-Volmer", butlervolmer, INPAR::S2I::kinetics_butlervolmer)));
            }

            {
              // Butler-Volmer-Peltier
              std::vector<Teuchos::RCP<ConditionComponent>> butlervolmerpeltier;
              butlervolmerpeltier.push_back(Teuchos::rcp(
                  new SeparatorConditionComponent("numscal")));  // total number of existing scalars
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp;
              intsepcomp.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("stoichiometries")));
              std::vector<Teuchos::RCP<IntVectorConditionComponent>>
                  intvectcomp;  // string separator in front of integer stoichiometry vector in
                                // input file line
              intvectcomp.push_back(Teuchos::rcp(new IntVectorConditionComponent(
                  "stoichiometries", 0)));  // integer vector of stoichiometric coefficients
              std::vector<Teuchos::RCP<SeparatorConditionComponent>>
                  realsepcomp;  // empty vector --> no separators for real vectors needed
              std::vector<Teuchos::RCP<RealVectorConditionComponent>>
                  realvectcomp;  // empty vector --> no real vectors needed
              butlervolmerpeltier.push_back(Teuchos::rcp(new IntRealBundle("stoichiometries",
                  Teuchos::rcp(new IntConditionComponent("numscal")), intsepcomp, intvectcomp,
                  realsepcomp, realvectcomp)));
              butlervolmerpeltier.push_back(Teuchos::rcp(new SeparatorConditionComponent("e-")));
              butlervolmerpeltier.push_back(Teuchos::rcp(new IntConditionComponent("e-")));
              butlervolmerpeltier.push_back(Teuchos::rcp(new SeparatorConditionComponent("k_r")));
              butlervolmerpeltier.push_back(Teuchos::rcp(new RealConditionComponent("k_r")));
              butlervolmerpeltier.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("alpha_a")));
              butlervolmerpeltier.push_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
              butlervolmerpeltier.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("alpha_c")));
              butlervolmerpeltier.push_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));
              butlervolmerpeltier.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("peltier")));
              butlervolmerpeltier.push_back(Teuchos::rcp(new RealConditionComponent("peltier")));

              kineticmodels.push_back(Teuchos::rcp(new CondCompBundle("Butler-Volmer-Peltier",
                  butlervolmerpeltier, INPAR::S2I::kinetics_butlervolmerpeltier)));
            }

            {
              // Butler-Volmer-reduced
              std::vector<Teuchos::RCP<ConditionComponent>> butlervolmerreduced;
              butlervolmerreduced.push_back(Teuchos::rcp(
                  new SeparatorConditionComponent("numscal")));  // total number of existing scalars
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp;
              intsepcomp.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("stoichiometries")));
              std::vector<Teuchos::RCP<IntVectorConditionComponent>>
                  intvectcomp;  // string separator in front of integer stoichiometry vector in
                                // input file line
              intvectcomp.push_back(Teuchos::rcp(new IntVectorConditionComponent(
                  "stoichiometries", 0)));  // integer vector of stoichiometric coefficients
              std::vector<Teuchos::RCP<SeparatorConditionComponent>>
                  realsepcomp;  // empty vector --> no separators for real vectors needed
              std::vector<Teuchos::RCP<RealVectorConditionComponent>>
                  realvectcomp;  // empty vector --> no real vectors needed
              butlervolmerreduced.push_back(Teuchos::rcp(new IntRealBundle("stoichiometries",
                  Teuchos::rcp(new IntConditionComponent("numscal")), intsepcomp, intvectcomp,
                  realsepcomp, realvectcomp)));
              butlervolmerreduced.push_back(Teuchos::rcp(new SeparatorConditionComponent("e-")));
              butlervolmerreduced.push_back(Teuchos::rcp(new IntConditionComponent("e-")));
              butlervolmerreduced.push_back(Teuchos::rcp(new SeparatorConditionComponent("k_r")));
              butlervolmerreduced.push_back(Teuchos::rcp(new RealConditionComponent("k_r")));
              butlervolmerreduced.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("alpha_a")));
              butlervolmerreduced.push_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
              butlervolmerreduced.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("alpha_c")));
              butlervolmerreduced.push_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));

              kineticmodels.push_back(Teuchos::rcp(new CondCompBundle("Butler-Volmer-reduced",
                  butlervolmerreduced, INPAR::S2I::kinetics_butlervolmerreduced)));
            }

            {
              // Butler-Volmer-Resistance
              std::vector<Teuchos::RCP<ConditionComponent>> butlervolmerresistance;
              butlervolmerresistance.push_back(Teuchos::rcp(
                  new SeparatorConditionComponent("numscal")));  // total number of existing scalars
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp;
              intsepcomp.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("stoichiometries")));
              std::vector<Teuchos::RCP<IntVectorConditionComponent>>
                  intvectcomp;  // string separator in front of integer stoichiometry vector in
                                // input file line
              intvectcomp.push_back(Teuchos::rcp(new IntVectorConditionComponent(
                  "stoichiometries", 0)));  // integer vector of stoichiometric coefficients
              std::vector<Teuchos::RCP<SeparatorConditionComponent>>
                  realsepcomp;  // empty vector --> no separators for real vectors needed
              std::vector<Teuchos::RCP<RealVectorConditionComponent>>
                  realvectcomp;  // empty vector --> no real vectors needed
              butlervolmerresistance.push_back(Teuchos::rcp(new IntRealBundle("stoichiometries",
                  Teuchos::rcp(new IntConditionComponent("numscal")), intsepcomp, intvectcomp,
                  realsepcomp, realvectcomp)));
              butlervolmerresistance.push_back(Teuchos::rcp(new SeparatorConditionComponent("e-")));
              butlervolmerresistance.push_back(Teuchos::rcp(new IntConditionComponent("e-")));
              butlervolmerresistance.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("k_r")));
              butlervolmerresistance.push_back(Teuchos::rcp(new RealConditionComponent("k_r")));
              butlervolmerresistance.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("alpha_a")));
              butlervolmerresistance.push_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
              butlervolmerresistance.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("alpha_c")));
              butlervolmerresistance.push_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));
              butlervolmerresistance.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("resistance")));
              butlervolmerresistance.push_back(
                  Teuchos::rcp(new RealConditionComponent("resistance")));
              butlervolmerresistance.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("CONVTOL_IMPLBUTLERVOLMER")));
              butlervolmerresistance.push_back(
                  Teuchos::rcp(new RealConditionComponent("CONVTOL_IMPLBUTLERVOLMER")));
              butlervolmerresistance.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("ITEMAX_IMPLBUTLERVOLMER")));
              butlervolmerresistance.push_back(
                  Teuchos::rcp(new IntConditionComponent("ITEMAX_IMPLBUTLERVOLMER")));

              kineticmodels.push_back(Teuchos::rcp(new CondCompBundle("Butler-Volmer-Resistance",
                  butlervolmerresistance, INPAR::S2I::kinetics_butlervolmerresistance)));
            }

            {
              // Butler-Volmer-Reduced with resistance
              std::vector<Teuchos::RCP<ConditionComponent>> butlervolmerreducedwithresistance;
              butlervolmerreducedwithresistance.push_back(Teuchos::rcp(
                  new SeparatorConditionComponent("numscal")));  // total number of existing scalars
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp;
              intsepcomp.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("stoichiometries")));
              std::vector<Teuchos::RCP<IntVectorConditionComponent>>
                  intvectcomp;  // string separator in front of integer stoichiometry vector in
                                // input file line
              intvectcomp.push_back(Teuchos::rcp(new IntVectorConditionComponent(
                  "stoichiometries", 0)));  // integer vector of stoichiometric coefficients
              std::vector<Teuchos::RCP<SeparatorConditionComponent>>
                  realsepcomp;  // empty vector --> no separators for real vectors needed
              std::vector<Teuchos::RCP<RealVectorConditionComponent>>
                  realvectcomp;  // empty vector --> no real vectors needed
              butlervolmerreducedwithresistance.push_back(Teuchos::rcp(new IntRealBundle(
                  "stoichiometries", Teuchos::rcp(new IntConditionComponent("numscal")), intsepcomp,
                  intvectcomp, realsepcomp, realvectcomp)));
              butlervolmerreducedwithresistance.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("e-")));
              butlervolmerreducedwithresistance.push_back(
                  Teuchos::rcp(new IntConditionComponent("e-")));
              butlervolmerreducedwithresistance.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("k_r")));
              butlervolmerreducedwithresistance.push_back(
                  Teuchos::rcp(new RealConditionComponent("k_r")));
              butlervolmerreducedwithresistance.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("alpha_a")));
              butlervolmerreducedwithresistance.push_back(
                  Teuchos::rcp(new RealConditionComponent("alpha_a")));
              butlervolmerreducedwithresistance.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("alpha_c")));
              butlervolmerreducedwithresistance.push_back(
                  Teuchos::rcp(new RealConditionComponent("alpha_c")));
              butlervolmerreducedwithresistance.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("resistance")));
              butlervolmerreducedwithresistance.push_back(
                  Teuchos::rcp(new RealConditionComponent("resistance")));
              butlervolmerreducedwithresistance.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("CONVTOL_IMPLBUTLERVOLMER")));
              butlervolmerreducedwithresistance.push_back(
                  Teuchos::rcp(new RealConditionComponent("CONVTOL_IMPLBUTLERVOLMER")));
              butlervolmerreducedwithresistance.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("ITEMAX_IMPLBUTLERVOLMER")));
              butlervolmerreducedwithresistance.push_back(
                  Teuchos::rcp(new IntConditionComponent("ITEMAX_IMPLBUTLERVOLMER")));

              kineticmodels.push_back(Teuchos::rcp(new CondCompBundle(
                  "Butler-Volmer-ReducedwithResistance", butlervolmerreducedwithresistance,
                  INPAR::S2I::kinetics_butlervolmerreducedresistance)));
            }

            {
              // Butler-Volmer-reduced-thermoresistance
              std::vector<Teuchos::RCP<ConditionComponent>> butlervolmerreducedthermo;
              butlervolmerreducedthermo.emplace_back(Teuchos::rcp(
                  new SeparatorConditionComponent("numscal")));  // total number of existing scalars
              std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp;
              intsepcomp.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("stoichiometries")));
              std::vector<Teuchos::RCP<IntVectorConditionComponent>>
                  intvectcomp;  // string separator in front of integer stoichiometry vector in
                                // input file line
              intvectcomp.emplace_back(Teuchos::rcp(new IntVectorConditionComponent(
                  "stoichiometries", 0)));  // integer vector of stoichiometric coefficients
              std::vector<Teuchos::RCP<SeparatorConditionComponent>>
                  realsepcomp;  // empty vector --> no separators for real vectors needed
              std::vector<Teuchos::RCP<RealVectorConditionComponent>>
                  realvectcomp;  // empty vector --> no real vectors needed
              butlervolmerreducedthermo.emplace_back(Teuchos::rcp(new IntRealBundle(
                  "stoichiometries", Teuchos::rcp(new IntConditionComponent("numscal")), intsepcomp,
                  intvectcomp, realsepcomp, realvectcomp)));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("e-")));
              butlervolmerreducedthermo.emplace_back(Teuchos::rcp(new IntConditionComponent("e-")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("k_r")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("k_r")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("alpha_a")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("alpha_a")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("alpha_c")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("alpha_c")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("thermoperm")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("thermoperm")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new SeparatorConditionComponent("molar_heat_capacity")));
              butlervolmerreducedthermo.emplace_back(
                  Teuchos::rcp(new RealConditionComponent("molar_heat_capacity")));

              kineticmodels.emplace_back(Teuchos::rcp(new CondCompBundle(
                  "Butler-Volmer-reduced-thermoresistance", butlervolmerreducedthermo,
                  INPAR::S2I::kinetics_butlervolmerreducedthermoresistance)));
            }

            {
              // constant interface resistance
              std::vector<Teuchos::RCP<ConditionComponent>> constantinterfaceresistance;
              constantinterfaceresistance.push_back(
                  Teuchos::rcp(new SeparatorConditionComponent("resistance")));
              constantinterfaceresistance.push_back(
                  Teuchos::rcp(new RealConditionComponent("resistance")));

              kineticmodels.push_back(Teuchos::rcp(new CondCompBundle("ConstantInterfaceResistance",
                  constantinterfaceresistance, INPAR::S2I::kinetics_constantinterfaceresistance)));
            }

            {
              // no interface flux
              std::vector<Teuchos::RCP<ConditionComponent>> nointerfaceflux;

              kineticmodels.push_back(Teuchos::rcp(new CondCompBundle(
                  "NoInterfaceFlux", nointerfaceflux, INPAR::S2I::kinetics_nointerfaceflux)));
            }
          }  // kinetic models for scatra-scatra interface coupling

          // insert kinetic models into vector with slave-side condition components
          slaveside.push_back(Teuchos::rcp(new SeparatorConditionComponent("KineticModel")));
          slaveside.push_back(Teuchos::rcp(new CondCompBundleSelector(
              "kinetic models for scatra-scatra interface coupling",
              Teuchos::rcp(new StringConditionComponent("kinetic model", "ConstantPermeability",
                  Teuchos::tuple<std::string>("ConstantPermeability", "Butler-Volmer",
                      "Butler-Volmer-Peltier", "Butler-VolmerReduced", "Butler-Volmer_Resistance",
                      "Butler-VolmerReduced_Resistance", "Butler-VolmerReduced_ThermoResistance",
                      "ConstantInterfaceResistance", "NoInterfaceFlux"),
                  Teuchos::tuple<int>(INPAR::S2I::kinetics_constperm,
                      INPAR::S2I::kinetics_butlervolmer, INPAR::S2I::kinetics_butlervolmerpeltier,
                      INPAR::S2I::kinetics_butlervolmerreduced,
                      INPAR::S2I::kinetics_butlervolmerresistance,
                      INPAR::S2I::kinetics_butlervolmerreducedresistance,
                      INPAR::S2I::kinetics_butlervolmerreducedthermoresistance,
                      INPAR::S2I::kinetics_constantinterfaceresistance,
                      INPAR::S2I::kinetics_nointerfaceflux))),
              kineticmodels)));

          // insert slave-side condition components into vector of interface sides
          interfacesides.push_back(
              Teuchos::rcp(new CondCompBundle("Slave", slaveside, INPAR::S2I::side_slave)));
        }

        {
          // master side
          std::vector<Teuchos::RCP<ConditionComponent>> masterside;

          // insert master-side condition components into vector of interface sides
          interfacesides.push_back(
              Teuchos::rcp(new CondCompBundle("Master", masterside, INPAR::S2I::side_master)));
        }
      }  // interface sides for scatra-scatra interface coupling

      // insert interface sides into vector with input file line components
      s2icomponents.push_back(Teuchos::rcp(
          new CondCompBundleSelector("interface sides for scatra-scatra interface coupling",
              Teuchos::rcp(new StringConditionComponent("interface side", "Undefined",
                  Teuchos::tuple<std::string>("Undefined", "Slave", "Master"),
                  Teuchos::tuple<int>(INPAR::S2I::side_undefined, INPAR::S2I::side_slave,
                      INPAR::S2I::side_master))),
              interfacesides)));
    }

    // insert input file line components into condition definitions
    for (unsigned i = 0; i < s2icomponents.size(); ++i)
    {
      s2iline->AddComponent(s2icomponents[i]);
      s2isurf->AddComponent(s2icomponents[i]);
    }

    // insert condition definitions into global list of valid condition definitions
    condlist.push_back(s2iline);
    condlist.push_back(s2isurf);
  }


  /*--------------------------------------------------------------------*/
  // scatra-scatra interface coupling involving interface layer growth
  {
    // definition of scatra-scatra interface coupling line condition involving interface layer
    // growth
    Teuchos::RCP<ConditionDefinition> s2igrowthline = Teuchos::rcp(
        new ConditionDefinition("DESIGN S2I COUPLING GROWTH LINE CONDITIONS", "S2ICouplingGrowth",
            "Scatra-scatra line interface coupling involving interface layer growth",
            DRT::Condition::S2ICouplingGrowth, true, DRT::Condition::Line));

    // definition of scatra-scatra interface coupling surface condition involving interface layer
    // growth
    Teuchos::RCP<ConditionDefinition> s2igrowthsurf = Teuchos::rcp(
        new ConditionDefinition("DESIGN S2I COUPLING GROWTH SURF CONDITIONS", "S2ICouplingGrowth",
            "Scatra-scatra surface interface coupling involving interface layer growth",
            DRT::Condition::S2ICouplingGrowth, true, DRT::Condition::Surface));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<ConditionComponent>> s2igrowthcomponents;
    {
      // interface ID
      s2igrowthcomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));

      // kinetic models for scatra-scatra interface coupling involving interface layer growth
      std::vector<Teuchos::RCP<CondCompBundle>> kineticmodels;
      {
        {
          // Butler-Volmer
          std::vector<Teuchos::RCP<ConditionComponent>> butlervolmer;
          butlervolmer.push_back(Teuchos::rcp(
              new SeparatorConditionComponent("numscal")));  // total number of existing scalars
          std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp;
          intsepcomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("stoichiometries")));
          std::vector<Teuchos::RCP<IntVectorConditionComponent>>
              intvectcomp;  // string separator in front of integer stoichiometry vector in input
                            // file line
          intvectcomp.push_back(Teuchos::rcp(new IntVectorConditionComponent(
              "stoichiometries", 0)));  // integer vector of stoichiometric coefficients
          std::vector<Teuchos::RCP<SeparatorConditionComponent>>
              realsepcomp;  // empty vector --> no separators for real vectors needed
          std::vector<Teuchos::RCP<RealVectorConditionComponent>>
              realvectcomp;  // empty vector --> no real vectors needed
          butlervolmer.push_back(Teuchos::rcp(new IntRealBundle("stoichiometries",
              Teuchos::rcp(new IntConditionComponent("numscal")), intsepcomp, intvectcomp,
              realsepcomp, realvectcomp)));
          butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("e-")));
          butlervolmer.push_back(Teuchos::rcp(new IntConditionComponent("e-")));
          butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("k_r")));
          butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("k_r")));
          butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_a")));
          butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
          butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_c")));
          butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));
          butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("molmass")));
          butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("molar mass")));
          butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("density")));
          butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("density")));
          butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("conductivity")));
          butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("conductivity")));
          butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("regtype")));
          butlervolmer.push_back(
              Teuchos::rcp(new StringConditionComponent("regularization type", "trigonometrical",
                  Teuchos::tuple<std::string>("none", "polynomial", "Hein", "trigonometrical"),
                  Teuchos::tuple<int>(INPAR::S2I::regularization_none,
                      INPAR::S2I::regularization_polynomial, INPAR::S2I::regularization_hein,
                      INPAR::S2I::regularization_trigonometrical))));
          butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("regpar")));
          butlervolmer.push_back(
              Teuchos::rcp(new RealConditionComponent("regularization parameter")));
          butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("initthickness")));
          butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("initial thickness")));

          kineticmodels.push_back(Teuchos::rcp(new CondCompBundle(
              "Butler-Volmer", butlervolmer, INPAR::S2I::growth_kinetics_butlervolmer)));
        }
      }  // kinetic models for scatra-scatra interface coupling involving interface layer growth

      // insert kinetic models into vector with input file line components
      s2igrowthcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("KineticModel")));
      s2igrowthcomponents.push_back(Teuchos::rcp(new CondCompBundleSelector(
          "kinetic models for scatra-scatra interface coupling involving interface layer growth",
          Teuchos::rcp(new StringConditionComponent("kinetic model", "Butler-Volmer",
              Teuchos::tuple<std::string>("Butler-Volmer"),
              Teuchos::tuple<int>(INPAR::S2I::growth_kinetics_butlervolmer))),
          kineticmodels)));
    }

    // insert input file line components into condition definitions
    for (unsigned i = 0; i < s2igrowthcomponents.size(); ++i)
    {
      s2igrowthline->AddComponent(s2igrowthcomponents[i]);
      s2igrowthsurf->AddComponent(s2igrowthcomponents[i]);
    }

    // insert condition definitions into global list of valid condition definitions
    condlist.push_back(s2igrowthline);
    condlist.push_back(s2igrowthsurf);
  }

  return;
}
