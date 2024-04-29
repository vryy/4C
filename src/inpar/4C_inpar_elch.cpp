/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for electrochemistry

\level 2


*/
/*----------------------------------------------------------------------*/
#include "4C_inpar_elch.hpp"

#include "4C_inpar_scatra.hpp"
#include "4C_lib_conditiondefinition.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

void INPAR::ELCH::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& elchcontrol =
      list->sublist("ELCH CONTROL", false, "control parameters for electrochemistry problems\n");

  CORE::UTILS::IntParameter("MOVBOUNDARYITEMAX", 10,
      "Maximum number of outer iterations in electrode shape change computations", &elchcontrol);
  CORE::UTILS::DoubleParameter("MOVBOUNDARYCONVTOL", 1e-6,
      "Convergence check tolerance for outer loop in electrode shape change computations",
      &elchcontrol);
  CORE::UTILS::DoubleParameter("TEMPERATURE", 298.0, "Constant temperature (Kelvin)", &elchcontrol);
  CORE::UTILS::IntParameter("TEMPERATURE_FROM_FUNCT", -1,
      "Homogeneous temperature within electrochemistry field that can be time dependent according "
      "to function definition",
      &elchcontrol);
  CORE::UTILS::DoubleParameter("FARADAY_CONSTANT", 9.64853399e4,
      "Faraday constant (in unit system as chosen in input file)", &elchcontrol);
  CORE::UTILS::DoubleParameter("GAS_CONSTANT", 8.314472,
      "(universal) gas constant (in unit system as chosen in input file)", &elchcontrol);
  // parameter for possible types of ELCH algorithms for deforming meshes
  setStringToIntegralParameter<int>("MOVINGBOUNDARY", "No", "ELCH algorithm for deforming meshes",
      tuple<std::string>("No", "pseudo-transient", "fully-transient"),
      tuple<std::string>("no moving boundary algorithm",
          "pseudo-transient moving boundary algorithm",
          "full moving boundary algorithm including fluid solve"),
      tuple<int>(
          elch_mov_bndry_no, elch_mov_bndry_pseudo_transient, elch_mov_bndry_fully_transient),
      &elchcontrol);
  CORE::UTILS::DoubleParameter(
      "MOLARVOLUME", 0.0, "Molar volume for electrode shape change computations", &elchcontrol);
  CORE::UTILS::DoubleParameter("MOVBOUNDARYTHETA", 0.0,
      "One-step-theta factor in electrode shape change computations", &elchcontrol);
  CORE::UTILS::BoolParameter("GALVANOSTATIC", "No", "flag for galvanostatic mode", &elchcontrol);
  setStringToIntegralParameter<int>("GSTAT_APPROX_ELECT_RESIST", "relation_pot_cur",
      "relation of potential and current flow",
      tuple<std::string>("relation_pot_cur", "effective_length_with_initial_cond",
          "effective_length_with_integrated_cond"),
      tuple<int>(approxelctresist_relpotcur, approxelctresist_effleninitcond,
          approxelctresist_efflenintegcond),
      &elchcontrol);
  CORE::UTILS::IntParameter(
      "GSTATCONDID_CATHODE", 0, "condition id of electrode kinetics for cathode", &elchcontrol);
  CORE::UTILS::IntParameter(
      "GSTATCONDID_ANODE", 1, "condition id of electrode kinetics for anode", &elchcontrol);
  CORE::UTILS::DoubleParameter(
      "GSTATCONVTOL", 1.e-5, "Convergence check tolerance for galvanostatic mode", &elchcontrol);
  CORE::UTILS::DoubleParameter("GSTATCURTOL", 1.e-15, "Current Tolerance", &elchcontrol);
  CORE::UTILS::IntParameter(
      "GSTATFUNCTNO", -1, "function number defining the imposed current curve", &elchcontrol);
  CORE::UTILS::IntParameter(
      "GSTATITEMAX", 10, "maximum number of iterations for galvanostatic mode", &elchcontrol);
  CORE::UTILS::DoubleParameter(
      "GSTAT_LENGTH_CURRENTPATH", 0.0, "average length of the current path", &elchcontrol);

  setStringToIntegralParameter<int>("EQUPOT", "Undefined",
      "type of closing equation for electric potential",
      tuple<std::string>(
          "Undefined", "ENC", "ENC_PDE", "ENC_PDE_ELIM", "Poisson", "Laplace", "divi"),
      tuple<int>(equpot_undefined, equpot_enc, equpot_enc_pde, equpot_enc_pde_elim, equpot_poisson,
          equpot_laplace, equpot_divi),
      &elchcontrol);
  CORE::UTILS::BoolParameter("BLOCKPRECOND", "NO",
      "Switch to block-preconditioned family of solvers, only works with block preconditioners "
      "like CheapSIMPLE!",
      &elchcontrol);
  CORE::UTILS::BoolParameter(
      "DIFFCOND_FORMULATION", "No", "Activation of diffusion-conduction formulation", &elchcontrol);
  CORE::UTILS::BoolParameter("INITPOTCALC", "No",
      "Automatically calculate initial field for electric potential", &elchcontrol);
  CORE::UTILS::BoolParameter("ONLYPOTENTIAL", "no",
      "Coupling of general ion transport equation with Laplace equation", &elchcontrol);
  CORE::UTILS::BoolParameter("COUPLE_BOUNDARY_FLUXES", "Yes",
      "Coupling of lithium-ion flux density and electric current density at Dirichlet and Neumann "
      "boundaries",
      &elchcontrol);
  CORE::UTILS::DoubleParameter(
      "CYCLING_TIMESTEP", -1., "modified time step size for CCCV cell cycling", &elchcontrol);
  CORE::UTILS::BoolParameter("ELECTRODE_INFO_EVERY_STEP", "No",
      "the cell voltage, SOC, and C-Rate will be written to the csv file every step, even if "
      "RESULTSEVRY is not 1",
      &elchcontrol);

  /*----------------------------------------------------------------------*/
  // attention: this list is a sublist of elchcontrol
  Teuchos::ParameterList& elchdiffcondcontrol = elchcontrol.sublist(
      "DIFFCOND", false, "control parameters for electrochemical diffusion conduction problems\n");

  CORE::UTILS::BoolParameter(
      "CURRENT_SOLUTION_VAR", "No", "Current as a solution variable", &elchdiffcondcontrol);
  CORE::UTILS::BoolParameter("MAT_DIFFCOND_DIFFBASED", "Yes",
      "Coupling terms of chemical diffusion for current equation are based on t and kappa",
      &elchdiffcondcontrol);

  /// dilute solution theory (diffusion potential in current equation):
  ///    A          B
  ///   |--|  |----------|
  ///   z_1 + (z_2 - z_1) t_1
  /// ------------------------ (RT/F kappa (1+f+-) 1/c_k grad c_k)
  ///      z_1 z_2
  ///     |________|
  ///         C
  //
  // default: concentrated solution theory according to Newman
  CORE::UTILS::DoubleParameter("MAT_NEWMAN_CONST_A", 2.0,
      "Constant A for the Newman model(term for the concentration overpotential)",
      &elchdiffcondcontrol);
  CORE::UTILS::DoubleParameter("MAT_NEWMAN_CONST_B", -2.0,
      "Constant B for the Newman model(term for the concentration overpotential)",
      &elchdiffcondcontrol);
  CORE::UTILS::DoubleParameter("MAT_NEWMAN_CONST_C", -1.0,
      "Constant C for the Newman model(term for the concentration overpotential)",
      &elchdiffcondcontrol);
  CORE::UTILS::DoubleParameter(
      "PERMITTIVITY_VACUUM", 8.8541878128e-12, "Vacuum permittivity", &elchdiffcondcontrol);

  /*----------------------------------------------------------------------*/
  // sublist for space-charge layers
  auto& sclcontrol = elchcontrol.sublist(
      "SCL", false, "control parameters for coupled probelms with space-charge layer formation\n");

  CORE::UTILS::BoolParameter(
      "ADD_MICRO_MACRO_COUPLING", "No", "flag for micro macro coupling with scls", &sclcontrol);
  CORE::UTILS::BoolParameter("COUPLING_OUTPUT", "No",
      "write coupled node gids and node coordinates to csv file", &sclcontrol);
  CORE::UTILS::BoolParameter(
      "INITPOTCALC", "No", "calculate initial potential field?", &sclcontrol);
  CORE::UTILS::IntParameter("SOLVER", -1, "solver for coupled SCL problem", &sclcontrol);
  setStringToIntegralParameter<CORE::LINALG::MatrixType>("MATRIXTYPE", "undefined",
      "type of global system matrix in global system of equations",
      tuple<std::string>("undefined", "block", "sparse"),
      tuple<CORE::LINALG::MatrixType>(CORE::LINALG::MatrixType::undefined,
          CORE::LINALG::MatrixType::block_field, CORE::LINALG::MatrixType::sparse),
      &sclcontrol);
  CORE::UTILS::IntParameter("ADAPT_TIME_STEP", -1,
      "time step when time step size should be updated to 'ADAPTED_TIME_STEP_SIZE'.", &sclcontrol);
  CORE::UTILS::DoubleParameter("ADAPTED_TIME_STEP_SIZE", -1.0, "new time step size.", &sclcontrol);

  setStringToIntegralParameter<int>("INITIALFIELD", "zero_field",
      "Initial Field for scalar transport problem",
      tuple<std::string>("zero_field", "field_by_function", "field_by_condition"),
      tuple<int>(SCATRA::initfield_zero_field, SCATRA::initfield_field_by_function,
          SCATRA::initfield_field_by_condition),
      &sclcontrol);

  CORE::UTILS::IntParameter(
      "INITFUNCNO", -1, "function number for scalar transport initial field", &sclcontrol);
}


void INPAR::ELCH::SetValidConditions(
    std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist)
{
  using namespace INPUT;

  /*--------------------------------------------------------------------*/
  // electrode state of charge
  {
    // definition of electrode state of charge surface and volume conditions
    auto electrodesocline =
        Teuchos::rcp(new ConditionDefinition("DESIGN ELECTRODE STATE OF CHARGE LINE CONDITIONS",
            "ElectrodeSOC", "electrode state of charge line condition",
            CORE::Conditions::ElectrodeSOC, true, CORE::Conditions::geometry_type_line));

    auto electrodesocsurf =
        Teuchos::rcp(new ConditionDefinition("DESIGN ELECTRODE STATE OF CHARGE SURF CONDITIONS",
            "ElectrodeSOC", "electrode state of charge surface condition",
            CORE::Conditions::ElectrodeSOC, true, CORE::Conditions::geometry_type_surface));
    auto electrodesocvol =
        Teuchos::rcp(new ConditionDefinition("DESIGN ELECTRODE STATE OF CHARGE VOL CONDITIONS",
            "ElectrodeSOC", "electrode state of charge volume condition",
            CORE::Conditions::ElectrodeSOC, true, CORE::Conditions::geometry_type_volume));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<INPUT::LineComponent>> electrodesoccomponents;

    {
      electrodesoccomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ID")));
      electrodesoccomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent("ConditionID")));
      electrodesoccomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("C_0%")));
      electrodesoccomponents.emplace_back(Teuchos::rcp(new INPUT::RealComponent("c_0%")));
      electrodesoccomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("C_100%")));
      electrodesoccomponents.emplace_back(Teuchos::rcp(new INPUT::RealComponent("c_100%")));
      electrodesoccomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ONE_HOUR")));
      electrodesoccomponents.emplace_back(Teuchos::rcp(new INPUT::RealComponent("one_hour")));
    }

    // insert input file line components into condition definitions
    for (auto& electrodesoccomponent : electrodesoccomponents)
    {
      electrodesocline->AddComponent(electrodesoccomponent);
      electrodesocsurf->AddComponent(electrodesoccomponent);
      electrodesocvol->AddComponent(electrodesoccomponent);
    }

    // insert condition definitions into global list of valid condition definitions
    condlist.emplace_back(electrodesocline);
    condlist.emplace_back(electrodesocsurf);
    condlist.emplace_back(electrodesocvol);
  }

  /*--------------------------------------------------------------------*/
  // cell voltage
  {
    // definition of cell voltage point, line, and surface conditions
    auto cellvoltagepoint = Teuchos::rcp(new ConditionDefinition(
        "DESIGN CELL VOLTAGE POINT CONDITIONS", "CellVoltagePoint", "cell voltage point condition",
        CORE::Conditions::CellVoltage, false, CORE::Conditions::geometry_type_point));

    auto cellvoltageline = Teuchos::rcp(new ConditionDefinition(
        "DESIGN CELL VOLTAGE LINE CONDITIONS", "CellVoltage", "cell voltage line condition",
        CORE::Conditions::CellVoltage, true, CORE::Conditions::geometry_type_line));

    auto cellvoltagesurf = Teuchos::rcp(new ConditionDefinition(
        "DESIGN CELL VOLTAGE SURF CONDITIONS", "CellVoltage", "cell voltage surface condition",
        CORE::Conditions::CellVoltage, true, CORE::Conditions::geometry_type_surface));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<INPUT::LineComponent>> cellvoltagecomponents;

    {
      cellvoltagecomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ID")));
      cellvoltagecomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent("ConditionID")));
    }

    // insert input file line components into condition definitions
    for (auto& cellvoltagecomponent : cellvoltagecomponents)
    {
      cellvoltagepoint->AddComponent(cellvoltagecomponent);
      cellvoltageline->AddComponent(cellvoltagecomponent);
      cellvoltagesurf->AddComponent(cellvoltagecomponent);
    }

    // insert condition definitions into global list of valid condition definitions
    condlist.emplace_back(cellvoltagepoint);
    condlist.emplace_back(cellvoltageline);
    condlist.emplace_back(cellvoltagesurf);
  }

  /*--------------------------------------------------------------------*/
  // electrode kinetics as boundary condition on electrolyte
  {
    std::map<int, std::pair<std::string, std::vector<Teuchos::RCP<INPUT::LineComponent>>>>
        reaction_model_choices;

    // Butler-Volmer
    std::vector<Teuchos::RCP<INPUT::LineComponent>> butlervolmer;
    butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_A")));
    butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("alpha_a")));
    butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_C")));
    butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("alpha_c")));
    butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("I0")));
    butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("i0")));
    butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("GAMMA")));
    butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("gamma")));
    butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("REFCON")));
    butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("refcon")));
    butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("DL_SPEC_CAP")));
    butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("dl_spec_cap")));
    reaction_model_choices.emplace(
        butler_volmer, std::make_pair("Butler-Volmer", std::move(butlervolmer)));

    // Butler-Volmer Yang
    // parameter are identical to Butler-Volmer
    std::vector<Teuchos::RCP<INPUT::LineComponent>> butlervolmeryang;
    butlervolmeryang.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_A")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new INPUT::RealComponent("alpha_a")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_C")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new INPUT::RealComponent("alpha_c")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("I0")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new INPUT::RealComponent("i0")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("GAMMA")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new INPUT::RealComponent("gamma")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("REFCON")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new INPUT::RealComponent("refcon")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("DL_SPEC_CAP")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new INPUT::RealComponent("dl_spec_cap")));
    reaction_model_choices.emplace(INPAR::ELCH::butler_volmer_yang1997,
        std::make_pair("Butler-Volmer-Yang1997", butlervolmeryang));

    // Tafel kinetics
    std::vector<Teuchos::RCP<INPUT::LineComponent>> tafel;
    tafel.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA")));
    tafel.emplace_back(Teuchos::rcp(new INPUT::RealComponent("alpha")));
    tafel.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("I0")));
    tafel.emplace_back(Teuchos::rcp(new INPUT::RealComponent("i0")));
    tafel.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("GAMMA")));
    tafel.emplace_back(Teuchos::rcp(new INPUT::RealComponent("gamma")));
    tafel.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("REFCON")));
    tafel.emplace_back(Teuchos::rcp(new INPUT::RealComponent("refcon")));
    tafel.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("DL_SPEC_CAP")));
    tafel.emplace_back(Teuchos::rcp(new INPUT::RealComponent("dl_spec_cap")));
    reaction_model_choices.emplace(INPAR::ELCH::tafel, std::make_pair("Tafel", tafel));

    // linear kinetics
    std::vector<Teuchos::RCP<INPUT::LineComponent>> linear;
    linear.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA")));
    linear.emplace_back(Teuchos::rcp(new INPUT::RealComponent("alpha")));
    linear.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("I0")));
    linear.emplace_back(Teuchos::rcp(new INPUT::RealComponent("i0")));
    linear.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("GAMMA")));
    linear.emplace_back(Teuchos::rcp(new INPUT::RealComponent("gamma")));
    linear.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("REFCON")));
    linear.emplace_back(Teuchos::rcp(new INPUT::RealComponent("refcon")));
    linear.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("DL_SPEC_CAP")));
    linear.emplace_back(Teuchos::rcp(new INPUT::RealComponent("dl_spec_cap")));
    reaction_model_choices.emplace(INPAR::ELCH::linear, std::make_pair("linear", linear));

    // Butler-Volmer-Newman: "Newman (book), 2004, p. 213, eq. 8.26"
    //                       "Wittmann (Bachelor thesis), 2011, p. 15, eq. 2.30"
    std::vector<Teuchos::RCP<INPUT::LineComponent>> bvnewman;
    bvnewman.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("K_A")));
    bvnewman.emplace_back(Teuchos::rcp(new INPUT::RealComponent("k_a")));
    bvnewman.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("K_C")));
    bvnewman.emplace_back(Teuchos::rcp(new INPUT::RealComponent("k_c")));
    bvnewman.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("BETA")));
    bvnewman.emplace_back(Teuchos::rcp(new INPUT::RealComponent("beta")));
    bvnewman.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("DL_SPEC_CAP")));
    bvnewman.emplace_back(Teuchos::rcp(new INPUT::RealComponent("dl_spec_cap")));
    reaction_model_choices.emplace(
        INPAR::ELCH::butler_volmer_newman, std::make_pair("Butler-Volmer-Newman", bvnewman));

    // Butler-Volmer-Newman: "Bard (book), 2001, p. 99, eq. 3.4.10"
    //                       "Wittmann (Bachelor thesis), 2011, p. 16, eq. 2.32"
    std::vector<Teuchos::RCP<INPUT::LineComponent>> bvbard;
    bvbard.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("E0")));
    bvbard.emplace_back(Teuchos::rcp(new INPUT::RealComponent("e0")));
    bvbard.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("K0")));
    bvbard.emplace_back(Teuchos::rcp(new INPUT::RealComponent("k0")));
    bvbard.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("BETA")));
    bvbard.emplace_back(Teuchos::rcp(new INPUT::RealComponent("beta")));
    bvbard.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("C_C0")));
    bvbard.emplace_back(Teuchos::rcp(new INPUT::RealComponent("c_c0")));
    bvbard.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("C_A0")));
    bvbard.emplace_back(Teuchos::rcp(new INPUT::RealComponent("c_a0")));
    bvbard.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("DL_SPEC_CAP")));
    bvbard.emplace_back(Teuchos::rcp(new INPUT::RealComponent("dl_spec_cap")));
    reaction_model_choices.emplace(
        INPAR::ELCH::butler_volmer_bard, std::make_pair("Butler-Volmer-Bard", bvbard));

    // Nernst equation:
    std::vector<Teuchos::RCP<INPUT::LineComponent>> nernst;
    nernst.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("E0")));
    nernst.emplace_back(Teuchos::rcp(new INPUT::RealComponent("e0")));
    nernst.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("C0")));
    nernst.emplace_back(Teuchos::rcp(new INPUT::RealComponent("c0")));
    nernst.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("DL_SPEC_CAP")));
    nernst.emplace_back(Teuchos::rcp(new INPUT::RealComponent("dl_spec_cap")));
    reaction_model_choices.emplace(INPAR::ELCH::nernst, std::make_pair("Nernst", nernst));

    std::vector<Teuchos::RCP<INPUT::LineComponent>> elechemcomponents;
    elechemcomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ID")));
    elechemcomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent("ConditionID")));
    elechemcomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("POT")));
    elechemcomponents.emplace_back(Teuchos::rcp(new INPUT::RealComponent("pot")));
    elechemcomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("FUNCT")));
    elechemcomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent("funct", {0, true, true})));
    elechemcomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("NUMSCAL")));
    elechemcomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent("numscal")));
    elechemcomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("STOICH")));
    elechemcomponents.emplace_back(
        Teuchos::rcp(new INPUT::IntVectorComponent("stoich", INPUT::LengthFromInt("numscal"))));

    elechemcomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("E-")));
    elechemcomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent("e-")));
    // porosity of electrode boundary, set to -1 if equal to porosity of electrolyte domain
    elechemcomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("EPSILON")));
    elechemcomponents.emplace_back(Teuchos::rcp(new INPUT::RealComponent("epsilon")));
    elechemcomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ZERO_CUR")));
    elechemcomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent("zero_cur")));
    elechemcomponents.emplace_back(Teuchos::rcp(
        new INPUT::SwitchComponent("kinetic model", butler_volmer, reaction_model_choices)));

    auto electrodeboundarykineticspoint =
        Teuchos::rcp(new ConditionDefinition("ELECTRODE BOUNDARY KINETICS POINT CONDITIONS",
            "ElchBoundaryKineticsPoint", "point electrode boundary kinetics",
            CORE::Conditions::ElchBoundaryKinetics, false, CORE::Conditions::geometry_type_point));

    auto electrodeboundarykineticsline =
        Teuchos::rcp(new ConditionDefinition("ELECTRODE BOUNDARY KINETICS LINE CONDITIONS",
            "ElchBoundaryKinetics", "line electrode boundary kinetics",
            CORE::Conditions::ElchBoundaryKinetics, true, CORE::Conditions::geometry_type_line));

    auto electrodeboundarykineticssurf =
        Teuchos::rcp(new ConditionDefinition("ELECTRODE BOUNDARY KINETICS SURF CONDITIONS",
            "ElchBoundaryKinetics", "surface electrode boundary kinetics",
            CORE::Conditions::ElchBoundaryKinetics, true, CORE::Conditions::geometry_type_surface));

    for (auto& elechemcomponent : elechemcomponents)
    {
      electrodeboundarykineticspoint->AddComponent(elechemcomponent);
      electrodeboundarykineticsline->AddComponent(elechemcomponent);
      electrodeboundarykineticssurf->AddComponent(elechemcomponent);
    }

    condlist.emplace_back(electrodeboundarykineticspoint);
    condlist.emplace_back(electrodeboundarykineticsline);
    condlist.emplace_back(electrodeboundarykineticssurf);
  }

  /*--------------------------------------------------------------------*/
  // electrode kinetics as domain condition within electrolyte
  {
    // definition of line, surface, and volume conditions for electrode domain kinetics
    auto electrodedomainkineticsline =
        Teuchos::rcp(new ConditionDefinition("ELECTRODE DOMAIN KINETICS LINE CONDITIONS",
            "ElchDomainKinetics", "line electrode domain kinetics",
            CORE::Conditions::ElchDomainKinetics, true, CORE::Conditions::geometry_type_line));

    auto electrodedomainkineticssurf =
        Teuchos::rcp(new ConditionDefinition("ELECTRODE DOMAIN KINETICS SURF CONDITIONS",
            "ElchDomainKinetics", "surface electrode domain kinetics",
            CORE::Conditions::ElchDomainKinetics, true, CORE::Conditions::geometry_type_surface));

    auto electrodedomainkineticsvol =
        Teuchos::rcp(new ConditionDefinition("ELECTRODE DOMAIN KINETICS VOL CONDITIONS",
            "ElchDomainKinetics", "volume electrode domain kinetics",
            CORE::Conditions::ElchDomainKinetics, true, CORE::Conditions::geometry_type_volume));

    // equip condition definition with input file line components
    std::vector<Teuchos::RCP<INPUT::LineComponent>> electrodedomainkineticscomponents;

    {
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("ID")));
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new INPUT::IntComponent("ConditionID")));
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("POT")));
      electrodedomainkineticscomponents.emplace_back(Teuchos::rcp(new INPUT::RealComponent("pot")));
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("FUNCT")));
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new INPUT::IntComponent("funct", {0, true, true})));
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("NUMSCAL")));

      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new INPUT::IntComponent("numscal")));
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("STOICH")));
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new INPUT::IntVectorComponent("stoich", INPUT::LengthFromInt("numscal"))));

      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("E-")));
      electrodedomainkineticscomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent("e-")));
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("ZERO_CUR")));
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new INPUT::IntComponent("zero_cur")));


      {
        // Butler-Volmer
        std::vector<Teuchos::RCP<INPUT::LineComponent>> butlervolmer;
        butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent(
            "A_S")));  // ratio of electrode-electrolyte interface area to total two-phase volume
        butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("A_s")));
        butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_A")));
        butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("alpha_a")));
        butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_C")));
        butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("alpha_c")));
        butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("I0")));
        butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("i0")));
        butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("GAMMA")));
        butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("gamma")));
        butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("REFCON")));
        butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("refcon")));
        butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("DL_SPEC_CAP")));
        butlervolmer.emplace_back(Teuchos::rcp(new INPUT::RealComponent("dl_spec_cap")));
        butlervolmer.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("END")));

        electrodedomainkineticscomponents.emplace_back(
            Teuchos::rcp(new INPUT::SwitchComponent("kinetic model", butler_volmer,
                {{butler_volmer, std::make_pair("Butler-Volmer", std::move(butlervolmer))}})));
      }
    }

    // insert input file line components into condition definitions
    for (auto& electrodedomainkineticscomponent : electrodedomainkineticscomponents)
    {
      electrodedomainkineticsline->AddComponent(electrodedomainkineticscomponent);
      electrodedomainkineticssurf->AddComponent(electrodedomainkineticscomponent);
      electrodedomainkineticsvol->AddComponent(electrodedomainkineticscomponent);
    }

    // insert condition definitions into global list of valid condition definitions
    condlist.emplace_back(electrodedomainkineticsline);
    condlist.emplace_back(electrodedomainkineticssurf);
    condlist.emplace_back(electrodedomainkineticsvol);
  }

  /*--------------------------------------------------------------------*/
  // boundary condition for constant-current constant-voltage (CCCV) cell cycling
  {
    // definition of point, line and surface conditions for CCCV cell cycling
    auto cccvcyclingpoint = Teuchos::rcp(
        new ConditionDefinition("DESIGN CCCV CELL CYCLING POINT CONDITIONS", "CCCVCycling",
            "line boundary condition for constant-current constant-voltage (CCCV) cell cycling",
            CORE::Conditions::CCCVCycling, true, CORE::Conditions::geometry_type_point));

    auto cccvcyclingline = Teuchos::rcp(
        new ConditionDefinition("DESIGN CCCV CELL CYCLING LINE CONDITIONS", "CCCVCycling",
            "line boundary condition for constant-current constant-voltage (CCCV) cell cycling",
            CORE::Conditions::CCCVCycling, true, CORE::Conditions::geometry_type_line));

    auto cccvcyclingsurf = Teuchos::rcp(
        new ConditionDefinition("DESIGN CCCV CELL CYCLING SURF CONDITIONS", "CCCVCycling",
            "surface boundary condition for constant-current constant-voltage (CCCV) cell cycling",
            CORE::Conditions::CCCVCycling, true, CORE::Conditions::geometry_type_surface));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<INPUT::LineComponent>> cccvcyclingcomponents;

    {
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("NUMBER_OF_HALF_CYCLES")));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new INPUT::IntComponent("NumberOfHalfCycles")));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("BEGIN_WITH_CHARGING")));
      cccvcyclingcomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent(
          "BeginWithCharging")));  // Boolean parameter represented by integer parameter
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("CONDITION_ID_FOR_CHARGE")));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new INPUT::IntComponent("ConditionIDForCharge", {0, false, true})));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("CONDITION_ID_FOR_DISCHARGE")));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new INPUT::IntComponent("ConditionIDForDischarge", {0, false, true})));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("INIT_RELAX_TIME")));
      cccvcyclingcomponents.emplace_back(Teuchos::rcp(new INPUT::RealComponent("InitRelaxTime")));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("ADAPTIVE_TIME_STEPPING_INIT_RELAX")));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new INPUT::IntComponent("AdaptiveTimeSteppingInitRelax")));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("NUM_ADD_ADAPT_TIME_STEPS")));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new INPUT::IntComponent("NumAddAdaptTimeSteps", {0, false, true})));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("MIN_TIME_STEPS_DURING_INIT_RELAX")));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new INPUT::IntComponent("MinTimeStepsDuringInitRelax", {0, false, true})));
    }

    // insert input file line components into condition definitions
    for (auto& cccvcyclingcomponent : cccvcyclingcomponents)
    {
      cccvcyclingpoint->AddComponent(cccvcyclingcomponent);
      cccvcyclingline->AddComponent(cccvcyclingcomponent);
      cccvcyclingsurf->AddComponent(cccvcyclingcomponent);
    }

    // insert condition definitions into global list of valid condition definitions
    condlist.emplace_back(cccvcyclingpoint);
    condlist.emplace_back(cccvcyclingline);
    condlist.emplace_back(cccvcyclingsurf);
  }

  /*--------------------------------------------------------------------*/
  // boundary condition for constant-current constant-voltage (CCCV) half-cycle
  {
    // definition of point, line and surface conditions for CCCV half-cycle
    auto cccvhalfcyclepoint = Teuchos::rcp(
        new ConditionDefinition("DESIGN CCCV HALF-CYCLE POINT CONDITIONS", "CCCVHalfCycle",
            "line boundary condition for constant-current constant-voltage (CCCV) half-cycle",
            CORE::Conditions::CCCVHalfCycle, true, CORE::Conditions::geometry_type_point));

    auto cccvhalfcycleline = Teuchos::rcp(
        new ConditionDefinition("DESIGN CCCV HALF-CYCLE LINE CONDITIONS", "CCCVHalfCycle",
            "line boundary condition for constant-current constant-voltage (CCCV) half-cycle",
            CORE::Conditions::CCCVHalfCycle, true, CORE::Conditions::geometry_type_line));

    auto cccvhalfcyclesurf = Teuchos::rcp(
        new ConditionDefinition("DESIGN CCCV HALF-CYCLE SURF CONDITIONS", "CCCVHalfCycle",
            "surface boundary condition for constant-current constant-voltage (CCCV) half-cycle",
            CORE::Conditions::CCCVHalfCycle, true, CORE::Conditions::geometry_type_surface));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<INPUT::LineComponent>> cccvhalfcyclecomponents;

    {
      cccvhalfcyclecomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ID")));
      cccvhalfcyclecomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent("ConditionID")));
      cccvhalfcyclecomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("CURRENT")));
      cccvhalfcyclecomponents.emplace_back(Teuchos::rcp(new INPUT::RealComponent("Current")));
      cccvhalfcyclecomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("CUT_OFF_VOLTAGE")));
      cccvhalfcyclecomponents.emplace_back(Teuchos::rcp(new INPUT::RealComponent("CutoffVoltage")));
      cccvhalfcyclecomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("CUT_OFF_C_RATE")));
      cccvhalfcyclecomponents.emplace_back(Teuchos::rcp(new INPUT::RealComponent("CutoffCRate")));
      cccvhalfcyclecomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("RELAX_TIME")));
      cccvhalfcyclecomponents.emplace_back(Teuchos::rcp(new INPUT::RealComponent("RelaxTime")));
      // switch adaptive time stepping on for different phases of half cycle: 1st: end of constant
      // current, 2nd: end of constant voltage, 3rd: end of relaxation
      cccvhalfcyclecomponents.emplace_back(
          Teuchos::rcp(new INPUT::SeparatorComponent("ADAPTIVE_TIME_STEPPING_PHASE_ON_OFF")));
      cccvhalfcyclecomponents.emplace_back(
          Teuchos::rcp(new INPUT::IntVectorComponent("AdaptiveTimeSteppingPhaseOnOff", 3)));
    }

    // insert input file line components into condition definitions
    for (auto& cccvhalfcyclecomponent : cccvhalfcyclecomponents)
    {
      cccvhalfcyclepoint->AddComponent(cccvhalfcyclecomponent);
      cccvhalfcycleline->AddComponent(cccvhalfcyclecomponent);
      cccvhalfcyclesurf->AddComponent(cccvhalfcyclecomponent);
    }

    // insert condition definitions into global list of valid condition definitions
    condlist.emplace_back(cccvhalfcyclepoint);
    condlist.emplace_back(cccvhalfcycleline);
    condlist.emplace_back(cccvhalfcyclesurf);
  }
}

FOUR_C_NAMESPACE_CLOSE
