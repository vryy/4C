/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for electrochemistry

\level 2


*/
/*----------------------------------------------------------------------*/
#include "inpar_elch.H"

#include "drt_validparameters.H"

#include "../drt_lib/drt_conditiondefinition.H"

void INPAR::ELCH::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& elchcontrol =
      list->sublist("ELCH CONTROL", false, "control parameters for electrochemistry problems\n");

  IntParameter("MOVBOUNDARYITEMAX", 10,
      "Maximum number of outer iterations in electrode shape change computations", &elchcontrol);
  DoubleParameter("MOVBOUNDARYCONVTOL", 1e-6,
      "Convergence check tolerance for outer loop in electrode shape change computations",
      &elchcontrol);
  DoubleParameter("TEMPERATURE", 298.0, "Constant temperature (Kelvin)", &elchcontrol);
  IntParameter("TEMPERATURE_FROM_FUNCT", -1,
      "Homogeneous temperature within electrochemistry field that can be time dependent according "
      "to function definition",
      &elchcontrol);
  DoubleParameter("FARADAY_CONSTANT", 9.64853399e4,
      "Faraday constant (in unit system as chosen in input file)", &elchcontrol);
  DoubleParameter("GAS_CONSTANT", 8.314472,
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
  DoubleParameter(
      "MOLARVOLUME", 0.0, "Molar volume for electrode shape change computations", &elchcontrol);
  DoubleParameter("MOVBOUNDARYTHETA", 0.0,
      "One-step-theta factor in electrode shape change computations", &elchcontrol);
  BoolParameter("GALVANOSTATIC", "No", "flag for galvanostatic mode", &elchcontrol);
  setStringToIntegralParameter<int>("GSTAT_APPROX_ELECT_RESIST", "relation_pot_cur",
      "relation of potential and current flow",
      tuple<std::string>("relation_pot_cur", "effective_length_with_initial_cond",
          "effective_length_with_integrated_cond"),
      tuple<int>(approxelctresist_relpotcur, approxelctresist_effleninitcond,
          approxelctresist_efflenintegcond),
      &elchcontrol);
  IntParameter(
      "GSTATCONDID_CATHODE", 0, "condition id of electrode kinetics for cathode", &elchcontrol);
  IntParameter(
      "GSTATCONDID_ANODE", 1, "condition id of electrode kinetics for anode", &elchcontrol);
  DoubleParameter(
      "GSTATCONVTOL", 1.e-5, "Convergence check tolerance for galvanostatic mode", &elchcontrol);
  DoubleParameter("GSTATCURTOL", 1.e-15, "Current Tolerance", &elchcontrol);
  IntParameter(
      "GSTATFUNCTNO", -1, "function number defining the imposed current curve", &elchcontrol);
  IntParameter(
      "GSTATITEMAX", 10, "maximum number of iterations for galvanostatic mode", &elchcontrol);
  DoubleParameter(
      "GSTAT_LENGTH_CURRENTPATH", 0.0, "average length of the current path", &elchcontrol);

  setStringToIntegralParameter<int>("EQUPOT", "Undefined",
      "type of closing equation for electric potential",
      tuple<std::string>(
          "Undefined", "ENC", "ENC_PDE", "ENC_PDE_ELIM", "Poisson", "Laplace", "divi"),
      tuple<int>(equpot_undefined, equpot_enc, equpot_enc_pde, equpot_enc_pde_elim, equpot_poisson,
          equpot_laplace, equpot_divi),
      &elchcontrol);
  BoolParameter("BLOCKPRECOND", "NO",
      "Switch to block-preconditioned family of solvers, only works with block preconditioners "
      "like CheapSIMPLE!",
      &elchcontrol);
  BoolParameter(
      "DIFFCOND_FORMULATION", "No", "Activation of diffusion-conduction formulation", &elchcontrol);
  BoolParameter("INITPOTCALC", "No", "Automatically calculate initial field for electric potential",
      &elchcontrol);
  BoolParameter("ONLYPOTENTIAL", "no",
      "Coupling of general ion transport equation with Laplace equation", &elchcontrol);
  BoolParameter("COUPLE_BOUNDARY_FLUXES", "Yes",
      "Coupling of lithium-ion flux density and electric current density at Dirichlet and Neumann "
      "boundaries",
      &elchcontrol);
  DoubleParameter(
      "CYCLING_TIMESTEP", -1., "modified time step size for CCCV cell cycling", &elchcontrol);

  /*----------------------------------------------------------------------*/
  // attention: this list is a sublist of elchcontrol
  Teuchos::ParameterList& elchdiffcondcontrol = elchcontrol.sublist(
      "DIFFCOND", false, "control parameters for electrochemical diffusion conduction problems\n");

  BoolParameter(
      "CURRENT_SOLUTION_VAR", "No", "Current as a solution variable", &elchdiffcondcontrol);
  BoolParameter("MAT_DIFFCOND_DIFFBASED", "Yes",
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
  DoubleParameter("MAT_NEWMAN_CONST_A", 2.0,
      "Constant A for the Newman model(term for the concentration overpotential)",
      &elchdiffcondcontrol);
  DoubleParameter("MAT_NEWMAN_CONST_B", -2.0,
      "Constant B for the Newman model(term for the concentration overpotential)",
      &elchdiffcondcontrol);
  DoubleParameter("MAT_NEWMAN_CONST_C", -1.0,
      "Constant C for the Newman model(term for the concentration overpotential)",
      &elchdiffcondcontrol);
}


void INPAR::ELCH::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // electrode state of charge
  {
    // definition of electrode state of charge surface and volume conditions
    Teuchos::RCP<ConditionDefinition> electrodesocline =
        Teuchos::rcp(new ConditionDefinition("DESIGN ELECTRODE STATE OF CHARGE LINE CONDITIONS",
            "ElectrodeSOC", "electrode state of charge line condition",
            DRT::Condition::ElectrodeSOC, true, DRT::Condition::Line));

    Teuchos::RCP<ConditionDefinition> electrodesocsurf =
        Teuchos::rcp(new ConditionDefinition("DESIGN ELECTRODE STATE OF CHARGE SURF CONDITIONS",
            "ElectrodeSOC", "electrode state of charge surface condition",
            DRT::Condition::ElectrodeSOC, true, DRT::Condition::Surface));
    Teuchos::RCP<ConditionDefinition> electrodesocvol =
        Teuchos::rcp(new ConditionDefinition("DESIGN ELECTRODE STATE OF CHARGE VOL CONDITIONS",
            "ElectrodeSOC", "electrode state of charge volume condition",
            DRT::Condition::ElectrodeSOC, true, DRT::Condition::Volume));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<ConditionComponent>> electrodesoccomponents;

    {
      electrodesoccomponents.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("ID")));
      electrodesoccomponents.emplace_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));
      electrodesoccomponents.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("c_0%")));
      electrodesoccomponents.emplace_back(Teuchos::rcp(new RealConditionComponent("c_0%")));
      electrodesoccomponents.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("c_100%")));
      electrodesoccomponents.emplace_back(Teuchos::rcp(new RealConditionComponent("c_100%")));
      electrodesoccomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("one_hour")));
      electrodesoccomponents.emplace_back(Teuchos::rcp(new RealConditionComponent("one_hour")));
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
    Teuchos::RCP<ConditionDefinition> cellvoltagepoint = Teuchos::rcp(new ConditionDefinition(
        "DESIGN CELL VOLTAGE POINT CONDITIONS", "CellVoltagePoint", "cell voltage point condition",
        DRT::Condition::CellVoltage, false, DRT::Condition::Point));

    Teuchos::RCP<ConditionDefinition> cellvoltageline = Teuchos::rcp(new ConditionDefinition(
        "DESIGN CELL VOLTAGE LINE CONDITIONS", "CellVoltage", "cell voltage line condition",
        DRT::Condition::CellVoltage, true, DRT::Condition::Line));

    Teuchos::RCP<ConditionDefinition> cellvoltagesurf = Teuchos::rcp(new ConditionDefinition(
        "DESIGN CELL VOLTAGE SURF CONDITIONS", "CellVoltage", "cell voltage surface condition",
        DRT::Condition::CellVoltage, true, DRT::Condition::Surface));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<ConditionComponent>> cellvoltagecomponents;

    {
      cellvoltagecomponents.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("ID")));
      cellvoltagecomponents.emplace_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));
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
    std::vector<Teuchos::RCP<CondCompBundle>> reactionmodel;

    // Butler-Volmer
    std::vector<Teuchos::RCP<ConditionComponent>> butlervolmer;
    butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_a")));
    butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
    butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_c")));
    butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));
    butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("i0")));
    butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("i0")));
    butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("gamma")));
    butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("gamma")));
    butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("refcon")));
    butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("refcon")));
    butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("dl_spec_cap")));
    butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("dl_spec_cap")));
    butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("END")));
    reactionmodel.emplace_back(Teuchos::rcp(
        new CondCompBundle("Butler-Volmer", butlervolmer, INPAR::ELCH::butler_volmer)));

    // Butler-Volmer Yang
    // parameter are identical to Butler-Volmer
    std::vector<Teuchos::RCP<ConditionComponent>> butlervolmeryang;
    butlervolmeryang.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_a")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_c")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("i0")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new RealConditionComponent("i0")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("gamma")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new RealConditionComponent("gamma")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("refcon")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new RealConditionComponent("refcon")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("dl_spec_cap")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new RealConditionComponent("dl_spec_cap")));
    butlervolmeryang.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("END")));
    reactionmodel.emplace_back(Teuchos::rcp(new CondCompBundle(
        "Butler-Volmer-Yang1997", butlervolmeryang, INPAR::ELCH::butler_volmer_yang1997)));

    // Tafel kinetics
    std::vector<Teuchos::RCP<ConditionComponent>> tafel;
    tafel.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("alpha")));
    tafel.emplace_back(Teuchos::rcp(new RealConditionComponent("alpha")));
    tafel.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("i0")));
    tafel.emplace_back(Teuchos::rcp(new RealConditionComponent("i0")));
    tafel.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("gamma")));
    tafel.emplace_back(Teuchos::rcp(new RealConditionComponent("gamma")));
    tafel.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("refcon")));
    tafel.emplace_back(Teuchos::rcp(new RealConditionComponent("refcon")));
    tafel.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("dl_spec_cap")));
    tafel.emplace_back(Teuchos::rcp(new RealConditionComponent("dl_spec_cap")));
    reactionmodel.emplace_back(
        Teuchos::rcp(new CondCompBundle("Tafel", tafel, INPAR::ELCH::tafel)));

    // linear kinetics
    std::vector<Teuchos::RCP<ConditionComponent>> linear;
    linear.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("alpha")));
    linear.emplace_back(Teuchos::rcp(new RealConditionComponent("alpha")));
    linear.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("i0")));
    linear.emplace_back(Teuchos::rcp(new RealConditionComponent("i0")));
    linear.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("gamma")));
    linear.emplace_back(Teuchos::rcp(new RealConditionComponent("gamma")));
    linear.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("refcon")));
    linear.emplace_back(Teuchos::rcp(new RealConditionComponent("refcon")));
    linear.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("dl_spec_cap")));
    linear.emplace_back(Teuchos::rcp(new RealConditionComponent("dl_spec_cap")));
    linear.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("END")));
    reactionmodel.emplace_back(
        Teuchos::rcp(new CondCompBundle("linear", linear, INPAR::ELCH::linear)));

    // Butler-Volmer-Newman: "Newman (book), 2004, p. 213, eq. 8.26"
    //                       "Wittmann (Bachelor thesis), 2011, p. 15, eq. 2.30"
    std::vector<Teuchos::RCP<ConditionComponent>> bvnewman;
    bvnewman.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("k_a")));
    bvnewman.emplace_back(Teuchos::rcp(new RealConditionComponent("k_a")));
    bvnewman.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("k_c")));
    bvnewman.emplace_back(Teuchos::rcp(new RealConditionComponent("k_c")));
    bvnewman.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("beta")));
    bvnewman.emplace_back(Teuchos::rcp(new RealConditionComponent("beta")));
    bvnewman.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("dl_spec_cap")));
    bvnewman.emplace_back(Teuchos::rcp(new RealConditionComponent("dl_spec_cap")));
    bvnewman.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("END")));
    reactionmodel.emplace_back(Teuchos::rcp(
        new CondCompBundle("Butler-Volmer-Newman", bvnewman, INPAR::ELCH::butler_volmer_newman)));

    // Butler-Volmer-Newman: "Bard (book), 2001, p. 99, eq. 3.4.10"
    //                       "Wittmann (Bachelor thesis), 2011, p. 16, eq. 2.32"
    std::vector<Teuchos::RCP<ConditionComponent>> bvbard;
    bvbard.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("e0")));
    bvbard.emplace_back(Teuchos::rcp(new RealConditionComponent("e0")));
    bvbard.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("k0")));
    bvbard.emplace_back(Teuchos::rcp(new RealConditionComponent("k0")));
    bvbard.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("beta")));
    bvbard.emplace_back(Teuchos::rcp(new RealConditionComponent("beta")));
    bvbard.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("c_c0")));
    bvbard.emplace_back(Teuchos::rcp(new RealConditionComponent("c_c0")));
    bvbard.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("c_a0")));
    bvbard.emplace_back(Teuchos::rcp(new RealConditionComponent("c_a0")));
    bvbard.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("dl_spec_cap")));
    bvbard.emplace_back(Teuchos::rcp(new RealConditionComponent("dl_spec_cap")));
    bvbard.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("END")));
    reactionmodel.emplace_back(Teuchos::rcp(
        new CondCompBundle("Butler-Volmer-Bard", bvbard, INPAR::ELCH::butler_volmer_bard)));

    // Nernst equation:
    std::vector<Teuchos::RCP<ConditionComponent>> nernst;
    nernst.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("e0")));
    nernst.emplace_back(Teuchos::rcp(new RealConditionComponent("e0")));
    nernst.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("c0")));
    nernst.emplace_back(Teuchos::rcp(new RealConditionComponent("c0")));
    nernst.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("dl_spec_cap")));
    nernst.emplace_back(Teuchos::rcp(new RealConditionComponent("dl_spec_cap")));
    reactionmodel.emplace_back(
        Teuchos::rcp(new CondCompBundle("Nernst", nernst, INPAR::ELCH::nernst)));

    // input: stoichiometry for reaction mechanism (IntRealBundle)
    // definition separator for int vectors
    std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepveccomp;
    intsepveccomp.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("stoich")));

    // definition int vectors
    std::vector<Teuchos::RCP<IntVectorConditionComponent>> intveccomp;
    intveccomp.emplace_back(Teuchos::rcp(new IntVectorConditionComponent("stoich", 2)));

    // definition separator for real vectors: length of the real vector is zero -> nothing is read
    std::vector<Teuchos::RCP<SeparatorConditionComponent>> realsepveccomp;

    // definition real vectors: length of the real vector is zero -> nothing is read
    std::vector<Teuchos::RCP<RealVectorConditionComponent>> realveccomp;


    std::vector<Teuchos::RCP<ConditionComponent>> elechemcomponents;
    elechemcomponents.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("ID")));
    elechemcomponents.emplace_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));
    elechemcomponents.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("pot")));
    elechemcomponents.emplace_back(Teuchos::rcp(new RealConditionComponent("pot")));
    elechemcomponents.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("funct")));
    elechemcomponents.emplace_back(Teuchos::rcp(new IntConditionComponent("funct", true, true)));
    elechemcomponents.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("numscal")));
    elechemcomponents.emplace_back(Teuchos::rcp(
        new IntRealBundle("intreal bundle", Teuchos::rcp(new IntConditionComponent("numscal")),
            intsepveccomp, intveccomp, realsepveccomp, realveccomp)));
    elechemcomponents.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("e-")));
    elechemcomponents.emplace_back(Teuchos::rcp(new IntConditionComponent("e-")));
    // porosity of electrode boundary, set to -1 if equal to porosity of electrolyte domain
    elechemcomponents.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("epsilon")));
    elechemcomponents.emplace_back(Teuchos::rcp(new RealConditionComponent("epsilon")));
    elechemcomponents.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("zero_cur")));
    elechemcomponents.emplace_back(Teuchos::rcp(new IntConditionComponent("zero_cur")));
    elechemcomponents.emplace_back(Teuchos::rcp(new CondCompBundleSelector("kinetic model bundle",
        Teuchos::rcp(new StringConditionComponent("kinetic model", "Butler-Volmer",
            Teuchos::tuple<std::string>("Butler-Volmer", "Butler-Volmer-Yang1997", "Tafel",
                "linear", "Butler-Volmer-Newman", "Butler-Volmer-Bard", "Nernst", "zero"),
            Teuchos::tuple<int>(INPAR::ELCH::butler_volmer, INPAR::ELCH::butler_volmer_yang1997,
                INPAR::ELCH::tafel, INPAR::ELCH::linear, INPAR::ELCH::butler_volmer_newman,
                INPAR::ELCH::butler_volmer_bard, INPAR::ELCH::nernst, INPAR::ELCH::zero))),
        reactionmodel)));

    Teuchos::RCP<ConditionDefinition> electrodeboundarykineticspoint =
        Teuchos::rcp(new ConditionDefinition("ELECTRODE BOUNDARY KINETICS POINT CONDITIONS",
            "ElchBoundaryKineticsPoint", "point electrode boundary kinetics",
            DRT::Condition::ElchBoundaryKinetics, false, DRT::Condition::Point));

    Teuchos::RCP<ConditionDefinition> electrodeboundarykineticsline =
        Teuchos::rcp(new ConditionDefinition("ELECTRODE BOUNDARY KINETICS LINE CONDITIONS",
            "ElchBoundaryKinetics", "line electrode boundary kinetics",
            DRT::Condition::ElchBoundaryKinetics, true, DRT::Condition::Line));

    Teuchos::RCP<ConditionDefinition> electrodeboundarykineticssurf =
        Teuchos::rcp(new ConditionDefinition("ELECTRODE BOUNDARY KINETICS SURF CONDITIONS",
            "ElchBoundaryKinetics", "surface electrode boundary kinetics",
            DRT::Condition::ElchBoundaryKinetics, true, DRT::Condition::Surface));

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
    Teuchos::RCP<ConditionDefinition> electrodedomainkineticsline =
        Teuchos::rcp(new ConditionDefinition("ELECTRODE DOMAIN KINETICS LINE CONDITIONS",
            "ElchDomainKinetics", "line electrode domain kinetics",
            DRT::Condition::ElchDomainKinetics, true, DRT::Condition::Line));

    Teuchos::RCP<ConditionDefinition> electrodedomainkineticssurf =
        Teuchos::rcp(new ConditionDefinition("ELECTRODE DOMAIN KINETICS SURF CONDITIONS",
            "ElchDomainKinetics", "surface electrode domain kinetics",
            DRT::Condition::ElchDomainKinetics, true, DRT::Condition::Surface));

    Teuchos::RCP<ConditionDefinition> electrodedomainkineticsvol =
        Teuchos::rcp(new ConditionDefinition("ELECTRODE DOMAIN KINETICS VOL CONDITIONS",
            "ElchDomainKinetics", "volume electrode domain kinetics",
            DRT::Condition::ElchDomainKinetics, true, DRT::Condition::Volume));

    // equip condition definition with input file line components
    std::vector<Teuchos::RCP<ConditionComponent>> electrodedomainkineticscomponents;

    {
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("ID")));
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new IntConditionComponent("ConditionID")));
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("pot")));
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new RealConditionComponent("pot")));
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("funct")));
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new IntConditionComponent("funct", true, true)));
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("numscal")));

      // input: stoichiometry for reaction mechanism (IntRealBundle)
      // definition separator for int vectors
      std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepveccomp;
      intsepveccomp.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("stoich")));

      // definition int vectors
      std::vector<Teuchos::RCP<IntVectorConditionComponent>> intveccomp;
      intveccomp.emplace_back(Teuchos::rcp(new IntVectorConditionComponent("stoich", 2)));

      // definition separator for real vectors: length of the real vector is zero -> nothing is read
      std::vector<Teuchos::RCP<SeparatorConditionComponent>> realsepveccomp;

      // definition real vectors: length of the real vector is zero -> nothing is read
      std::vector<Teuchos::RCP<RealVectorConditionComponent>> realveccomp;

      electrodedomainkineticscomponents.emplace_back(Teuchos::rcp(
          new IntRealBundle("intreal bundle", Teuchos::rcp(new IntConditionComponent("numscal")),
              intsepveccomp, intveccomp, realsepveccomp, realveccomp)));

      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("e-")));
      electrodedomainkineticscomponents.emplace_back(Teuchos::rcp(new IntConditionComponent("e-")));
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("zero_cur")));
      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new IntConditionComponent("zero_cur")));

      // kinetic models
      std::vector<Teuchos::RCP<CondCompBundle>> kineticmodels;

      {
        // Butler-Volmer
        std::vector<Teuchos::RCP<ConditionComponent>> butlervolmer;
        butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent(
            "A_s")));  // ratio of electrode-electrolyte interface area to total two-phase volume
        butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("A_s")));
        butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_a")));
        butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
        butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_c")));
        butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));
        butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("i0")));
        butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("i0")));
        butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("gamma")));
        butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("gamma")));
        butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("refcon")));
        butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("refcon")));
        butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("dl_spec_cap")));
        butlervolmer.emplace_back(Teuchos::rcp(new RealConditionComponent("dl_spec_cap")));
        butlervolmer.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("END")));
        kineticmodels.emplace_back(Teuchos::rcp(
            new CondCompBundle("Butler-Volmer", butlervolmer, INPAR::ELCH::butler_volmer)));
      }

      electrodedomainkineticscomponents.emplace_back(
          Teuchos::rcp(new CondCompBundleSelector("kinetic model bundle",
              Teuchos::rcp(new StringConditionComponent("kinetic model", "Butler-Volmer",
                  Teuchos::tuple<std::string>("Butler-Volmer"),
                  Teuchos::tuple<int>(INPAR::ELCH::butler_volmer))),
              kineticmodels)));
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
    // definition of line and surface conditions for CCCV cell cycling
    Teuchos::RCP<ConditionDefinition> cccvcyclingline = Teuchos::rcp(
        new ConditionDefinition("DESIGN CCCV CELL CYCLING LINE CONDITIONS", "CCCVCycling",
            "line boundary condition for constant-current constant-voltage (CCCV) cell cycling",
            DRT::Condition::CCCVCycling, true, DRT::Condition::Line));

    Teuchos::RCP<ConditionDefinition> cccvcyclingsurf = Teuchos::rcp(
        new ConditionDefinition("DESIGN CCCV CELL CYCLING SURF CONDITIONS", "CCCVCycling",
            "surface boundary condition for constant-current constant-voltage (CCCV) cell cycling",
            DRT::Condition::CCCVCycling, true, DRT::Condition::Surface));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<ConditionComponent>> cccvcyclingcomponents;

    {
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("NumberOfHalfCycles")));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new IntConditionComponent("NumberOfHalfCycles")));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("BeginWithCharging")));
      cccvcyclingcomponents.emplace_back(Teuchos::rcp(new IntConditionComponent(
          "BeginWithCharging")));  // Boolean parameter represented by integer parameter
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("ConditionIDForCharge")));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new IntConditionComponent("ConditionIDForCharge", false, true)));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("ConditionIDForDischarge")));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new IntConditionComponent("ConditionIDForDischarge", false, true)));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("InitRelaxTime")));
      cccvcyclingcomponents.emplace_back(Teuchos::rcp(new RealConditionComponent("InitRelaxTime")));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("AdaptiveTimeSteppingInitRelax")));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new IntConditionComponent("AdaptiveTimeSteppingInitRelax")));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("NumAddAdaptTimeSteps")));
      cccvcyclingcomponents.emplace_back(
          Teuchos::rcp(new IntConditionComponent("NumAddAdaptTimeSteps", false, true)));
    }

    // insert input file line components into condition definitions
    for (auto& cccvcyclingcomponent : cccvcyclingcomponents)
    {
      cccvcyclingline->AddComponent(cccvcyclingcomponent);
      cccvcyclingsurf->AddComponent(cccvcyclingcomponent);
    }

    // insert condition definitions into global list of valid condition definitions
    condlist.emplace_back(cccvcyclingline);
    condlist.emplace_back(cccvcyclingsurf);
  }

  /*--------------------------------------------------------------------*/
  // boundary condition for constant-current constant-voltage (CCCV) half-cycle
  {
    // definition of line and surface conditions for CCCV half-cycle
    Teuchos::RCP<ConditionDefinition> cccvhalfcycleline = Teuchos::rcp(
        new ConditionDefinition("DESIGN CCCV HALF-CYCLE LINE CONDITIONS", "CCCVHalfCycle",
            "line boundary condition for constant-current constant-voltage (CCCV) half-cycle",
            DRT::Condition::CCCVHalfCycle, true, DRT::Condition::Line));

    Teuchos::RCP<ConditionDefinition> cccvhalfcyclesurf = Teuchos::rcp(
        new ConditionDefinition("DESIGN CCCV HALF-CYCLE SURF CONDITIONS", "CCCVHalfCycle",
            "surface boundary condition for constant-current constant-voltage (CCCV) half-cycle",
            DRT::Condition::CCCVHalfCycle, true, DRT::Condition::Surface));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<ConditionComponent>> cccvhalfcyclecomponents;

    {
      cccvhalfcyclecomponents.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("ID")));
      cccvhalfcyclecomponents.emplace_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));
      cccvhalfcyclecomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("Current")));
      cccvhalfcyclecomponents.emplace_back(Teuchos::rcp(new RealConditionComponent("Current")));
      cccvhalfcyclecomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("CutoffVoltage")));
      cccvhalfcyclecomponents.emplace_back(
          Teuchos::rcp(new RealConditionComponent("CutoffVoltage")));
      cccvhalfcyclecomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("CutoffCRate")));
      cccvhalfcyclecomponents.emplace_back(Teuchos::rcp(new RealConditionComponent("CutoffCRate")));
      cccvhalfcyclecomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("RelaxTime")));
      cccvhalfcyclecomponents.emplace_back(Teuchos::rcp(new RealConditionComponent("RelaxTime")));
      // switch adaptive time stepping on for different phases of half cycle: 1st: end of constant
      // current, 2nd: end of constant voltage, 3rd: end of relaxation
      cccvhalfcyclecomponents.emplace_back(
          Teuchos::rcp(new SeparatorConditionComponent("AdaptiveTimeSteppingPhaseOnOff")));
      cccvhalfcyclecomponents.emplace_back(
          Teuchos::rcp(new IntVectorConditionComponent("AdaptiveTimeSteppingPhaseOnOff", 3)));
    }

    // insert input file line components into condition definitions
    for (auto& cccvhalfcyclecomponent : cccvhalfcyclecomponents)
    {
      cccvhalfcycleline->AddComponent(cccvhalfcyclecomponent);
      cccvhalfcyclesurf->AddComponent(cccvhalfcyclecomponent);
    }

    // insert condition definitions into global list of valid condition definitions
    condlist.emplace_back(cccvhalfcycleline);
    condlist.emplace_back(cccvhalfcyclesurf);
  }
}
