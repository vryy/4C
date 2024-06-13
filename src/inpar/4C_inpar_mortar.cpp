/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameters for mortar coupling

\level 1

*/
/*----------------------------------------------------------------------*/



#include "4C_inpar_mortar.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::Mortar::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  /* parameters for mortar coupling */
  Teuchos::ParameterList& mortar = list->sublist("MORTAR COUPLING", false, "");

  setStringToIntegralParameter<int>("LM_SHAPEFCN", "Dual",
      "Type of employed set of shape functions",
      tuple<std::string>(
          "Dual", "dual", "Standard", "standard", "std", "PetrovGalerkin", "petrovgalerkin", "pg"),
      tuple<int>(shape_dual, shape_dual, shape_standard, shape_standard, shape_standard,
          shape_petrovgalerkin, shape_petrovgalerkin, shape_petrovgalerkin),
      &mortar);

  setStringToIntegralParameter<int>("SEARCH_ALGORITHM", "Binarytree", "Type of contact search",
      tuple<std::string>("BruteForce", "bruteforce", "BruteForceEleBased", "bruteforceelebased",
          "BinaryTree", "Binarytree", "binarytree"),
      tuple<int>(search_bfele, search_bfele, search_bfele, search_bfele, search_binarytree,
          search_binarytree, search_binarytree),
      &mortar);

  setStringToIntegralParameter<int>("BINARYTREE_UPDATETYPE", "BottomUp",
      "Type of binary tree update, which is either a bottom up or a top down approach.",
      tuple<std::string>("BottomUp", "TopDown"),
      tuple<int>(binarytree_bottom_up, binarytree_top_down), &mortar);

  Core::UTILS::DoubleParameter(
      "SEARCH_PARAM", 0.3, "Radius / Bounding volume inflation for contact search", &mortar);

  Core::UTILS::BoolParameter("SEARCH_USE_AUX_POS", "Yes",
      "If chosen auxiliary position is used for computing dops", &mortar);

  setStringToIntegralParameter<int>("LM_QUAD", "undefined",
      "Type of LM interpolation for quadratic FE",
      tuple<std::string>(
          "undefined", "quad", "quadratic", "pwlin", "piecewiselinear", "lin", "linear", "const"),
      tuple<int>(lagmult_undefined, lagmult_quad, lagmult_quad, lagmult_pwlin, lagmult_pwlin,
          lagmult_lin, lagmult_lin, lagmult_const),
      &mortar);

  Core::UTILS::BoolParameter("CROSSPOINTS", "No",
      "If chosen, multipliers are removed from crosspoints / edge nodes", &mortar);

  setStringToIntegralParameter<int>("LM_DUAL_CONSISTENT", "boundary",
      "For which elements should the dual basis be calculated on EXACTLY the same GPs as the "
      "contact terms",
      tuple<std::string>("none", "boundary", "all"),
      tuple<int>(consistent_none, consistent_boundary, consistent_all), &mortar);

  setStringToIntegralParameter<int>("MESH_RELOCATION", "Initial", "Type of mesh relocation",
      tuple<std::string>("Initial", "initial", "Every_Timestep", "every_timestep", "No", "no"),
      tuple<int>(relocation_initial, relocation_initial, relocation_timestep, relocation_timestep,
          relocation_none, relocation_none),
      &mortar);

  setStringToIntegralParameter<int>("ALGORITHM", "Mortar", "Type of meshtying/contact algorithm",
      tuple<std::string>("mortar", "Mortar", "nts", "NTS", "gpts", "GPTS", "lts", "LTS", "ltl",
          "LTL", "stl", "STL"),
      tuple<int>(algorithm_mortar, algorithm_mortar, algorithm_nts, algorithm_nts, algorithm_gpts,
          algorithm_gpts, algorithm_lts, algorithm_lts, algorithm_ltl, algorithm_ltl, algorithm_stl,
          algorithm_stl),
      &mortar);

  setStringToIntegralParameter<int>("INTTYPE", "Segments", "Type of numerical integration scheme",
      tuple<std::string>(
          "Segments", "segments", "Elements", "elements", "Elements_BS", "elements_BS"),
      tuple<int>(inttype_segments, inttype_segments, inttype_elements, inttype_elements,
          inttype_elements_BS, inttype_elements_BS),
      &mortar);

  Core::UTILS::IntParameter(
      "NUMGP_PER_DIM", 0, "Number of employed integration points per dimension", &mortar);

  setStringToIntegralParameter<int>("TRIANGULATION", "Delaunay",
      "Type of triangulation for segment-based integration",
      tuple<std::string>("Delaunay", "delaunay", "Center", "center"),
      tuple<int>(triangulation_delaunay, triangulation_delaunay, triangulation_center,
          triangulation_center),
      &mortar);

  Core::UTILS::BoolParameter("RESTART_WITH_MESHTYING", "No",
      "Must be chosen if a non-meshtying simulation is to be restarted with meshtying", &mortar);

  Core::UTILS::BoolParameter("OUTPUT_INTERFACES", "No",
      "Write output for each mortar interface separately.\nThis is an additional feature, purely "
      "to enhance visualization. Currently, this is limited to solid meshtying and contact w/o "
      "friction.",
      &mortar);

  /*--------------------------------------------------------------------*/
  // parameters for parallel redistribution of mortar interfaces
  Teuchos::ParameterList& parallelRedist = mortar.sublist("PARALLEL REDISTRIBUTION", false,
      "Parameters to control parallel redistribution of mortar interfaces");

  Core::UTILS::BoolParameter("EXPLOIT_PROXIMITY", "Yes",
      "Exploit information on geometric proximity to split slave interface into close and "
      "non-close parts and redistribute them independently. [Contact only]",
      &parallelRedist);

  setStringToIntegralParameter<ExtendGhosting>("GHOSTING_STRATEGY", "redundant_master",
      "Type of interface ghosting and ghosting extension algorithm",
      tuple<std::string>("redundant_all", "redundant_master", "round_robin", "binning"),
      tuple<ExtendGhosting>(ExtendGhosting::redundant_all, ExtendGhosting::redundant_master,
          ExtendGhosting::roundrobin, ExtendGhosting::binning),
      &parallelRedist);

  Core::UTILS::DoubleParameter("IMBALANCE_TOL", 1.1,
      "Max. relative imbalance of subdomain size after redistribution", &parallelRedist);

  Core::UTILS::DoubleParameter("MAX_BALANCE_EVAL_TIME", 2.0,
      "Max-to-min ratio of contact evalation time per processor to triggger parallel "
      "redistribution",
      &parallelRedist);

  Core::UTILS::DoubleParameter("MAX_BALANCE_SLAVE_ELES", 0.5,
      "Max-to-min ratio of mortar slave elements per processor to triggger parallel "
      "redistribution",
      &parallelRedist);

  Core::UTILS::IntParameter("MIN_ELEPROC", 0,
      "Minimum no. of elements per processor for parallel redistribution", &parallelRedist);

  setStringToIntegralParameter<ParallelRedist>("PARALLEL_REDIST", "Static",
      "Type of redistribution algorithm",
      tuple<std::string>("None", "none", "No", "no", "Static", "static", "Dynamic", "dynamic"),
      tuple<ParallelRedist>(ParallelRedist::redist_none, ParallelRedist::redist_none,
          ParallelRedist::redist_none, ParallelRedist::redist_none, ParallelRedist::redist_static,
          ParallelRedist::redist_static, ParallelRedist::redist_dynamic,
          ParallelRedist::redist_dynamic),
      &parallelRedist);

  Core::UTILS::BoolParameter("PRINT_DISTRIBUTION", "Yes",
      "Print details of the parallel distribution, i.e. number of nodes/elements for each rank.",
      &parallelRedist);
}

void Inpar::Mortar::SetValidConditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;
  /*--------------------------------------------------------------------*/
  // mortar contact

  std::vector<Teuchos::RCP<Input::LineComponent>> contactcomponents;

  contactcomponents.push_back(Teuchos::rcp(new Input::IntComponent("Interface ID")));
  contactcomponents.push_back(Teuchos::rcp(new Input::SelectionComponent("Side", "Master",
      Teuchos::tuple<std::string>("Master", "Slave", "Selfcontact"),
      Teuchos::tuple<std::string>("Master", "Slave", "Selfcontact"))));
  contactcomponents.push_back(Teuchos::rcp(new Input::SelectionComponent("Initialization",
      "Inactive", Teuchos::tuple<std::string>("Inactive", "Active"),
      Teuchos::tuple<std::string>("Inactive", "Active"), true)));

  contactcomponents.push_back(
      Teuchos::rcp(new Input::SeparatorComponent("FrCoeffOrBound", "", true)));
  contactcomponents.push_back(Teuchos::rcp(new Input::RealComponent("FrCoeffOrBound", {0., true})));

  contactcomponents.push_back(
      Teuchos::rcp(new Input::SeparatorComponent("AdhesionBound", "", true)));
  contactcomponents.push_back(Teuchos::rcp(new Input::RealComponent("AdhesionBound")));

  contactcomponents.push_back(
      Teuchos::rcp(new Input::SelectionComponent("Application", "Solidcontact",
          Teuchos::tuple<std::string>("Solidcontact", "Beamtosolidcontact", "Beamtosolidmeshtying"),
          Teuchos::tuple<std::string>("Solidcontact", "Beamtosolidcontact", "Beamtosolidmeshtying"),
          true)));

  // optional DBC handling
  contactcomponents.push_back(Teuchos::rcp(new Input::SelectionComponent("dbc_handling",
      "DoNothing", Teuchos::tuple<std::string>("DoNothing", "RemoveDBCSlaveNodes"),
      Teuchos::tuple<int>(static_cast<int>(DBCHandling::do_nothing),
          static_cast<int>(DBCHandling::remove_dbc_nodes_from_slave_side)),
      true)));

  // optional two half pass approach
  contactcomponents.push_back(Teuchos::rcp(new Input::SeparatorComponent("TwoHalfPass", "", true)));
  contactcomponents.push_back(Teuchos::rcp(new Input::RealComponent("TwoHalfPass")));

  // optional reference configuration check for non-smooth self contact surfaces
  contactcomponents.push_back(Teuchos::rcp(
      new Input::SeparatorComponent("RefConfCheckNonSmoothSelfContactSurface", "", true)));
  contactcomponents.push_back(
      Teuchos::rcp(new Input::RealComponent("RefConfCheckNonSmoothSelfContactSurface")));

  contactcomponents.push_back(
      Teuchos::rcp(new Input::SeparatorComponent("ConstitutiveLawID", "", true)));
  contactcomponents.push_back(
      Teuchos::rcp(new Input::IntComponent("ConstitutiveLawID", {0, false, true, true})));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> linecontact =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE MORTAR CONTACT CONDITIONS 2D", "Contact", "Line Contact Coupling",
          Core::Conditions::Contact, true, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfcontact =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN SURF MORTAR CONTACT CONDITIONS 3D", "Contact", "Surface Contact Coupling",
          Core::Conditions::Contact, true, Core::Conditions::geometry_type_surface));

  for (unsigned i = 0; i < contactcomponents.size(); ++i)
  {
    linecontact->add_component(contactcomponents[i]);
    surfcontact->add_component(contactcomponents[i]);
  }

  condlist.push_back(linecontact);
  condlist.push_back(surfcontact);

  /*--------------------------------------------------------------------*/
  // mortar coupling (for ALL kinds of interface problems except contact)

  std::vector<Teuchos::RCP<Input::LineComponent>> mortarcomponents;

  mortarcomponents.push_back(Teuchos::rcp(new Input::IntComponent("Interface ID")));
  mortarcomponents.push_back(Teuchos::rcp(new Input::SelectionComponent("Side", "Master",
      Teuchos::tuple<std::string>("Master", "Slave"),
      Teuchos::tuple<std::string>("Master", "Slave"))));
  mortarcomponents.push_back(Teuchos::rcp(new Input::SelectionComponent("Initialization",
      "Inactive", Teuchos::tuple<std::string>("Inactive", "Active"),
      Teuchos::tuple<std::string>("Inactive", "Active"), true)));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> linemortar =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE MORTAR COUPLING CONDITIONS 2D", "Mortar", "Line Mortar Coupling",
          Core::Conditions::Mortar, true, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfmortar =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN SURF MORTAR COUPLING CONDITIONS 3D", "Mortar", "Surface Mortar Coupling",
          Core::Conditions::Mortar, true, Core::Conditions::geometry_type_surface));

  for (unsigned i = 0; i < mortarcomponents.size(); ++i)
  {
    linemortar->add_component(mortarcomponents[i]);
    surfmortar->add_component(mortarcomponents[i]);
  }

  condlist.push_back(linemortar);
  condlist.push_back(surfmortar);

  /*--------------------------------------------------------------------*/
  // mortar coupling symmetry condition

  std::vector<Teuchos::RCP<Input::LineComponent>> mrtrsymcomponents;
  mrtrsymcomponents.push_back(Teuchos::rcp(new Input::SeparatorComponent("ONOFF")));
  mrtrsymcomponents.push_back(Teuchos::rcp(new Input::IntVectorComponent("onoff", 3)));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> linemrtrsym = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN LINE MORTAR SYMMETRY CONDITIONS 3D",
          "mrtrsym", "Symmetry plane normal for 3D contact", Core::Conditions::LineMrtrSym, true,
          Core::Conditions::geometry_type_line));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointmrtrsym = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN POINT MORTAR SYMMETRY CONDITIONS 2D/3D",
          "mrtrsym", "Symmetry plane normal for 2D/3D contact", Core::Conditions::PointMrtrSym,
          true, Core::Conditions::geometry_type_point));

  for (unsigned i = 0; i < mrtrsymcomponents.size(); ++i)
  {
    linemrtrsym->add_component(mrtrsymcomponents[i]);
    pointmrtrsym->add_component(mrtrsymcomponents[i]);
  }

  condlist.push_back(linemrtrsym);
  condlist.push_back(pointmrtrsym);

  /*--------------------------------------------------------------------*/
  // mortar edge/corner condition

  Teuchos::RCP<Core::Conditions::ConditionDefinition> edgemrtr =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE MORTAR EDGE CONDITIONS 3D", "mrtredge", "Geometrical edge for 3D contact",
          Core::Conditions::EdgeMrtr, true, Core::Conditions::geometry_type_line));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> cornermrtr = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN POINT MORTAR CORNER CONDITIONS 2D/3D",
          "mrtrcorner", "Geometrical corner for 2D/3D contact", Core::Conditions::CornerMrtr, true,
          Core::Conditions::geometry_type_point));

  condlist.push_back(edgemrtr);
  condlist.push_back(cornermrtr);


  {
    /*--------------------------------------------------------------------*/
    // mortar coupling (for ALL kinds of interface problems except contact)
    std::vector<Teuchos::RCP<Input::LineComponent>> mortarcomponents;

    mortarcomponents.push_back(Teuchos::rcp(new Input::IntComponent("Interface ID")));
    mortarcomponents.push_back(Teuchos::rcp(new Input::SelectionComponent("Side", "Master",
        Teuchos::tuple<std::string>("Master", "Slave"),
        Teuchos::tuple<std::string>("Master", "Slave"))));
    mortarcomponents.push_back(Teuchos::rcp(new Input::SelectionComponent("Initialization",
        "Inactive", Teuchos::tuple<std::string>("Inactive", "Active"),
        Teuchos::tuple<std::string>("Inactive", "Active"), true)));

    Teuchos::RCP<Core::Conditions::ConditionDefinition> linemortar = Teuchos::rcp(
        new Core::Conditions::ConditionDefinition("DESIGN LINE MORTAR MULTI-COUPLING CONDITIONS 2D",
            "MortarMulti", "Line Mortar Multi-Coupling", Core::Conditions::MortarMulti, true,
            Core::Conditions::geometry_type_line));
    Teuchos::RCP<Core::Conditions::ConditionDefinition> surfmortar = Teuchos::rcp(
        new Core::Conditions::ConditionDefinition("DESIGN SURF MORTAR MULTI-COUPLING CONDITIONS 3D",
            "MortarMulti", "Surface Mortar Multi-Coupling", Core::Conditions::MortarMulti, true,
            Core::Conditions::geometry_type_surface));

    for (unsigned i = 0; i < mortarcomponents.size(); ++i)
    {
      linemortar->add_component(mortarcomponents[i]);
      surfmortar->add_component(mortarcomponents[i]);
    }

    condlist.push_back(linemortar);
    condlist.push_back(surfmortar);
  }
}

FOUR_C_NAMESPACE_CLOSE
