/*----------------------------------------------------------------------*/
/*!
\file inpar_mortar.cpp

\brief Input parameters for mortar coupling

\level 1

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_mortar.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::MORTAR::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::Array<std::string> yesnotuple =
      tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

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

  DoubleParameter(
      "SEARCH_PARAM", 0.3, "Radius / Bounding volume inflation for contact search", &mortar);

  setStringToIntegralParameter<int>("SEARCH_USE_AUX_POS", "Yes",
      "If chosen auxiliary position is used for computing dops", yesnotuple, yesnovalue, &mortar);

  setStringToIntegralParameter<int>("LM_QUAD", "undefined",
      "Type of LM interpolation for quadratic FE",
      tuple<std::string>(
          "undefined", "quad", "quadratic", "pwlin", "piecewiselinear", "lin", "linear", "const"),
      tuple<int>(lagmult_undefined, lagmult_quad, lagmult_quad, lagmult_pwlin, lagmult_pwlin,
          lagmult_lin, lagmult_lin, lagmult_const),
      &mortar);

  setStringToIntegralParameter<int>("CROSSPOINTS", "No",
      "If chosen, multipliers are removed from crosspoints / edge nodes", yesnotuple, yesnovalue,
      &mortar);

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

  setStringToIntegralParameter<int>("REDUNDANT_STORAGE", "Master",
      "Type of redundancy in interface storage",
      tuple<std::string>("All", "all", "Master", "master", "None", "none"),
      tuple<int>(redundant_all, redundant_all, redundant_master, redundant_master, redundant_none,
          redundant_none),
      &mortar);

  setStringToIntegralParameter<int>("PARALLEL_STRATEGY", "redundant_ghosting",
      "Type of parallel interface evaluation",
      tuple<std::string>("rg", "redundant_ghosting", "ghosting", "rre", "roundrobinevaluate",
          "RoundRobinEvaluate", "rrg", "roundrobinghost", "RoundRobinGhost", "bs",
          "binningstrategy", "binning"),
      tuple<int>(ghosting_redundant, ghosting_redundant, ghosting_redundant, roundrobinevaluate,
          roundrobinevaluate, roundrobinevaluate, roundrobinghost, roundrobinghost, roundrobinghost,
          binningstrategy, binningstrategy, binningstrategy),
      &mortar);

  setStringToIntegralParameter<int>("PARALLEL_REDIST", "Static", "Type of redistribution algorithm",
      tuple<std::string>("None", "none", "No", "no", "Static", "static", "Dynamic", "dynamic"),
      tuple<int>(parredist_none, parredist_none, parredist_none, parredist_none, parredist_static,
          parredist_static, parredist_dynamic, parredist_dynamic),
      &mortar);

  setStringToIntegralParameter<int>("ALGORITHM", "Mortar", "Type of meshtying/contact algorithm",
      tuple<std::string>("mortar", "Mortar", "nts", "NTS", "gpts", "GPTS", "lts", "LTS", "ltl",
          "LTL", "stl", "STL"),
      tuple<int>(algorithm_mortar, algorithm_mortar, algorithm_nts, algorithm_nts, algorithm_gpts,
          algorithm_gpts, algorithm_lts, algorithm_lts, algorithm_ltl, algorithm_ltl, algorithm_stl,
          algorithm_stl),
      &mortar);

  DoubleParameter("MAX_BALANCE", 2.0,
      "Maximum value of load balance measure before parallel redistribution", &mortar);
  IntParameter("MIN_ELEPROC", 0,
      "Minimum no. of elements per processor for parallel redistribution", &mortar);

  setStringToIntegralParameter<int>("INTTYPE", "Segments", "Type of numerical integration scheme",
      tuple<std::string>(
          "Segments", "segments", "Elements", "elements", "Elements_BS", "elements_BS"),
      tuple<int>(inttype_segments, inttype_segments, inttype_elements, inttype_elements,
          inttype_elements_BS, inttype_elements_BS),
      &mortar);

  IntParameter("NUMGP_PER_DIM", 0, "Number of employed integration points per dimension", &mortar);

  setStringToIntegralParameter<int>("TRIANGULATION", "Delaunay",
      "Type of triangulation for segment-based integration",
      tuple<std::string>("Delaunay", "delaunay", "Center", "center"),
      tuple<int>(triangulation_delaunay, triangulation_delaunay, triangulation_center,
          triangulation_center),
      &mortar);

  setStringToIntegralParameter<int>("RESTART_WITH_MESHTYING", "No",
      "Must be chosen if a non-meshtying simulation is to be restarted with meshtying", yesnotuple,
      yesnovalue, &mortar);
}



void INPAR::MORTAR::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;
  /*--------------------------------------------------------------------*/
  // mortar contact

  std::vector<Teuchos::RCP<ConditionComponent>> contactcomponents;

  contactcomponents.push_back(Teuchos::rcp(new IntConditionComponent("Interface ID")));
  contactcomponents.push_back(Teuchos::rcp(new StringConditionComponent("Side", "Master",
      Teuchos::tuple<std::string>("Master", "Slave", "Selfcontact"),
      Teuchos::tuple<std::string>("Master", "Slave", "Selfcontact"))));
  contactcomponents.push_back(Teuchos::rcp(new StringConditionComponent("Initialization",
      "Inactive", Teuchos::tuple<std::string>("Inactive", "Active"),
      Teuchos::tuple<std::string>("Inactive", "Active"), true)));

  contactcomponents.push_back(
      Teuchos::rcp(new SeparatorConditionComponent("FrCoeffOrBound", true)));
  contactcomponents.push_back(Teuchos::rcp(new RealConditionComponent("FrCoeffOrBound")));

  contactcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("AdhesionBound", true)));
  contactcomponents.push_back(Teuchos::rcp(new RealConditionComponent("AdhesionBound")));

  contactcomponents.push_back(
      Teuchos::rcp(new StringConditionComponent("Application", "Solidcontact",
          Teuchos::tuple<std::string>("Solidcontact", "Beamtosolidcontact", "Beamtosolidmeshtying"),
          Teuchos::tuple<std::string>("Solidcontact", "Beamtosolidcontact", "Beamtosolidmeshtying"),
          true)));

  // optional DBC handling
  contactcomponents.push_back(Teuchos::rcp(new StringConditionComponent("dbc_handling", "DoNothing",
      Teuchos::tuple<std::string>("DoNothing", "RemoveDBCSlaveNodes"),
      Teuchos::tuple<int>(static_cast<int>(DBCHandling::do_nothing),
          static_cast<int>(DBCHandling::remove_dbc_nodes_from_slave_side)),
      true)));

  Teuchos::RCP<ConditionDefinition> linecontact =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE MORTAR CONTACT CONDITIONS 2D", "Contact",
          "Line Contact Coupling", DRT::Condition::Contact, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfcontact =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURF MORTAR CONTACT CONDITIONS 3D", "Contact",
          "Surface Contact Coupling", DRT::Condition::Contact, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < contactcomponents.size(); ++i)
  {
    linecontact->AddComponent(contactcomponents[i]);
    surfcontact->AddComponent(contactcomponents[i]);
  }

  condlist.push_back(linecontact);
  condlist.push_back(surfcontact);

  /*--------------------------------------------------------------------*/
  // mortar coupling (for ALL kinds of interface problems except contact)

  std::vector<Teuchos::RCP<ConditionComponent>> mortarcomponents;

  mortarcomponents.push_back(Teuchos::rcp(new IntConditionComponent("Interface ID")));
  mortarcomponents.push_back(Teuchos::rcp(
      new StringConditionComponent("Side", "Master", Teuchos::tuple<std::string>("Master", "Slave"),
          Teuchos::tuple<std::string>("Master", "Slave"))));
  mortarcomponents.push_back(Teuchos::rcp(new StringConditionComponent("Initialization", "Inactive",
      Teuchos::tuple<std::string>("Inactive", "Active"),
      Teuchos::tuple<std::string>("Inactive", "Active"), true)));

  Teuchos::RCP<ConditionDefinition> linemortar =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE MORTAR COUPLING CONDITIONS 2D", "Mortar",
          "Line Mortar Coupling", DRT::Condition::Mortar, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfmortar =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURF MORTAR COUPLING CONDITIONS 3D", "Mortar",
          "Surface Mortar Coupling", DRT::Condition::Mortar, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < mortarcomponents.size(); ++i)
  {
    linemortar->AddComponent(mortarcomponents[i]);
    surfmortar->AddComponent(mortarcomponents[i]);
  }

  condlist.push_back(linemortar);
  condlist.push_back(surfmortar);

  /*--------------------------------------------------------------------*/
  // mortar coupling symmetry condition

  std::vector<Teuchos::RCP<ConditionComponent>> mrtrsymcomponents;
  mrtrsymcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
  mrtrsymcomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("onoff", 3)));

  Teuchos::RCP<ConditionDefinition> linemrtrsym =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE MORTAR SYMMETRY CONDITIONS 3D", "mrtrsym",
          "Symmetry plane normal for 3D contact", DRT::Condition::LineMrtrSym, true,
          DRT::Condition::Line));

  Teuchos::RCP<ConditionDefinition> pointmrtrsym =
      Teuchos::rcp(new ConditionDefinition("DESIGN POINT MORTAR SYMMETRY CONDITIONS 2D/3D",
          "mrtrsym", "Symmetry plane normal for 2D/3D contact", DRT::Condition::PointMrtrSym, true,
          DRT::Condition::Point));

  for (unsigned i = 0; i < mrtrsymcomponents.size(); ++i)
  {
    linemrtrsym->AddComponent(mrtrsymcomponents[i]);
    pointmrtrsym->AddComponent(mrtrsymcomponents[i]);
  }

  condlist.push_back(linemrtrsym);
  condlist.push_back(pointmrtrsym);

  /*--------------------------------------------------------------------*/
  // mortar edge/corner condition

  //  std::vector<Teuchos::RCP<ConditionComponent> > nscomponents;
  //  nscomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
  //  nscomponents.push_back(Teuchos::rcp(new IntVectorConditionComponent("onoff",3)));

  Teuchos::RCP<ConditionDefinition> edgemrtr =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE MORTAR EDGE CONDITIONS 3D", "mrtredge",
          "Geometrical edge for 3D contact", DRT::Condition::EdgeMrtr, true, DRT::Condition::Line));

  Teuchos::RCP<ConditionDefinition> cornermrtr =
      Teuchos::rcp(new ConditionDefinition("DESIGN POINT MORTAR CORNER CONDITIONS 2D/3D",
          "mrtrcorner", "Geometrical corner for 2D/3D contact", DRT::Condition::CornerMrtr, true,
          DRT::Condition::Point));

  //  for (unsigned i=0; i<nscomponents.size(); ++i)
  //  {
  //    edgemrtr->AddComponent(nscomponents[i]);
  //    cornermrtr->AddComponent(nscomponents[i]);
  //  }

  condlist.push_back(edgemrtr);
  condlist.push_back(cornermrtr);


  {
    /*--------------------------------------------------------------------*/
    // mortar coupling (for ALL kinds of interface problems except contact)
    std::vector<Teuchos::RCP<ConditionComponent>> mortarcomponents;

    mortarcomponents.push_back(Teuchos::rcp(new IntConditionComponent("Interface ID")));
    mortarcomponents.push_back(Teuchos::rcp(new StringConditionComponent("Side", "Master",
        Teuchos::tuple<std::string>("Master", "Slave"),
        Teuchos::tuple<std::string>("Master", "Slave"))));
    mortarcomponents.push_back(Teuchos::rcp(new StringConditionComponent("Initialization",
        "Inactive", Teuchos::tuple<std::string>("Inactive", "Active"),
        Teuchos::tuple<std::string>("Inactive", "Active"), true)));

    Teuchos::RCP<ConditionDefinition> linemortar = Teuchos::rcp(
        new ConditionDefinition("DESIGN LINE MORTAR MULTI-COUPLING CONDITIONS 2D", "MortarMulti",
            "Line Mortar Multi-Coupling", DRT::Condition::MortarMulti, true, DRT::Condition::Line));
    Teuchos::RCP<ConditionDefinition> surfmortar =
        Teuchos::rcp(new ConditionDefinition("DESIGN SURF MORTAR MULTI-COUPLING CONDITIONS 3D",
            "MortarMulti", "Surface Mortar Multi-Coupling", DRT::Condition::MortarMulti, true,
            DRT::Condition::Surface));

    for (unsigned i = 0; i < mortarcomponents.size(); ++i)
    {
      linemortar->AddComponent(mortarcomponents[i]);
      surfmortar->AddComponent(mortarcomponents[i]);
    }

    condlist.push_back(linemortar);
    condlist.push_back(surfmortar);
  }
}
