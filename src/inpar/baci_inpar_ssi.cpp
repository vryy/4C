/*----------------------------------------------------------------------*/
/*! \file
\brief input parameters for solid-scatra-interaction

\level 2


*/

/*----------------------------------------------------------------------*/

#include "baci_inpar_ssi.hpp"

#include "baci_inpar_s2i.hpp"
#include "baci_inpar_scatra.hpp"
#include "baci_inpar_validparameters.hpp"
#include "baci_lib_conditiondefinition.hpp"
#include "baci_linalg_equilibrate.hpp"
#include "baci_linalg_sparseoperator.hpp"

BACI_NAMESPACE_OPEN

void INPAR::SSI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& ssidyn =
      list->sublist("SSI CONTROL", false, "Control paramters for scatra structure interaction");

  // Output type
  DoubleParameter(
      "RESTARTEVRYTIME", 0, "write restart possibility every RESTARTEVRY steps", &ssidyn);
  IntParameter("RESTARTEVRY", 1, "write restart possibility every RESTARTEVRY steps", &ssidyn);
  // Time loop control
  IntParameter("NUMSTEP", 200, "maximum number of Timesteps", &ssidyn);
  DoubleParameter("MAXTIME", 1000.0, "total simulation time", &ssidyn);
  DoubleParameter("TIMESTEP", -1, "time step size dt", &ssidyn);
  BoolParameter("DIFFTIMESTEPSIZE", "No", "use different step size for scatra and solid", &ssidyn);
  DoubleParameter("RESULTSEVRYTIME", 0, "increment for writing solution", &ssidyn);
  IntParameter("RESULTSEVRY", 1, "increment for writing solution", &ssidyn);
  IntParameter("ITEMAX", 10, "maximum number of iterations over fields", &ssidyn);
  BoolParameter("SCATRA_FROM_RESTART_FILE", "No",
      "read scatra result from restart files (use option 'restartfromfile' during execution of "
      "baci)",
      &ssidyn);
  StringParameter(
      "SCATRA_FILENAME", "nil", "Control-file name for reading scatra results in SSI", &ssidyn);

  // Type of coupling strategy between the two fields
  setStringToIntegralParameter<FieldCoupling>("FIELDCOUPLING", "volume_matching",
      "Type of coupling strategy between fields",
      tuple<std::string>("volume_matching", "volume_nonmatching", "boundary_nonmatching",
          "volumeboundary_matching"),
      tuple<FieldCoupling>(FieldCoupling::volume_match, FieldCoupling::volume_nonmatch,
          FieldCoupling::boundary_nonmatch, FieldCoupling::volumeboundary_match),
      &ssidyn);

  // Coupling strategy for SSI solvers
  setStringToIntegralParameter<SolutionSchemeOverFields>("COUPALGO", "ssi_IterStagg",
      "Coupling strategies for SSI solvers",
      tuple<std::string>("ssi_OneWay_ScatraToSolid", "ssi_OneWay_SolidToScatra",
          //                                "ssi_SequStagg_ScatraToSolid",
          //                                "ssi_SequStagg_SolidToScatra",
          "ssi_IterStagg", "ssi_IterStaggFixedRel_ScatraToSolid",
          "ssi_IterStaggFixedRel_SolidToScatra", "ssi_IterStaggAitken_ScatraToSolid",
          "ssi_IterStaggAitken_SolidToScatra", "ssi_Monolithic"),
      tuple<SolutionSchemeOverFields>(SolutionSchemeOverFields::ssi_OneWay_ScatraToSolid,
          SolutionSchemeOverFields::ssi_OneWay_SolidToScatra,
          //                                ssi_SequStagg_ScatraToSolid,
          //                                ssi_SequStagg_SolidToScatra,
          SolutionSchemeOverFields::ssi_IterStagg,
          SolutionSchemeOverFields::ssi_IterStaggFixedRel_ScatraToSolid,
          SolutionSchemeOverFields::ssi_IterStaggFixedRel_SolidToScatra,
          SolutionSchemeOverFields::ssi_IterStaggAitken_ScatraToSolid,
          SolutionSchemeOverFields::ssi_IterStaggAitken_SolidToScatra,
          SolutionSchemeOverFields::ssi_Monolithic),
      &ssidyn);

  // type of scalar transport time integration
  setStringToIntegralParameter<ScaTraTimIntType>("SCATRATIMINTTYPE", "Standard",
      "scalar transport time integration type is needed to instantiate correct scalar transport "
      "time integration scheme for ssi problems",
      tuple<std::string>("Standard", "Cardiac_Monodomain", "Elch"),
      tuple<ScaTraTimIntType>(
          ScaTraTimIntType::standard, ScaTraTimIntType::cardiac_monodomain, ScaTraTimIntType::elch),
      &ssidyn);

  // Restart from Structure problem instead of SSI
  BoolParameter("RESTART_FROM_STRUCTURE", "no",
      "restart from structure problem (e.g. from prestress calculations) instead of ssi", &ssidyn);

  // Adaptive time stepping
  BoolParameter("ADAPTIVE_TIMESTEPPING", "no", "flag for adaptive time stepping", &ssidyn);

  // do redistribution by binning of solid mechanics discretization (scatra dis is cloned from solid
  // dis for volume_matching and volumeboundary_matching)
  BoolParameter("REDISTRIBUTE_SOLID", "No",
      "redistribution by binning of solid mechanics discretization", &ssidyn);

  /*----------------------------------------------------------------------*/
  /* parameters for partitioned SSI */
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& ssidynpart = ssidyn.sublist("PARTITIONED", false,
      "Partitioned Structure Scalar Interaction\n"
      "Control section for partitioned SSI");

  // Solver parameter for relaxation of iterative staggered partitioned SSI
  DoubleParameter("MAXOMEGA", 10.0, "largest omega allowed for Aitken relaxation", &ssidynpart);
  DoubleParameter("MINOMEGA", 0.1, "smallest omega allowed for Aitken relaxation", &ssidynpart);
  DoubleParameter("STARTOMEGA", 1.0, "fixed relaxation parameter", &ssidynpart);

  // convergence tolerance of outer iteration loop
  DoubleParameter("CONVTOL", 1e-6,
      "tolerance for convergence check of outer iteration within partitioned SSI", &ssidynpart);

  /*----------------------------------------------------------------------*/
  /* parameters for monolithic SSI */
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& ssidynmono = ssidyn.sublist("MONOLITHIC", false,
      "Monolithic Structure Scalar Interaction\n"
      "Control section for monolithic SSI");

  // convergence tolerances of Newton-Raphson iteration loop
  DoubleParameter("ABSTOLRES", 1.e-14,
      "absolute tolerance for deciding if global residual of nonlinear problem is already zero",
      &ssidynmono);
  DoubleParameter("CONVTOL", 1.e-6,
      "tolerance for convergence check of Newton-Raphson iteration within monolithic SSI",
      &ssidynmono);

  // ID of linear solver for global system of equations
  IntParameter(
      "LINEAR_SOLVER", -1, "ID of linear solver for global system of equations", &ssidynmono);

  // type of global system matrix in global system of equations
  setStringToIntegralParameter<CORE::LINALG::MatrixType>("MATRIXTYPE", "undefined",
      "type of global system matrix in global system of equations",
      tuple<std::string>("undefined", "block", "sparse"),
      tuple<CORE::LINALG::MatrixType>(CORE::LINALG::MatrixType::undefined,
          CORE::LINALG::MatrixType::block_field, CORE::LINALG::MatrixType::sparse),
      &ssidynmono);

  setStringToIntegralParameter<CORE::LINALG::EquilibrationMethod>("EQUILIBRATION", "none",
      "flag for equilibration of global system of equations",
      tuple<std::string>("none", "rows_full", "rows_maindiag", "columns_full", "columns_maindiag",
          "rowsandcolumns_full", "rowsandcolumns_maindiag", "local"),
      tuple<CORE::LINALG::EquilibrationMethod>(CORE::LINALG::EquilibrationMethod::none,
          CORE::LINALG::EquilibrationMethod::rows_full,
          CORE::LINALG::EquilibrationMethod::rows_maindiag,
          CORE::LINALG::EquilibrationMethod::columns_full,
          CORE::LINALG::EquilibrationMethod::columns_maindiag,
          CORE::LINALG::EquilibrationMethod::rowsandcolumns_full,
          CORE::LINALG::EquilibrationMethod::rowsandcolumns_maindiag,
          CORE::LINALG::EquilibrationMethod::local),
      &ssidynmono);

  setStringToIntegralParameter<CORE::LINALG::EquilibrationMethod>("EQUILIBRATION_STRUCTURE", "none",
      "flag for equilibration of structural equations",
      tuple<std::string>(
          "none", "rows_maindiag", "columns_maindiag", "rowsandcolumns_maindiag", "symmetry"),
      tuple<CORE::LINALG::EquilibrationMethod>(CORE::LINALG::EquilibrationMethod::none,
          CORE::LINALG::EquilibrationMethod::rows_maindiag,
          CORE::LINALG::EquilibrationMethod::columns_maindiag,
          CORE::LINALG::EquilibrationMethod::rowsandcolumns_maindiag,
          CORE::LINALG::EquilibrationMethod::symmetry),
      &ssidynmono);

  setStringToIntegralParameter<CORE::LINALG::EquilibrationMethod>("EQUILIBRATION_SCATRA", "none",
      "flag for equilibration of scatra equations",
      tuple<std::string>(
          "none", "rows_maindiag", "columns_maindiag", "rowsandcolumns_maindiag", "symmetry"),
      tuple<CORE::LINALG::EquilibrationMethod>(CORE::LINALG::EquilibrationMethod::none,
          CORE::LINALG::EquilibrationMethod::rows_maindiag,
          CORE::LINALG::EquilibrationMethod::columns_maindiag,
          CORE::LINALG::EquilibrationMethod::rowsandcolumns_maindiag,
          CORE::LINALG::EquilibrationMethod::symmetry),
      &ssidynmono);

  BoolParameter("PRINT_MAT_RHS_MAP_MATLAB", "no",
      "print system matrix, rhs vector, and full map to matlab readable file after solution of "
      "time step",
      &ssidynmono);

  DoubleParameter("RELAX_LIN_SOLVER_TOLERANCE", 1.0,
      "relax the tolerance of the linear solver in case it is an iterative solver by scaling the "
      "convergence tolerance with factor RELAX_LIN_SOLVER_TOLERANCE",
      &ssidynmono);

  IntParameter("RELAX_LIN_SOLVER_STEP", -1,
      "relax the tolerance of the linear solver within the first RELAX_LIN_SOLVER_STEP steps",
      &ssidynmono);

  /*----------------------------------------------------------------------*/
  /* parameters for SSI with manifold */
  /*----------------------------------------------------------------------*/

  Teuchos::ParameterList& ssidynmanifold = ssidyn.sublist("MANIFOLD", false,
      "Monolithic Structure Scalar Interaction with additional scalar transport on manifold");

  BoolParameter("ADD_MANIFOLD", "no", "activate additional manifold?", &ssidynmanifold);

  BoolParameter("MESHTYING_MANIFOLD", "no",
      "activate meshtying betweeen all manifold fields in case they intersect?", &ssidynmanifold);

  setStringToIntegralParameter<int>("INITIALFIELD", "zero_field",
      "Initial field for scalar transport on manifold",
      tuple<std::string>("zero_field", "field_by_function", "field_by_condition"),
      tuple<int>(INPAR::SCATRA::initfield_zero_field, INPAR::SCATRA::initfield_field_by_function,
          INPAR::SCATRA::initfield_field_by_condition),
      &ssidynmanifold);

  IntParameter("INITFUNCNO", -1, "function number for scalar transport on manifold initial field",
      &ssidynmanifold);

  IntParameter(
      "LINEAR_SOLVER", -1, "linear solver for scalar transport on manifold", &ssidynmanifold);

  BoolParameter("OUTPUT_INFLOW", "no",
      "write output of inflow of scatra manifold - scatra coupling into scatra manifold to csv "
      "file",
      &ssidynmanifold);

  /*----------------------------------------------------------------------*/
  /* parameters for SSI with elch */
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& ssidynelch = ssidyn.sublist(
      "ELCH", false, "Monolithic Structure Scalar Interaction with Elch as SCATRATIMINTTYPE");
  BoolParameter("INITPOTCALC", "No", "Automatically calculate initial field for electric potential",
      &ssidynelch);
}

/*--------------------------------------------------------------------
--------------------------------------------------------------------*/
void INPAR::SSI::SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist)
{
  using namespace INPUT;


  /*--------------------------------------------------------------------*/
  auto linessiplain = Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING LINE CONDITIONS",
      "SSICoupling", "SSI Coupling", DRT::Condition::SSICoupling, true, DRT::Condition::Line));
  auto surfssiplain = Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SURF CONDITIONS",
      "SSICoupling", "SSI Coupling", DRT::Condition::SSICoupling, true, DRT::Condition::Surface));
  auto volssiplain = Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING VOL CONDITIONS",
      "SSICoupling", "SSI Coupling", DRT::Condition::SSICoupling, true, DRT::Condition::Volume));

  // equip condition definitions with input file line components
  std::vector<Teuchos::RCP<INPUT::LineComponent>> ssicoupcomponentsplain;
  ssicoupcomponentsplain.emplace_back(Teuchos::rcp(new INPUT::IntComponent("coupling id")));

  // insert input file line components into condition definitions
  for (auto& ssicoupcomponentplain : ssicoupcomponentsplain)
  {
    linessiplain->AddComponent(ssicoupcomponentplain);
    surfssiplain->AddComponent(ssicoupcomponentplain);
    volssiplain->AddComponent(ssicoupcomponentplain);
  }

  condlist.push_back(linessiplain);
  condlist.push_back(surfssiplain);
  condlist.push_back(volssiplain);

  /*--------------------------------------------------------------------*/
  //! set solid dofset on scatra discretization
  auto linessi =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SOLIDTOSCATRA LINE CONDITIONS",
          "SSICouplingSolidToScatra", "SSI Coupling SolidToScatra",
          DRT::Condition::SSICouplingSolidToScatra, true, DRT::Condition::Line));
  auto surfssi =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SOLIDTOSCATRA SURF CONDITIONS",
          "SSICouplingSolidToScatra", "SSI Coupling SolidToScatra",
          DRT::Condition::SSICouplingSolidToScatra, true, DRT::Condition::Surface));
  auto volssi =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SOLIDTOSCATRA VOL CONDITIONS",
          "SSICouplingSolidToScatra", "SSI Coupling SolidToScatra",
          DRT::Condition::SSICouplingSolidToScatra, true, DRT::Condition::Volume));

  // equip condition definitions with input file line components
  std::vector<Teuchos::RCP<INPUT::LineComponent>> ssicoupcomponents;
  ssicoupcomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent("coupling id")));

  // insert input file line components into condition definitions
  for (auto& ssicoupcomponent : ssicoupcomponents)
  {
    linessi->AddComponent(ssicoupcomponent);
    surfssi->AddComponent(ssicoupcomponent);
    volssi->AddComponent(ssicoupcomponent);
  }

  condlist.push_back(linessi);
  condlist.push_back(surfssi);
  condlist.push_back(volssi);

  /*--------------------------------------------------------------------*/
  //! set scatra dofset on solid discretization
  auto linessi2 =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SCATRATOSOLID LINE CONDITIONS",
          "SSICouplingScatraToSolid", "SSI Coupling ScatraToSolid",
          DRT::Condition::SSICouplingScatraToSolid, true, DRT::Condition::Line));
  auto surfssi2 =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SCATRATOSOLID SURF CONDITIONS",
          "SSICouplingScatraToSolid", "SSI Coupling ScatraToSolid",
          DRT::Condition::SSICouplingScatraToSolid, true, DRT::Condition::Surface));
  auto volssi2 =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SCATRATOSOLID VOL CONDITIONS",
          "SSICouplingScatraToSolid", "SSI Coupling ScatraToSolid",
          DRT::Condition::SSICouplingScatraToSolid, true, DRT::Condition::Volume));

  // equip condition definitions with input file line components
  std::vector<Teuchos::RCP<INPUT::LineComponent>> ssicoupcomponents2;
  ssicoupcomponents2.emplace_back(Teuchos::rcp(new INPUT::IntComponent("coupling id")));

  // insert input file line components into condition definitions
  for (auto& ssicoupcomponent2 : ssicoupcomponents2)
  {
    linessi2->AddComponent(ssicoupcomponent2);
    surfssi2->AddComponent(ssicoupcomponent2);
    volssi2->AddComponent(ssicoupcomponent2);
  }

  condlist.push_back(linessi2);
  condlist.push_back(surfssi2);
  condlist.push_back(volssi2);

  /*--------------------------------------------------------------------*/
  // set ScaTra-Structure interaction interface meshtying condition
  auto pointssiinterfacemeshtying =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI INTERFACE MESHTYING POINT CONDITIONS",
          "SSIInterfaceMeshtying", "SSI Interface Meshtying", DRT::Condition::SSIInterfaceMeshtying,
          true, DRT::Condition::Point));
  auto linessiinterfacemeshtying =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI INTERFACE MESHTYING LINE CONDITIONS",
          "SSIInterfaceMeshtying", "SSI Interface Meshtying", DRT::Condition::SSIInterfaceMeshtying,
          true, DRT::Condition::Line));
  auto surfssiinterfacemeshtying =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI INTERFACE MESHTYING SURF CONDITIONS",
          "SSIInterfaceMeshtying", "SSI Interface Meshtying", DRT::Condition::SSIInterfaceMeshtying,
          true, DRT::Condition::Surface));

  // equip condition definitions with input file line components
  //
  // REMARK: it would be cleaner to also set a reference to the structural meshtying condition here
  // and not only to the S2ICoupling condition. Of course, then also the structural meshtying should
  // be used which could/should be the long-term goal. However, to date, a simple structural
  // meshtying version for matching node is implemented within the SSI framework and therefore no
  // reference is necessary.
  std::vector<Teuchos::RCP<INPUT::LineComponent>> ssiinterfacemeshtying;
  ssiinterfacemeshtying.emplace_back(Teuchos::rcp(new INPUT::IntComponent("ConditionID")));
  ssiinterfacemeshtying.emplace_back(Teuchos::rcp(new INPUT::SelectionComponent("interface side",
      "Undefined", Teuchos::tuple<std::string>("Undefined", "Slave", "Master"),
      Teuchos::tuple<int>(
          INPAR::S2I::side_undefined, INPAR::S2I::side_slave, INPAR::S2I::side_master))));
  ssiinterfacemeshtying.emplace_back(
      Teuchos::rcp(new INPUT::SeparatorComponent("S2I_KINETICS_ID")));
  ssiinterfacemeshtying.emplace_back(Teuchos::rcp(new INPUT::IntComponent("S2IKineticsID")));

  // insert input file line components into condition definitions
  for (auto& conditioncomponent : ssiinterfacemeshtying)
  {
    pointssiinterfacemeshtying->AddComponent(conditioncomponent);
    linessiinterfacemeshtying->AddComponent(conditioncomponent);
    surfssiinterfacemeshtying->AddComponent(conditioncomponent);
  }

  condlist.push_back(pointssiinterfacemeshtying);
  condlist.push_back(linessiinterfacemeshtying);
  condlist.push_back(surfssiinterfacemeshtying);

  /*--------------------------------------------------------------------*/
  // condition, where additional scatra field on manifold is created
  auto ssisurfacemanifold = Teuchos::rcp(new ConditionDefinition(
      "DESIGN SSI MANIFOLD SURF CONDITIONS", "SSISurfaceManifold", "scalar transport on manifold",
      DRT::Condition::SSISurfaceManifold, true, DRT::Condition::Surface));

  ssisurfacemanifold->AddComponent(Teuchos::rcp(new INPUT::IntComponent("ConditionID")));

  ssisurfacemanifold->AddComponent(Teuchos::rcp(new INPUT::SeparatorComponent("ImplType")));
  ssisurfacemanifold->AddComponent(
      Teuchos::rcp(new INPUT::SelectionComponent("ImplType", "Undefined",
          Teuchos::tuple<std::string>("Undefined", "Standard", "ElchElectrode", "ElchDiffCond"),
          Teuchos::tuple<int>(INPAR::SCATRA::impltype_undefined, INPAR::SCATRA::impltype_std,
              INPAR::SCATRA::impltype_elch_electrode, INPAR::SCATRA::impltype_elch_diffcond))));
  ssisurfacemanifold->AddComponent(Teuchos::rcp(new INPUT::SeparatorComponent("thickness")));
  ssisurfacemanifold->AddComponent(Teuchos::rcp(new INPUT::RealComponent("thickness")));

  condlist.emplace_back(ssisurfacemanifold);

  /*--------------------------------------------------------------------*/
  // initial field by condition for scatra on manifold
  auto surfmanifoldinitfields =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURF SCATRA MANIFOLD INITIAL FIELD CONDITIONS",
          "ScaTraManifoldInitfield", "Surface ScaTra Manifold Initfield",
          DRT::Condition::SurfaceInitfield, false, DRT::Condition::Surface));

  surfmanifoldinitfields->AddComponent(Teuchos::rcp(new INPUT::SelectionComponent("Field", "ScaTra",
      Teuchos::tuple<std::string>("ScaTra"), Teuchos::tuple<std::string>("ScaTra"))));

  surfmanifoldinitfields->AddComponent(Teuchos::rcp(new INPUT::IntVectorComponent("funct", 1)));

  condlist.emplace_back(surfmanifoldinitfields);

  /*--------------------------------------------------------------------*/
  // kinetics condition for flux scatra <-> scatra on manifold
  auto surfmanifoldkinetics =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI MANIFOLD KINETICS SURF CONDITIONS",
          "SSISurfaceManifoldKinetics", "kinetics model for coupling scatra <-> scatra on manifold",
          DRT::Condition::SSISurfaceManifoldKinetics, true, DRT::Condition::Surface));

  {
    surfmanifoldkinetics->AddComponent(Teuchos::rcp(new INPUT::IntComponent("ConditionID")));

    surfmanifoldkinetics->AddComponent(
        Teuchos::rcp(new INPUT::SeparatorComponent("ManifoldConditionID")));
    surfmanifoldkinetics->AddComponent(
        Teuchos::rcp(new INPUT::IntComponent("ManifoldConditionID")));

    std::map<int, std::pair<std::string, std::vector<Teuchos::RCP<INPUT::LineComponent>>>>
        kinetic_model_choices;
    {
      {
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

        kinetic_model_choices.emplace(INPAR::S2I::kinetics_constantinterfaceresistance,
            std::make_pair("ConstantInterfaceResistance", constantinterfaceresistance));
      }

      {
        // Butler-Volmer-reduced
        std::vector<Teuchos::RCP<INPUT::LineComponent>> butlervolmerreduced;
        // total number of existing scalars
        butlervolmerreduced.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("NUMSCAL")));
        butlervolmerreduced.emplace_back(Teuchos::rcp(new INPUT::IntComponent("numscal")));
        butlervolmerreduced.emplace_back(
            Teuchos::rcp(new INPUT::SeparatorComponent("STOICHIOMETRIES")));
        butlervolmerreduced.emplace_back(Teuchos::rcp(
            new INPUT::IntVectorComponent("stoichiometries", INPUT::LengthFromInt("numscal"))));

        butlervolmerreduced.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("E-")));
        butlervolmerreduced.emplace_back(Teuchos::rcp(new INPUT::IntComponent("e-")));
        butlervolmerreduced.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("K_R")));
        butlervolmerreduced.emplace_back(Teuchos::rcp(new INPUT::RealComponent("k_r")));
        butlervolmerreduced.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_A")));
        butlervolmerreduced.emplace_back(Teuchos::rcp(new INPUT::RealComponent("alpha_a")));
        butlervolmerreduced.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ALPHA_C")));
        butlervolmerreduced.emplace_back(Teuchos::rcp(new INPUT::RealComponent("alpha_c")));

        kinetic_model_choices.emplace(INPAR::S2I::kinetics_butlervolmerreduced,
            std::make_pair("Butler-VolmerReduced", butlervolmerreduced));
      }

      {
        std::vector<Teuchos::RCP<INPUT::LineComponent>> noflux;

        kinetic_model_choices.emplace(
            INPAR::S2I::kinetics_nointerfaceflux, std::make_pair("NoInterfaceFlux", noflux));
      }
    }

    surfmanifoldkinetics->AddComponent(
        Teuchos::rcp(new INPUT::SeparatorComponent("KINETIC_MODEL")));
    surfmanifoldkinetics->AddComponent(Teuchos::rcp(new INPUT::SwitchComponent(
        "kinetic model", INPAR::S2I::kinetics_constantinterfaceresistance, kinetic_model_choices)));
  }

  condlist.emplace_back(surfmanifoldkinetics);

  /*--------------------------------------------------------------------*/
  // Dirichlet conditions for scatra on manifold
  auto pointmanifolddirichlet = Teuchos::rcp(
      new ConditionDefinition("DESIGN POINT MANIFOLD DIRICH CONDITIONS", "ManifoldDirichlet",
          "Point Dirichlet", DRT::Condition::PointDirichlet, false, DRT::Condition::Point));
  auto linemanifolddirichlet = Teuchos::rcp(
      new ConditionDefinition("DESIGN LINE MANIFOLD DIRICH CONDITIONS", "ManifoldDirichlet",
          "Line Dirichlet", DRT::Condition::LineDirichlet, false, DRT::Condition::Line));
  auto surfmanifolddirichlet = Teuchos::rcp(
      new ConditionDefinition("DESIGN SURF MANIFOLD DIRICH CONDITIONS", "ManifoldDirichlet",
          "Surface Dirichlet", DRT::Condition::SurfaceDirichlet, false, DRT::Condition::Surface));

  const auto add_dirichlet_manifold_components = [](ConditionDefinition& definition)
  {
    definition.AddComponent(Teuchos::rcp(new INPUT::SeparatorComponent("NUMDOF")));
    definition.AddComponent(Teuchos::rcp(new INPUT::IntComponent("numdof")));

    definition.AddComponent(Teuchos::rcp(new INPUT::SeparatorComponent("ONOFF")));
    definition.AddComponent(
        Teuchos::rcp(new INPUT::IntVectorComponent("onoff", INPUT::LengthFromInt("numdof"))));

    definition.AddComponent(Teuchos::rcp(new INPUT::SeparatorComponent("VAL")));
    definition.AddComponent(
        Teuchos::rcp(new INPUT::RealVectorComponent("val", INPUT::LengthFromInt("numdof"))));

    definition.AddComponent(Teuchos::rcp(new INPUT::SeparatorComponent("FUNCT")));
    definition.AddComponent(
        Teuchos::rcp(new INPUT::IntVectorComponent("funct", INPUT::LengthFromInt("numdof"),
            {/*default=*/0, /*fortranstyle=*/false, /*noneallowed=*/true, /*optional=*/false})));
  };

  {
    add_dirichlet_manifold_components(*pointmanifolddirichlet);
    add_dirichlet_manifold_components(*linemanifolddirichlet);
    add_dirichlet_manifold_components(*surfmanifolddirichlet);
  }

  condlist.push_back(pointmanifolddirichlet);
  condlist.push_back(linemanifolddirichlet);
  condlist.push_back(surfmanifolddirichlet);

  /*--------------------------------------------------------------------*/
  // set ScaTra-Structure Interaction interface contact condition
  auto linessiinterfacecontact = Teuchos::rcp(new ConditionDefinition(
      "DESIGN SSI INTERFACE CONTACT LINE CONDITIONS", "SSIInterfaceContact",
      "SSI Interface Contact", DRT::Condition::SSIInterfaceContact, true, DRT::Condition::Line));
  auto surfssiinterfacecontact = Teuchos::rcp(new ConditionDefinition(
      "DESIGN SSI INTERFACE CONTACT SURF CONDITIONS", "SSIInterfaceContact",
      "SSI Interface Contact", DRT::Condition::SSIInterfaceContact, true, DRT::Condition::Surface));

  // equip condition definitions with input file line components
  std::vector<Teuchos::RCP<INPUT::LineComponent>> ssiinterfacecontact;
  ssiinterfacecontact.emplace_back(Teuchos::rcp(new INPUT::IntComponent("ConditionID")));
  ssiinterfacecontact.emplace_back(Teuchos::rcp(new INPUT::SelectionComponent("interface side",
      "Undefined", Teuchos::tuple<std::string>("Undefined", "Slave", "Master"),
      Teuchos::tuple<int>(
          INPAR::S2I::side_undefined, INPAR::S2I::side_slave, INPAR::S2I::side_master))));
  ssiinterfacecontact.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("S2I_KINETICS_ID")));
  ssiinterfacecontact.emplace_back(Teuchos::rcp(new INPUT::IntComponent("S2IKineticsID")));
  ssiinterfacecontact.emplace_back(
      Teuchos::rcp(new INPUT::SeparatorComponent("CONTACT_CONDITION_ID")));
  ssiinterfacecontact.emplace_back(Teuchos::rcp(new INPUT::IntComponent("ContactConditionID")));

  // insert input file line components into condition definitions
  for (const auto& ssiinterfacecontactcomponent : ssiinterfacecontact)
  {
    linessiinterfacecontact->AddComponent(ssiinterfacecontactcomponent);
    surfssiinterfacecontact->AddComponent(ssiinterfacecontactcomponent);
  }

  condlist.push_back(linessiinterfacecontact);
  condlist.push_back(surfssiinterfacecontact);
}

BACI_NAMESPACE_CLOSE
