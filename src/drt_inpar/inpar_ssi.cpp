/*----------------------------------------------------------------------*/
/*! \file
\brief input parameters for solid-scatra-interaction

\level 2


*/

/*----------------------------------------------------------------------*/

#include "drt_validparameters.H"
#include "inpar_ssi.H"
#include "inpar_s2i.H"
#include "inpar_scatra.H"
#include "../drt_lib/drt_conditiondefinition.H"
#include "../linalg/linalg_equilibrate.H"
#include "../linalg/linalg_sparseoperator.H"

void INPAR::SSI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
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
  setStringToIntegralParameter<LINALG::MatrixType>("MATRIXTYPE", "undefined",
      "type of global system matrix in global system of equations",
      tuple<std::string>("undefined", "block", "sparse"),
      tuple<LINALG::MatrixType>(LINALG::MatrixType::undefined, LINALG::MatrixType::block_field,
          LINALG::MatrixType::sparse),
      &ssidynmono);

  setStringToIntegralParameter<LINALG::EquilibrationMethod>("EQUILIBRATION", "none",
      "flag for equilibration of global system of equations",
      tuple<std::string>("none", "rows_full", "rows_maindiag", "columns_full", "columns_maindiag",
          "rowsandcolumns_full", "rowsandcolumns_maindiag", "local"),
      tuple<LINALG::EquilibrationMethod>(LINALG::EquilibrationMethod::none,
          LINALG::EquilibrationMethod::rows_full, LINALG::EquilibrationMethod::rows_maindiag,
          LINALG::EquilibrationMethod::columns_full, LINALG::EquilibrationMethod::columns_maindiag,
          LINALG::EquilibrationMethod::rowsandcolumns_full,
          LINALG::EquilibrationMethod::rowsandcolumns_maindiag, LINALG::EquilibrationMethod::local),
      &ssidynmono);

  setStringToIntegralParameter<LINALG::EquilibrationMethod>("EQUILIBRATION_STRUCTURE", "none",
      "flag for equilibration of structural equations",
      tuple<std::string>(
          "none", "rows_maindiag", "columns_maindiag", "rowsandcolumns_maindiag", "symmetry"),
      tuple<LINALG::EquilibrationMethod>(LINALG::EquilibrationMethod::none,
          LINALG::EquilibrationMethod::rows_maindiag, LINALG::EquilibrationMethod::columns_maindiag,
          LINALG::EquilibrationMethod::rowsandcolumns_maindiag,
          LINALG::EquilibrationMethod::symmetry),
      &ssidynmono);

  setStringToIntegralParameter<LINALG::EquilibrationMethod>("EQUILIBRATION_SCATRA", "none",
      "flag for equilibration of scatra equations",
      tuple<std::string>(
          "none", "rows_maindiag", "columns_maindiag", "rowsandcolumns_maindiag", "symmetry"),
      tuple<LINALG::EquilibrationMethod>(LINALG::EquilibrationMethod::none,
          LINALG::EquilibrationMethod::rows_maindiag, LINALG::EquilibrationMethod::columns_maindiag,
          LINALG::EquilibrationMethod::rowsandcolumns_maindiag,
          LINALG::EquilibrationMethod::symmetry),
      &ssidynmono);

  BoolParameter("PRINT_MAT_RHS_MAP_MATLAB", "no",
      "print system matrix, rhs vector, and full map to matlab readable file after solution of "
      "time step",
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
void INPAR::SSI::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;


  /*--------------------------------------------------------------------*/
  auto linessiplain = Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING LINE CONDITIONS",
      "SSICoupling", "SSI Coupling", DRT::Condition::SSICoupling, true, DRT::Condition::Line));
  auto surfssiplain = Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SURF CONDITIONS",
      "SSICoupling", "SSI Coupling", DRT::Condition::SSICoupling, true, DRT::Condition::Surface));
  auto volssiplain = Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING VOL CONDITIONS",
      "SSICoupling", "SSI Coupling", DRT::Condition::SSICoupling, true, DRT::Condition::Volume));

  // equip condition definitions with input file line components
  std::vector<Teuchos::RCP<ConditionComponent>> ssicoupcomponentsplain;
  ssicoupcomponentsplain.emplace_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

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
  std::vector<Teuchos::RCP<ConditionComponent>> ssicoupcomponents;
  ssicoupcomponents.emplace_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

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
  std::vector<Teuchos::RCP<ConditionComponent>> ssicoupcomponents2;
  ssicoupcomponents2.emplace_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

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
  std::vector<Teuchos::RCP<ConditionComponent>> ssiinterfacemeshtying;
  ssiinterfacemeshtying.emplace_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  ssiinterfacemeshtying.emplace_back(Teuchos::rcp(new StringConditionComponent("interface side",
      "Undefined", Teuchos::tuple<std::string>("Undefined", "Slave", "Master"),
      Teuchos::tuple<int>(
          INPAR::S2I::side_undefined, INPAR::S2I::side_slave, INPAR::S2I::side_master))));
  ssiinterfacemeshtying.emplace_back(
      Teuchos::rcp(new SeparatorConditionComponent("S2IKineticsID")));
  ssiinterfacemeshtying.emplace_back(Teuchos::rcp(new IntConditionComponent("S2IKineticsID")));

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

  ssisurfacemanifold->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  ssisurfacemanifold->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("ImplType")));
  ssisurfacemanifold->AddComponent(
      Teuchos::rcp(new StringConditionComponent("ImplType", "Undefined",
          Teuchos::tuple<std::string>("Undefined", "Standard", "ElchElectrode", "ElchDiffCond"),
          Teuchos::tuple<int>(INPAR::SCATRA::impltype_undefined, INPAR::SCATRA::impltype_std,
              INPAR::SCATRA::impltype_elch_electrode, INPAR::SCATRA::impltype_elch_diffcond))));

  condlist.emplace_back(ssisurfacemanifold);

  /*--------------------------------------------------------------------*/
  // initial field by condition for scatra on manifold
  auto surfmanifoldinitfields =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURF SCATRA MANIFOLD INITIAL FIELD CONDITIONS",
          "ScaTraManifoldInitfield", "Surface ScaTra Manifold Initfield",
          DRT::Condition::SurfaceInitfield, false, DRT::Condition::Surface));

  surfmanifoldinitfields->AddComponent(Teuchos::rcp(new StringConditionComponent("Field", "ScaTra",
      Teuchos::tuple<std::string>("ScaTra"), Teuchos::tuple<std::string>("ScaTra"))));

  surfmanifoldinitfields->AddComponent(Teuchos::rcp(new IntVectorConditionComponent("funct", 1)));

  condlist.emplace_back(surfmanifoldinitfields);

  /*--------------------------------------------------------------------*/
  // kinetics condition for flux scatra <-> scatra on manifold
  auto surfmanifoldkinetics =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI MANIFOLD KINETICS SURF CONDITIONS",
          "SSISurfaceManifoldKinetics", "kinetics model for coupling scatra <-> scatra on manifold",
          DRT::Condition::SSISurfaceManifoldKinetics, true, DRT::Condition::Surface));

  {
    surfmanifoldkinetics->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));

    surfmanifoldkinetics->AddComponent(
        Teuchos::rcp(new SeparatorConditionComponent("ManifoldConditionID")));
    surfmanifoldkinetics->AddComponent(
        Teuchos::rcp(new IntConditionComponent("ManifoldConditionID")));

    std::vector<Teuchos::RCP<CondCompBundle>> kineticmodels;
    {
      {
        std::vector<Teuchos::RCP<ConditionComponent>> constantinterfaceresistance;
        constantinterfaceresistance.emplace_back(
            Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
        constantinterfaceresistance.emplace_back(
            Teuchos::rcp(new IntVectorConditionComponent("onoff", 2)));

        constantinterfaceresistance.emplace_back(
            Teuchos::rcp(new SeparatorConditionComponent("resistance")));
        constantinterfaceresistance.emplace_back(
            Teuchos::rcp(new RealConditionComponent("resistance")));
        constantinterfaceresistance.emplace_back(new SeparatorConditionComponent("e-"));
        constantinterfaceresistance.emplace_back(new IntConditionComponent("e-"));

        kineticmodels.emplace_back(Teuchos::rcp(new CondCompBundle("ConstantInterfaceResistance",
            constantinterfaceresistance, INPAR::S2I::kinetics_constantinterfaceresistance)));
      }

      {
        // Butler-Volmer-reduced
        std::vector<Teuchos::RCP<ConditionComponent>> butlervolmerreduced;
        // total number of existing scalars
        butlervolmerreduced.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("numscal")));
        // string separator in front of integer stoichiometry vector in input file line
        std::vector<Teuchos::RCP<SeparatorConditionComponent>> intsepcomp;
        intsepcomp.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("stoichiometries")));
        // integer vector of stoichiometric coefficients
        std::vector<Teuchos::RCP<IntVectorConditionComponent>> intvectcomp;
        intvectcomp.emplace_back(
            Teuchos::rcp(new IntVectorConditionComponent("stoichiometries", 0)));
        // empty vector --> no separators for real vectors needed
        std::vector<Teuchos::RCP<SeparatorConditionComponent>> realsepcomp;
        // empty vector --> no real vectors needed
        std::vector<Teuchos::RCP<RealVectorConditionComponent>> realvectcomp;
        butlervolmerreduced.emplace_back(Teuchos::rcp(
            new IntRealBundle("stoichiometries", Teuchos::rcp(new IntConditionComponent("numscal")),
                intsepcomp, intvectcomp, realsepcomp, realvectcomp)));
        butlervolmerreduced.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("e-")));
        butlervolmerreduced.emplace_back(Teuchos::rcp(new IntConditionComponent("e-")));
        butlervolmerreduced.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("k_r")));
        butlervolmerreduced.emplace_back(Teuchos::rcp(new RealConditionComponent("k_r")));
        butlervolmerreduced.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_a")));
        butlervolmerreduced.emplace_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
        butlervolmerreduced.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_c")));
        butlervolmerreduced.emplace_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));

        kineticmodels.emplace_back(Teuchos::rcp(new CondCompBundle("Butler-VolmerReduced",
            butlervolmerreduced, INPAR::S2I::kinetics_butlervolmerreduced)));
      }

      {
        std::vector<Teuchos::RCP<ConditionComponent>> noflux;

        kineticmodels.emplace_back(Teuchos::rcp(
            new CondCompBundle("NoInterfaceFlux", noflux, INPAR::S2I::kinetics_nointerfaceflux)));
      }
    }

    surfmanifoldkinetics->AddComponent(
        Teuchos::rcp(new SeparatorConditionComponent("KineticModel")));
    surfmanifoldkinetics->AddComponent(
        Teuchos::rcp(new CondCompBundleSelector("kinetic model", kineticmodels)));
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

  std::vector<Teuchos::RCP<SeparatorConditionComponent>> dirichletintsepveccomponents;
  std::vector<Teuchos::RCP<IntVectorConditionComponent>> dirichletintveccomponents;
  std::vector<Teuchos::RCP<SeparatorConditionComponent>> dirichletrealsepveccomponents;
  std::vector<Teuchos::RCP<RealVectorConditionComponent>> dirichletrealveccomponents;
  std::vector<Teuchos::RCP<ConditionComponent>> dirichletbundcomponents;

  dirichletintsepveccomponents.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
  dirichletintveccomponents.emplace_back(Teuchos::rcp(new IntVectorConditionComponent("onoff", 1)));
  dirichletrealsepveccomponents.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("VAL")));
  dirichletrealveccomponents.emplace_back(Teuchos::rcp(new RealVectorConditionComponent("val", 1)));
  dirichletintsepveccomponents.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("FUNCT")));
  dirichletintveccomponents.emplace_back(
      Teuchos::rcp(new IntVectorConditionComponent("funct", 1, false, true, false)));

  dirichletbundcomponents.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("NUMDOF")));
  dirichletbundcomponents.emplace_back(Teuchos::rcp(new DirichletNeumannBundle("dirichbund",
      Teuchos::rcp(new IntConditionComponent("numdof")), dirichletintsepveccomponents,
      dirichletintveccomponents, dirichletrealsepveccomponents, dirichletrealveccomponents)));

  for (auto& dirichletbundcomponent : dirichletbundcomponents)
  {
    pointmanifolddirichlet->AddComponent(dirichletbundcomponent);
    linemanifolddirichlet->AddComponent(dirichletbundcomponent);
    surfmanifolddirichlet->AddComponent(dirichletbundcomponent);
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
  std::vector<Teuchos::RCP<ConditionComponent>> ssiinterfacecontact;
  ssiinterfacecontact.emplace_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  ssiinterfacecontact.emplace_back(Teuchos::rcp(new SeparatorConditionComponent("S2IKineticsID")));
  ssiinterfacecontact.emplace_back(Teuchos::rcp(new IntConditionComponent("S2IKineticsID")));
  ssiinterfacecontact.emplace_back(
      Teuchos::rcp(new SeparatorConditionComponent("ContactConditionID")));
  ssiinterfacecontact.emplace_back(Teuchos::rcp(new IntConditionComponent("ContactConditionID")));

  // insert input file line components into condition definitions
  for (const auto& ssiinterfacecontactcomponent : ssiinterfacecontact)
  {
    linessiinterfacecontact->AddComponent(ssiinterfacecontactcomponent);
    surfssiinterfacecontact->AddComponent(ssiinterfacecontactcomponent);
  }

  condlist.push_back(linessiinterfacecontact);
  condlist.push_back(surfssiinterfacecontact);
}
