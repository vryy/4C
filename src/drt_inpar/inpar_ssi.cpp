/*----------------------------------------------------------------------*/
/*! \file
\brief input parameters for solid-scatra-interaction

\level 2


*/

/*----------------------------------------------------------------------*/

#include "drt_validparameters.H"
#include "inpar_ssi.H"
#include "../drt_lib/drt_conditiondefinition.H"
#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_equilibrate.H"

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
  setStringToIntegralParameter<int>("FIELDCOUPLING", "volume_matching",
      "Type of coupling strategy between fields",
      tuple<std::string>("volume_matching", "volume_nonmatching", "boundary_nonmatching",
          "volumeboundary_matching"),
      tuple<int>(coupling_volume_match, coupling_volume_nonmatch, coupling_boundary_nonmatch,
          coupling_volumeboundary_match),
      &ssidyn);

  // Coupling strategy for SSI solvers
  setStringToIntegralParameter<int>("COUPALGO", "ssi_IterStagg",
      "Coupling strategies for SSI solvers",
      tuple<std::string>("ssi_OneWay_ScatraToSolid", "ssi_OneWay_SolidToScatra",
          //                                "ssi_SequStagg_ScatraToSolid",
          //                                "ssi_SequStagg_SolidToScatra",
          "ssi_IterStagg", "ssi_IterStaggFixedRel_ScatraToSolid",
          "ssi_IterStaggFixedRel_SolidToScatra", "ssi_IterStaggAitken_ScatraToSolid",
          "ssi_IterStaggAitken_SolidToScatra", "ssi_Monolithic"),
      tuple<int>(ssi_OneWay_ScatraToSolid, ssi_OneWay_SolidToScatra,
          //                                ssi_SequStagg_ScatraToSolid,
          //                                ssi_SequStagg_SolidToScatra,
          ssi_IterStagg, ssi_IterStaggFixedRel_ScatraToSolid, ssi_IterStaggFixedRel_SolidToScatra,
          ssi_IterStaggAitken_ScatraToSolid, ssi_IterStaggAitken_SolidToScatra, ssi_Monolithic),
      &ssidyn);

  // type of scalar transport time integration
  setStringToIntegralParameter<int>("SCATRATIMINTTYPE", "Standard",
      "scalar transport time integration type is needed to instantiate correct scalar transport "
      "time integration scheme for ssi problems",
      tuple<std::string>("Standard", "Cardiac_Monodomain", "Elch"),
      tuple<int>(INPAR::SSI::scatratiminttype_standard,
          INPAR::SSI::scatratiminttype_cardiac_monodomain, INPAR::SSI::scatratiminttype_elch),
      &ssidyn);

  // Restart from Structure problem instead of SSI
  BoolParameter("RESTART_FROM_STRUCTURE", "no",
      "restart from structure problem (e.g. from prestress calculations) instead of ssi", &ssidyn);

  // Adaptive time stepping
  BoolParameter("ADAPTIVE_TIMESTEPPING", "no", "flag for adaptive time stepping", &ssidyn);

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
          "rowsandcolumns_full", "rowsandcolumns_maindiag"),
      tuple<LINALG::EquilibrationMethod>(LINALG::EquilibrationMethod::none,
          LINALG::EquilibrationMethod::rows_full, LINALG::EquilibrationMethod::rows_maindiag,
          LINALG::EquilibrationMethod::columns_full, LINALG::EquilibrationMethod::columns_maindiag,
          LINALG::EquilibrationMethod::rowsandcolumns_full,
          LINALG::EquilibrationMethod::rowsandcolumns_maindiag),
      &ssidynmono);
}



void INPAR::SSI::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;


  /*--------------------------------------------------------------------*/
  Teuchos::RCP<ConditionDefinition> linessiplain =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING LINE CONDITIONS", "SSICoupling",
          "SSI Coupling", DRT::Condition::SSICoupling, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfssiplain =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SURF CONDITIONS", "SSICoupling",
          "SSI Coupling", DRT::Condition::SSICoupling, true, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volssiplain =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING VOL CONDITIONS", "SSICoupling",
          "SSI Coupling", DRT::Condition::SSICoupling, true, DRT::Condition::Volume));

  // equip condition definitions with input file line components
  std::vector<Teuchos::RCP<ConditionComponent>> ssicoupcomponentsplain;
  ssicoupcomponentsplain.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

  // insert input file line components into condition definitions
  for (unsigned i = 0; i < ssicoupcomponentsplain.size(); ++i)
  {
    linessiplain->AddComponent(ssicoupcomponentsplain[i]);
    surfssiplain->AddComponent(ssicoupcomponentsplain[i]);
    volssiplain->AddComponent(ssicoupcomponentsplain[i]);
  }

  condlist.push_back(linessiplain);
  condlist.push_back(surfssiplain);
  condlist.push_back(volssiplain);

  /*--------------------------------------------------------------------*/
  //! set solid dofset on scatra discretization
  Teuchos::RCP<ConditionDefinition> linessi =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SOLIDTOSCATRA LINE CONDITIONS",
          "SSICouplingSolidToScatra", "SSI Coupling SolidToScatra",
          DRT::Condition::SSICouplingSolidToScatra, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfssi =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SOLIDTOSCATRA SURF CONDITIONS",
          "SSICouplingSolidToScatra", "SSI Coupling SolidToScatra",
          DRT::Condition::SSICouplingSolidToScatra, true, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volssi =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SOLIDTOSCATRA VOL CONDITIONS",
          "SSICouplingSolidToScatra", "SSI Coupling SolidToScatra",
          DRT::Condition::SSICouplingSolidToScatra, true, DRT::Condition::Volume));

  // equip condition definitions with input file line components
  std::vector<Teuchos::RCP<ConditionComponent>> ssicoupcomponents;
  ssicoupcomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

  // insert input file line components into condition definitions
  for (unsigned i = 0; i < ssicoupcomponents.size(); ++i)
  {
    linessi->AddComponent(ssicoupcomponents[i]);
    surfssi->AddComponent(ssicoupcomponents[i]);
    volssi->AddComponent(ssicoupcomponents[i]);
  }

  condlist.push_back(linessi);
  condlist.push_back(surfssi);
  condlist.push_back(volssi);

  /*--------------------------------------------------------------------*/
  //! set scatra dofset on solid discretization
  Teuchos::RCP<ConditionDefinition> linessi2 =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SCATRATOSOLID LINE CONDITIONS",
          "SSICouplingScatraToSolid", "SSI Coupling ScatraToSolid",
          DRT::Condition::SSICouplingScatraToSolid, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfssi2 =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SCATRATOSOLID SURF CONDITIONS",
          "SSICouplingScatraToSolid", "SSI Coupling ScatraToSolid",
          DRT::Condition::SSICouplingScatraToSolid, true, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volssi2 =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI COUPLING SCATRATOSOLID VOL CONDITIONS",
          "SSICouplingScatraToSolid", "SSI Coupling ScatraToSolid",
          DRT::Condition::SSICouplingScatraToSolid, true, DRT::Condition::Volume));

  // equip condition definitions with input file line components
  std::vector<Teuchos::RCP<ConditionComponent>> ssicoupcomponents2;
  ssicoupcomponents2.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

  // insert input file line components into condition definitions
  for (unsigned i = 0; i < ssicoupcomponents2.size(); ++i)
  {
    linessi2->AddComponent(ssicoupcomponents2[i]);
    surfssi2->AddComponent(ssicoupcomponents2[i]);
    volssi2->AddComponent(ssicoupcomponents2[i]);
  }

  condlist.push_back(linessi2);
  condlist.push_back(surfssi2);
  condlist.push_back(volssi2);

  /*--------------------------------------------------------------------*/
  // set Scalar-Structure interaction interface meshtying condition
  Teuchos::RCP<ConditionDefinition> linessiinterfacemeshtying =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI INTERFACE MESHTYING LINE CONDITIONS",
          "SSIInterfaceMeshtying", "SSI Interface Meshtying", DRT::Condition::SSIInterfaceMeshtying,
          true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfssiinterfacemeshtying =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSI INTERFACE MESHTYING SURF CONDITIONS",
          "SSIInterfaceMeshtying", "SSI Interface Meshtying", DRT::Condition::SSIInterfaceMeshtying,
          true, DRT::Condition::Surface));

  // equip condition definitions with input file line components
  //
  // REMARK: it would be cleaner to also set a reference to the structural meshtying condition here
  // and not only to the S2ICoupling condition. Of course, then also the structural meshtying should
  // be used which could/should be the long-term goal. However, to date, a simple structural
  // meshtying version for matching node is implemented within the SSI framework and therefore no
  // reference is neccessary.
  std::vector<Teuchos::RCP<ConditionComponent>> ssiinterfacemeshtying;
  ssiinterfacemeshtying.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  ssiinterfacemeshtying.push_back(Teuchos::rcp(
      new StringConditionComponent("Side", "Master", Teuchos::tuple<std::string>("Master", "Slave"),
          Teuchos::tuple<std::string>("Master", "Slave"))));
  ssiinterfacemeshtying.push_back(Teuchos::rcp(new SeparatorConditionComponent("S2ICouplingID")));
  ssiinterfacemeshtying.push_back(Teuchos::rcp(new IntConditionComponent("S2ICouplingID")));

  // insert input file line components into condition definitions
  for (unsigned i = 0; i < ssiinterfacemeshtying.size(); ++i)
  {
    linessiinterfacemeshtying->AddComponent(ssiinterfacemeshtying[i]);
    surfssiinterfacemeshtying->AddComponent(ssiinterfacemeshtying[i]);
  }

  condlist.push_back(linessiinterfacemeshtying);
  condlist.push_back(surfssiinterfacemeshtying);
}
