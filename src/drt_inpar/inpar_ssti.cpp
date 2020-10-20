/*----------------------------------------------------------------------*/
/*! \file
\brief input parameters for solid-scatra-thermo-interaction

\level 2

*/

/*----------------------------------------------------------------------*/

#include "drt_validparameters.H"
#include "inpar_ssti.H"
#include "inpar_scatra.H"
#include "../drt_lib/drt_conditiondefinition.H"
#include "../linalg/linalg_equilibrate.H"
#include "../linalg/linalg_sparseoperator.H"

void INPAR::SSTI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& sstidyn = list->sublist(
      "SSTI CONTROL", false, "Control paramters for scatra structure thermo interaction");

  DoubleParameter(
      "RESTARTEVRYTIME", 0, "write restart possibility every RESTARTEVRY steps", &sstidyn);
  IntParameter("RESTARTEVRY", 1, "write restart possibility every RESTARTEVRY steps", &sstidyn);
  IntParameter("NUMSTEP", 200, "maximum number of Timesteps", &sstidyn);
  DoubleParameter("MAXTIME", 1000.0, "total simulation time", &sstidyn);
  DoubleParameter("TIMESTEP", -1, "time step size dt", &sstidyn);
  DoubleParameter("RESULTSEVRYTIME", 0, "increment for writing solution", &sstidyn);
  IntParameter("RESULTSEVRY", 1, "increment for writing solution", &sstidyn);
  IntParameter("ITEMAX", 10, "maximum number of iterations over fields", &sstidyn);
  BoolParameter("SCATRA_FROM_RESTART_FILE", "No",
      "read scatra result from restart files (use option 'restartfromfile' during execution of "
      "baci)",
      &sstidyn);
  StringParameter(
      "SCATRA_FILENAME", "nil", "Control-file name for reading scatra results in SSTI", &sstidyn);
  setStringToIntegralParameter<SolutionScheme>("COUPALGO", "ssti_Monolithic",
      "Coupling strategies for SSTI solvers", tuple<std::string>("ssti_Monolithic"),
      tuple<INPAR::SSTI::SolutionScheme>(SolutionScheme::monolithic), &sstidyn);
  setStringToIntegralParameter<ScaTraTimIntType>("SCATRATIMINTTYPE", "Elch",
      "scalar transport time integration type is needed to instantiate correct scalar transport "
      "time integration scheme for ssi problems",
      tuple<std::string>("Elch"), tuple<ScaTraTimIntType>(ScaTraTimIntType::elch), &sstidyn);
  BoolParameter("ADAPTIVE_TIMESTEPPING", "no", "flag for adaptive time stepping", &sstidyn);

  /*----------------------------------------------------------------------*/
  /* parameters for monolithic SSTI                                       */
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& sstidynmono = sstidyn.sublist("MONOLITHIC", false,
      "Monolithic Structure Scalar Interaction\n"
      "Control section for monolithic SSI");
  DoubleParameter("ABSTOLRES", 1.0e-14,
      "absolute tolerance for deciding if global residual of nonlinear problem is already zero",
      &sstidynmono);
  DoubleParameter("CONVTOL", 1.0e-6,
      "tolerance for convergence check of Newton-Raphson iteration within monolithic SSI",
      &sstidynmono);
  IntParameter(
      "LINEAR_SOLVER", -1, "ID of linear solver for global system of equations", &sstidynmono);
  setStringToIntegralParameter<LINALG::MatrixType>("MATRIXTYPE", "undefined",
      "type of global system matrix in global system of equations",
      tuple<std::string>("undefined", "block", "sparse"),
      tuple<LINALG::MatrixType>(LINALG::MatrixType::undefined, LINALG::MatrixType::block_field,
          LINALG::MatrixType::sparse),
      &sstidynmono);
  setStringToIntegralParameter<LINALG::EquilibrationMethod>("EQUILIBRATION", "none",
      "flag for equilibration of global system of equations",
      tuple<std::string>(
          "none", "rows_full", "rows_maindiag", "rowsandcolumns_full", "rowsandcolumns_maindiag"),
      tuple<LINALG::EquilibrationMethod>(LINALG::EquilibrationMethod::none,
          LINALG::EquilibrationMethod::rows_full, LINALG::EquilibrationMethod::rows_maindiag,
          LINALG::EquilibrationMethod::rowsandcolumns_full,
          LINALG::EquilibrationMethod::rowsandcolumns_maindiag),
      &sstidynmono);

  /*----------------------------------------------------------------------*/
  /* parameters for thermo                                                */
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& thermodyn =
      sstidyn.sublist("THERMO", false, "Parameters for thermo subproblem");
  IntParameter("INITTHERMOFUNCT", -1, "initial function for thermo field", &thermodyn);
  IntParameter("LINEAR_SOLVER", -1, "linear solver for thermo field", &thermodyn);
  setStringToIntegralParameter<INPAR::SCATRA::InitialField>("INITIALFIELD", "field_by_function",
      "defines, how to set the initial field",
      tuple<std::string>("field_by_function", "field_by_condition"),
      tuple<INPAR::SCATRA::InitialField>(INPAR::SCATRA::InitialField::initfield_field_by_function,
          INPAR::SCATRA::InitialField::initfield_field_by_condition),
      &thermodyn);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void INPAR::SSTI::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // set Scalar-Structure-Thermo interaction interface meshtying condition
  Teuchos::RCP<ConditionDefinition> linesstiinterfacemeshtying =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSTI INTERFACE MESHTYING LINE CONDITIONS",
          "SSTIInterfaceMeshtying", "SSTI Interface Meshtying",
          DRT::Condition::SSTIInterfaceMeshtying, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfsstiinterfacemeshtying =
      Teuchos::rcp(new ConditionDefinition("DESIGN SSTI INTERFACE MESHTYING SURF CONDITIONS",
          "SSTIInterfaceMeshtying", "SSTI Interface Meshtying",
          DRT::Condition::SSTIInterfaceMeshtying, true, DRT::Condition::Surface));

  // equip condition definitions with input file line components
  std::vector<Teuchos::RCP<ConditionComponent>> sstiinterfacemeshtyingcomponents;
  sstiinterfacemeshtyingcomponents.emplace_back(
      Teuchos::rcp(new IntConditionComponent("ConditionID")));
  sstiinterfacemeshtyingcomponents.emplace_back(Teuchos::rcp(
      new StringConditionComponent("Side", "Master", Teuchos::tuple<std::string>("Master", "Slave"),
          Teuchos::tuple<std::string>("Master", "Slave"))));
  sstiinterfacemeshtyingcomponents.emplace_back(
      Teuchos::rcp(new SeparatorConditionComponent("S2ICouplingID")));
  sstiinterfacemeshtyingcomponents.emplace_back(
      Teuchos::rcp(new IntConditionComponent("S2ICouplingID")));

  // insert input file line components into condition definitions
  for (const auto& sstiinterfacemeshtyingcomponent : sstiinterfacemeshtyingcomponents)
  {
    linesstiinterfacemeshtying->AddComponent(sstiinterfacemeshtyingcomponent);
    surfsstiinterfacemeshtying->AddComponent(sstiinterfacemeshtyingcomponent);
  }

  condlist.emplace_back(linesstiinterfacemeshtying);
  condlist.emplace_back(surfsstiinterfacemeshtying);
}
