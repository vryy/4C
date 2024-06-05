/*----------------------------------------------------------------------*/
/*! \file
\brief input parameters for solid-scatra-thermo-interaction

\level 2

*/

/*----------------------------------------------------------------------*/

#include "4C_inpar_ssti.hpp"

#include "4C_discretization_condition_definition.hpp"
#include "4C_inpar_s2i.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_linalg_equilibrate.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

void INPAR::SSTI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& sstidyn = list->sublist(
      "SSTI CONTROL", false, "Control paramters for scatra structure thermo interaction");

  CORE::UTILS::DoubleParameter(
      "RESTARTEVRYTIME", 0, "write restart possibility every RESTARTEVRY steps", &sstidyn);
  CORE::UTILS::IntParameter(
      "RESTARTEVRY", 1, "write restart possibility every RESTARTEVRY steps", &sstidyn);
  CORE::UTILS::IntParameter("NUMSTEP", 200, "maximum number of Timesteps", &sstidyn);
  CORE::UTILS::DoubleParameter("MAXTIME", 1000.0, "total simulation time", &sstidyn);
  CORE::UTILS::DoubleParameter("TIMESTEP", -1, "time step size dt", &sstidyn);
  CORE::UTILS::DoubleParameter("RESULTSEVRYTIME", 0, "increment for writing solution", &sstidyn);
  CORE::UTILS::IntParameter("RESULTSEVRY", 1, "increment for writing solution", &sstidyn);
  CORE::UTILS::IntParameter("ITEMAX", 10, "maximum number of iterations over fields", &sstidyn);
  CORE::UTILS::BoolParameter("SCATRA_FROM_RESTART_FILE", "No",
      "read scatra result from restart files (use option 'restartfromfile' during execution of "
      "4C)",
      &sstidyn);
  CORE::UTILS::StringParameter(
      "SCATRA_FILENAME", "nil", "Control-file name for reading scatra results in SSTI", &sstidyn);
  setStringToIntegralParameter<SolutionScheme>("COUPALGO", "ssti_Monolithic",
      "Coupling strategies for SSTI solvers", tuple<std::string>("ssti_Monolithic"),
      tuple<INPAR::SSTI::SolutionScheme>(SolutionScheme::monolithic), &sstidyn);
  setStringToIntegralParameter<ScaTraTimIntType>("SCATRATIMINTTYPE", "Elch",
      "scalar transport time integration type is needed to instantiate correct scalar transport "
      "time integration scheme for ssi problems",
      tuple<std::string>("Elch"), tuple<ScaTraTimIntType>(ScaTraTimIntType::elch), &sstidyn);
  CORE::UTILS::BoolParameter(
      "ADAPTIVE_TIMESTEPPING", "no", "flag for adaptive time stepping", &sstidyn);

  /*----------------------------------------------------------------------*/
  /* parameters for monolithic SSTI                                       */
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& sstidynmono = sstidyn.sublist("MONOLITHIC", false,
      "Monolithic Structure Scalar Interaction\n"
      "Control section for monolithic SSI");
  CORE::UTILS::DoubleParameter("ABSTOLRES", 1.0e-14,
      "absolute tolerance for deciding if global residual of nonlinear problem is already zero",
      &sstidynmono);
  CORE::UTILS::DoubleParameter("CONVTOL", 1.0e-6,
      "tolerance for convergence check of Newton-Raphson iteration within monolithic SSI",
      &sstidynmono);
  CORE::UTILS::IntParameter(
      "LINEAR_SOLVER", -1, "ID of linear solver for global system of equations", &sstidynmono);
  setStringToIntegralParameter<CORE::LINALG::MatrixType>("MATRIXTYPE", "undefined",
      "type of global system matrix in global system of equations",
      tuple<std::string>("undefined", "block", "sparse"),
      tuple<CORE::LINALG::MatrixType>(CORE::LINALG::MatrixType::undefined,
          CORE::LINALG::MatrixType::block_field, CORE::LINALG::MatrixType::sparse),
      &sstidynmono);
  setStringToIntegralParameter<CORE::LINALG::EquilibrationMethod>("EQUILIBRATION", "none",
      "flag for equilibration of global system of equations",
      tuple<std::string>("none", "rows_full", "rows_maindiag", "rowsandcolumns_full",
          "rowsandcolumns_maindiag", "local"),
      tuple<CORE::LINALG::EquilibrationMethod>(CORE::LINALG::EquilibrationMethod::none,
          CORE::LINALG::EquilibrationMethod::rows_full,
          CORE::LINALG::EquilibrationMethod::rows_maindiag,
          CORE::LINALG::EquilibrationMethod::rowsandcolumns_full,
          CORE::LINALG::EquilibrationMethod::rowsandcolumns_maindiag,
          CORE::LINALG::EquilibrationMethod::local),
      &sstidynmono);
  setStringToIntegralParameter<CORE::LINALG::EquilibrationMethod>("EQUILIBRATION_STRUCTURE", "none",
      "flag for equilibration of structural equations",
      tuple<std::string>(
          "none", "rows_maindiag", "columns_maindiag", "rowsandcolumns_maindiag", "symmetry"),
      tuple<CORE::LINALG::EquilibrationMethod>(CORE::LINALG::EquilibrationMethod::none,
          CORE::LINALG::EquilibrationMethod::rows_maindiag,
          CORE::LINALG::EquilibrationMethod::columns_maindiag,
          CORE::LINALG::EquilibrationMethod::rowsandcolumns_maindiag,
          CORE::LINALG::EquilibrationMethod::symmetry),
      &sstidynmono);
  setStringToIntegralParameter<CORE::LINALG::EquilibrationMethod>("EQUILIBRATION_SCATRA", "none",
      "flag for equilibration of scatra equations",
      tuple<std::string>(
          "none", "rows_maindiag", "columns_maindiag", "rowsandcolumns_maindiag", "symmetry"),
      tuple<CORE::LINALG::EquilibrationMethod>(CORE::LINALG::EquilibrationMethod::none,
          CORE::LINALG::EquilibrationMethod::rows_maindiag,
          CORE::LINALG::EquilibrationMethod::columns_maindiag,
          CORE::LINALG::EquilibrationMethod::rowsandcolumns_maindiag,
          CORE::LINALG::EquilibrationMethod::symmetry),
      &sstidynmono);
  setStringToIntegralParameter<CORE::LINALG::EquilibrationMethod>("EQUILIBRATION_THERMO", "none",
      "flag for equilibration of scatra equations",
      tuple<std::string>(
          "none", "rows_maindiag", "columns_maindiag", "rowsandcolumns_maindiag", "symmetry"),
      tuple<CORE::LINALG::EquilibrationMethod>(CORE::LINALG::EquilibrationMethod::none,
          CORE::LINALG::EquilibrationMethod::rows_maindiag,
          CORE::LINALG::EquilibrationMethod::columns_maindiag,
          CORE::LINALG::EquilibrationMethod::rowsandcolumns_maindiag,
          CORE::LINALG::EquilibrationMethod::symmetry),
      &sstidynmono);
  CORE::UTILS::BoolParameter("EQUILIBRATION_INIT_SCATRA", "no",
      "use equilibration method of ScaTra to equilibrate initial calculation of potential",
      &sstidynmono);

  /*----------------------------------------------------------------------*/
  /* parameters for thermo                                                */
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& thermodyn =
      sstidyn.sublist("THERMO", false, "Parameters for thermo subproblem");
  CORE::UTILS::IntParameter("INITTHERMOFUNCT", -1, "initial function for thermo field", &thermodyn);
  CORE::UTILS::IntParameter("LINEAR_SOLVER", -1, "linear solver for thermo field", &thermodyn);
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
    std::vector<Teuchos::RCP<CORE::Conditions::ConditionDefinition>>& condlist)
{
  using namespace INPUT;

  /*--------------------------------------------------------------------*/
  // set Scalar-Structure-Thermo interaction interface meshtying condition
  Teuchos::RCP<CORE::Conditions::ConditionDefinition> linesstiinterfacemeshtying = Teuchos::rcp(
      new CORE::Conditions::ConditionDefinition("DESIGN SSTI INTERFACE MESHTYING LINE CONDITIONS",
          "SSTIInterfaceMeshtying", "SSTI Interface Meshtying",
          CORE::Conditions::SSTIInterfaceMeshtying, true, CORE::Conditions::geometry_type_line));
  Teuchos::RCP<CORE::Conditions::ConditionDefinition> surfsstiinterfacemeshtying = Teuchos::rcp(
      new CORE::Conditions::ConditionDefinition("DESIGN SSTI INTERFACE MESHTYING SURF CONDITIONS",
          "SSTIInterfaceMeshtying", "SSTI Interface Meshtying",
          CORE::Conditions::SSTIInterfaceMeshtying, true, CORE::Conditions::geometry_type_surface));

  // equip condition definitions with input file line components
  std::vector<Teuchos::RCP<INPUT::LineComponent>> sstiinterfacemeshtyingcomponents;
  sstiinterfacemeshtyingcomponents.emplace_back(
      Teuchos::rcp(new INPUT::IntComponent("ConditionID")));
  sstiinterfacemeshtyingcomponents.emplace_back(Teuchos::rcp(new INPUT::SelectionComponent(
      "interface side", "Undefined", Teuchos::tuple<std::string>("Undefined", "Slave", "Master"),
      Teuchos::tuple<int>(
          INPAR::S2I::side_undefined, INPAR::S2I::side_slave, INPAR::S2I::side_master))));
  sstiinterfacemeshtyingcomponents.emplace_back(
      Teuchos::rcp(new INPUT::SeparatorComponent("S2I_KINETICS_ID")));
  sstiinterfacemeshtyingcomponents.emplace_back(
      Teuchos::rcp(new INPUT::IntComponent("S2IKineticsID")));

  // insert input file line components into condition definitions
  for (const auto& sstiinterfacemeshtyingcomponent : sstiinterfacemeshtyingcomponents)
  {
    linesstiinterfacemeshtying->AddComponent(sstiinterfacemeshtyingcomponent);
    surfsstiinterfacemeshtying->AddComponent(sstiinterfacemeshtyingcomponent);
  }

  condlist.emplace_back(linesstiinterfacemeshtying);
  condlist.emplace_back(surfsstiinterfacemeshtying);
}

FOUR_C_NAMESPACE_CLOSE
