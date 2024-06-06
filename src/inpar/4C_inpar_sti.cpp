/*----------------------------------------------------------------------*/
/*! \file
\brief input quantities and globally accessible enumerations for scatra-thermo interaction

\level 2


*/
/*----------------------------------------------------------------------*/
#include "4C_inpar_sti.hpp"

#include "4C_discretization_condition_definition.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_linalg_sparseoperator.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*------------------------------------------------------------------------*
 | set valid parameters for scatra-thermo interaction          fang 10/16 |
 *------------------------------------------------------------------------*/
void Inpar::STI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& stidyn = list->sublist(
      "STI DYNAMIC", false, "general control parameters for scatra-thermo interaction problems");

  // type of scalar transport time integration
  setStringToIntegralParameter<ScaTraTimIntType>("SCATRATIMINTTYPE", "Standard",
      "scalar transport time integration type is needed to instantiate correct scalar transport "
      "time integration scheme for scatra-thermo interaction problems",
      tuple<std::string>("Standard", "Elch"),
      tuple<ScaTraTimIntType>(ScaTraTimIntType::standard, ScaTraTimIntType::elch), &stidyn);

  // type of coupling between scatra and thermo fields
  setStringToIntegralParameter<CouplingType>("COUPLINGTYPE", "Undefined",
      "type of coupling between scatra and thermo fields",
      tuple<std::string>("Undefined", "Monolithic", "OneWay_ScatraToThermo",
          "OneWay_ThermoToScatra", "TwoWay_ScatraToThermo", "TwoWay_ScatraToThermo_Aitken",
          "TwoWay_ScatraToThermo_Aitken_Dofsplit", "TwoWay_ThermoToScatra",
          "TwoWay_ThermoToScatra_Aitken"),
      tuple<CouplingType>(CouplingType::undefined, CouplingType::monolithic,
          CouplingType::oneway_scatratothermo, CouplingType::oneway_thermotoscatra,
          CouplingType::twoway_scatratothermo, CouplingType::twoway_scatratothermo_aitken,
          CouplingType::twoway_scatratothermo_aitken_dofsplit, CouplingType::twoway_thermotoscatra,
          CouplingType::twoway_thermotoscatra_aitken),
      &stidyn);

  // specification of initial temperature field
  setStringToIntegralParameter<int>("THERMO_INITIALFIELD", "zero_field",
      "initial temperature field for scatra-thermo interaction problems",
      tuple<std::string>("zero_field", "field_by_function", "field_by_condition"),
      tuple<int>(Inpar::ScaTra::initfield_zero_field, Inpar::ScaTra::initfield_field_by_function,
          Inpar::ScaTra::initfield_field_by_condition),
      &stidyn);

  // function number for initial temperature field
  Core::UTILS::IntParameter("THERMO_INITFUNCNO", -1,
      "function number for initial temperature field for scatra-thermo interaction problems",
      &stidyn);

  // ID of linear solver for temperature field
  Core::UTILS::IntParameter(
      "THERMO_LINEAR_SOLVER", -1, "ID of linear solver for temperature field", &stidyn);

  // flag for double condensation of linear equations associated with temperature field
  Core::UTILS::BoolParameter("THERMO_CONDENSATION", "No",
      "flag for double condensation of linear equations associated with temperature field",
      &stidyn);

  /*----------------------------------------------------------------------*/
  // valid parameters for monolithic scatra-thermo interaction
  Teuchos::ParameterList& stidyn_monolithic = stidyn.sublist(
      "MONOLITHIC", false, "control parameters for monolithic scatra-thermo interaction problems");

  // ID of linear solver for global system of equations
  Core::UTILS::IntParameter("LINEAR_SOLVER", -1,
      "ID of linear solver for global system of equations", &stidyn_monolithic);

  // type of global system matrix in global system of equations
  setStringToIntegralParameter<Core::LinAlg::MatrixType>("MATRIXTYPE", "block",
      "type of global system matrix in global system of equations",
      tuple<std::string>("block", "sparse"),
      tuple<Core::LinAlg::MatrixType>(
          Core::LinAlg::MatrixType::block_condition, Core::LinAlg::MatrixType::sparse),
      &stidyn_monolithic);

  /*----------------------------------------------------------------------*/
  // valid parameters for partitioned scatra-thermo interaction
  Teuchos::ParameterList& stidyn_partitioned = stidyn.sublist("PARTITIONED", false,
      "control parameters for partitioned scatra-thermo interaction problems");

  // relaxation parameter
  Core::UTILS::DoubleParameter("OMEGA", 1., "relaxation parameter", &stidyn_partitioned);

  // maximum value of Aitken relaxation parameter
  Core::UTILS::DoubleParameter("OMEGAMAX", 0.,
      "maximum value of Aitken relaxation parameter (0.0 = no constraint)", &stidyn_partitioned);

  return;
}


/*------------------------------------------------------------------------*
 | set valid conditions for scatra-thermo interaction          fang 10/16 |
 *------------------------------------------------------------------------*/
void Inpar::STI::SetValidConditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  return;
}

FOUR_C_NAMESPACE_CLOSE
