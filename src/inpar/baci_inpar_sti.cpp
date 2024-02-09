/*----------------------------------------------------------------------*/
/*! \file
\brief input quantities and globally accessible enumerations for scatra-thermo interaction

\level 2


*/
/*----------------------------------------------------------------------*/
#include "baci_inpar_sti.hpp"

#include "baci_inpar_scatra.hpp"
#include "baci_inpar_validparameters.hpp"
#include "baci_lib_conditiondefinition.hpp"
#include "baci_linalg_sparseoperator.hpp"

BACI_NAMESPACE_OPEN

/*------------------------------------------------------------------------*
 | set valid parameters for scatra-thermo interaction          fang 10/16 |
 *------------------------------------------------------------------------*/
void INPAR::STI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace INPUT;
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
      tuple<int>(INPAR::SCATRA::initfield_zero_field, INPAR::SCATRA::initfield_field_by_function,
          INPAR::SCATRA::initfield_field_by_condition),
      &stidyn);

  // function number for initial temperature field
  IntParameter("THERMO_INITFUNCNO", -1,
      "function number for initial temperature field for scatra-thermo interaction problems",
      &stidyn);

  // ID of linear solver for temperature field
  IntParameter("THERMO_LINEAR_SOLVER", -1, "ID of linear solver for temperature field", &stidyn);

  // flag for double condensation of linear equations associated with temperature field
  BoolParameter("THERMO_CONDENSATION", "No",
      "flag for double condensation of linear equations associated with temperature field",
      &stidyn);

  /*----------------------------------------------------------------------*/
  // valid parameters for monolithic scatra-thermo interaction
  Teuchos::ParameterList& stidyn_monolithic = stidyn.sublist(
      "MONOLITHIC", false, "control parameters for monolithic scatra-thermo interaction problems");

  // ID of linear solver for global system of equations
  IntParameter("LINEAR_SOLVER", -1, "ID of linear solver for global system of equations",
      &stidyn_monolithic);

  // type of global system matrix in global system of equations
  setStringToIntegralParameter<CORE::LINALG::MatrixType>("MATRIXTYPE", "block",
      "type of global system matrix in global system of equations",
      tuple<std::string>("block", "sparse"),
      tuple<CORE::LINALG::MatrixType>(
          CORE::LINALG::MatrixType::block_condition, CORE::LINALG::MatrixType::sparse),
      &stidyn_monolithic);

  /*----------------------------------------------------------------------*/
  // valid parameters for partitioned scatra-thermo interaction
  Teuchos::ParameterList& stidyn_partitioned = stidyn.sublist("PARTITIONED", false,
      "control parameters for partitioned scatra-thermo interaction problems");

  // relaxation parameter
  DoubleParameter("OMEGA", 1., "relaxation parameter", &stidyn_partitioned);

  // maximum value of Aitken relaxation parameter
  DoubleParameter("OMEGAMAX", 0.,
      "maximum value of Aitken relaxation parameter (0.0 = no constraint)", &stidyn_partitioned);

  return;
}


/*------------------------------------------------------------------------*
 | set valid conditions for scatra-thermo interaction          fang 10/16 |
 *------------------------------------------------------------------------*/
void INPAR::STI::SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist)
{
  using namespace INPUT;

  return;
}

BACI_NAMESPACE_CLOSE
