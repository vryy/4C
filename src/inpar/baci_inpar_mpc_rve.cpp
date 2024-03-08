/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for multi point constraints used for periodic boundary conditions
 for representative volume elements (RVEs)
\level 3
*/
/*----------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_inpar_mpc_rve.hpp"

#include "baci_inpar_validparameters.hpp"
#include "baci_lib_conditiondefinition.hpp"

BACI_NAMESPACE_OPEN
// set the mpc specific parameters
void INPAR::RVE_MPC::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace INPUT;
  using Teuchos::setStringToIntegralParameter;


  Teuchos::ParameterList& mpc = list->sublist("MULTI POINT CONSTRAINTS", false, "");

  Teuchos::setStringToIntegralParameter<int>("RVE_REFERENCE_POINTS", "automatic",
      "Method of definition of the reference points of an RVE",
      Teuchos::tuple<std::string>("automatic", "manual"), Teuchos::tuple<int>(automatic, manual),
      &mpc);

  Teuchos::setStringToIntegralParameter<int>("ENFORCEMENT", "penalty_method",
      "Method to enforce the multi point constraint",
      Teuchos::tuple<std::string>("penalty_method", "lagrange_multiplier_method"),
      Teuchos::tuple<int>(penalty, lagrangeMultiplier), &mpc);

  Teuchos::setDoubleParameter("PENALTY_PARAM", 1e5, "Value of the penalty parameter", &mpc);
}

// set mpc specific conditions
void INPAR::RVE_MPC::SetValidConditions(
    std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist)
{
  using namespace INPUT;

  // ================================================================================================
  Teuchos::RCP<ConditionDefinition> rve_lineperiodic_condition =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE PERIODIC RVE 2D BOUNDARY CONDITIONS",
          "LinePeriodicRve", "definition of edges forming 2D periodic boundary conditions",
          DRT::Condition::LineRvePeriodic, false, DRT::Condition::Line));

  rve_lineperiodic_condition->AddComponent(Teuchos::rcp(new INPUT::SeparatorComponent("EDGE")));

  rve_lineperiodic_condition->AddComponent(Teuchos::rcp(new INPUT::SelectionComponent("EdgeLineId",
      "undefined", Teuchos::tuple<std::string>("x+", "x-", "y+", "y-", "undefined"),
      Teuchos::tuple<std::string>("x+", "x-", "y+", "y-", "undefined"), true)));

  condlist.push_back(rve_lineperiodic_condition);

  // ================================================================================================
  Teuchos::RCP<ConditionDefinition> rve_surfperiodic_condition =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURF PERIODIC RVE 3D BOUNDARY CONDITIONS",
          "SurfacePeriodicRve", "definition of surfaces forming 3D periodic boundary conditions",
          DRT::Condition::SurfaceRvePeriodic, false, DRT::Condition::Surface));

  rve_surfperiodic_condition->AddComponent(Teuchos::rcp(new INPUT::SeparatorComponent("SURF")));

  rve_surfperiodic_condition->AddComponent(Teuchos::rcp(new INPUT::SelectionComponent("SurfId",
      "undefined", Teuchos::tuple<std::string>("x+", "x-", "y+", "y-", "z+", "z-", "undefined"),
      Teuchos::tuple<std::string>("x+", "x-", "y+", "y-", "z+", "z-", "undefined"), true)));

  condlist.push_back(rve_surfperiodic_condition);

  // ================================================================================================
  Teuchos::RCP<ConditionDefinition> rve_cornerpoint_condition = Teuchos::rcp(
      new ConditionDefinition("DESIGN POINT PERIODIC RVE 2D BOUNDARY REFERENCE CONDITIONS",
          "PointPeriodicRveReferenceNode",
          "definition of reference points defining the reference vector of the periodic boundary"
          "condition -  only required if RVE_REFERENCE_POINTS = automatic",
          DRT::Condition::PointRvePeriodicReference, false, DRT::Condition::Point));

  rve_cornerpoint_condition->AddComponent(Teuchos::rcp(new INPUT::SeparatorComponent("POSITION")));

  rve_cornerpoint_condition->AddComponent(
      Teuchos::rcp(new INPUT::SelectionComponent("referenceNode", "undefined",
          Teuchos::tuple<std::string>("N1L", "N1B", "N2", "N4", "N1", "N3", "undefined"),
          Teuchos::tuple<std::string>("N1L", "N1B", "N2", "N4", "N1", "N3", "undefined"), true)));

  condlist.push_back(rve_cornerpoint_condition);

  // ================================================================================================
  Teuchos::RCP<ConditionDefinition> linear_ce = Teuchos::rcp(new ConditionDefinition(
      "DESIGN POINT COUPLED DOF EQUATION CONDITIONS", "PointLinearCoupledEquation",
      "definition of the term of a linear couple equation coupling different degrees of freedom in "
      "2d",
      DRT::Condition::PointLinearCoupledEquation, false, DRT::Condition::Point));

  linear_ce->AddComponent(Teuchos::rcp(new INPUT::SeparatorComponent("EQUATION")));

  linear_ce->AddComponent(Teuchos::rcp(new INPUT::IntComponent("EQUATION_ID")));

  linear_ce->AddComponent(Teuchos::rcp(new INPUT::SeparatorComponent("ADD")));

  linear_ce->AddComponent(Teuchos::rcp(new INPUT::SelectionComponent("DOF", "undefined",
      Teuchos::tuple<std::string>("dispx", "dispy", "undefined"),
      Teuchos::tuple<std::string>("dispx", "dispy", "undefined"), true)));

  linear_ce->AddComponent(Teuchos::rcp(new INPUT::SeparatorComponent("COEFFICIENT")));

  linear_ce->AddComponent(Teuchos::rcp(new INPUT::RealComponent("COEFFICIENT")));

  condlist.push_back(linear_ce);
  /*--------------------------------------------------------------------*/
}
BACI_NAMESPACE_CLOSE