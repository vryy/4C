/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for multi point constraints used for periodic boundary conditions
 for representative volume elements (RVEs)
\level 3
*/
/*----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_inpar_mpc_rve.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_inpar_validparameters.hpp"

FOUR_C_NAMESPACE_OPEN
// set the mpc specific parameters
void Inpar::RveMpc::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
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
void Inpar::RveMpc::SetValidConditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  // ================================================================================================
  Teuchos::RCP<Core::Conditions::ConditionDefinition> rve_lineperiodic_condition = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN LINE PERIODIC RVE 2D BOUNDARY CONDITIONS",
          "LinePeriodicRve", "definition of edges forming 2D periodic boundary conditions",
          Core::Conditions::LineRvePeriodic, false, Core::Conditions::geometry_type_line));

  rve_lineperiodic_condition->AddComponent(Teuchos::rcp(new Input::SeparatorComponent("EDGE")));

  rve_lineperiodic_condition->AddComponent(Teuchos::rcp(new Input::SelectionComponent("EdgeLineId",
      "undefined", Teuchos::tuple<std::string>("x+", "x-", "y+", "y-", "undefined"),
      Teuchos::tuple<std::string>("x+", "x-", "y+", "y-", "undefined"), true)));

  condlist.push_back(rve_lineperiodic_condition);

  // ================================================================================================
  Teuchos::RCP<Core::Conditions::ConditionDefinition> rve_surfperiodic_condition = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN SURF PERIODIC RVE 3D BOUNDARY CONDITIONS",
          "SurfacePeriodicRve", "definition of surfaces forming 3D periodic boundary conditions",
          Core::Conditions::SurfaceRvePeriodic, false, Core::Conditions::geometry_type_surface));

  rve_surfperiodic_condition->AddComponent(Teuchos::rcp(new Input::SeparatorComponent("SURF")));

  rve_surfperiodic_condition->AddComponent(Teuchos::rcp(new Input::SelectionComponent("SurfId",
      "undefined", Teuchos::tuple<std::string>("x+", "x-", "y+", "y-", "z+", "z-", "undefined"),
      Teuchos::tuple<std::string>("x+", "x-", "y+", "y-", "z+", "z-", "undefined"), true)));

  condlist.push_back(rve_surfperiodic_condition);

  // ================================================================================================
  Teuchos::RCP<Core::Conditions::ConditionDefinition> rve_cornerpoint_condition =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT PERIODIC RVE 2D BOUNDARY REFERENCE CONDITIONS",
          "PointPeriodicRveReferenceNode",
          "definition of reference points defining the reference vector of the periodic boundary"
          "condition -  only required if RVE_REFERENCE_POINTS = automatic",
          Core::Conditions::PointRvePeriodicReference, false,
          Core::Conditions::geometry_type_point));

  rve_cornerpoint_condition->AddComponent(Teuchos::rcp(new Input::SeparatorComponent("POSITION")));

  rve_cornerpoint_condition->AddComponent(
      Teuchos::rcp(new Input::SelectionComponent("referenceNode", "undefined",
          Teuchos::tuple<std::string>("N1L", "N1B", "N2", "N4", "N1", "N3", "undefined"),
          Teuchos::tuple<std::string>("N1L", "N1B", "N2", "N4", "N1", "N3", "undefined"), true)));

  condlist.push_back(rve_cornerpoint_condition);

  // ================================================================================================
  Teuchos::RCP<Core::Conditions::ConditionDefinition> linear_ce =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN POINT COUPLED DOF EQUATION CONDITIONS", "PointLinearCoupledEquation",
          "definition of the term of a linear couple equation coupling different degrees of "
          "freedom in "
          "2d",
          Core::Conditions::PointLinearCoupledEquation, false,
          Core::Conditions::geometry_type_point));

  linear_ce->AddComponent(Teuchos::rcp(new Input::SeparatorComponent("EQUATION")));

  linear_ce->AddComponent(Teuchos::rcp(new Input::IntComponent("EQUATION_ID")));

  linear_ce->AddComponent(Teuchos::rcp(new Input::SeparatorComponent("ADD")));

  linear_ce->AddComponent(Teuchos::rcp(new Input::SelectionComponent("DOF", "undefined",
      Teuchos::tuple<std::string>("dispx", "dispy", "undefined"),
      Teuchos::tuple<std::string>("dispx", "dispy", "undefined"), true)));

  linear_ce->AddComponent(Teuchos::rcp(new Input::SeparatorComponent("COEFFICIENT")));

  linear_ce->AddComponent(Teuchos::rcp(new Input::RealComponent("COEFFICIENT")));

  condlist.push_back(linear_ce);
  /*--------------------------------------------------------------------*/
}
FOUR_C_NAMESPACE_CLOSE