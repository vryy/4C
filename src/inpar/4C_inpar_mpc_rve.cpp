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
#include "4C_io_linecomponent.hpp"

FOUR_C_NAMESPACE_OPEN
// set the mpc specific parameters
void Inpar::RveMpc::set_valid_parameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;


  Teuchos::ParameterList& mpc = list->sublist("MULTI POINT CONSTRAINTS", false, "");

  Teuchos::setStringToIntegralParameter<Inpar::RveMpc::RveReferenceDeformationDefinition>(
      "RVE_REFERENCE_POINTS", "automatic", "Method of definition of the reference points of an RVE",
      Teuchos::tuple<std::string>("automatic", "manual"),
      Teuchos::tuple<Inpar::RveMpc::RveReferenceDeformationDefinition>(automatic, manual), &mpc);

  Teuchos::setStringToIntegralParameter<Inpar::RveMpc::EnforcementStrategy>("ENFORCEMENT",
      "penalty_method", "Method to enforce the multi point constraint",
      Teuchos::tuple<std::string>("penalty_method", "lagrange_multiplier_method"),
      Teuchos::tuple<Inpar::RveMpc::EnforcementStrategy>(penalty, lagrangeMultiplier), &mpc);

  Teuchos::setDoubleParameter("PENALTY_PARAM", 1e5, "Value of the penalty parameter", &mpc);
}

// set mpc specific conditions
void Inpar::RveMpc::set_valid_conditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  // ================================================================================================
  Teuchos::RCP<Core::Conditions::ConditionDefinition> rve_lineperiodic_condition = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN LINE PERIODIC RVE 2D BOUNDARY CONDITIONS",
          "LinePeriodicRve", "definition of edges forming 2D periodic boundary conditions",
          Core::Conditions::LineRvePeriodic, false, Core::Conditions::geometry_type_line));

  add_named_selection_component(rve_lineperiodic_condition, "EDGE", "edge line id", "undefined",
      Teuchos::tuple<std::string>("x+", "x-", "y+", "y-", "undefined"),
      Teuchos::tuple<std::string>("x+", "x-", "y+", "y-", "undefined"), true);

  condlist.push_back(rve_lineperiodic_condition);

  // ================================================================================================
  Teuchos::RCP<Core::Conditions::ConditionDefinition> rve_surfperiodic_condition = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN SURF PERIODIC RVE 3D BOUNDARY CONDITIONS",
          "SurfacePeriodicRve", "definition of surfaces forming 3D periodic boundary conditions",
          Core::Conditions::SurfaceRvePeriodic, false, Core::Conditions::geometry_type_surface));

  add_named_selection_component(rve_surfperiodic_condition, "SURF", "surface id", "undefined",
      Teuchos::tuple<std::string>("x+", "x-", "y+", "y-", "z+", "z-", "undefined"),
      Teuchos::tuple<std::string>("x+", "x-", "y+", "y-", "z+", "z-", "undefined"), true);

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

  add_named_selection_component(rve_cornerpoint_condition, "POSITION", "position of reference node",
      "undefined", Teuchos::tuple<std::string>("N1L", "N1B", "N2", "N4", "N1", "N3", "undefined"),
      Teuchos::tuple<std::string>("N1L", "N1B", "N2", "N4", "N1", "N3", "undefined"), true);

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

  add_named_int(linear_ce, "EQUATION", "EQUATION");
  add_named_selection_component(linear_ce, "ADD", "degrees of freedom", "undefined",
      Teuchos::tuple<std::string>("dispx", "dispy", "undefined"),
      Teuchos::tuple<std::string>("dispx", "dispy", "undefined"), true);
  add_named_real(linear_ce, "COEFFICIENT");

  condlist.push_back(linear_ce);
  /*--------------------------------------------------------------------*/
}
FOUR_C_NAMESPACE_CLOSE