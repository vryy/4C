/*----------------------------------------------------------------------*/
/*! \file


\brief Input parameters for ALE mesh motion

\level 2
*/

/*----------------------------------------------------------------------*/

#include "4C_inpar_ale.hpp"

#include "4C_discretization_condition_definition.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::ALE::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& adyn = list->sublist("ALE DYNAMIC", false, "");

  Core::UTILS::DoubleParameter("TIMESTEP", 0.1, "time step size", &adyn);
  Core::UTILS::IntParameter("NUMSTEP", 41, "max number of time steps", &adyn);
  Core::UTILS::DoubleParameter("MAXTIME", 4.0, "max simulation time", &adyn);

  setStringToIntegralParameter<int>("ALE_TYPE", "solid", "ale mesh movement algorithm",
      tuple<std::string>("solid", "solid_linear", "laplace_material", "laplace_spatial",
          "springs_material", "springs_spatial"),
      tuple<int>(solid, solid_linear, laplace_material, laplace_spatial, springs_material,
          springs_spatial),
      &adyn);

  Core::UTILS::BoolParameter("ASSESSMESHQUALITY", "no",
      "Evaluate element quality measure according to [Oddy et al. 1988]", &adyn);

  Core::UTILS::BoolParameter("UPDATEMATRIX", "no",
      "Update stiffness matrix in every time step (only for linear/material strategies)", &adyn);

  Core::UTILS::IntParameter("MAXITER", 1, "Maximum number of newton iterations.", &adyn);
  Core::UTILS::DoubleParameter(
      "TOLRES", 1.0e-06, "Absolute tolerance for length scaled L2 residual norm ", &adyn);
  Core::UTILS::DoubleParameter(
      "TOLDISP", 1.0e-06, "Absolute tolerance for length scaled L2 increment norm ", &adyn);

  Core::UTILS::IntParameter("NUM_INITSTEP", 0, "", &adyn);
  Core::UTILS::IntParameter("RESTARTEVRY", 1, "write restart data every RESTARTEVRY steps", &adyn);
  Core::UTILS::IntParameter("RESULTSEVRY", 0, "write results every RESULTSTEVRY steps", &adyn);
  setStringToIntegralParameter<int>("DIVERCONT", "continue",
      "What to do if nonlinear solver does not converge?", tuple<std::string>("stop", "continue"),
      tuple<int>(divcont_stop, divcont_continue), &adyn);

  setStringToIntegralParameter<int>("MESHTYING", "no",
      "Flag to (de)activate mesh tying and mesh sliding algorithm",
      tuple<std::string>("no", "meshtying", "meshsliding"),
      tuple<int>(no_meshtying, meshtying, meshsliding), &adyn);

  // Initial displacement
  setStringToIntegralParameter<int>("INITIALDISP", "zero_displacement",
      "Initial displacement for structure problem",
      tuple<std::string>("zero_displacement", "displacement_by_function"),
      tuple<int>(initdisp_zero_disp, initdisp_disp_by_function), &adyn);

  // Function to evaluate initial displacement
  Core::UTILS::IntParameter("STARTFUNCNO", -1, "Function for Initial displacement", &adyn);

  // linear solver id used for scalar ale problems
  Core::UTILS::IntParameter(
      "LINEAR_SOLVER", -1, "number of linear solver used for ale problems...", &adyn);
}



void Inpar::ALE::SetValidConditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  /*--------------------------------------------------------------------*/
  // Ale update boundary condition

  std::vector<Teuchos::RCP<Input::LineComponent>> aleupdatecomponents;

  aleupdatecomponents.push_back(Teuchos::rcp(new Input::SeparatorComponent("COUPLING")));
  aleupdatecomponents.push_back(Teuchos::rcp(new Input::SelectionComponent("coupling", "lagrange",
      Teuchos::tuple<std::string>("lagrange", "heightfunction", "sphereHeightFunction",
          "meantangentialvelocity", "meantangentialvelocityscaled"),
      Teuchos::tuple<std::string>("lagrange", "heightfunction", "sphereHeightFunction",
          "meantangentialvelocity", "meantangentialvelocityscaled"),
      true)));

  aleupdatecomponents.push_back(Teuchos::rcp(new Input::SeparatorComponent("VAL")));
  aleupdatecomponents.push_back(Teuchos::rcp(new Input::RealComponent("val")));

  aleupdatecomponents.push_back(Teuchos::rcp(new Input::SeparatorComponent("NODENORMALFUNCT")));
  aleupdatecomponents.push_back(Teuchos::rcp(new Input::IntComponent("nodenormalfunct")));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> linealeupdate =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN ALE UPDATE LINE CONDITIONS",
          "ALEUPDATECoupling", "ALEUPDATE Coupling", Core::Conditions::ALEUPDATECoupling, true,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfaleupdate =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN ALE UPDATE SURF CONDITIONS",
          "ALEUPDATECoupling", "ALEUPDATE Coupling", Core::Conditions::ALEUPDATECoupling, true,
          Core::Conditions::geometry_type_surface));

  for (unsigned i = 0; i < aleupdatecomponents.size(); ++i)
  {
    linealeupdate->AddComponent(aleupdatecomponents[i]);
    surfaleupdate->AddComponent(aleupdatecomponents[i]);
  }

  condlist.push_back(linealeupdate);
  condlist.push_back(surfaleupdate);
}

FOUR_C_NAMESPACE_CLOSE
