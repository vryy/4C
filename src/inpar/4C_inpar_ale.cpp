/*----------------------------------------------------------------------*/
/*! \file


\brief Input parameters for ALE mesh motion

\level 2
*/

/*----------------------------------------------------------------------*/

#include "4C_inpar_ale.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_io_linecomponent.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::ALE::set_valid_parameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& adyn = list->sublist("ALE DYNAMIC", false, "");

  Core::UTILS::double_parameter("TIMESTEP", 0.1, "time step size", &adyn);
  Core::UTILS::int_parameter("NUMSTEP", 41, "max number of time steps", &adyn);
  Core::UTILS::double_parameter("MAXTIME", 4.0, "max simulation time", &adyn);

  setStringToIntegralParameter<int>("ALE_TYPE", "solid", "ale mesh movement algorithm",
      tuple<std::string>("solid", "solid_linear", "laplace_material", "laplace_spatial",
          "springs_material", "springs_spatial"),
      tuple<int>(solid, solid_linear, laplace_material, laplace_spatial, springs_material,
          springs_spatial),
      &adyn);

  Core::UTILS::bool_parameter("ASSESSMESHQUALITY", "no",
      "Evaluate element quality measure according to [Oddy et al. 1988]", &adyn);

  Core::UTILS::bool_parameter("UPDATEMATRIX", "no",
      "Update stiffness matrix in every time step (only for linear/material strategies)", &adyn);

  Core::UTILS::int_parameter("MAXITER", 1, "Maximum number of newton iterations.", &adyn);
  Core::UTILS::double_parameter(
      "TOLRES", 1.0e-06, "Absolute tolerance for length scaled L2 residual norm ", &adyn);
  Core::UTILS::double_parameter(
      "TOLDISP", 1.0e-06, "Absolute tolerance for length scaled L2 increment norm ", &adyn);

  Core::UTILS::int_parameter("NUM_INITSTEP", 0, "", &adyn);
  Core::UTILS::int_parameter("RESTARTEVRY", 1, "write restart data every RESTARTEVRY steps", &adyn);
  Core::UTILS::int_parameter("RESULTSEVRY", 0, "write results every RESULTSTEVRY steps", &adyn);
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
  Core::UTILS::int_parameter("STARTFUNCNO", -1, "Function for Initial displacement", &adyn);

  // linear solver id used for scalar ale problems
  Core::UTILS::int_parameter(
      "LINEAR_SOLVER", -1, "number of linear solver used for ale problems...", &adyn);
}



void Inpar::ALE::set_valid_conditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  /*--------------------------------------------------------------------*/
  // Ale update boundary condition

  Teuchos::RCP<Core::Conditions::ConditionDefinition> linealeupdate =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN ALE UPDATE LINE CONDITIONS",
          "ALEUPDATECoupling", "ALEUPDATE Coupling", Core::Conditions::ALEUPDATECoupling, true,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfaleupdate =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN ALE UPDATE SURF CONDITIONS",
          "ALEUPDATECoupling", "ALEUPDATE Coupling", Core::Conditions::ALEUPDATECoupling, true,
          Core::Conditions::geometry_type_surface));

  for (const auto& cond : {linealeupdate, surfaleupdate})
  {
    add_named_selection_component(cond, "COUPLING", "", "lagrange",
        Teuchos::tuple<std::string>("lagrange", "heightfunction", "sphereHeightFunction",
            "meantangentialvelocity", "meantangentialvelocityscaled"),
        Teuchos::tuple<std::string>("lagrange", "heightfunction", "sphereHeightFunction",
            "meantangentialvelocity", "meantangentialvelocityscaled"),
        true);
    add_named_real(cond, "VAL");
    add_named_int(cond, "NODENORMALFUNCT");

    condlist.emplace_back(cond);
  }
}

FOUR_C_NAMESPACE_CLOSE
