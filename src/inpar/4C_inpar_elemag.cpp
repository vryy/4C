/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for electromagnetic simulations

\level 3


*/
/*----------------------------------------------------------------------*/

#include "4C_inpar_elemag.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_linalg_equilibrate.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

void Inpar::EleMag::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& electromagneticdyn = list->sublist(
      "ELECTROMAGNETIC DYNAMIC", false, "control parameters for electromagnetic problems\n");

  // general settings for time-integration scheme
  Core::UTILS::DoubleParameter("TIMESTEP", 0.01, "Time-step length dt", &electromagneticdyn);
  Core::UTILS::DoubleParameter("TAU", 1, "Stabilization parameter", &electromagneticdyn);
  Core::UTILS::IntParameter("NUMSTEP", 100, "Number of time steps", &electromagneticdyn);
  Core::UTILS::DoubleParameter("MAXTIME", 1.0, "Total simulation time", &electromagneticdyn);

  // additional parameters
  Core::UTILS::IntParameter(
      "RESULTSEVRY", 1, "Increment for writing solution", &electromagneticdyn);
  Core::UTILS::IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &electromagneticdyn);
  Core::UTILS::IntParameter("LINEAR_SOLVER", -1,
      "Number of linear solver used for electromagnetic problem", &electromagneticdyn);
  Core::UTILS::IntParameter("STARTFUNCNO", -1, "Function for initial field", &electromagneticdyn);
  Core::UTILS::IntParameter(
      "SOURCEFUNCNO", -1, "Function for source term in volume", &electromagneticdyn);
  // Core::UTILS::BoolParameter("DOUBLEORFLOAT","Yes","Yes, if evaluation with double, no if with
  // float",&electromagneticdyn); Core::UTILS::BoolParameter("ALLELESEQUAL","No","Yes, if all
  // elements have same shape and material",&electromagneticdyn);

  {
    // time integration

    Teuchos::Tuple<std::string, 8> name;
    Teuchos::Tuple<int, 8> label;
    name[0] = "One_Step_Theta";
    label[0] = elemag_ost;
    name[1] = "BDF1";
    label[1] = elemag_bdf1;
    name[2] = "BDF2";
    label[2] = elemag_bdf2;
    name[3] = "BDF4";
    label[3] = elemag_bdf4;
    name[4] = "GenAlpha";
    label[4] = elemag_genAlpha;
    name[5] = "Explicit_Euler";
    label[5] = elemag_explicit_euler;
    name[6] = "Runge_Kutta";
    label[6] = elemag_rk;
    name[7] = "Crank_Nicolson";
    label[7] = elemag_cn;

    setStringToIntegralParameter<int>("TIMEINT", "One_Step_Theta",
        "Type of time integration scheme", name, label, &electromagneticdyn);
  }

  {
    // a standard Teuchos::tuple can have at maximum 10 entries! We have to circumvent this here.
    Teuchos::Tuple<std::string, 4> name;
    Teuchos::Tuple<int, 4> label;
    name[0] = "zero_field";
    label[0] = initfield_zero_field;
    name[1] = "field_by_function";
    label[1] = initfield_field_by_function;
    name[2] = "field_by_steady_state";
    label[2] = initfield_scatra;
    name[3] = "field_by_steady_state_hdg";
    label[3] = initfield_scatra_hdg;

    setStringToIntegralParameter<int>("INITIALFIELD", "zero_field", "Initial field for ele problem",
        name, label, &electromagneticdyn);

    // Error calculation
    Core::UTILS::BoolParameter(
        "CALCERR", "No", "Calc the error wrt ERRORFUNCNO?", &electromagneticdyn);

    // Post process solution?
    Core::UTILS::BoolParameter(
        "POSTPROCESS", "No", "Postprocess solution? (very slow)", &electromagneticdyn);
  }

  Core::UTILS::IntParameter(
      "ERRORFUNCNO", -1, "Function for error calculation", &electromagneticdyn);

  // flag for equilibration of global system of equations
  setStringToIntegralParameter<Core::LinAlg::EquilibrationMethod>("EQUILIBRATION", "none",
      "flag for equilibration of global system of equations",
      tuple<std::string>("none", "rows_full", "rows_maindiag", "columns_full", "columns_maindiag",
          "rowsandcolumns_full", "rowsandcolumns_maindiag"),
      tuple<Core::LinAlg::EquilibrationMethod>(Core::LinAlg::EquilibrationMethod::none,
          Core::LinAlg::EquilibrationMethod::rows_full,
          Core::LinAlg::EquilibrationMethod::rows_maindiag,
          Core::LinAlg::EquilibrationMethod::columns_full,
          Core::LinAlg::EquilibrationMethod::columns_maindiag,
          Core::LinAlg::EquilibrationMethod::rowsandcolumns_full,
          Core::LinAlg::EquilibrationMethod::rowsandcolumns_maindiag),
      &electromagneticdyn);

  // PML
  // Core::UTILS::StringParameter("PML_DEFINITION_FILE","none.txt","Filename of file containing the
  // pml definition",&electromagneticdyn);
}

/// set specific electromagnetic conditions
void Inpar::EleMag::SetValidConditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  std::vector<Teuchos::RCP<Input::LineComponent>> abcbundcomponents;

  abcbundcomponents.emplace_back(Teuchos::rcp(new Input::SeparatorComponent("NUMDOF")));
  abcbundcomponents.emplace_back(Teuchos::rcp(new Input::IntComponent("numdof")));
  abcbundcomponents.emplace_back(Teuchos::rcp(new Input::SeparatorComponent("ONOFF")));
  abcbundcomponents.emplace_back(
      Teuchos::rcp(new Input::IntVectorComponent("onoff", Input::LengthFromInt("numdof"))));
  abcbundcomponents.emplace_back(Teuchos::rcp(new Input::SeparatorComponent("FUNCT")));
  abcbundcomponents.emplace_back(Teuchos::rcp(new Input::IntVectorComponent(
      "funct", Input::LengthFromInt("numdof"), {0, false, true, false})));
  abcbundcomponents.emplace_back(Teuchos::rcp(new Input::SeparatorComponent("VAL")));
  abcbundcomponents.emplace_back(
      Teuchos::rcp(new Input::RealVectorComponent("val", Input::LengthFromInt("numdof"))));

  //*--------------------------------------------------------------------* /
  // absorbing boundary condition for electromagnetic problems
  // line
  Teuchos::RCP<Core::Conditions::ConditionDefinition> silvermueller_line = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN LINE SILVER-MUELLER CONDITIONS",
          "Silver-Mueller", "Absorbing-emitting line for electromagnetics",
          Core::Conditions::SilverMueller, true, Core::Conditions::geometry_type_line));

  // surface
  Teuchos::RCP<Core::Conditions::ConditionDefinition> silvermueller_surface = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN SURF SILVER-MUELLER CONDITIONS",
          "Silver-Mueller", "Absorbing-emitting surface for electromagnetics",
          Core::Conditions::SilverMueller, true, Core::Conditions::geometry_type_surface));

  for (unsigned i = 0; i < abcbundcomponents.size(); ++i)
  {
    silvermueller_line->add_component(abcbundcomponents[i]);
    silvermueller_surface->add_component(abcbundcomponents[i]);
  }

  condlist.push_back(silvermueller_line);
  condlist.push_back(silvermueller_surface);
}

FOUR_C_NAMESPACE_CLOSE
