/*----------------------------------------------------------------------*/
/*! \file
\brief fpsi parameters
\level 1
*/
/*----------------------------------------------------------------------*/



#include "4C_inpar_fpsi.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::FPSI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& fpsidyn = list->sublist("FPSI DYNAMIC", false,
      "Fluid Porous Structure Interaction\n"
      "FPSI solver with various coupling methods");

  Teuchos::Tuple<std::string, 1> fpsiname;
  Teuchos::Tuple<int, 1> fpsilabel;

  Teuchos::Array<std::string> name(1);
  Teuchos::Array<int> label(1);
  name[0] = "fpsi_monolithic_plain";
  label[0] = fpsi_monolithic_plain;
  setStringToIntegralParameter<int>("COUPALGO", "fpsi_monolithic_plain",
      "Iteration Scheme over the fields", name, label, &fpsidyn);

  Core::UTILS::BoolParameter("SHAPEDERIVATIVES", "No",
      "Include linearization with respect to mesh movement in Navier Stokes equation.\n"
      "Supported in monolithic FPSI for now.",
      &fpsidyn);

  Core::UTILS::BoolParameter("USESHAPEDERIVATIVES", "No",
      "Add linearization with respect to mesh movement in Navier Stokes equation to stiffness "
      "matrix.\n"
      "Supported in monolithic FPSI for now.",
      &fpsidyn);

  setStringToIntegralParameter<int>("PARTITIONED", "RobinNeumann",
      "Coupling strategies for partitioned FPSI solvers.",
      tuple<std::string>("RobinNeumann", "monolithic", "nocoupling"),
      tuple<int>(RobinNeumann, monolithic, nocoupling), &fpsidyn);

  Core::UTILS::BoolParameter(
      "SECONDORDER", "No", "Second order coupling at the interface.", &fpsidyn);

  // Iterationparameters
  Core::UTILS::StringParameter("RESTOL", "1e-8 1e-8 1e-8 1e-8 1e-8 1e-8",
      "Tolerances for single fields in the residual norm for the Newton iteration. \n"
      "For NORM_RESF != Abs_sys_split only the first value is used for all fields. \n"
      "Order of fields: porofluidvelocity, porofluidpressure, porostructure, fluidvelocity, "
      "fluidpressure, ale",
      &fpsidyn);

  Core::UTILS::StringParameter("INCTOL", "1e-8 1e-8 1e-8 1e-8 1e-8 1e-8",
      "Tolerance in the increment norm for the Newton iteration. \n"
      "For NORM_INC != \\*_split only the first value is used for all fields. \n"
      "Order of fields: porofluidvelocity, porofluidpressure, porostructure, fluidvelocity, "
      "fluidpressure, ale",
      &fpsidyn);

  setStringToIntegralParameter<int>("NORM_INC", "Abs",
      "Type of norm for primary variables convergence check.  \n"
      "Abs: absolute values, Abs_sys_split: absolute values with correction of systemsize for "
      "every field seperate, Rel_sys: relative values with correction of systemsize.",
      tuple<std::string>("Abs", "Abs_sys_split", "Rel_sys"),
      tuple<int>(
          absoluteconvergencenorm, absoluteconvergencenorm_sys_split, relativconvergencenorm_sys),
      &fpsidyn);

  setStringToIntegralParameter<int>("NORM_RESF", "Abs",
      "Type of norm for primary variables convergence check. \n"
      "Abs: absolute values, Abs_sys_split: absolute values with correction of systemsize for "
      "every field seperate, Rel_sys: relative values with correction of systemsize.",
      tuple<std::string>("Abs", "Abs_sys_split", "Rel_sys"),
      tuple<int>(
          absoluteconvergencenorm, absoluteconvergencenorm_sys_split, relativconvergencenorm_sys),
      &fpsidyn);

  setStringToIntegralParameter<int>("NORMCOMBI_RESFINC", "And",
      "binary operator to combine primary variables and residual force values",
      tuple<std::string>("And", "Or"), tuple<int>(bop_and, bop_or), &fpsidyn);

  setStringToIntegralParameter<int>("LineSearch", "No",
      "adapt increment in case of non-monotonic residual convergence or residual oscillations",
      tuple<std::string>("Yes", "No"), tuple<int>(1, 0), &fpsidyn);

  setStringToIntegralParameter<int>("FDCheck", "No",
      "perform FPSIFDCheck() finite difference check", tuple<std::string>("Yes", "No"),
      tuple<int>(1, 0), &fpsidyn);

  // number of linear solver used for poroelasticity
  Core::UTILS::IntParameter(
      "LINEAR_SOLVER", -1, "number of linear solver used for FPSI problems", &fpsidyn);

  Core::UTILS::IntParameter("ITEMAX", 10, "maximum number of iterations over fields", &fpsidyn);
  Core::UTILS::IntParameter("ITEMIN", 1, "minimal number of iterations over fields", &fpsidyn);
  Core::UTILS::IntParameter("NUMSTEP", 200, "Total number of Timesteps", &fpsidyn);
  Core::UTILS::IntParameter("ITEMAX", 100, "Maximum number of iterations over fields", &fpsidyn);
  Core::UTILS::IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &fpsidyn);
  Core::UTILS::IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &fpsidyn);

  Core::UTILS::IntParameter("FDCheck_row", 0, "print row value during fd_check", &fpsidyn);
  Core::UTILS::IntParameter("FDCheck_column", 0, "print column value during fd_check", &fpsidyn);

  Core::UTILS::DoubleParameter("TIMESTEP", 0.1, "Time increment dt", &fpsidyn);
  Core::UTILS::DoubleParameter("MAXTIME", 1000.0, "Total simulation time", &fpsidyn);
  Core::UTILS::DoubleParameter("CONVTOL", 1e-6, "Tolerance for iteration over fields", &fpsidyn);
  Core::UTILS::DoubleParameter("ALPHABJ", 1.0,
      "Beavers-Joseph-Coefficient for Slip-Boundary-Condition at Fluid-Porous-Interface (0.1-4)",
      &fpsidyn);
}



void Inpar::FPSI::SetValidConditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  /*--------------------------------------------------------------------*/
  // FPSI

  std::vector<Teuchos::RCP<Input::LineComponent>> fpsicomponents;

  fpsicomponents.push_back(Teuchos::rcp(new Input::IntComponent("coupling id")));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> linefpsi =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN FPSI COUPLING LINE CONDITIONS",
          "fpsi_coupling", "FPSI Coupling", Core::Conditions::fpsi_coupling, true,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surffpsi =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN FPSI COUPLING SURF CONDITIONS",
          "fpsi_coupling", "FPSI Coupling", Core::Conditions::fpsi_coupling, true,
          Core::Conditions::geometry_type_surface));

  for (unsigned i = 0; i < fpsicomponents.size(); ++i)
  {
    linefpsi->add_component(fpsicomponents[i]);
    surffpsi->add_component(fpsicomponents[i]);
  }

  condlist.push_back(linefpsi);
  condlist.push_back(surffpsi);


  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in fpsi problems
  // necessary where neumann term needs to be integrated in interface
  // elements which share a node with the fpsi interface. Tangential
  // Beaver-Joseph-Condition must not be overwritten by prescribed value!

  Teuchos::RCP<Core::Conditions::ConditionDefinition> neumannintegration_surf =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN SURFACE NEUMANN INTEGRATION",
          "NeumannIntegration", "Neumann Integration", Core::Conditions::NeumannIntegration, true,
          Core::Conditions::geometry_type_surface));

  condlist.push_back(neumannintegration_surf);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in fpsi problems

  Teuchos::RCP<Core::Conditions::ConditionDefinition> neumannintegration_line =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN LINE NEUMANN INTEGRATION",
          "NeumannIntegration", "Neumann Integration", Core::Conditions::NeumannIntegration, true,
          Core::Conditions::geometry_type_line));

  condlist.push_back(neumannintegration_line);
}

FOUR_C_NAMESPACE_CLOSE
