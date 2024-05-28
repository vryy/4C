/*----------------------------------------------------------------------*/
/*! \file
\brief fpsi parameters
\level 1
*/
/*----------------------------------------------------------------------*/



#include "4C_inpar_fpsi.hpp"

#include "4C_io_condition_definition.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void INPAR::FPSI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace INPUT;
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

  CORE::UTILS::BoolParameter("SHAPEDERIVATIVES", "No",
      "Include linearization with respect to mesh movement in Navier Stokes equation.\n"
      "Supported in monolithic FPSI for now.",
      &fpsidyn);

  CORE::UTILS::BoolParameter("USESHAPEDERIVATIVES", "No",
      "Add linearization with respect to mesh movement in Navier Stokes equation to stiffness "
      "matrix.\n"
      "Supported in monolithic FPSI for now.",
      &fpsidyn);

  setStringToIntegralParameter<int>("PARTITIONED", "RobinNeumann",
      "Coupling strategies for partitioned FPSI solvers.",
      tuple<std::string>("RobinNeumann", "monolithic", "nocoupling"),
      tuple<int>(RobinNeumann, monolithic, nocoupling), &fpsidyn);

  CORE::UTILS::BoolParameter(
      "SECONDORDER", "No", "Second order coupling at the interface.", &fpsidyn);

  // Iterationparameters
  CORE::UTILS::StringParameter("RESTOL", "1e-8 1e-8 1e-8 1e-8 1e-8 1e-8",
      "Tolerances for single fields in the residual norm for the Newton iteration. \n"
      "For NORM_RESF != Abs_sys_split only the first value is used for all fields. \n"
      "Order of fields: porofluidvelocity, porofluidpressure, porostructure, fluidvelocity, "
      "fluidpressure, ale",
      &fpsidyn);

  CORE::UTILS::StringParameter("INCTOL", "1e-8 1e-8 1e-8 1e-8 1e-8 1e-8",
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
  CORE::UTILS::IntParameter(
      "LINEAR_SOLVER", -1, "number of linear solver used for FPSI problems", &fpsidyn);

  CORE::UTILS::IntParameter("ITEMAX", 10, "maximum number of iterations over fields", &fpsidyn);
  CORE::UTILS::IntParameter("ITEMIN", 1, "minimal number of iterations over fields", &fpsidyn);
  CORE::UTILS::IntParameter("NUMSTEP", 200, "Total number of Timesteps", &fpsidyn);
  CORE::UTILS::IntParameter("ITEMAX", 100, "Maximum number of iterations over fields", &fpsidyn);
  CORE::UTILS::IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &fpsidyn);
  CORE::UTILS::IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &fpsidyn);

  CORE::UTILS::IntParameter("FDCheck_row", 0, "print row value during fd_check", &fpsidyn);
  CORE::UTILS::IntParameter("FDCheck_column", 0, "print column value during fd_check", &fpsidyn);

  CORE::UTILS::DoubleParameter("TIMESTEP", 0.1, "Time increment dt", &fpsidyn);
  CORE::UTILS::DoubleParameter("MAXTIME", 1000.0, "Total simulation time", &fpsidyn);
  CORE::UTILS::DoubleParameter("CONVTOL", 1e-6, "Tolerance for iteration over fields", &fpsidyn);
  CORE::UTILS::DoubleParameter("ALPHABJ", 1.0,
      "Beavers-Joseph-Coefficient for Slip-Boundary-Condition at Fluid-Porous-Interface (0.1-4)",
      &fpsidyn);
}



void INPAR::FPSI::SetValidConditions(
    std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist)
{
  using namespace INPUT;

  /*--------------------------------------------------------------------*/
  // FPSI

  std::vector<Teuchos::RCP<INPUT::LineComponent>> fpsicomponents;

  fpsicomponents.push_back(Teuchos::rcp(new INPUT::IntComponent("coupling id")));

  Teuchos::RCP<ConditionDefinition> linefpsi = Teuchos::rcp(new ConditionDefinition(
      "DESIGN FPSI COUPLING LINE CONDITIONS", "fpsi_coupling", "FPSI Coupling",
      CORE::Conditions::fpsi_coupling, true, CORE::Conditions::geometry_type_line));
  Teuchos::RCP<ConditionDefinition> surffpsi = Teuchos::rcp(new ConditionDefinition(
      "DESIGN FPSI COUPLING SURF CONDITIONS", "fpsi_coupling", "FPSI Coupling",
      CORE::Conditions::fpsi_coupling, true, CORE::Conditions::geometry_type_surface));

  for (unsigned i = 0; i < fpsicomponents.size(); ++i)
  {
    linefpsi->AddComponent(fpsicomponents[i]);
    surffpsi->AddComponent(fpsicomponents[i]);
  }

  condlist.push_back(linefpsi);
  condlist.push_back(surffpsi);


  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in fpsi problems
  // necessary where neumann term needs to be integrated in interface
  // elements which share a node with the fpsi interface. Tangential
  // Beaver-Joseph-Condition must not be overwritten by prescribed value!

  Teuchos::RCP<ConditionDefinition> neumannintegration_surf = Teuchos::rcp(new ConditionDefinition(
      "DESIGN SURFACE NEUMANN INTEGRATION", "NeumannIntegration", "Neumann Integration",
      CORE::Conditions::NeumannIntegration, true, CORE::Conditions::geometry_type_surface));

  condlist.push_back(neumannintegration_surf);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in fpsi problems

  Teuchos::RCP<ConditionDefinition> neumannintegration_line = Teuchos::rcp(new ConditionDefinition(
      "DESIGN LINE NEUMANN INTEGRATION", "NeumannIntegration", "Neumann Integration",
      CORE::Conditions::NeumannIntegration, true, CORE::Conditions::geometry_type_line));

  condlist.push_back(neumannintegration_line);
}

FOUR_C_NAMESPACE_CLOSE
