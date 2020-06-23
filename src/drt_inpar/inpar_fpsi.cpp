/*----------------------------------------------------------------------*/
/*! \file
\brief fpsi parameters
\level 1
\maintainer  Johannes Kremheller
*/
/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_fpsi.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::FPSI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::Array<std::string> yesnotuple =
      tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

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

  setStringToIntegralParameter<int>("SHAPEDERIVATIVES", "No",
      "Include linearization with respect to mesh movement in Navier Stokes equation.\n"
      "Supported in monolithic FPSI for now.",
      yesnotuple, yesnovalue, &fpsidyn);

  setStringToIntegralParameter<int>("USESHAPEDERIVATIVES", "No",
      "Add linearization with respect to mesh movement in Navier Stokes equation to stiffness "
      "matrix.\n"
      "Supported in monolithic FPSI for now.",
      yesnotuple, yesnovalue, &fpsidyn);

  setStringToIntegralParameter<int>("PARTITIONED", "RobinNeumann",
      "Coupling strategies for partitioned FPSI solvers.",
      tuple<std::string>("RobinNeumann", "monolithic", "nocoupling"),
      tuple<int>(RobinNeumann, monolithic, nocoupling), &fpsidyn);

  setStringToIntegralParameter<int>("SECONDORDER", "No", "Second order coupling at the interface.",
      yesnotuple, yesnovalue, &fpsidyn);

  // Iterationparameters
  StringParameter("RESTOL", "1e-8 1e-8 1e-8 1e-8 1e-8 1e-8",
      "tolerances for single fields in the residual norm for the Newton iteration \n"
      "for NORM_RESF != *_split only the first value is used for all fields \n"
      "order of fields: porofluidvelocity, porofluidpressure, porostructure, fluidvelocity, "
      "fluidpressure, ale",
      &fpsidyn);

  StringParameter("INCTOL", "1e-8 1e-8 1e-8 1e-8 1e-8 1e-8",
      "tolerance in the increment norm for the Newton iteration \n"
      "for NORM_INC != *_split only the first value is used for all fields \n"
      "order of fields: porofluidvelocity, porofluidpressure, porostructure, fluidvelocity, "
      "fluidpressure, ale",
      &fpsidyn);

  setStringToIntegralParameter<int>("NORM_INC", "Abs",
      "type of norm for primary variables convergence check \n"
      "Abs: absolute values, Abs_sys_split: absolute values with correction of systemsize for "
      "every field seperate, Rel_sys: relative values with correction of systemsize",
      tuple<std::string>("Abs", "Abs_sys_split", "Rel_sys"),
      tuple<int>(
          absoluteconvergencenorm, absoluteconvergencenorm_sys_split, relativconvergencenorm_sys),
      &fpsidyn);

  setStringToIntegralParameter<int>("NORM_RESF", "Abs",
      "type of norm for primary variables convergence check \n"
      "Abs: absolute values, Abs_sys_split: absolute values with correction of systemsize for "
      "every field seperate, Rel_sys: relative values with correction of systemsize",
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
  IntParameter("LINEAR_SOLVER", -1, "number of linear solver used for FPSI problems", &fpsidyn);

  IntParameter("ITEMAX", 10, "maximum number of iterations over fields", &fpsidyn);
  IntParameter("ITEMIN", 1, "minimal number of iterations over fields", &fpsidyn);
  IntParameter("NUMSTEP", 200, "Total number of Timesteps", &fpsidyn);
  IntParameter("ITEMAX", 100, "Maximum number of iterations over fields", &fpsidyn);
  IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &fpsidyn);
  IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &fpsidyn);

  IntParameter("FDCheck_row", 0, "print row value during FDCheck", &fpsidyn);
  IntParameter("FDCheck_column", 0, "print column value during FDCheck", &fpsidyn);

  DoubleParameter("TIMESTEP", 0.1, "Time increment dt", &fpsidyn);
  DoubleParameter("MAXTIME", 1000.0, "Total simulation time", &fpsidyn);
  DoubleParameter("CONVTOL", 1e-6, "Tolerance for iteration over fields", &fpsidyn);
  DoubleParameter("ALPHABJ", 1.0,
      "Beavers-Joseph-Coefficient for Slip-Boundary-Condition at Fluid-Porous-Interface (0.1-4)",
      &fpsidyn);
}



void INPAR::FPSI::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // FPSI

  std::vector<Teuchos::RCP<ConditionComponent>> fpsicomponents;

  fpsicomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

  Teuchos::RCP<ConditionDefinition> linefpsi =
      Teuchos::rcp(new ConditionDefinition("DESIGN FPSI COUPLING LINE CONDITIONS", "FPSICoupling",
          "FPSI Coupling", DRT::Condition::FPSICoupling, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surffpsi =
      Teuchos::rcp(new ConditionDefinition("DESIGN FPSI COUPLING SURF CONDITIONS", "FPSICoupling",
          "FPSI Coupling", DRT::Condition::FPSICoupling, true, DRT::Condition::Surface));

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
      DRT::Condition::NeumannIntegration, true, DRT::Condition::Surface));

  condlist.push_back(neumannintegration_surf);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in fpsi problems

  Teuchos::RCP<ConditionDefinition> neumannintegration_line =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE NEUMANN INTEGRATION", "NeumannIntegration",
          "Neumann Integration", DRT::Condition::NeumannIntegration, true, DRT::Condition::Line));

  condlist.push_back(neumannintegration_line);
}
