/*--------------------------------------------------------------------------*/
/*! \file

\brief Elastohydrodynamic lubrication (lubrication structure interaction) parameters

\level 3

*/
/*--------------------------------------------------------------------------*/


#include "4C_inpar_ehl.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::EHL::set_valid_parameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& ehldyn = list->sublist("ELASTO HYDRO DYNAMIC", false,
      "Elastohydrodynamic paramters for elastohydrodynamic lubrication (lubrication structure "
      "interaction)");

  // Output type
  Core::UTILS::double_parameter(
      "RESTARTEVRYTIME", 0, "write restart possibility every RESTARTEVRY steps", &ehldyn);
  Core::UTILS::int_parameter(
      "RESTARTEVRY", 1, "write restart possibility every RESTARTEVRY steps", &ehldyn);
  // Time loop control
  Core::UTILS::int_parameter("NUMSTEP", 200, "maximum number of Timesteps", &ehldyn);
  Core::UTILS::double_parameter("MAXTIME", 1000.0, "total simulation time", &ehldyn);
  Core::UTILS::double_parameter("TIMESTEP", -1, "time step size dt", &ehldyn);
  Core::UTILS::bool_parameter(
      "DIFFTIMESTEPSIZE", "No", "use different step size for lubrication and solid", &ehldyn);
  Core::UTILS::double_parameter("RESULTSEVRYTIME", 0, "increment for writing solution", &ehldyn);
  Core::UTILS::int_parameter("RESULTSEVRY", 1, "increment for writing solution", &ehldyn);
  Core::UTILS::int_parameter("ITEMAX", 10, "maximum number of iterations over fields", &ehldyn);
  Core::UTILS::int_parameter("ITEMIN", 1, "minimal number of iterations over fields", &ehldyn);

  // Type of coupling strategy between the two fields
  setStringToIntegralParameter<FieldCoupling>("FIELDCOUPLING", "none",
      "Type of coupling strategy between fields", tuple<std::string>("none", "matching"),
      tuple<FieldCoupling>(coupling_none, coupling_matching), &ehldyn);

  // Coupling strategy for EHL solvers
  setStringToIntegralParameter<SolutionSchemeOverFields>("COUPALGO", "ehl_Monolithic",
      "Coupling strategies for EHL solvers", tuple<std::string>("ehl_IterStagg", "ehl_Monolithic"),
      tuple<SolutionSchemeOverFields>(ehl_IterStagg, ehl_Monolithic), &ehldyn);

  /*----------------------------------------------------------------------*/
  /* parameters for monolithic EHL */
  Teuchos::ParameterList& ehldynmono = ehldyn.sublist("MONOLITHIC", false,
      "Monolithic Thermo Structure Interaction\n"
      "Dynamic section for monolithic EHL");

  // convergence tolerance of EHL residual
  Core::UTILS::double_parameter(
      "CONVTOL", 1e-6, "tolerance for convergence check of EHL", &ehldynmono);
  // Iterationparameters
  Core::UTILS::double_parameter("TOLINC", 1.0e-6,
      "tolerance for convergence check of EHL-increment in monolithic EHL", &ehldynmono);

  setStringToIntegralParameter<ConvNorm>("NORM_RESF", "Abs",
      "type of norm for residual convergence check", tuple<std::string>("Abs", "Rel", "Mix"),
      tuple<ConvNorm>(convnorm_abs, convnorm_rel, convnorm_mix), &ehldynmono);

  setStringToIntegralParameter<ConvNorm>("NORM_INC", "Abs",
      "type of norm for convergence check of primary variables in EHL",
      tuple<std::string>("Abs", "Rel", "Mix"),
      tuple<ConvNorm>(convnorm_abs, convnorm_rel, convnorm_mix), &ehldynmono);


  setStringToIntegralParameter<BinaryOp>("NORMCOMBI_RESFINC", "Coupl_And_Singl",
      "binary operator to combine primary variables and residual force values",
      tuple<std::string>("And", "Or", "Coupl_Or_Singl", "Coupl_And_Singl", "And_Singl", "Or_Singl"),
      tuple<BinaryOp>(
          bop_and, bop_or, bop_coupl_or_singl, bop_coupl_and_singl, bop_and_singl, bop_or_singl),
      &ehldynmono);

  setStringToIntegralParameter<VectorNorm>("ITERNORM", "Rms",
      "type of norm to be applied to residuals",
      tuple<std::string>("L1", "L1_Scaled", "L2", "Rms", "Inf"),
      tuple<VectorNorm>(norm_l1, norm_l1_scaled, norm_l2, norm_rms, norm_inf), &ehldynmono);

  Core::UTILS::double_parameter("PTCDT", 0.1,
      "pseudo time step for pseudo-transient continuation (PTC) stabilised Newton procedure",
      &ehldynmono);

  // number of linear solver used for monolithic EHL
  Core::UTILS::int_parameter(
      "LINEAR_SOLVER", -1, "number of linear solver used for monolithic EHL problems", &ehldynmono);

  // convergence criteria adaptivity of monolithic EHL solver
  Core::UTILS::bool_parameter("ADAPTCONV", "No",
      "Switch on adaptive control of linear solver tolerance for nonlinear solution", &ehldynmono);
  Core::UTILS::double_parameter("ADAPTCONV_BETTER", 0.1,
      "The linear solver shall be this much better than the current nonlinear residual in the "
      "nonlinear convergence limit",
      &ehldynmono);

  Core::UTILS::bool_parameter(
      "INFNORMSCALING", "yes", "Scale blocks of matrix with row infnorm?", &ehldynmono);

  /*----------------------------------------------------------------------*/
  /* parameters for partitioned EHL */
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& ehldynpart = ehldyn.sublist("PARTITIONED", false,
      "Partitioned Structure Scalar Interaction\n"
      "Control section for partitioned EHL");

  // Solver parameter for relaxation of iterative staggered partitioned EHL
  Core::UTILS::double_parameter(
      "MAXOMEGA", 10.0, "largest omega allowed for Aitken relaxation", &ehldynpart);
  Core::UTILS::double_parameter(
      "MINOMEGA", 0.1, "smallest omega allowed for Aitken relaxation", &ehldynpart);
  Core::UTILS::double_parameter("STARTOMEGA", 1.0, "fixed relaxation parameter", &ehldynpart);

  // convergence tolerance of outer iteration loop
  Core::UTILS::double_parameter("CONVTOL", 1e-6,
      "tolerance for convergence check of outer iteration within partitioned EHL", &ehldynpart);

  // set unprojectable nodes to zero pressure via Dirichlet condition
  Core::UTILS::bool_parameter("UNPROJ_ZERO_DBC", "No",
      "set unprojectable nodes to zero pressure via Dirichlet condition", &ehldyn);

  // use dry contact model
  Core::UTILS::bool_parameter("DRY_CONTACT_MODEL", "No",
      "set unprojectable nodes to zero pressure via Dirichlet condition", &ehldyn);
}


void Inpar::EHL::set_valid_conditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;
  /*--------------------------------------------------------------------*/
  // ehl mortar coupling

  Teuchos::RCP<Core::Conditions::ConditionDefinition> lineehl =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE EHL MORTAR COUPLING CONDITIONS 2D", "EHLCoupling", "Line EHL Coupling",
          Core::Conditions::EHLCoupling, true, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfehl =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN SURF EHL MORTAR COUPLING CONDITIONS 3D", "EHLCoupling", "Surface EHL Coupling",
          Core::Conditions::EHLCoupling, true, Core::Conditions::geometry_type_surface));

  for (const auto& cond : {lineehl, surfehl})
  {
    cond->add_component(Teuchos::rcp(new Input::IntComponent("Interface ID")));
    cond->add_component(Teuchos::rcp(new Input::SelectionComponent("Side", "Master",
        Teuchos::tuple<std::string>("Master", "Slave"),
        Teuchos::tuple<std::string>("Master", "Slave"))));
    cond->add_component(Teuchos::rcp(new Input::SelectionComponent("Initialization", "Active",
        Teuchos::tuple<std::string>("Inactive", "Active"),
        Teuchos::tuple<std::string>("Inactive", "Active"), true)));
    add_named_real(cond, "FrCoeffOrBound", "", 0.0, true);

    condlist.push_back(cond);
  }
}

FOUR_C_NAMESPACE_CLOSE
