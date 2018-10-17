/*--------------------------------------------------------------------------*/
/*!
\file inpar_ehl.cpp

\brief Elastohydrodynamic lubrication (lubrication structure interaction) parameters

\level 3

\maintainer Mostafa Faraji
*/
/*--------------------------------------------------------------------------*/


#include "drt_validparameters.H"
#include "inpar_ehl.H"
#include "inpar_scatra.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::EHL::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::Array<std::string> yesnotuple =
      tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

  Teuchos::ParameterList& ehldyn = list->sublist("ELASTO HYDRO DYNAMIC", false,
      "Elastohydrodynamic paramters for elastohydrodynamic lubrication (lubrication structure "
      "interaction)");

  // Output type
  DoubleParameter(
      "RESTARTEVRYTIME", 0, "write restart possibility every RESTARTEVRY steps", &ehldyn);
  IntParameter("RESTARTEVRY", 1, "write restart possibility every RESTARTEVRY steps", &ehldyn);
  // Time loop control
  IntParameter("NUMSTEP", 200, "maximum number of Timesteps", &ehldyn);
  DoubleParameter("MAXTIME", 1000.0, "total simulation time", &ehldyn);
  DoubleParameter("TIMESTEP", -1, "time step size dt", &ehldyn);
  BoolParameter(
      "DIFFTIMESTEPSIZE", "No", "use different step size for lubrication and solid", &ehldyn);
  DoubleParameter("RESULTSEVRYTIME", 0, "increment for writing solution", &ehldyn);
  IntParameter("RESULTSEVRY", 1, "increment for writing solution", &ehldyn);
  IntParameter("ITEMAX", 10, "maximum number of iterations over fields", &ehldyn);
  IntParameter("ITEMIN", 1, "minimal number of iterations over fields", &ehldyn);

  // Type of coupling strategy between the two fields
  setStringToIntegralParameter<int>("FIELDCOUPLING", "none",
      "Type of coupling strategy between fields", tuple<std::string>("none", "matching"),
      tuple<int>(coupling_none, coupling_matching), &ehldyn);

  // Coupling strategy for EHL solvers
  setStringToIntegralParameter<int>("COUPALGO", "ehl_Monolithic",
      "Coupling strategies for EHL solvers", tuple<std::string>("ehl_IterStagg", "ehl_Monolithic"),
      tuple<int>(ehl_IterStagg, ehl_Monolithic), &ehldyn);

  /*----------------------------------------------------------------------*/
  /* parameters for monolithic EHL */
  Teuchos::ParameterList& ehldynmono = ehldyn.sublist("MONOLITHIC", false,
      "Monolithic Thermo Structure Interaction\n"
      "Dynamic section for monolithic EHL");

  // convergence tolerance of EHL residual
  DoubleParameter("CONVTOL", 1e-6, "tolerance for convergence check of EHL", &ehldynmono);
  // Iterationparameters
  DoubleParameter("TOLINC", 1.0e-6,
      "tolerance for convergence check of EHL-increment in monolithic EHL", &ehldynmono);

  setStringToIntegralParameter<int>("NORM_RESF", "Abs",
      "type of norm for residual convergence check", tuple<std::string>("Abs", "Rel", "Mix"),
      tuple<int>(convnorm_abs, convnorm_rel, convnorm_mix), &ehldynmono);

  setStringToIntegralParameter<int>("NORM_INC", "Abs",
      "type of norm for convergence check of primary variables in EHL",
      tuple<std::string>("Abs", "Rel", "Mix"), tuple<int>(convnorm_abs, convnorm_rel, convnorm_mix),
      &ehldynmono);


  setStringToIntegralParameter<int>("NORMCOMBI_RESFINC", "Coupl_And_Singl",
      "binary operator to combine primary variables and residual force values",
      tuple<std::string>("And", "Or", "Coupl_Or_Singl", "Coupl_And_Singl", "And_Singl", "Or_Singl"),
      tuple<int>(
          bop_and, bop_or, bop_coupl_or_singl, bop_coupl_and_singl, bop_and_singl, bop_or_singl),
      &ehldynmono);

  setStringToIntegralParameter<int>("ITERNORM", "Rms", "type of norm to be applied to residuals",
      tuple<std::string>("L1", "L1_Scaled", "L2", "Rms", "Inf"),
      tuple<int>(norm_l1, norm_l1_scaled, norm_l2, norm_rms, norm_inf), &ehldynmono);

  setStringToIntegralParameter<int>("NLNSOL", "fullnewton", "Nonlinear solution technique",
      tuple<std::string>("fullnewton"), tuple<int>(soltech_newtonfull), &ehldynmono);

  DoubleParameter("PTCDT", 0.1,
      "pseudo time step for pseudo-transient continuation (PTC) stabilised Newton procedure",
      &ehldynmono);

  // number of linear solver used for monolithic EHL
  IntParameter(
      "LINEAR_SOLVER", -1, "number of linear solver used for monolithic EHL problems", &ehldynmono);

  // convergence criteria adaptivity of monolithic EHL solver
  setStringToIntegralParameter<int>("ADAPTCONV", "No",
      "Switch on adaptive control of linear solver tolerance for nonlinear solution", yesnotuple,
      yesnovalue, &ehldynmono);
  DoubleParameter("ADAPTCONV_BETTER", 0.1,
      "The linear solver shall be this much better than the current nonlinear residual in the "
      "nonlinear convergence limit",
      &ehldynmono);

  setStringToIntegralParameter<int>("INFNORMSCALING", "yes",
      "Scale blocks of matrix with row infnorm?", yesnotuple, yesnovalue, &ehldynmono);

  // merge EHL block matrix to enable use of direct solver in monolithic EHL
  // default: "No", i.e. use block matrix
  BoolParameter("MERGE_EHL_BLOCK_MATRIX", "No", "Merge EHL block matrix", &ehldynmono);

  /*----------------------------------------------------------------------*/
  /* parameters for partitioned EHL */
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& ehldynpart = ehldyn.sublist("PARTITIONED", false,
      "Partitioned Structure Scalar Interaction\n"
      "Control section for partitioned EHL");

  // Solver parameter for relaxation of iterative staggered partitioned EHL
  DoubleParameter("MAXOMEGA", 10.0, "largest omega allowed for Aitken relaxation", &ehldynpart);
  DoubleParameter("MINOMEGA", 0.1, "smallest omega allowed for Aitken relaxation", &ehldynpart);
  DoubleParameter("STARTOMEGA", 1.0, "fixed relaxation parameter", &ehldynpart);

  // convergence tolerance of outer iteration loop
  DoubleParameter("CONVTOL", 1e-6,
      "tolerance for convergence check of outer iteration within partitioned EHL", &ehldynpart);

  // set unprojectable nodes to zero pressure via Dirichlet condition
  setStringToIntegralParameter<int>("UNPROJ_ZERO_DBC", "No",
      "set unprojectable nodes to zero pressure via Dirichlet condition", yesnotuple, yesnovalue,
      &ehldyn);

  // use dry contact model
  setStringToIntegralParameter<int>("DRY_CONTACT_MODEL", "No",
      "set unprojectable nodes to zero pressure via Dirichlet condition", yesnotuple, yesnovalue,
      &ehldyn);
}


void INPAR::EHL::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;
  /*--------------------------------------------------------------------*/
  // ehl mortar coupling

  std::vector<Teuchos::RCP<ConditionComponent>> ehlcomponents;

  ehlcomponents.push_back(Teuchos::rcp(new IntConditionComponent("Interface ID")));
  ehlcomponents.push_back(Teuchos::rcp(
      new StringConditionComponent("Side", "Master", Teuchos::tuple<std::string>("Master", "Slave"),
          Teuchos::tuple<std::string>("Master", "Slave"))));
  ehlcomponents.push_back(Teuchos::rcp(new StringConditionComponent("Initialization", "Active",
      Teuchos::tuple<std::string>("Inactive", "Active"),
      Teuchos::tuple<std::string>("Inactive", "Active"), true)));

  ehlcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("FrCoeffOrBound", true)));
  ehlcomponents.push_back(Teuchos::rcp(new RealConditionComponent("FrCoeffOrBound")));

  Teuchos::RCP<ConditionDefinition> lineehl = Teuchos::rcp(
      new ConditionDefinition("DESIGN LINE EHL MORTAR COUPLING CONDITIONS 2D", "EHLCoupling",
          "Line EHL Coupling", DRT::Condition::EHLCoupling, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfehl = Teuchos::rcp(
      new ConditionDefinition("DESIGN SURF EHL MORTAR COUPLING CONDITIONS 3D", "EHLCoupling",
          "Surface EHL Coupling", DRT::Condition::EHLCoupling, true, DRT::Condition::Surface));


  for (unsigned i = 0; i < ehlcomponents.size(); ++i)
  {
    lineehl->AddComponent(ehlcomponents[i]);
    surfehl->AddComponent(ehlcomponents[i]);
  }

  condlist.push_back(lineehl);
  condlist.push_back(surfehl);
}
