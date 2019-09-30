/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for level set

\level 2

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_levelset.H"
#include "inpar_scatra.H"
#include "../drt_lib/drt_conditiondefinition.H"


void INPAR::LEVELSET::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& levelsetcontrol =
      list->sublist("LEVEL-SET CONTROL", false, "control parameters for level-set problems\n");

  IntParameter("NUMSTEP", 24, "Total number of time steps", &levelsetcontrol);
  DoubleParameter("TIMESTEP", 0.1, "Time increment dt", &levelsetcontrol);
  DoubleParameter("MAXTIME", 1000.0, "Total simulation time", &levelsetcontrol);
  IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &levelsetcontrol);
  IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &levelsetcontrol);

  setStringToIntegralParameter<int>("CALCERROR", "No",
      "compute error compared to analytical solution", tuple<std::string>("No", "InitialField"),
      tuple<int>(INPAR::SCATRA::calcerror_no_ls, INPAR::SCATRA::calcerror_initial_field),
      &levelsetcontrol);

  BoolParameter("EXTRACT_INTERFACE_VEL", "No",
      "replace computed velocity at nodes of given distance of interface by approximated interface "
      "velocity",
      &levelsetcontrol);
  IntParameter("NUM_CONVEL_LAYERS", -1,
      "number of layers around the interface which keep their computed convective velocity",
      &levelsetcontrol);


  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& ls_reinit = levelsetcontrol.sublist("REINITIALIZATION", false, "");

  setStringToIntegralParameter<int>("REINITIALIZATION", "None",
      "Type of reinitialization strategy for level set function",
      tuple<std::string>("None", "Signed_Distance_Function", "Sussman", "EllipticEq"),
      tuple<int>(INPAR::SCATRA::reinitaction_none,
          INPAR::SCATRA::reinitaction_signeddistancefunction, INPAR::SCATRA::reinitaction_sussman,
          INPAR::SCATRA::reinitaction_ellipticeq),
      &ls_reinit);

  BoolParameter("REINIT_INITIAL", "No",
      "Has level set field to be reinitialized before first time step?", &ls_reinit);
  IntParameter("REINITINTERVAL", 1, "reinitialization interval", &ls_reinit);

  // parameters for signed distance reinitialization
  BoolParameter("REINITBAND", "No",
      "reinitialization only within a band around the interface, or entire domain?", &ls_reinit);
  DoubleParameter("REINITBANDWIDTH", 1.0,
      "level-set value defining band width for reinitialization", &ls_reinit);

  // parameters for reinitialization equation
  IntParameter("NUMSTEPSREINIT", 1, "(maximal) number of pseudo-time steps", &ls_reinit);
  // this parameter selects the tau definition applied
  setStringToIntegralParameter<int>("LINEARIZATIONREINIT", "fixed_point",
      "linearization of reinitialization equation", tuple<std::string>("newton", "fixed_point"),
      tuple<int>(INPAR::SCATRA::newton, INPAR::SCATRA::fixed_point), &ls_reinit);
  DoubleParameter("TIMESTEPREINIT", 1.0,
      "pseudo-time step length (usually a * characteristic element length of discretization with "
      "a>0)",
      &ls_reinit);
  DoubleParameter(
      "THETAREINIT", 1.0, "theta for time discretization of reinitialization equation", &ls_reinit);
  setStringToIntegralParameter<int>("STABTYPEREINIT", "SUPG", "type of stabilization (if any)",
      tuple<std::string>("no_stabilization", "SUPG", "GLS", "USFEM"),
      tuple<std::string>(
          "Do not use any stabilization -> only reasonable for low-Peclet-number flows", "Use SUPG",
          "Use GLS", "Use USFEM"),
      tuple<int>(INPAR::SCATRA::stabtype_no_stabilization, INPAR::SCATRA::stabtype_SUPG,
          INPAR::SCATRA::stabtype_GLS, INPAR::SCATRA::stabtype_USFEM),
      &ls_reinit);
  // this parameter selects the tau definition applied
  setStringToIntegralParameter<int>("DEFINITION_TAU_REINIT", "Taylor_Hughes_Zarins",
      "Definition of tau",
      tuple<std::string>("Taylor_Hughes_Zarins", "Taylor_Hughes_Zarins_wo_dt", "Franca_Valentin",
          "Franca_Valentin_wo_dt", "Shakib_Hughes_Codina", "Shakib_Hughes_Codina_wo_dt", "Codina",
          "Codina_wo_dt", "Exact_1D", "Zero"),
      tuple<int>(INPAR::SCATRA::tau_taylor_hughes_zarins,
          INPAR::SCATRA::tau_taylor_hughes_zarins_wo_dt, INPAR::SCATRA::tau_franca_valentin,
          INPAR::SCATRA::tau_franca_valentin_wo_dt, INPAR::SCATRA::tau_shakib_hughes_codina,
          INPAR::SCATRA::tau_shakib_hughes_codina_wo_dt, INPAR::SCATRA::tau_codina,
          INPAR::SCATRA::tau_codina_wo_dt, INPAR::SCATRA::tau_exact_1d, INPAR::SCATRA::tau_zero),
      &ls_reinit);
  // this parameter governs whether all-scale subgrid diffusivity is included
  setStringToIntegralParameter<int>("ARTDIFFREINIT", "no",
      "potential incorporation of all-scale subgrid diffusivity (a.k.a. discontinuity-capturing) "
      "term",
      tuple<std::string>("no", "isotropic", "crosswind"),
      tuple<std::string>("no artificial diffusion", "homogeneous artificial diffusion",
          "artificial diffusion in crosswind directions only"),
      tuple<int>(INPAR::SCATRA::artdiff_none, INPAR::SCATRA::artdiff_isotropic,
          INPAR::SCATRA::artdiff_crosswind),
      &ls_reinit);
  // this parameter selects the all-scale subgrid-diffusivity definition applied
  setStringToIntegralParameter<int>("DEFINITION_ARTDIFFREINIT", "artificial_linear",
      "Definition of (all-scale) subgrid diffusivity",
      tuple<std::string>("artificial_linear", "artificial_linear_reinit",
          "Hughes_etal_86_nonlinear", "Tezduyar_Park_86_nonlinear",
          "Tezduyar_Park_86_nonlinear_wo_phizero", "doCarmo_Galeao_91_nonlinear",
          "Almeida_Silva_97_nonlinear", "YZbeta_nonlinear", "Codina_nonlinear"),
      tuple<std::string>("simple linear artificial diffusion",
          "simple linear artificial diffusion const*h",
          "nonlinear isotropic according to Hughes et al. (1986)",
          "nonlinear isotropic according to Tezduyar and Park (1986)",
          "nonlinear isotropic according to Tezduyar and Park (1986) without user parameter "
          "phi_zero",
          "nonlinear isotropic according to doCarmo and Galeao (1991)",
          "nonlinear isotropic according to Almeida and Silva (1997)", "nonlinear YZ beta model",
          "nonlinear isotropic according to Codina"),
      tuple<int>(INPAR::SCATRA::assgd_artificial, INPAR::SCATRA::assgd_lin_reinit,
          INPAR::SCATRA::assgd_hughes, INPAR::SCATRA::assgd_tezduyar,
          INPAR::SCATRA::assgd_tezduyar_wo_phizero, INPAR::SCATRA::assgd_docarmo,
          INPAR::SCATRA::assgd_almeida, INPAR::SCATRA::assgd_yzbeta, INPAR::SCATRA::assgd_codina),
      &ls_reinit);

  setStringToIntegralParameter<int>("SMOOTHED_SIGN_TYPE", "SussmanSmerekaOsher1994",
      "sign function for reinitialization equation",
      tuple<std::string>("NonSmoothed",
          "SussmanFatemi1999",  // smeared-out Heaviside function
          "SussmanSmerekaOsher1994", "PengEtAl1999"),
      tuple<int>(INPAR::SCATRA::signtype_nonsmoothed, INPAR::SCATRA::signtype_SussmanFatemi1999,
          INPAR::SCATRA::signtype_SussmanSmerekaOsher1994, INPAR::SCATRA::signtype_PengEtAl1999),
      &ls_reinit);
  setStringToIntegralParameter<int>("CHARELELENGTHREINIT", "root_of_volume",
      "characteristic element length for sign function",
      tuple<std::string>("root_of_volume", "streamlength"),
      tuple<int>(INPAR::SCATRA::root_of_volume_reinit, INPAR::SCATRA::streamlength_reinit),
      &ls_reinit);
  DoubleParameter("INTERFACE_THICKNESS", 1.0,
      "factor for interface thickness (multiplied by element length)", &ls_reinit);
  setStringToIntegralParameter<int>("VELREINIT", "integration_point_based",
      "evaluate velocity at integration point or compute node-based velocity",
      tuple<std::string>("integration_point_based", "node_based"),
      tuple<int>(
          INPAR::SCATRA::vel_reinit_integration_point_based, INPAR::SCATRA::vel_reinit_node_based),
      &ls_reinit);
  setStringToIntegralParameter<int>("LINEARIZATIONREINIT", "newton",
      "linearization scheme for nonlinear convective term of reinitialization equation",
      tuple<std::string>("newton", "fixed_point"),
      tuple<int>(INPAR::SCATRA::newton, INPAR::SCATRA::fixed_point), &ls_reinit);
  BoolParameter("CORRECTOR_STEP", "yes",
      "correction of interface position via volume constraint according to Sussman & Fatemi",
      &ls_reinit);
  DoubleParameter("CONVTOL_REINIT", -1.0,
      "tolerance for convergence check according to Sussman et al. 1994 (turned off negative)",
      &ls_reinit);

  BoolParameter(
      "REINITVOLCORRECTION", "No", "volume correction after reinitialization", &ls_reinit);

  DoubleParameter(
      "PENALTY_PARA", -1.0, "penalty parameter for elliptic reinitialization", &ls_reinit);

  setStringToIntegralParameter<int>("DIMENSION", "3D",
      "number of space dimensions for handling of quasi-2D problems with 3D elements",
      tuple<std::string>("3D", "2Dx", "2Dy", "2Dz"),
      tuple<int>(INPAR::SCATRA::ls_3D, INPAR::SCATRA::ls_2Dx, INPAR::SCATRA::ls_2Dy,
          INPAR::SCATRA::ls_2Dz),
      &ls_reinit);

  BoolParameter(
      "PROJECTION", "yes", "use L2-projection for grad phi and related quantities", &ls_reinit);
  DoubleParameter("PROJECTION_DIFF", 0.0, "use diffusive term for L2-projection", &ls_reinit);
  BoolParameter("LUMPING", "no", "use lumped mass matrix for L2-projection", &ls_reinit);

  setStringToIntegralParameter<int>("DIFF_FUNC", "hyperbolic", "function for diffusivity",
      tuple<std::string>("hyperbolic", "hyperbolic_smoothed_positive", "hyperbolic_clipped_05",
          "hyperbolic_clipped_1"),
      tuple<int>(INPAR::SCATRA::hyperbolic, INPAR::SCATRA::hyperbolic_smoothed_positive,
          INPAR::SCATRA::hyperbolic_clipped_05, INPAR::SCATRA::hyperbolic_clipped_1),
      &ls_reinit);
}



void INPAR::LEVELSET::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // Taylor Galerkin outflow Boundaries for level set transport equation

  Teuchos::RCP<ConditionDefinition> surfOutflowTaylorGalerkin =
      Teuchos::rcp(new ConditionDefinition("TAYLOR GALERKIN OUTFLOW SURF CONDITIONS",
          "TaylorGalerkinOutflow", "Surface Taylor Galerkin Outflow",
          DRT::Condition::TaylorGalerkinOutflow, true, DRT::Condition::Surface));

  condlist.push_back(surfOutflowTaylorGalerkin);

  /*--------------------------------------------------------------------*/

  Teuchos::RCP<ConditionDefinition> surfneumanninflowTaylorGalerkin =
      Teuchos::rcp(new ConditionDefinition("TAYLOR GALERKIN NEUMANN INFLOW SURF CONDITIONS",
          "TaylorGalerkinNeumannInflow", "Surface Taylor Galerkin Neumann Inflow",
          DRT::Condition::TaylorGalerkinNeumannInflow, true, DRT::Condition::Surface));

  condlist.push_back(surfneumanninflowTaylorGalerkin);


  /*--------------------------------------------------------------------*/
  // Characteristic Galerkin Boundaries for LevelSet-Reinitialization

  Teuchos::RCP<ConditionDefinition> surfreinitializationtaylorgalerkin =
      Teuchos::rcp(new ConditionDefinition("REINITIALIZATION TAYLOR GALERKIN SURF CONDITIONS",
          "ReinitializationTaylorGalerkin", "Surface Reinitialization Taylor Galerkin",
          DRT::Condition::ReinitializationTaylorGalerkin, true, DRT::Condition::Surface));

  condlist.push_back(surfreinitializationtaylorgalerkin);

  /*--------------------------------------------------------------------*/
  // level-set condition for contact points

  Teuchos::RCP<ConditionDefinition> linelscontact =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE LEVEL SET CONTACT CONDITION", "LsContact",
          "level-set condition for contact points", DRT::Condition::LsContact, false,
          DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> pointlscontact =
      Teuchos::rcp(new ConditionDefinition("DESIGN POINT LEVEL SET CONTACT CONDITION", "LsContact",
          "level-set condition for contact points", DRT::Condition::LsContact, false,
          DRT::Condition::Point));

  condlist.push_back(linelscontact);
  condlist.push_back(pointlscontact);
}
