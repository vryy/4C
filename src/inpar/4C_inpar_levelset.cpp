/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for level set

\level 2


*/
/*----------------------------------------------------------------------*/



#include "4C_inpar_levelset.hpp"

#include "4C_discretization_condition_definition.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


void Inpar::LevelSet::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& levelsetcontrol =
      list->sublist("LEVEL-SET CONTROL", false, "control parameters for level-set problems\n");

  Core::UTILS::IntParameter("NUMSTEP", 24, "Total number of time steps", &levelsetcontrol);
  Core::UTILS::DoubleParameter("TIMESTEP", 0.1, "Time increment dt", &levelsetcontrol);
  Core::UTILS::DoubleParameter("MAXTIME", 1000.0, "Total simulation time", &levelsetcontrol);
  Core::UTILS::IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &levelsetcontrol);
  Core::UTILS::IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &levelsetcontrol);

  setStringToIntegralParameter<int>("CALCERROR", "No",
      "compute error compared to analytical solution", tuple<std::string>("No", "InitialField"),
      tuple<int>(Inpar::ScaTra::calcerror_no_ls, Inpar::ScaTra::calcerror_initial_field),
      &levelsetcontrol);

  Core::UTILS::BoolParameter("EXTRACT_INTERFACE_VEL", "No",
      "replace computed velocity at nodes of given distance of interface by approximated interface "
      "velocity",
      &levelsetcontrol);
  Core::UTILS::IntParameter("NUM_CONVEL_LAYERS", -1,
      "number of layers around the interface which keep their computed convective velocity",
      &levelsetcontrol);


  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& ls_reinit = levelsetcontrol.sublist("REINITIALIZATION", false, "");

  setStringToIntegralParameter<int>("REINITIALIZATION", "None",
      "Type of reinitialization strategy for level set function",
      tuple<std::string>("None", "Signed_Distance_Function", "Sussman", "EllipticEq"),
      tuple<int>(Inpar::ScaTra::reinitaction_none,
          Inpar::ScaTra::reinitaction_signeddistancefunction, Inpar::ScaTra::reinitaction_sussman,
          Inpar::ScaTra::reinitaction_ellipticeq),
      &ls_reinit);

  Core::UTILS::BoolParameter("REINIT_INITIAL", "No",
      "Has level set field to be reinitialized before first time step?", &ls_reinit);
  Core::UTILS::IntParameter("REINITINTERVAL", 1, "reinitialization interval", &ls_reinit);

  // parameters for signed distance reinitialization
  Core::UTILS::BoolParameter("REINITBAND", "No",
      "reinitialization only within a band around the interface, or entire domain?", &ls_reinit);
  Core::UTILS::DoubleParameter("REINITBANDWIDTH", 1.0,
      "level-set value defining band width for reinitialization", &ls_reinit);

  // parameters for reinitialization equation
  Core::UTILS::IntParameter(
      "NUMSTEPSREINIT", 1, "(maximal) number of pseudo-time steps", &ls_reinit);
  // this parameter selects the tau definition applied
  setStringToIntegralParameter<int>("LINEARIZATIONREINIT", "fixed_point",
      "linearization of reinitialization equation", tuple<std::string>("newton", "fixed_point"),
      tuple<int>(Inpar::ScaTra::newton, Inpar::ScaTra::fixed_point), &ls_reinit);
  Core::UTILS::DoubleParameter("TIMESTEPREINIT", 1.0,
      "pseudo-time step length (usually a * characteristic element length of discretization with "
      "a>0)",
      &ls_reinit);
  Core::UTILS::DoubleParameter(
      "THETAREINIT", 1.0, "theta for time discretization of reinitialization equation", &ls_reinit);
  setStringToIntegralParameter<int>("STABTYPEREINIT", "SUPG", "type of stabilization (if any)",
      tuple<std::string>("no_stabilization", "SUPG", "GLS", "USFEM"),
      tuple<std::string>(
          "Do not use any stabilization -> only reasonable for low-Peclet-number flows", "Use SUPG",
          "Use GLS", "Use USFEM"),
      tuple<int>(Inpar::ScaTra::stabtype_no_stabilization, Inpar::ScaTra::stabtype_SUPG,
          Inpar::ScaTra::stabtype_GLS, Inpar::ScaTra::stabtype_USFEM),
      &ls_reinit);
  // this parameter selects the tau definition applied
  setStringToIntegralParameter<int>("DEFINITION_TAU_REINIT", "Taylor_Hughes_Zarins",
      "Definition of tau",
      tuple<std::string>("Taylor_Hughes_Zarins", "Taylor_Hughes_Zarins_wo_dt", "Franca_Valentin",
          "Franca_Valentin_wo_dt", "Shakib_Hughes_Codina", "Shakib_Hughes_Codina_wo_dt", "Codina",
          "Codina_wo_dt", "Exact_1D", "Zero"),
      tuple<int>(Inpar::ScaTra::tau_taylor_hughes_zarins,
          Inpar::ScaTra::tau_taylor_hughes_zarins_wo_dt, Inpar::ScaTra::tau_franca_valentin,
          Inpar::ScaTra::tau_franca_valentin_wo_dt, Inpar::ScaTra::tau_shakib_hughes_codina,
          Inpar::ScaTra::tau_shakib_hughes_codina_wo_dt, Inpar::ScaTra::tau_codina,
          Inpar::ScaTra::tau_codina_wo_dt, Inpar::ScaTra::tau_exact_1d, Inpar::ScaTra::tau_zero),
      &ls_reinit);
  // this parameter governs whether all-scale subgrid diffusivity is included
  setStringToIntegralParameter<int>("ARTDIFFREINIT", "no",
      "potential incorporation of all-scale subgrid diffusivity (a.k.a. discontinuity-capturing) "
      "term",
      tuple<std::string>("no", "isotropic", "crosswind"),
      tuple<std::string>("no artificial diffusion", "homogeneous artificial diffusion",
          "artificial diffusion in crosswind directions only"),
      tuple<int>(Inpar::ScaTra::artdiff_none, Inpar::ScaTra::artdiff_isotropic,
          Inpar::ScaTra::artdiff_crosswind),
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
      tuple<int>(Inpar::ScaTra::assgd_artificial, Inpar::ScaTra::assgd_lin_reinit,
          Inpar::ScaTra::assgd_hughes, Inpar::ScaTra::assgd_tezduyar,
          Inpar::ScaTra::assgd_tezduyar_wo_phizero, Inpar::ScaTra::assgd_docarmo,
          Inpar::ScaTra::assgd_almeida, Inpar::ScaTra::assgd_yzbeta, Inpar::ScaTra::assgd_codina),
      &ls_reinit);

  setStringToIntegralParameter<int>("SMOOTHED_SIGN_TYPE", "SussmanSmerekaOsher1994",
      "sign function for reinitialization equation",
      tuple<std::string>("NonSmoothed",
          "SussmanFatemi1999",  // smeared-out Heaviside function
          "SussmanSmerekaOsher1994", "PengEtAl1999"),
      tuple<int>(Inpar::ScaTra::signtype_nonsmoothed, Inpar::ScaTra::signtype_SussmanFatemi1999,
          Inpar::ScaTra::signtype_SussmanSmerekaOsher1994, Inpar::ScaTra::signtype_PengEtAl1999),
      &ls_reinit);
  setStringToIntegralParameter<int>("CHARELELENGTHREINIT", "root_of_volume",
      "characteristic element length for sign function",
      tuple<std::string>("root_of_volume", "streamlength"),
      tuple<int>(Inpar::ScaTra::root_of_volume_reinit, Inpar::ScaTra::streamlength_reinit),
      &ls_reinit);
  Core::UTILS::DoubleParameter("INTERFACE_THICKNESS", 1.0,
      "factor for interface thickness (multiplied by element length)", &ls_reinit);
  setStringToIntegralParameter<int>("VELREINIT", "integration_point_based",
      "evaluate velocity at integration point or compute node-based velocity",
      tuple<std::string>("integration_point_based", "node_based"),
      tuple<int>(
          Inpar::ScaTra::vel_reinit_integration_point_based, Inpar::ScaTra::vel_reinit_node_based),
      &ls_reinit);
  setStringToIntegralParameter<int>("LINEARIZATIONREINIT", "newton",
      "linearization scheme for nonlinear convective term of reinitialization equation",
      tuple<std::string>("newton", "fixed_point"),
      tuple<int>(Inpar::ScaTra::newton, Inpar::ScaTra::fixed_point), &ls_reinit);
  Core::UTILS::BoolParameter("CORRECTOR_STEP", "yes",
      "correction of interface position via volume constraint according to Sussman & Fatemi",
      &ls_reinit);
  Core::UTILS::DoubleParameter("CONVTOL_REINIT", -1.0,
      "tolerance for convergence check according to Sussman et al. 1994 (turned off negative)",
      &ls_reinit);

  Core::UTILS::BoolParameter(
      "REINITVOLCORRECTION", "No", "volume correction after reinitialization", &ls_reinit);

  Core::UTILS::DoubleParameter(
      "PENALTY_PARA", -1.0, "penalty parameter for elliptic reinitialization", &ls_reinit);

  setStringToIntegralParameter<int>("DIMENSION", "3D",
      "number of space dimensions for handling of quasi-2D problems with 3D elements",
      tuple<std::string>("3D", "2Dx", "2Dy", "2Dz"),
      tuple<int>(Inpar::ScaTra::ls_3D, Inpar::ScaTra::ls_2Dx, Inpar::ScaTra::ls_2Dy,
          Inpar::ScaTra::ls_2Dz),
      &ls_reinit);

  Core::UTILS::BoolParameter(
      "PROJECTION", "yes", "use L2-projection for grad phi and related quantities", &ls_reinit);
  Core::UTILS::DoubleParameter(
      "PROJECTION_DIFF", 0.0, "use diffusive term for L2-projection", &ls_reinit);
  Core::UTILS::BoolParameter(
      "LUMPING", "no", "use lumped mass matrix for L2-projection", &ls_reinit);

  setStringToIntegralParameter<int>("DIFF_FUNC", "hyperbolic", "function for diffusivity",
      tuple<std::string>("hyperbolic", "hyperbolic_smoothed_positive", "hyperbolic_clipped_05",
          "hyperbolic_clipped_1"),
      tuple<int>(Inpar::ScaTra::hyperbolic, Inpar::ScaTra::hyperbolic_smoothed_positive,
          Inpar::ScaTra::hyperbolic_clipped_05, Inpar::ScaTra::hyperbolic_clipped_1),
      &ls_reinit);
}



void Inpar::LevelSet::SetValidConditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  /*--------------------------------------------------------------------*/
  // Taylor Galerkin outflow Boundaries for level set transport equation

  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfOutflowTaylorGalerkin = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("TAYLOR GALERKIN OUTFLOW SURF CONDITIONS",
          "TaylorGalerkinOutflow", "Surface Taylor Galerkin Outflow",
          Core::Conditions::TaylorGalerkinOutflow, true, Core::Conditions::geometry_type_surface));

  condlist.push_back(surfOutflowTaylorGalerkin);

  /*--------------------------------------------------------------------*/

  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfneumanninflowTaylorGalerkin =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "TAYLOR GALERKIN NEUMANN INFLOW SURF CONDITIONS", "TaylorGalerkinNeumannInflow",
          "Surface Taylor Galerkin Neumann Inflow", Core::Conditions::TaylorGalerkinNeumannInflow,
          true, Core::Conditions::geometry_type_surface));

  condlist.push_back(surfneumanninflowTaylorGalerkin);


  /*--------------------------------------------------------------------*/
  // Characteristic Galerkin Boundaries for LevelSet-reinitialization

  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfreinitializationtaylorgalerkin =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "REINITIALIZATION TAYLOR GALERKIN SURF CONDITIONS", "ReinitializationTaylorGalerkin",
          "Surface reinitialization Taylor Galerkin",
          Core::Conditions::ReinitializationTaylorGalerkin, true,
          Core::Conditions::geometry_type_surface));

  condlist.push_back(surfreinitializationtaylorgalerkin);

  /*--------------------------------------------------------------------*/
  // level-set condition for contact points

  Teuchos::RCP<Core::Conditions::ConditionDefinition> linelscontact = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN LINE LEVEL SET CONTACT CONDITION",
          "LsContact", "level-set condition for contact points", Core::Conditions::LsContact, false,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> pointlscontact = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN POINT LEVEL SET CONTACT CONDITION",
          "LsContact", "level-set condition for contact points", Core::Conditions::LsContact, false,
          Core::Conditions::geometry_type_point));

  condlist.push_back(linelscontact);
  condlist.push_back(pointlscontact);
}

FOUR_C_NAMESPACE_CLOSE
