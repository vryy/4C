/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for fluid and related problems

\maintainer Martin Kronbichler

\level 1
*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_fluid.H"
#include "inpar_turbulence.H"
#include "inpar.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::FLUID::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::Array<std::string> yesnotuple =
      tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

  Teuchos::ParameterList& fdyn = list->sublist("FLUID DYNAMIC", false, "");

  // physical type of fluid flow (incompressible, varying density, loma, Boussinesq approximation,
  // temperature-dependent water)
  setStringToIntegralParameter<int>("PHYSICAL_TYPE", "Incompressible", "Physical Type",
      tuple<std::string>("Incompressible", "Weakly_compressible", "Weakly_compressible_stokes",
          "Weakly_compressible_dens_mom", "Weakly_compressible_stokes_dens_mom",
          "Artificial_compressibility", "Varying_density", "Loma", "Temp_dep_water", "Boussinesq",
          "Topology_optimization", "Stokes", "Oseen"),
      tuple<int>(incompressible, weakly_compressible, weakly_compressible_stokes,
          weakly_compressible_dens_mom, weakly_compressible_stokes_dens_mom, artcomp,
          varying_density, loma, tempdepwater, boussinesq, topopt, stokes, oseen),
      &fdyn);

  // number of linear solver used for fluid problem
  IntParameter("LINEAR_SOLVER", -1, "number of linear solver used for fluid dynamics", &fdyn);

  // number of linear solver used for fluid problem (former fluid pressure solver for SIMPLER
  // preconditioning with fluid)
  IntParameter("SIMPLER_SOLVER", -1,
      "number of linear solver used for fluid dynamics (ONLY NECESSARY FOR BlockGaussSeidel solver "
      "block within fluid mehstying case any more!!!!)",
      &fdyn);

  // Flag to define the way of calculating stresses and wss
  setStringToIntegralParameter<int>("WSS_TYPE", "Standard", "which type of stresses and wss",
      tuple<std::string>("Standard", "Aggregation", "Mean"),
      tuple<std::string>(
          "calculate 'normal' wss", "calculate aggregated wss", "calculate mean wss"),
      tuple<int>(wss_standard, wss_aggregation, wss_mean), &fdyn);

  // Set ML-solver number for smooting of residual-based calculated wallshearstress via plain
  // aggregation.
  IntParameter("WSS_ML_AGR_SOLVER", -1,
      "Set ML-solver number for smoothing of residual-based calculated wallshearstress via plain "
      "aggregation.",
      &fdyn);

  setStringToIntegralParameter<int>("TIMEINTEGR", "One_Step_Theta", "Time Integration Scheme",
      tuple<std::string>("Stationary", "Np_Gen_Alpha", "Af_Gen_Alpha", "One_Step_Theta", "BDF2"),
      tuple<int>(timeint_stationary, timeint_npgenalpha, timeint_afgenalpha, timeint_one_step_theta,
          timeint_bdf2),
      &fdyn);

  setStringToIntegralParameter<int>("OST_CONT_PRESS", "Cont_normal_Press_normal",
      "One step theta option for time discretization of continuity eq. and pressure",
      tuple<std::string>(
          "Cont_normal_Press_normal", "Cont_impl_Press_normal", "Cont_impl_Press_impl"),
      tuple<int>(Cont_normal_Press_normal, Cont_impl_Press_normal, Cont_impl_Press_impl), &fdyn);

  setStringToIntegralParameter<int>("GEOMETRY", "full", "How the geometry is specified",
      tuple<std::string>("full", "box", "file"),
      tuple<int>(INPAR::geometry_full, INPAR::geometry_box, INPAR::geometry_file), &fdyn);

  setStringToIntegralParameter<int>("NONLINITER", "fixed_point_like", "Nonlinear iteration scheme",
      tuple<std::string>("fixed_point_like", "Newton"), tuple<int>(fixed_point_like, Newton),
      &fdyn);

  setStringToIntegralParameter<int>("PREDICTOR", "steady_state",
      "Predictor for first guess in nonlinear iteration",
      tuple<std::string>("steady_state", "zero_acceleration", "constant_acceleration",
          "constant_increment", "explicit_second_order_midpoint", "TangVel"),
      tuple<int>(1, 2, 3, 4, 5, 6), &fdyn);

  setStringToIntegralParameter<int>("CONVCHECK", "L_2_norm", "norm for convergence check",
      tuple<std::string>("L_2_norm"),
      tuple<std::string>("compute L2 errors of increments (relative) and residuals (absolute)"),
      tuple<int>(fncc_L2), &fdyn);

  BoolParameter("INCONSISTENT_RESIDUAL", "No",
      "do not evaluate residual after solution has converged (->faster)", &fdyn);

  {
    // a standard Teuchos::tuple can have at maximum 10 entries! We have to circumvent this here.
    Teuchos::Tuple<std::string, 13> name;
    Teuchos::Tuple<int, 13> label;
    name[0] = "zero_field";
    label[0] = initfield_zero_field;
    name[1] = "field_by_function";
    label[1] = initfield_field_by_function;
    name[2] = "disturbed_field_from_function";
    label[2] = initfield_disturbed_field_from_function;
    name[3] = "FLAME_VORTEX_INTERACTION";
    label[3] = initfield_flame_vortex_interaction;
    name[4] = "BELTRAMI-FLOW";
    label[4] = initfield_beltrami_flow;
    name[5] = "KIM-MOIN-FLOW";
    label[5] = initfield_kim_moin_flow;
    name[6] = "BOCHEV-TEST";
    label[6] = initfield_bochev_test;
    name[7] = "hit_comte_bellot_corrsin_initial_field";
    label[7] = initfield_hit_comte_bellot_corrsin;
    name[8] = "forced_hit_simple_algebraic_spectrum";
    label[8] = initfield_forced_hit_simple_algebraic_spectrum;
    name[9] = "forced_hit_numeric_spectrum";
    label[9] = initfield_forced_hit_numeric_spectrum;
    name[10] = "forced_hit_passive";
    label[10] = initfield_passive_hit_const_input;
    name[11] = "channel_weakly_compressible";
    label[11] = initfield_channel_weakly_compressible;
    name[12] = "channel_weakly_compressible_fourier_3";
    label[12] = initfield_channel_weakly_compressible_fourier_3;

    setStringToIntegralParameter<int>(
        "INITIALFIELD", "zero_field", "Initial field for fluid problem", name, label, &fdyn);
  }

  IntParameter("OSEENFIELDFUNCNO", -1, "function number of Oseen advective field", &fdyn);

  BoolParameter("LIFTDRAG", "No", "Calculate lift and drag forces along specified boundary", &fdyn);

  setStringToIntegralParameter<int>("CONVFORM", "convective", "form of convective term",
      tuple<std::string>("convective", "conservative"), tuple<int>(0, 1), &fdyn);

  setStringToIntegralParameter<int>("NONLINEARBC", "no",
      "Flag to activate check for potential nonlinear boundary conditions",
      tuple<std::string>("no", "yes"),
      tuple<std::string>(
          "no nonlinear boundary conditions", "nonlinear boundary conditions might be present"),
      tuple<int>(0, 1), &fdyn);

  setStringToIntegralParameter<int>("MESHTYING", "no", "Flag to (de)activate mesh tying algorithm",
      tuple<std::string>("no", "Condensed_Smat", "Condensed_Bmat", "Condensed_Bmat_merged"),
      tuple<int>(no_meshtying, condensed_smat, condensed_bmat, condensed_bmat_merged), &fdyn);

  setStringToIntegralParameter<int>("GRIDVEL", "BE",
      "scheme for determination of gridvelocity from displacements",
      tuple<std::string>("BE", "BDF2", "OST"), tuple<int>(BE, BDF2, OST), &fdyn);

  BoolParameter("ALLDOFCOUPLED", "Yes", "all dof (incl. pressure) are coupled", &fdyn);

  {
    Teuchos::Tuple<std::string, 19> name;
    Teuchos::Tuple<int, 19> label;

    name[0] = "no";
    label[0] = no_error_calculation;
    name[1] = "beltrami_flow";
    label[1] = beltrami_flow;
    name[2] = "channel2D";
    label[2] = channel2D;
    name[3] = "gravitation";
    label[3] = gravitation;
    name[4] = "shear_flow";
    label[4] = shear_flow;
    name[5] = "jeffery_hamel_flow";
    label[5] = jeffery_hamel_flow;
    name[6] = "byfunct";
    label[6] = byfunct;
    name[7] = "beltrami_stat_stokes";
    label[7] = beltrami_stat_stokes;
    name[8] = "beltrami_stat_navier_stokes";
    label[8] = beltrami_stat_navier_stokes;
    name[9] = "beltrami_instat_stokes";
    label[9] = beltrami_instat_stokes;
    name[10] = "beltrami_instat_navier_stokes";
    label[10] = beltrami_instat_navier_stokes;
    name[11] = "kimmoin_stat_stokes";
    label[11] = kimmoin_stat_stokes;
    name[12] = "kimmoin_stat_navier_stokes";
    label[12] = kimmoin_stat_navier_stokes;
    name[13] = "kimmoin_instat_stokes";
    label[13] = kimmoin_instat_stokes;
    name[14] = "kimmoin_instat_navier_stokes";
    label[14] = kimmoin_instat_navier_stokes;
    name[15] = "fsi_fluid_pusher";
    label[15] = fsi_fluid_pusher;
    name[16] = "topopt_channel";
    label[16] = topoptchannel;
    name[17] = "channel_weakly_compressible";
    label[17] = channel_weakly_compressible;
    name[18] = "channel_weakly_compressible_fourier_3";
    label[18] = channel_weakly_compressible_fourier_3;

    setStringToIntegralParameter<int>(
        "CALCERROR", "no", "Flag to (de)activate error calculations", name, label, &fdyn);
  }
  IntParameter("CALCERRORFUNCNO", -1, "Function for Error Calculation", &fdyn);

  IntParameter("CORRTERMFUNCNO", -1,
      "Function for calculation of the correction term for the weakly compressible problem", &fdyn);

  IntParameter("BODYFORCEFUNCNO", -1,
      "Function for calculation of the body force for the weakly compressible problem", &fdyn);

  DoubleParameter("STAB_DEN_REF", 0.0,
      "Reference stabilization parameter for the density for the HDG weakly compressible "
      "formulation",
      &fdyn);

  DoubleParameter("STAB_MOM_REF", 0.0,
      "Reference stabilization parameter for the momentum for the HDG weakly compressible "
      "formulation",
      &fdyn);

  IntParameter("VARVISCFUNCNO", -1,
      "Function for calculation of a variable viscosity for the weakly compressible problem",
      &fdyn);

  {
    Teuchos::Tuple<std::string, 2> name;
    Teuchos::Tuple<int, 2> label;

    name[0] = "no";
    label[0] = no_pressure_average_bc;
    name[1] = "yes";
    label[1] = yes_pressure_average_bc;

    setStringToIntegralParameter<int>("PRESSAVGBC", "no",
        "Flag to (de)activate imposition of boundary condition for the considered element average "
        "pressure",
        name, label, &fdyn);
  }

  DoubleParameter("REFMACH", 1.0, "Reference Mach number", &fdyn);

  setStringToIntegralParameter<int>("SIMPLER", "no",
      "Switch on SIMPLE family of solvers, only works with block preconditioners like CheapSIMPLE!",
      yesnotuple, yesnovalue, &fdyn);

  setStringToIntegralParameter<int>("ADAPTCONV", "yes",
      "Switch on adaptive control of linear solver tolerance for nonlinear solution", yesnotuple,
      yesnovalue, &fdyn);
  DoubleParameter("ADAPTCONV_BETTER", 0.1,
      "The linear solver shall be this much better than the current nonlinear residual in the "
      "nonlinear convergence limit",
      &fdyn);

  setStringToIntegralParameter<int>("INFNORMSCALING", "no",
      "Scale blocks of matrix with row infnorm?", yesnotuple, yesnovalue, &fdyn);

  BoolParameter("GMSH_OUTPUT", "No", "write output to gmsh files", &fdyn);
  BoolParameter(
      "COMPUTE_DIVU", "No", "Compute divergence of velocity field at the element center", &fdyn);
  BoolParameter("COMPUTE_EKIN", "No",
      "Compute kinetic energy at the end of each time step and write it to file.", &fdyn);
  BoolParameter("NEW_OST", "No",
      "Solve the Navier-Stokes equation with the new One Step Theta algorithm",
      &fdyn);  // TODO: To be removed.
  IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &fdyn);
  IntParameter("RESTARTEVRY", 20, "Increment for writing restart", &fdyn);
  IntParameter("NUMSTEP", 1, "Total number of Timesteps", &fdyn);
  IntParameter("STEADYSTEP", -1, "steady state check every step", &fdyn);
  IntParameter("NUMSTASTEPS", 0, "Number of Steps for Starting Scheme", &fdyn);
  IntParameter("STARTFUNCNO", -1, "Function for Initial Starting Field", &fdyn);
  IntParameter("ITEMAX", 10, "max. number of nonlin. iterations", &fdyn);
  IntParameter("INITSTATITEMAX", 5,
      "max number of nonlinear iterations for initial stationary solution", &fdyn);
  DoubleParameter("TIMESTEP", 0.01, "Time increment dt", &fdyn);
  DoubleParameter("MAXTIME", 1000.0, "Total simulation time", &fdyn);
  DoubleParameter("ALPHA_M", 1.0, "Time integration factor", &fdyn);
  DoubleParameter("ALPHA_F", 1.0, "Time integration factor", &fdyn);
  DoubleParameter("GAMMA", 1.0, "Time integration factor", &fdyn);
  DoubleParameter("THETA", 0.66, "Time integration factor", &fdyn);

  DoubleParameter("START_THETA", 1.0, "Time integration factor for starting scheme", &fdyn);

  setStringToIntegralParameter<int>("STRONG_REDD_3D_COUPLING_TYPE", "no",
      "Flag to (de)activate potential Strong 3D redD coupling", tuple<std::string>("no", "yes"),
      tuple<std::string>("Weak coupling", "Strong coupling"), tuple<int>(0, 1), &fdyn);

  IntParameter("VELGRAD_PROJ_SOLVER", -1, "Number of linear solver used for L2 projection", &fdyn);

  setStringToIntegralParameter<int>("VELGRAD_PROJ_METHOD", "none",
      "Flag to (de)activate gradient reconstruction.",
      tuple<std::string>("none", "superconvergent_patch_recovery", "L2_projection"),
      tuple<std::string>("no gradient reconstruction",
          "gradient reconstruction via superconvergent patch recovery",
          "gracient reconstruction via l2-projection"),
      tuple<int>(gradreco_none,  // no convective streamline edge-based stabilization
          gradreco_spr,  // convective streamline edge-based stabilization on the entire domain
          gradreco_l2    // pressure edge-based stabilization as ghost penalty around cut elements
          ),
      &fdyn);

  BoolParameter("OFF_PROC_ASSEMBLY", "No",
      "Do not evaluate ghosted elements but communicate them --> faster if element call is "
      "expensive",
      &fdyn);
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fdyn_nln = fdyn.sublist("NONLINEAR SOLVER TOLERANCES", false, "");

  DoubleParameter(
      "TOL_VEL_RES", 1e-6, "Tolerance for convergence check of velocity residual", &fdyn_nln);

  DoubleParameter(
      "TOL_VEL_INC", 1e-6, "Tolerance for convergence check of velocity increment", &fdyn_nln);

  DoubleParameter(
      "TOL_PRES_RES", 1e-6, "Tolerance for convergence check of pressure residual", &fdyn_nln);

  DoubleParameter(
      "TOL_PRES_INC", 1e-6, "Tolerance for convergence check of pressure increment", &fdyn_nln);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fdyn_stab = fdyn.sublist("RESIDUAL-BASED STABILIZATION", false, "");

  // this parameter defines various stabilized methods
  setStringToIntegralParameter<int>("STABTYPE", "residual_based",
      "Apply (un)stabilized fluid formulation",
      tuple<std::string>("no_stabilization", "residual_based", "edge_based", "pressure_projection"),
      tuple<std::string>("Do not use any stabilization -> inf-sup stable elements required!",
          "Use a residual-based stabilization or, more generally, a stabilization \nbased on the "
          "concept of the residual-based variational multiscale method...\nExpecting additional "
          "input",
          "Use an edge-based stabilization, especially for XFEM",
          "Element/cell based polynomial pressure projection, see Dohrmann/Bochev 2004, IJNMF"),
      tuple<int>(
          stabtype_nostab, stabtype_residualbased, stabtype_edgebased, stabtype_pressureprojection),
      &fdyn_stab);

  BoolParameter("INCONSISTENT", "No",
      "residual based without second derivatives (i.e. only consistent for tau->0, but faster)",
      &fdyn_stab);

  BoolParameter("Reconstruct_Sec_Der", "No",
      "residual computed with a reconstruction of the second derivatives via projection or "
      "superconvergent patch recovery",
      &fdyn_stab);

  // the following parameters are necessary only if a residual based stabilized method is applied
  setStringToIntegralParameter<int>("TDS", "quasistatic",
      "Flag to allow time dependency of subscales for residual-based stabilization.",
      tuple<std::string>("quasistatic", "time_dependent"),
      tuple<std::string>("Use a quasi-static residual-based stabilization (standard case)",
          "Residual-based stabilization including time evolution equations for subscales"),
      tuple<int>(subscales_quasistatic, subscales_time_dependent), &fdyn_stab);

  setStringToIntegralParameter<int>("TRANSIENT", "no_transient",
      "Specify how to treat the transient term.",
      tuple<std::string>("no_transient", "yes_transient", "transient_complete"),
      tuple<std::string>(
          "Do not use transient term (currently only opportunity for quasistatic stabilization)",
          "Use transient term (recommended for time dependent subscales)",
          "Use transient term including a linearisation of 1/tau"),
      tuple<int>(inertia_stab_drop, inertia_stab_keep, inertia_stab_keep_complete), &fdyn_stab);

  BoolParameter("PSPG", "Yes", "Flag to (de)activate PSPG stabilization.", &fdyn_stab);
  BoolParameter("SUPG", "Yes", "Flag to (de)activate SUPG stabilization.", &fdyn_stab);
  BoolParameter("GRAD_DIV", "Yes", "Flag to (de)activate grad-div term.", &fdyn_stab);

  setStringToIntegralParameter<int>("VSTAB", "no_vstab",
      "Flag to (de)activate viscous term in residual-based stabilization.",
      tuple<std::string>(
          "no_vstab", "vstab_gls", "vstab_gls_rhs", "vstab_usfem", "vstab_usfem_rhs"),
      tuple<std::string>("No viscous term in stabilization", "Viscous stabilization of GLS type",
          "Viscous stabilization of GLS type, included only on the right hand side",
          "Viscous stabilization of USFEM type",
          "Viscous stabilization of USFEM type, included only on the right hand side"),
      tuple<int>(viscous_stab_none, viscous_stab_gls, viscous_stab_gls_only_rhs, viscous_stab_usfem,
          viscous_stab_usfem_only_rhs),
      &fdyn_stab);

  setStringToIntegralParameter<int>("RSTAB", "no_rstab",
      "Flag to (de)activate reactive term in residual-based stabilization.",
      tuple<std::string>("no_rstab", "rstab_gls", "rstab_usfem"),
      tuple<std::string>("no reactive term in stabilization", "reactive stabilization of GLS type",
          "reactive stabilization of USFEM type"),
      tuple<int>(reactive_stab_none, reactive_stab_gls, reactive_stab_usfem), &fdyn_stab);

  setStringToIntegralParameter<int>("CROSS-STRESS", "no_cross",
      "Flag to (de)activate cross-stress term -> residual-based VMM.",
      tuple<std::string>("no_cross", "yes_cross", "cross_rhs"
          //"cross_complete"
          ),
      tuple<std::string>("No cross-stress term",
          "Include the cross-stress term with a linearization of the convective part",
          "Include cross-stress term, but only explicitly on right hand side"
          //""
          ),
      tuple<int>(cross_stress_stab_none, cross_stress_stab, cross_stress_stab_only_rhs),
      &fdyn_stab);

  setStringToIntegralParameter<int>("REYNOLDS-STRESS", "no_reynolds",
      "Flag to (de)activate Reynolds-stress term -> residual-based VMM.",
      tuple<std::string>("no_reynolds", "yes_reynolds", "reynolds_rhs"
          //"reynolds_complete"
          ),
      tuple<std::string>(
          "No Reynolds-stress term", "Include Reynolds-stress term with linearisation",
          "Include Reynolds-stress term explicitly on right hand side"
          //""
          ),
      tuple<int>(reynolds_stress_stab_none, reynolds_stress_stab, reynolds_stress_stab_only_rhs),
      &fdyn_stab);

  {
    // this parameter selects the tau definition applied
    // a standard Teuchos::tuple can have at maximum 10 entries! We have to circumvent this here.
    Teuchos::Tuple<std::string, 16> name;
    Teuchos::Tuple<int, 16> label;
    name[0] = "Taylor_Hughes_Zarins";
    label[0] = tau_taylor_hughes_zarins;
    name[1] = "Taylor_Hughes_Zarins_wo_dt";
    label[1] = tau_taylor_hughes_zarins_wo_dt;
    name[2] = "Taylor_Hughes_Zarins_Whiting_Jansen";
    label[2] = tau_taylor_hughes_zarins_whiting_jansen;
    name[3] = "Taylor_Hughes_Zarins_Whiting_Jansen_wo_dt";
    label[3] = tau_taylor_hughes_zarins_whiting_jansen_wo_dt;
    name[4] = "Taylor_Hughes_Zarins_scaled";
    label[4] = tau_taylor_hughes_zarins_scaled;
    name[5] = "Taylor_Hughes_Zarins_scaled_wo_dt";
    label[5] = tau_taylor_hughes_zarins_scaled_wo_dt;
    name[6] = "Franca_Barrenechea_Valentin_Frey_Wall";
    label[6] = tau_franca_barrenechea_valentin_frey_wall;
    name[7] = "Franca_Barrenechea_Valentin_Frey_Wall_wo_dt";
    label[7] = tau_franca_barrenechea_valentin_frey_wall_wo_dt;
    name[8] = "Shakib_Hughes_Codina";
    label[8] = tau_shakib_hughes_codina;
    name[9] = "Shakib_Hughes_Codina_wo_dt";
    label[9] = tau_shakib_hughes_codina_wo_dt;
    name[10] = "Codina";
    label[10] = tau_codina;
    name[11] = "Codina_wo_dt";
    label[11] = tau_codina_wo_dt;
    name[12] = "Codina_convscaled";
    label[12] = tau_codina_convscaled;
    name[13] = "Franca_Madureira_Valentin_Badia_Codina";
    label[13] = tau_franca_madureira_valentin_badia_codina;
    name[14] = "Franca_Madureira_Valentin_Badia_Codina_wo_dt";
    label[14] = tau_franca_madureira_valentin_badia_codina_wo_dt;
    name[15] = "Hughes_Franca_Balestra_wo_dt";
    label[15] = tau_hughes_franca_balestra_wo_dt;

    setStringToIntegralParameter<int>("DEFINITION_TAU", "Franca_Barrenechea_Valentin_Frey_Wall",
        "Definition of tau_M and Tau_C", name, label, &fdyn_stab);
  }

  // this parameter selects the characteristic element length for tau_Mu for all
  // stabilization parameter definitions requiring such a length
  setStringToIntegralParameter<int>("CHARELELENGTH_U", "streamlength",
      "Characteristic element length for tau_Mu",
      tuple<std::string>("streamlength", "volume_equivalent_diameter", "root_of_volume"),
      tuple<int>(streamlength_u, volume_equivalent_diameter_u, root_of_volume_u), &fdyn_stab);

  // this parameter selects the characteristic element length for tau_Mp and tau_C for
  // all stabilization parameter definitions requiring such a length
  setStringToIntegralParameter<int>("CHARELELENGTH_PC", "volume_equivalent_diameter",
      "Characteristic element length for tau_Mp/tau_C",
      tuple<std::string>("streamlength", "volume_equivalent_diameter", "root_of_volume"),
      tuple<int>(streamlength_pc, volume_equivalent_diameter_pc, root_of_volume_pc), &fdyn_stab);

  // this parameter selects the location where tau is evaluated
  setStringToIntegralParameter<int>("EVALUATION_TAU", "element_center",
      "Location where tau is evaluated", tuple<std::string>("element_center", "integration_point"),
      tuple<std::string>("evaluate tau at element center", "evaluate tau at integration point"),
      tuple<int>(0, 1), &fdyn_stab);

  // this parameter selects the location where the material law is evaluated
  // (does not fit here very well, but parameter transfer is easier)
  setStringToIntegralParameter<int>("EVALUATION_MAT", "element_center",
      "Location where material law is evaluated",
      tuple<std::string>("element_center", "integration_point"),
      tuple<std::string>(
          "evaluate material law at element center", "evaluate material law at integration point"),
      tuple<int>(0, 1), &fdyn_stab);

  // these parameters active additional terms in loma continuity equation
  // which might be identified as SUPG-/cross- and Reynolds-stress term
  BoolParameter("LOMA_CONTI_SUPG", "No",
      "Flag to (de)activate SUPG stabilization in loma continuity equation.", &fdyn_stab);

  setStringToIntegralParameter<int>("LOMA_CONTI_CROSS_STRESS", "no_cross",
      "Flag to (de)activate cross-stress term loma continuity equation-> residual-based VMM.",
      tuple<std::string>("no_cross", "yes_cross", "cross_rhs"
          //"cross_complete"
          ),
      tuple<std::string>("No cross-stress term",
          "Include the cross-stress term with a linearization of the convective part",
          "Include cross-stress term, but only explicitly on right hand side"
          //""
          ),
      tuple<int>(cross_stress_stab_none, cross_stress_stab, cross_stress_stab_only_rhs),
      &fdyn_stab);

  setStringToIntegralParameter<int>("LOMA_CONTI_REYNOLDS_STRESS", "no_reynolds",
      "Flag to (de)activate Reynolds-stress term loma continuity equation-> residual-based VMM.",
      tuple<std::string>("no_reynolds", "yes_reynolds", "reynolds_rhs"),
      tuple<std::string>(
          "No Reynolds-stress term", "Include Reynolds-stress term with linearisation",
          "Include Reynolds-stress term explicitly on right hand side"
          //""
          ),
      tuple<int>(reynolds_stress_stab_none, reynolds_stress_stab, reynolds_stress_stab_only_rhs),
      &fdyn_stab);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fdyn_edge_based_stab =
      fdyn.sublist("EDGE-BASED STABILIZATION", false, "");

  //! Flag to (de)activate edge-based (EOS) pressure stabilization
  setStringToIntegralParameter<int>("EOS_PRES", "none",
      "Flag to (de)activate pressure edge-based stabilization.",
      tuple<std::string>("none", "std_eos", "xfem_gp"),
      tuple<std::string>("do not use pressure edge-based stabilization",
          "use pressure edge-based stabilization as standard edge-based stabilization on the "
          "entire domain",
          "use pressure edge-based stabilization as xfem ghost-penalty stabilization just around "
          "cut elements"),
      tuple<int>(EOS_PRES_none,  // no pressure edge-based stabilization
          EOS_PRES_std_eos,      // pressure edge-based stabilization on the entire domain
          EOS_PRES_xfem_gp       // pressure edge-based stabilization as ghost penalty around cut
                                 // elements
          ),
      &fdyn_edge_based_stab);

  //! Flag to (de)activate edge-based (EOS) convective streamline stabilization
  setStringToIntegralParameter<int>("EOS_CONV_STREAM", "none",
      "Flag to (de)activate convective streamline edge-based stabilization.",
      tuple<std::string>("none", "std_eos", "xfem_gp"),
      tuple<std::string>("do not use convective streamline edge-based stabilization",
          "use convective streamline edge-based stabilization as standard edge-based stabilization "
          "on the entire domain",
          "use convective streamline edge-based stabilization as xfem ghost-penalty stabilization "
          "just around cut elements"),
      tuple<int>(EOS_CONV_STREAM_none,  // no convective streamline edge-based stabilization
          EOS_CONV_STREAM_std_eos,  // convective streamline edge-based stabilization on the entire
                                    // domain
          EOS_CONV_STREAM_xfem_gp   // pressure edge-based stabilization as ghost penalty around cut
                                    // elements
          ),
      &fdyn_edge_based_stab);

  //! Flag to (de)activate edge-based (EOS) convective crosswind stabilization
  setStringToIntegralParameter<int>("EOS_CONV_CROSS", "none",
      "Flag to (de)activate convective crosswind edge-based stabilization.",
      tuple<std::string>("none", "std_eos", "xfem_gp"),
      tuple<std::string>("do not use convective crosswind edge-based stabilization",
          "use convective crosswind edge-based stabilization as standard edge-based stabilization "
          "on the entire domain",
          "use convective crosswind edge-based stabilization as xfem ghost-penalty stabilization "
          "just around cut elements"),
      tuple<int>(EOS_CONV_CROSS_none,  // no convective crosswind edge-based stabilization
          EOS_CONV_CROSS_std_eos,  // convective crosswind edge-based stabilization on the entire
                                   // domain
          EOS_CONV_CROSS_xfem_gp   // convective crosswind edge-based stabilization as ghost penalty
                                   // around cut elements
          ),
      &fdyn_edge_based_stab);

  //! Flag to (de)activate edge-based (EOS) divergence stabilization
  setStringToIntegralParameter<int>("EOS_DIV", "none",
      "Flag to (de)activate divergence edge-based stabilization.",
      tuple<std::string>(
          "none", "vel_jump_std_eos", "vel_jump_xfem_gp", "div_jump_std_eos", "div_jump_xfem_gp"),
      tuple<std::string>("do not use divergence edge-based stabilization",
          "divergence edge-based stabilization based on velocity jump on the entire domain",
          "divergence edge-based stabilization based on divergence jump just around cut elements",
          "divergence edge-based stabilization based on velocity jump on the entire domain",
          "divergence edge-based stabilization based on divergence jump just around cut elements"),
      tuple<int>(EOS_DIV_none,       // no convective edge-based stabilization
          EOS_DIV_vel_jump_std_eos,  // streamline convective edge-based stabilization
          EOS_DIV_vel_jump_xfem_gp,  // streamline convective edge-based stabilization
          EOS_DIV_div_jump_std_eos,  // crosswind convective edge-based stabilization
          EOS_DIV_div_jump_xfem_gp   // crosswind convective edge-based stabilization
          ),
      &fdyn_edge_based_stab);

  //! special least-squares condition for pseudo 2D examples where pressure level is determined via
  //! Krylov-projection
  BoolParameter("PRES_KRYLOV_2Dz", "No",
      "residual based without second derivatives (i.e. only consistent for tau->0, but faster)",
      &fdyn_edge_based_stab);

  //! this parameter selects the definition of Edge-based stabilization parameter
  setStringToIntegralParameter<int>("EOS_DEFINITION_TAU", "Burman_Hansbo_DAngelo_Zunino",
      "Definition of stabilization parameter for edge-based stabilization",
      tuple<std::string>("Burman_Fernandez_Hansbo", "Burman_Fernandez_Hansbo_wo_dt",
          "Braack_Burman_John_Lube", "Braack_Burman_John_Lube_wo_divjump",
          "Franca_Barrenechea_Valentin_Wall", "Burman_Fernandez", "Burman_Hansbo_DAngelo_Zunino",
          "Burman_Hansbo_DAngelo_Zunino_wo_dt", "Schott_Massing_Burman_DAngelo_Zunino",
          "Schott_Massing_Burman_DAngelo_Zunino_wo_dt", "Burman",
          "Taylor_Hughes_Zarins_Whiting_Jansen_Codina_scaling", "tau_not_defined"),
      tuple<std::string>("definition of burman_fernandez_hansbo",
          "definition of burman_fernandez_hansbo for stationary problems",
          "definition of braack_burman_john_lube",
          "definition of braack_burman_john_lube without explicit inclusion of divergence jump",
          "definition of tau_franca_barrenechea_valentin_wall",
          "definition of EOS_tau_burman_fernandez",
          "definition of EOS_tau_burman_hansbo_dangelo_zunino",
          "definition of EOS_tau_burman_hansbo_dangelo_zunino for stationary problems",
          "definition of EOS_tau_schott_massing_burman_dangelo_zunino",
          "definition of EOS_tau_schott_massing_burman_dangelo_zunino for stationary problems",
          "definition of EOS_tau_burman",
          "definition of EOS_tau related to residual-based stabilization", "no chosen definition"),
      tuple<int>(INPAR::FLUID::EOS_tau_burman_fernandez_hansbo,
          INPAR::FLUID::EOS_tau_burman_fernandez_hansbo_wo_dt,
          INPAR::FLUID::EOS_tau_braack_burman_john_lube,
          INPAR::FLUID::EOS_tau_braack_burman_john_lube_wo_divjump,
          INPAR::FLUID::EOS_tau_franca_barrenechea_valentin_wall,
          INPAR::FLUID::EOS_tau_burman_fernandez,
          INPAR::FLUID::EOS_tau_burman_hansbo_dangelo_zunino,
          INPAR::FLUID::EOS_tau_burman_hansbo_dangelo_zunino_wo_dt,
          INPAR::FLUID::EOS_tau_schott_massing_burman_dangelo_zunino,
          INPAR::FLUID::EOS_tau_schott_massing_burman_dangelo_zunino_wo_dt,
          INPAR::FLUID::EOS_tau_burman,
          INPAR::FLUID::EOS_tau_Taylor_Hughes_Zarins_Whiting_Jansen_Codina_scaling,
          INPAR::FLUID::EOS_tau_not_defined),
      &fdyn_edge_based_stab);

  //! this parameter selects how the element length of Edge-based stabilization is defined
  setStringToIntegralParameter<int>("EOS_H_DEFINITION", "EOS_he_max_diameter_to_opp_surf",
      "Definition of element length for edge-based stabilization",
      tuple<std::string>("EOS_he_max_diameter_to_opp_surf", "EOS_he_max_dist_to_opp_surf",
          "EOS_he_surf_with_max_diameter", "EOS_hk_max_diameter", "EOS_he_surf_diameter",
          "EOS_he_vol_eq_diameter"),
      tuple<std::string>("take the maximal (nsd-1)D diameter of faces that connect the internal "
                         "face to its opposite faces",
          "take the maximal 1D distance along 1D edge to opposite surface for both parent elements",
          "take the maximal (nsd-1)D face diameter of all faces for both parent elements",
          "maximal nD diameter of the neighboring elements",
          "maximal (n-1)D diameter of the internal face/edge",
          "take the maximal volume eqivalent diameter of adjecent elements"),
      tuple<int>(EOS_he_max_diameter_to_opp_surf, EOS_he_max_dist_to_opp_surf,
          EOS_he_surf_with_max_diameter, EOS_hk_max_diameter, EOS_he_surf_diameter,
          EOS_he_vol_eq_diameter),
      &fdyn_edge_based_stab);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fdyn_porostab = fdyn.sublist("POROUS-FLOW STABILIZATION", false, "");

  BoolParameter("STAB_BIOT", "No", "Flag to (de)activate BIOT stabilization.", &fdyn_porostab);
  DoubleParameter("STAB_BIOT_SCALING", 1.0,
      "Scaling factor for stabilization parameter for biot stabilization of porous flow.",
      &fdyn_porostab);

  // this parameter defines various stabilized methods
  setStringToIntegralParameter<int>("STABTYPE", "residual_based",
      "Apply (un)stabilized fluid formulation",
      tuple<std::string>("no_stabilization", "residual_based", "edge_based"),
      tuple<std::string>("Do not use any stabilization -> inf-sup stable elements required!",
          "Use a residual-based stabilization or, more generally, a stabilization \nbased on the "
          "concept of the residual-based variational multiscale method...\nExpecting additional "
          "input",
          "Use an edge-based stabilization, especially for XFEM"),
      tuple<int>(stabtype_nostab, stabtype_residualbased, stabtype_edgebased), &fdyn_porostab);

  BoolParameter("INCONSISTENT", "No",
      "residual based without second derivatives (i.e. only consistent for tau->0, but faster)",
      &fdyn_porostab);

  BoolParameter("Reconstruct_Sec_Der", "No",
      "residual computed with a reconstruction of the second derivatives via projection or "
      "superconvergent patch recovery",
      &fdyn_porostab);

  // the following parameters are necessary only if a residual based stabilized method is applied
  setStringToIntegralParameter<int>("TDS", "quasistatic",
      "Flag to allow time dependency of subscales for residual-based stabilization.",
      tuple<std::string>("quasistatic", "time_dependent"),
      tuple<std::string>("Use a quasi-static residual-based stabilization (standard case)",
          "Residual-based stabilization including time evolution equations for subscales"),
      tuple<int>(subscales_quasistatic, subscales_time_dependent), &fdyn_porostab);

  setStringToIntegralParameter<int>("TRANSIENT", "no_transient",
      "Specify how to treat the transient term.",
      tuple<std::string>("no_transient", "yes_transient", "transient_complete"),
      tuple<std::string>(
          "Do not use transient term (currently only opportunity for quasistatic stabilization)",
          "Use transient term (recommended for time dependent subscales)",
          "Use transient term including a linearisation of 1/tau"),
      tuple<int>(inertia_stab_drop, inertia_stab_keep, inertia_stab_keep_complete), &fdyn_porostab);

  BoolParameter("PSPG", "Yes", "Flag to (de)activate PSPG stabilization.", &fdyn_porostab);
  BoolParameter("SUPG", "Yes", "Flag to (de)activate SUPG stabilization.", &fdyn_porostab);
  BoolParameter("GRAD_DIV", "Yes", "Flag to (de)activate grad-div term.", &fdyn_porostab);

  setStringToIntegralParameter<int>("VSTAB", "no_vstab",
      "Flag to (de)activate viscous term in residual-based stabilization.",
      tuple<std::string>(
          "no_vstab", "vstab_gls", "vstab_gls_rhs", "vstab_usfem", "vstab_usfem_rhs"),
      tuple<std::string>("No viscous term in stabilization", "Viscous stabilization of GLS type",
          "Viscous stabilization of GLS type, included only on the right hand side",
          "Viscous stabilization of USFEM type",
          "Viscous stabilization of USFEM type, included only on the right hand side"),
      tuple<int>(viscous_stab_none, viscous_stab_gls, viscous_stab_gls_only_rhs, viscous_stab_usfem,
          viscous_stab_usfem_only_rhs),
      &fdyn_porostab);

  setStringToIntegralParameter<int>("RSTAB", "no_rstab",
      "Flag to (de)activate reactive term in residual-based stabilization.",
      tuple<std::string>("no_rstab", "rstab_gls", "rstab_usfem"),
      tuple<std::string>("no reactive term in stabilization", "reactive stabilization of GLS type",
          "reactive stabilization of USFEM type"),
      tuple<int>(reactive_stab_none, reactive_stab_gls, reactive_stab_usfem), &fdyn_porostab);

  setStringToIntegralParameter<int>("CROSS-STRESS", "no_cross",
      "Flag to (de)activate cross-stress term -> residual-based VMM.",
      tuple<std::string>("no_cross", "yes_cross", "cross_rhs"
          //"cross_complete"
          ),
      tuple<std::string>("No cross-stress term",
          "Include the cross-stress term with a linearization of the convective part",
          "Include cross-stress term, but only explicitly on right hand side"
          //""
          ),
      tuple<int>(cross_stress_stab_none, cross_stress_stab, cross_stress_stab_only_rhs),
      &fdyn_porostab);

  setStringToIntegralParameter<int>("REYNOLDS-STRESS", "no_reynolds",
      "Flag to (de)activate Reynolds-stress term -> residual-based VMM.",
      tuple<std::string>("no_reynolds", "yes_reynolds", "reynolds_rhs"
          //"reynolds_complete"
          ),
      tuple<std::string>(
          "No Reynolds-stress term", "Include Reynolds-stress term with linearisation",
          "Include Reynolds-stress term explicitly on right hand side"
          //""
          ),
      tuple<int>(reynolds_stress_stab_none, reynolds_stress_stab, reynolds_stress_stab_only_rhs),
      &fdyn_porostab);

  // this parameter selects the tau definition applied
  setStringToIntegralParameter<int>("DEFINITION_TAU", "Franca_Barrenechea_Valentin_Frey_Wall",
      "Definition of tau_M and Tau_C",
      tuple<std::string>("Taylor_Hughes_Zarins", "Taylor_Hughes_Zarins_wo_dt",
          "Taylor_Hughes_Zarins_Whiting_Jansen", "Taylor_Hughes_Zarins_Whiting_Jansen_wo_dt",
          "Taylor_Hughes_Zarins_scaled", "Taylor_Hughes_Zarins_scaled_wo_dt",
          "Franca_Barrenechea_Valentin_Frey_Wall", "Franca_Barrenechea_Valentin_Frey_Wall_wo_dt",
          "Shakib_Hughes_Codina", "Shakib_Hughes_Codina_wo_dt", "Codina", "Codina_wo_dt",
          "Franca_Madureira_Valentin_Badia_Codina", "Franca_Madureira_Valentin_Badia_Codina_wo_dt"),
      tuple<int>(tau_taylor_hughes_zarins, tau_taylor_hughes_zarins_wo_dt,
          tau_taylor_hughes_zarins_whiting_jansen, tau_taylor_hughes_zarins_whiting_jansen_wo_dt,
          tau_taylor_hughes_zarins_scaled, tau_taylor_hughes_zarins_scaled_wo_dt,
          tau_franca_barrenechea_valentin_frey_wall,
          tau_franca_barrenechea_valentin_frey_wall_wo_dt, tau_shakib_hughes_codina,
          tau_shakib_hughes_codina_wo_dt, tau_codina, tau_codina_wo_dt,
          tau_franca_madureira_valentin_badia_codina,
          tau_franca_madureira_valentin_badia_codina_wo_dt),
      &fdyn_porostab);

  // this parameter selects the characteristic element length for tau_Mu for all
  // stabilization parameter definitions requiring such a length
  setStringToIntegralParameter<int>("CHARELELENGTH_U", "streamlength",
      "Characteristic element length for tau_Mu",
      tuple<std::string>("streamlength", "volume_equivalent_diameter", "root_of_volume"),
      tuple<int>(streamlength_u, volume_equivalent_diameter_u, root_of_volume_u), &fdyn_porostab);

  // this parameter selects the characteristic element length for tau_Mp and tau_C for
  // all stabilization parameter definitions requiring such a length
  setStringToIntegralParameter<int>("CHARELELENGTH_PC", "volume_equivalent_diameter",
      "Characteristic element length for tau_Mp/tau_C",
      tuple<std::string>("streamlength", "volume_equivalent_diameter", "root_of_volume"),
      tuple<int>(streamlength_pc, volume_equivalent_diameter_pc, root_of_volume_pc),
      &fdyn_porostab);

  // this parameter selects the location where tau is evaluated
  setStringToIntegralParameter<int>("EVALUATION_TAU", "element_center",
      "Location where tau is evaluated", tuple<std::string>("element_center", "integration_point"),
      tuple<std::string>("evaluate tau at element center", "evaluate tau at integration point"),
      tuple<int>(0, 1), &fdyn_porostab);

  // this parameter selects the location where the material law is evaluated
  // (does not fit here very well, but parameter transfer is easier)
  setStringToIntegralParameter<int>("EVALUATION_MAT", "element_center",
      "Location where material law is evaluated",
      tuple<std::string>("element_center", "integration_point"),
      tuple<std::string>(
          "evaluate material law at element center", "evaluate material law at integration point"),
      tuple<int>(0, 1), &fdyn_porostab);

  // these parameters active additional terms in loma continuity equation
  // which might be identified as SUPG-/cross- and Reynolds-stress term
  BoolParameter("LOMA_CONTI_SUPG", "No",
      "Flag to (de)activate SUPG stabilization in loma continuity equation.", &fdyn_porostab);

  setStringToIntegralParameter<int>("LOMA_CONTI_CROSS_STRESS", "no_cross",
      "Flag to (de)activate cross-stress term loma continuity equation-> residual-based VMM.",
      tuple<std::string>("no_cross", "yes_cross", "cross_rhs"
          //"cross_complete"
          ),
      tuple<std::string>("No cross-stress term",
          "Include the cross-stress term with a linearization of the convective part",
          "Include cross-stress term, but only explicitly on right hand side"
          //""
          ),
      tuple<int>(cross_stress_stab_none, cross_stress_stab, cross_stress_stab_only_rhs),
      &fdyn_porostab);

  setStringToIntegralParameter<int>("LOMA_CONTI_REYNOLDS_STRESS", "no_reynolds",
      "Flag to (de)activate Reynolds-stress term loma continuity equation-> residual-based VMM.",
      tuple<std::string>("no_reynolds", "yes_reynolds", "reynolds_rhs"),
      tuple<std::string>(
          "No Reynolds-stress term", "Include Reynolds-stress term with linearisation",
          "Include Reynolds-stress term explicitly on right hand side"
          //""
          ),
      tuple<int>(reynolds_stress_stab_none, reynolds_stress_stab, reynolds_stress_stab_only_rhs),
      &fdyn_porostab);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fdyn_turbu = fdyn.sublist("TURBULENCE MODEL", false, "");

  //----------------------------------------------------------------------
  // modeling strategies
  //----------------------------------------------------------------------

  setStringToIntegralParameter<int>("TURBULENCE_APPROACH", "DNS_OR_RESVMM_LES",
      "There are several options to deal with turbulent flows.",
      tuple<std::string>("DNS_OR_RESVMM_LES", "CLASSICAL_LES"),
      tuple<std::string>("Try to solve flow as an underresolved DNS.\nMind that your stabilisation "
                         "already acts as a kind of turbulence model!",
          "Perform a classical Large Eddy Simulation adding \naddititional turbulent viscosity. "
          "This may be based on various physical models."),
      tuple<int>(0, 1), &fdyn_turbu);

  setStringToIntegralParameter<int>("PHYSICAL_MODEL", "no_model",
      "Classical LES approaches require an additional model for\nthe turbulent viscosity.",
      tuple<std::string>("no_model", "Smagorinsky", "Smagorinsky_with_van_Driest_damping",
          "Dynamic_Smagorinsky", "Multifractal_Subgrid_Scales", "Vreman", "Dynamic_Vreman"),
      tuple<std::string>("If classical LES is our turbulence approach, this is a contradiction and "
                         "should cause a dserror.",
          "Classical constant coefficient Smagorinsky model. Be careful if you \nhave a wall "
          "bounded flow domain!",
          "Use an exponential damping function for the turbulent viscosity \nclose to the wall. "
          "This is only implemented for a channel geometry of \nheight 2 in y direction. The "
          "viscous lengthscale l_tau is \nrequired as additional input.",
          "The solution is filtered and by comparison of the filtered \nvelocity field with the "
          "real solution, the Smagorinsky constant is \nestimated in each step --- mind that this "
          "procedure includes \nan averaging in the xz plane, hence this implementation will only "
          "work \nfor a channel flow.",
          "Multifractal Subgrid-Scale Modeling based on the work of burton",
          "Vremans constant model", "Dynamic Vreman model according to You and Moin (2007)"),
      tuple<int>(0, 1, 2, 3, 4, 5, 6), &fdyn_turbu);

  setStringToIntegralParameter<int>("FSSUGRVISC", "No", "fine-scale subgrid viscosity",
      tuple<std::string>("No", "Smagorinsky_all", "Smagorinsky_small"),
      tuple<int>(no_fssgv, smagorinsky_all, smagorinsky_small), &fdyn_turbu);

  //----------------------------------------------------------------------
  // turbulence specific output and statistics
  //----------------------------------------------------------------------

  IntParameter(
      "SAMPLING_START", 10000000, "Time step after when sampling shall be started", &fdyn_turbu);
  IntParameter("SAMPLING_STOP", 1, "Time step when sampling shall be stopped", &fdyn_turbu);
  IntParameter("DUMPING_PERIOD", 1,
      "Period of time steps after which statistical data shall be dumped", &fdyn_turbu);

  BoolParameter("SUBGRID_DISSIPATION", "No",
      "Flag to (de)activate estimation of subgrid-scale dissipation (only for seclected flows).",
      &fdyn_turbu);

  BoolParameter("OUTMEAN", "No", "Flag to (de)activate averaged paraview output", &fdyn_turbu);

  BoolParameter("TURBMODEL_LS", "Yes",
      "Flag to (de)activate turbulence model in level-set equation", &fdyn_turbu);

  //----------------------------------------------------------------------
  // turbulent flow problem and general characteristics
  //----------------------------------------------------------------------

  {
    // a standard Teuchos::tuple can have at maximum 10 entries! We have to circumvent this here.
    // Otherwise BACI DEBUG version will crash during runtime!
    Teuchos::Tuple<std::string, 22> name;
    Teuchos::Tuple<int, 22> label;
    name[0] = "no";
    label[0] = 0;
    name[1] = "time_averaging";
    label[1] = 1;
    name[2] = "channel_flow_of_height_2";
    label[2] = 2;
    name[3] = "lid_driven_cavity";
    label[3] = 3;
    name[4] = "backward_facing_step";
    label[4] = 4;
    name[5] = "square_cylinder";
    label[5] = 5;
    name[6] = "square_cylinder_nurbs";
    label[6] = 6;
    name[7] = "rotating_circular_cylinder_nurbs";
    label[7] = 7;
    name[8] = "rotating_circular_cylinder_nurbs_scatra";
    label[8] = 8;
    name[9] = "loma_channel_flow_of_height_2";
    label[9] = 9;
    name[10] = "loma_lid_driven_cavity";
    label[10] = 10;
    name[11] = "loma_backward_facing_step";
    label[11] = 11;
    name[12] = "combust_oracles";
    label[12] = 12;
    name[13] = "bubbly_channel_flow";
    label[13] = 13;
    name[14] = "scatra_channel_flow_of_height_2";
    label[14] = 14;
    name[15] = "decaying_homogeneous_isotropic_turbulence";
    label[15] = 15;
    name[16] = "forced_homogeneous_isotropic_turbulence";
    label[16] = 16;
    name[17] = "scatra_forced_homogeneous_isotropic_turbulence";
    label[17] = 17;
    name[18] = "taylor_green_vortex";
    label[18] = 18;
    name[19] = "periodic_hill";
    label[19] = 18;
    name[20] = "blood_fda_flow";
    label[20] = 20;
    name[21] = "backward_facing_step2";
    label[21] = 21;

    Teuchos::Tuple<std::string, 22> description;
    description[0] =
        "The flow is not further specified, so spatial averaging \nand hence the standard sampling "
        "procedure is not possible";
    description[1] =
        "The flow is not further specified, but time averaging of velocity and pressure field is "
        "performed";
    description[2] =
        "For this flow, all statistical data could be averaged in \nthe homogenous planes --- it "
        "is essentially a statistically one dimensional flow.";
    description[3] =
        "For this flow, all statistical data are evaluated on the center lines of the xy-midplane, "
        "averaged only over time.";
    description[4] =
        "For this flow, statistical data are evaluated on various lines, averaged over time and z.";
    description[5] =
        "For this flow, statistical data are evaluated on various lines of the xy-midplane, "
        "averaged only over time.";
    description[6] =
        "For this flow, statistical data are evaluated on various lines of the xy-midplane, "
        "averaged over time and eventually in one hom.direction.";
    description[7] =
        "For this flow, statistical data is computed in concentric surfaces and averaged. in time "
        "and in one hom. direction";
    description[8] =
        "For this flow with mass transport, statistical data is computed in concentric surfaces "
        "and averaged. in time and in one hom. direction";
    description[9] =
        "For this low-Mach-number flow, all statistical data could be averaged in \nthe homogenous "
        "planes --- it is essentially a statistically one dimensional flow.";
    description[10] =
        "For this low-Mach-number flow, all statistical data are evaluated on the center lines of "
        "the xy-midplane, averaged only over time.";
    description[11] =
        "For this low-Mach-number flow, statistical data are evaluated on various lines, averaged "
        "over time and z.";
    description[12] = "ORACLES test rig for turbulent premixed combustion.";
    description[13] =
        "Turbulent two-phase flow: bubbly channel flow, statistical data are averaged in "
        "homogeneous planse and over time.";
    description[14] =
        "For this flow, all statistical data could be averaged in \nthe homogenous planes --- it "
        "is essentially a statistically one dimensional flow.";
    description[15] =
        "For this flow, all statistical data could be averaged in \nthe in all homogeneous "
        "directions  --- it is essentially a statistically zero dimensional flow.";
    description[16] =
        "For this flow, all statistical data could be averaged in \nthe in all homogeneous "
        "directions  --- it is essentially a statistically zero dimensional flow.";
    description[17] =
        "For this flow, all statistical data could be averaged in \nthe in all homogeneous "
        "directions  --- it is essentially a statistically zero dimensional flow.";
    description[18] =
        "For this flow, dissipation rate could be averaged in \nthe in all homogeneous directions  "
        "--- it is essentially a statistically zero dimensional flow.";
    description[19] =
        "For this flow, statistical data is evaluated on various lines, averaged over time and z.";
    description[20] = "For this flow, statistical data is evaluated on various planes.";
    description[21] = "For this flow, statistical data is evaluated on various planes.";

    setStringToIntegralParameter<int>("CANONICAL_FLOW", "no",
        "Sampling is different for different canonical flows \n--- so specify what kind of flow "
        "you've got",
        name, description, label, &fdyn_turbu);
  }

  setStringToIntegralParameter<int>("HOMDIR", "not_specified",
      "Specify the homogenous direction(s) of a flow",
      tuple<std::string>("not_specified", "x", "y", "z", "xy", "xz", "yz", "xyz"),
      tuple<std::string>(
          "no homogeneous directions available, averaging is restricted to time averaging",
          "average along x-direction", "average along y-direction", "average along z-direction",
          "Wall normal direction is z, average in x and y direction",
          "Wall normal direction is y, average in x and z direction (standard case)",
          "Wall normal direction is x, average in y and z direction",
          "averageing in all directions"),
      tuple<int>(0, 1, 2, 3, 4, 5, 6, 7), &fdyn_turbu);

  //---------------------------------------
  // further problem-specific parameters

  // CHANNEL FLOW
  //--------------

  DoubleParameter("CHAN_AMPL_INIT_DIST", 0.1,
      "Max. amplitude of the random disturbance in percent of the initial value in mean flow "
      "direction.",
      &fdyn_turbu);

  setStringToIntegralParameter<int>("FORCING_TYPE",
      "linear_compensation_from_intermediate_spectrum", "forcing strategy",
      tuple<std::string>("linear_compensation_from_intermediate_spectrum", "fixed_power_input"),
      tuple<int>(linear_compensation_from_intermediate_spectrum, fixed_power_input), &fdyn_turbu);

  IntParameter(
      "CHA_NUMSUBDIVISIONS", 5, "Number of homogenious sampling planes in element", &fdyn_turbu);

  // HIT
  //--------------

  IntParameter("FORCING_TIME_STEPS", 0,
      "Number of time steps during which forcing is applied. Decaying homogeneous isotropic "
      "turbulence only.",
      &fdyn_turbu);

  DoubleParameter("THRESHOLD_WAVENUMBER", 0.0,
      "Forcing is only applied to wave numbers lower or equal than the given threshold wave "
      "number.",
      &fdyn_turbu);

  DoubleParameter("POWER_INPUT", 0.0, "power of forcing", &fdyn_turbu);

  setStringToIntegralParameter<int>("SCALAR_FORCING", "no", "Define forcing for scalar field.",
      tuple<std::string>("no", "isotropic", "mean_scalar_gradient"),
      tuple<std::string>("Do not force the scalar field",
          "Force scalar field isotropically such as the fluid field.",
          "Force scalar field by imposed mean-scalar gradient."),
      tuple<int>(0, 1, 2), &fdyn_turbu);

  DoubleParameter("MEAN_SCALAR_GRADIENT", 0.0,
      "Value of imposed mean-scalar gradient to force scalar field.", &fdyn_turbu);

  // filtering with xfem
  //--------------

  BoolParameter("EXCLUDE_XFEM", "No",
      "Flag to (de)activate XFEM dofs in calculation of fine-scale velocity.", &fdyn_turbu);

  /*----------------------------------------------------------------------*/
  // sublist with additional input parameters for Smagorinsky model
  Teuchos::ParameterList& fdyn_turbsgv = fdyn.sublist("SUBGRID VISCOSITY", false, "");

  DoubleParameter("C_SMAGORINSKY", 0.0,
      "Constant for the Smagorinsky model. Something between 0.1 to 0.24. Vreman constant if the "
      "constant vreman model is applied (something between 0.07 and 0.01).",
      &fdyn_turbsgv);
  DoubleParameter("C_YOSHIZAWA", -1.0,
      "Constant for the compressible Smagorinsky model: isotropic part of subgrid-stress tensor. "
      "About 0.09 or 0.0066. Ci will not be squared!",
      &fdyn_turbsgv);
  BoolParameter("C_SMAGORINSKY_AVERAGED", "No",
      "Flag to (de)activate averaged Smagorinksy constant", &fdyn_turbsgv);
  BoolParameter("C_INCLUDE_CI", "No", "Flag to (de)inclusion of Yoshizawa model", &fdyn_turbsgv);
  // remark: following Moin et al. 1991, the extension of the dynamic Smagorinsky model to
  // compressibel flow
  //        also contains a model for the isotropic part of the subgrid-stress tensor according to
  //        Yoshizawa 1989 although used in literature for turbulent variable-density flow at low
  //        Mach number, this additional term seems to destabilize the simulation when the flow is
  //        only weakly compressible therefore C_INCLUDE_CI allows to exclude this term if
  //        C_SMAGORINSKY_AVERAGED == true
  //           if C_INCLUDE_CI==true and C_YOSHIZAWA>=0.0 then the given value C_YOSHIZAWA is used
  //           if C_INCLUDE_CI==true and C_YOSHIZAWA<0.0 then C_YOSHIZAWA is determined dynamically
  //        else all values are taken from input

  DoubleParameter("CHANNEL_L_TAU", 0.0,
      "Used for normalisation of the wall normal distance in the Van \nDriest Damping function. "
      "May be taken from the output of \nthe apply_mesh_stretching.pl preprocessing script.",
      &fdyn_turbsgv);

  DoubleParameter("C_TURBPRANDTL", 1.0,
      "(Constant) turbulent Prandtl number for the Smagorinsky model in scalar transport.",
      &fdyn_turbsgv);

  setStringToIntegralParameter<int>("FILTER_WIDTH", "CubeRootVol",
      "The Vreman model requires a filter width.",
      tuple<std::string>("CubeRootVol", "Direction_dependent", "Minimum_length"),
      tuple<int>(cuberootvol, dir_dep, min_len), &fdyn_turbsgv);


  /*----------------------------------------------------------------------*/
  // sublist with additional input parameters for Smagorinsky model
  Teuchos::ParameterList& fdyn_wallmodel = fdyn.sublist("WALL MODEL", false, "");

  BoolParameter("X_WALL", "No", "Flag to switch on the xwall model", &fdyn_wallmodel);

  setStringToIntegralParameter<int>("Tauw_Type", "constant",
      "Methods for calculating/updating the wall shear stress necessary for Spalding's law.",
      tuple<std::string>("constant", "between_steps"),
      tuple<std::string>(
          "Use the constant wall shear stress given in the input file for the whole simulation.",
          "Calculate wall shear stress in between time steps."),
      tuple<int>(0, 1), &fdyn_wallmodel);

  setStringToIntegralParameter<int>("Tauw_Calc_Type", "residual",
      "Methods for calculating the wall shear stress necessary for Spalding's law.",
      tuple<std::string>("residual", "gradient", "gradient_to_residual"),
      tuple<std::string>("Residual (force) devided by area.",
          "Gradient via shape functions and nodal values.", "First gradient, then residual."),
      tuple<int>(0, 1, 3), &fdyn_wallmodel);

  IntParameter("Switch_Step", -1, "Switch from gradient to residual based tauw.", &fdyn_wallmodel);

  setStringToIntegralParameter<int>("Projection", "No",
      "Flag to switch projection of the enriched dofs after updating tauw, alternatively with or "
      "without continuity constraint.",
      tuple<std::string>("No", "onlyl2projection", "l2projectionwithcontinuityconstraint"),
      tuple<std::string>("Switch off projection.", "Only l2 projection.",
          "L2 projection with continuity constraint."),
      tuple<int>(0, 1, 2), &fdyn_wallmodel);

  DoubleParameter("C_Tauw", 1.0, "Constant wall shear stress for Spalding's law, if applicable",
      &fdyn_wallmodel);

  DoubleParameter("Min_Tauw", 2.0e-9,
      "Minimum wall shear stress preventing system to become singular", &fdyn_wallmodel);

  DoubleParameter(
      "Inc_Tauw", 1.0, "Increment of Tauw of full step, between 0.0 and 1.0", &fdyn_wallmodel);

  setStringToIntegralParameter<int>("Blending_Type", "none",
      "Methods for blending the enrichment space.", tuple<std::string>("none", "ramp_function"),
      tuple<std::string>("No ramp function, does not converge!",
          "Enrichment is multiplied with linear ramp function resulting in zero enrichment at the "
          "interface"),
      tuple<int>(0, 1), &fdyn_wallmodel);

  IntParameter("GP_Wall_Normal", 3, "Gauss points in wall normal direction", &fdyn_wallmodel);
  IntParameter("GP_Wall_Normal_Off_Wall", 3,
      "Gauss points in wall normal direction, off-wall elements", &fdyn_wallmodel);
  IntParameter("GP_Wall_Parallel", 3, "Gauss points in wall parallel direction", &fdyn_wallmodel);

  BoolParameter("Treat_Tauw_on_Dirichlet_Inflow", "No",
      "Flag to treat residual on Dirichlet inflow nodes for calculation of wall shear stress",
      &fdyn_wallmodel);

  IntParameter("PROJECTION_SOLVER", -1, "Set solver number for l2-projection.", &fdyn_wallmodel);

  /*----------------------------------------------------------------------*/
  // sublist with additional input parameters for multifractal subgrid-scales
  Teuchos::ParameterList& fdyn_turbmfs = fdyn.sublist("MULTIFRACTAL SUBGRID SCALES", false, "");

  DoubleParameter("CSGS", 0.0, "Modelparameter of multifractal subgrid-scales.", &fdyn_turbmfs);

  setStringToIntegralParameter<int>("SCALE_SEPARATION", "no_scale_sep",
      "Specify the filter type for scale separation in LES",
      tuple<std::string>("no_scale_sep", "box_filter", "algebraic_multigrid_operator"),
      tuple<std::string>("no scale separation", "classical box filter",
          "scale separation by algebraic multigrid operator"),
      tuple<int>(0, 1, 2), &fdyn_turbmfs);

  IntParameter("ML_SOLVER", -1,
      "Set solver number for scale separation via level set transfer operators from plain "
      "aggregation.",
      &fdyn_turbmfs);

  BoolParameter("CALC_N", "No", "Flag to (de)activate calculation of N from the Reynolds number.",
      &fdyn_turbmfs);

  DoubleParameter("N", 1.0, "Set grid to viscous scale ratio.", &fdyn_turbmfs);

  setStringToIntegralParameter<int>("REF_LENGTH", "cube_edge",
      "Specify the reference length for Re-dependent N.",
      tuple<std::string>(
          "cube_edge", "sphere_diameter", "streamlength", "gradient_based", "metric_tensor"),
      tuple<std::string>("edge length of volume equivalent cube",
          "diameter of volume equivalent sphere", "streamlength taken from stabilization",
          "gradient based length taken from stabilization",
          "metric tensor taken from stabilization"),
      tuple<int>(0, 1, 2, 3, 4), &fdyn_turbmfs);

  setStringToIntegralParameter<int>("REF_VELOCITY", "strainrate",
      "Specify the reference velocity for Re-dependent N.",
      tuple<std::string>("strainrate", "resolved", "fine_scale"),
      tuple<std::string>("norm of strain rate", "resolved velocity", "fine-scale velocity"),
      tuple<int>(0, 1, 2), &fdyn_turbmfs);

  DoubleParameter("C_NU", 1.0,
      "Proportionality constant between Re and ratio viscous scale to element length.",
      &fdyn_turbmfs);

  BoolParameter("NEAR_WALL_LIMIT", "No", "Flag to (de)activate near-wall limit.", &fdyn_turbmfs);

  setStringToIntegralParameter<int>("EVALUATION_B", "element_center",
      "Location where B is evaluated", tuple<std::string>("element_center", "integration_point"),
      tuple<std::string>("evaluate B at element center", "evaluate B at integration point"),
      tuple<int>(0, 1), &fdyn_turbmfs);

  DoubleParameter(
      "BETA", 0.0, "Cross- and Reynolds-stress terms only on right-hand-side.", &fdyn_turbmfs);

  setStringToIntegralParameter<int>("CONVFORM", "convective", "form of convective term",
      tuple<std::string>("convective", "conservative"), tuple<int>(0, 1), &fdyn_turbmfs);

  DoubleParameter("CSGS_PHI", 0.0,
      "Modelparameter of multifractal subgrid-scales for scalar transport.", &fdyn_turbmfs);

  BoolParameter(
      "ADAPT_CSGS_PHI", "No", "Flag to (de)activate adaption of CsgsD to CsgsB.", &fdyn_turbmfs);

  BoolParameter("NEAR_WALL_LIMIT_CSGS_PHI", "No",
      "Flag to (de)activate near-wall limit for scalar field.", &fdyn_turbmfs);

  BoolParameter("CONSISTENT_FLUID_RESIDUAL", "No",
      "Flag to (de)activate the consistency term for residual-based stabilization.", &fdyn_turbmfs);

  DoubleParameter("C_DIFF", 1.0,
      "Proportionality constant between Re*Pr and ratio dissipative scale to element length. "
      "Usually equal cnu.",
      &fdyn_turbmfs);

  BoolParameter("SET_FINE_SCALE_VEL", "No",
      "Flag to set fine-scale velocity for parallel nightly tests.", &fdyn_turbmfs);

  // activate cross- and Reynolds-stress terms in loma continuity equation
  BoolParameter("LOMA_CONTI", "No",
      "Flag to (de)activate cross- and Reynolds-stress terms in loma continuity equation.",
      &fdyn_turbmfs);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& fdyn_turbinf = fdyn.sublist("TURBULENT INFLOW", false, "");

  BoolParameter("TURBULENTINFLOW", "No",
      "Flag to (de)activate potential separate turbulent inflow section", &fdyn_turbinf);

  setStringToIntegralParameter<int>("INITIALINFLOWFIELD", "zero_field",
      "Initial field for inflow section",
      tuple<std::string>("zero_field", "field_by_function", "disturbed_field_from_function"),
      tuple<int>(initfield_zero_field, initfield_field_by_function,
          initfield_disturbed_field_from_function),
      &fdyn_turbinf);

  IntParameter(
      "INFLOWFUNC", -1, "Function number for initial flow field in inflow section", &fdyn_turbinf);

  DoubleParameter("INFLOW_INIT_DIST", 0.1,
      "Max. amplitude of the random disturbance in percent of the initial value in mean flow "
      "direction.",
      &fdyn_turbinf);

  IntParameter("NUMINFLOWSTEP", 1, "Total number of time steps for development of turbulent flow",
      &fdyn_turbinf);

  setStringToIntegralParameter<int>("CANONICAL_INFLOW", "no",
      "Sampling is different for different canonical flows \n--- so specify what kind of flow "
      "you've got",
      tuple<std::string>("no", "time_averaging", "channel_flow_of_height_2",
          "loma_channel_flow_of_height_2", "scatra_channel_flow_of_height_2"),
      tuple<std::string>("The flow is not further specified, so spatial averaging \nand hence the "
                         "standard sampling procedure is not possible",
          "The flow is not further specified, but time averaging of velocity and pressure field is "
          "performed",
          "For this flow, all statistical data could be averaged in \nthe homogenous planes --- it "
          "is essentially a statistically one dimensional flow.",
          "For this low-Mach-number flow, all statistical data could be averaged in \nthe "
          "homogenous planes --- it is essentially a statistically one dimensional flow.",
          "For this flow, all statistical data could be averaged in \nthe homogenous planes --- it "
          "is essentially a statistically one dimensional flow."),
      tuple<int>(0, 1, 2, 3, 4), &fdyn_turbinf);

  DoubleParameter("INFLOW_CHA_SIDE", 0.0,
      "Most right side of inflow channel. Necessary to define sampling domain.", &fdyn_turbinf);

  setStringToIntegralParameter<int>("INFLOW_HOMDIR", "not_specified",
      "Specify the homogenous direction(s) of a flow",
      tuple<std::string>("not_specified", "x", "y", "z", "xy", "xz", "yz"),
      tuple<std::string>(
          "no homogeneous directions available, averaging is restricted to time averaging",
          "average along x-direction", "average along y-direction", "average along z-direction",
          "Wall normal direction is z, average in x and y direction",
          "Wall normal direction is y, average in x and z direction (standard case)",
          "Wall normal direction is x, average in y and z direction"),
      tuple<int>(0, 1, 2, 3, 4, 5, 6), &fdyn_turbinf);

  IntParameter("INFLOW_SAMPLING_START", 10000000, "Time step after when sampling shall be started",
      &fdyn_turbinf);
  IntParameter(
      "INFLOW_SAMPLING_STOP", 1, "Time step when sampling shall be stopped", &fdyn_turbinf);
  IntParameter("INFLOW_DUMPING_PERIOD", 1,
      "Period of time steps after which statistical data shall be dumped", &fdyn_turbinf);

  /*----------------------------------------------------------------------*/
  // sublist with additional input parameters for time adaptivity in fluid/ coupled problems
  Teuchos::ParameterList& fdyn_timintada = fdyn.sublist("TIMEADAPTIVITY", false, "");
  setStringToIntegralParameter<int>("ADAPTIVE_TIME_STEP_ESTIMATOR", "none",
      "Method used to determine adaptive time step size.",
      tuple<std::string>("none", "cfl_number", "only_print_cfl_number"),
      tuple<std::string>(
          "constant time step", "evaluated via CFL number", "CFL number evaluated and printed"
          //""
          ),
      tuple<int>(const_dt, cfl_number, only_print_cfl_number), &fdyn_timintada);

  DoubleParameter("CFL_NUMBER", -1.0, "CFL number for adaptive time step", &fdyn_timintada);
  IntParameter("FREEZE_ADAPTIVE_DT_AT", 1000000,
      "keep time step constant after this step, otherwise turbulence statistics sampling is not "
      "consistent",
      &fdyn_timintada);
  DoubleParameter(
      "ADAPTIVE_DT_INC", 0.8, "Increment of whole step for adaptive dt via CFL", &fdyn_timintada);


  /*----------------------------------------------------------------------*/
  // TODO: Is this used anywhere?
  Teuchos::ParameterList& flucthydro = list->sublist("FLUCTUATING HYDRODYNAMICS", false, "");
  DoubleParameter("TEMPERATURE", 300, "Temperature in K", &flucthydro);
  DoubleParameter("BOLTZMANNCONST", 1.380650424e-23, "Boltzmann constant", &flucthydro);
  setStringToIntegralParameter<int>("SEEDCONTROL", "No",
      "control seeding with given unsigned integer", yesnotuple, yesnovalue, &flucthydro);
  IntParameter("SEEDVARIABLE", 0, "seed variable", &flucthydro);
  IntParameter("SAMPLEPERIOD", 1, "sample period", &flucthydro);
}



void INPAR::LOMA::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& lomacontrol = list->sublist(
      "LOMA CONTROL", false, "control parameters for low-Mach-number flow problems\n");

  BoolParameter("MONOLITHIC", "no", "monolithic solver", &lomacontrol);
  IntParameter("NUMSTEP", 24, "Total number of time steps", &lomacontrol);
  DoubleParameter("TIMESTEP", 0.1, "Time increment dt", &lomacontrol);
  DoubleParameter("MAXTIME", 1000.0, "Total simulation time", &lomacontrol);
  IntParameter("ITEMAX", 10, "Maximum number of outer iterations", &lomacontrol);
  IntParameter("ITEMAX_BEFORE_SAMPLING", 1,
      "Maximum number of outer iterations before sampling (for turbulent flows only)",
      &lomacontrol);
  DoubleParameter("CONVTOL", 1e-6, "Tolerance for convergence check", &lomacontrol);
  IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &lomacontrol);
  IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &lomacontrol);
  setStringToIntegralParameter<int>("CONSTHERMPRESS", "Yes",
      "treatment of thermodynamic pressure in time",
      tuple<std::string>("No_energy", "No_mass", "Yes"), tuple<int>(0, 1, 2), &lomacontrol);
  BoolParameter("SGS_MATERIAL_UPDATE", "no", "update material by adding subgrid-scale scalar field",
      &lomacontrol);

  // number of linear solver used for LOMA solver
  IntParameter("LINEAR_SOLVER", -1, "number of linear solver used for LOMA problem", &lomacontrol);
}



void INPAR::FLUID::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // transfer boundary condition for turbulent inflow

  std::vector<Teuchos::RCP<ConditionComponent>> tbc_turb_inflow_components;

  tbc_turb_inflow_components.push_back(Teuchos::rcp(new SeparatorConditionComponent("ID")));
  tbc_turb_inflow_components.push_back(Teuchos::rcp(new IntConditionComponent("id", true)));
  tbc_turb_inflow_components.push_back(Teuchos::rcp(new StringConditionComponent("toggle", "master",
      Teuchos::tuple<std::string>("master", "slave"),
      Teuchos::tuple<std::string>("master", "slave"))));
  tbc_turb_inflow_components.push_back(Teuchos::rcp(new SeparatorConditionComponent("DIRECTION")));
  tbc_turb_inflow_components.push_back(
      Teuchos::rcp(new StringConditionComponent("transfer direction", "x",
          Teuchos::tuple<std::string>("x", "y", "z"), Teuchos::tuple<std::string>("x", "y", "z"))));
  tbc_turb_inflow_components.push_back(
      Teuchos::rcp(new IntVectorConditionComponent("curve", 1, true, true)));

  Teuchos::RCP<ConditionDefinition> tbc_turb_inflow = Teuchos::rcp(new ConditionDefinition(
      "DESIGN SURF TURBULENT INFLOW TRANSFER", "TransferTurbulentInflow", "TransferTurbulentInflow",
      DRT::Condition::TransferTurbulentInflow, true, DRT::Condition::Surface));

  // we attach all the components of this condition to this weak line DBC
  for (unsigned i = 0; i < tbc_turb_inflow_components.size(); ++i)
  {
    tbc_turb_inflow->AddComponent(tbc_turb_inflow_components[i]);
  }

  // and append it to the list of all conditions
  condlist.push_back(tbc_turb_inflow);

  /*--------------------------------------------------------------------*/
  // separate domain for turbulent inflow generation

  Teuchos::RCP<ConditionDefinition> turbulentinflowgeneration =
      Teuchos::rcp(new ConditionDefinition("FLUID TURBULENT INFLOW VOLUME",
          "TurbulentInflowSection", "TurbulentInflowSection",
          DRT::Condition::TurbulentInflowSection, true, DRT::Condition::Volume));

  condlist.push_back(turbulentinflowgeneration);


  /*--------------------------------------------------------------------*/
  // flow-dependent pressure conditions

  std::vector<Teuchos::RCP<ConditionComponent>> flowdeppressurecomponents;

  // flow-dependent pressure conditions can be imposed either based on
  // (out)flow rate or (out)flow volume (e.g., for air-cushion condition)
  flowdeppressurecomponents.push_back(
      Teuchos::rcp(new StringConditionComponent("type of flow dependence", "flow_rate",
          Teuchos::tuple<std::string>("flow_rate", "flow_volume", "fixed_pressure"),
          Teuchos::tuple<std::string>("flow_rate", "flow_volume", "fixed_pressure"))));

  // constant coefficient for (linear) flow-rate-based condition
  // and constant fixed pressure
  flowdeppressurecomponents.push_back(Teuchos::rcp(new RealConditionComponent("ConstCoeff")));

  // linear coefficient for (linear) flow-rate-based condition
  flowdeppressurecomponents.push_back(Teuchos::rcp(new RealConditionComponent("LinCoeff")));

  // initial (air-cushion) volume outside of boundary
  flowdeppressurecomponents.push_back(Teuchos::rcp(new RealConditionComponent("InitialVolume")));

  // reference pressure outside of boundary
  flowdeppressurecomponents.push_back(
      Teuchos::rcp(new RealConditionComponent("ReferencePressure")));

  // adiabatic exponent
  flowdeppressurecomponents.push_back(
      Teuchos::rcp(new RealConditionComponent("AdiabaticExponent")));

  // values for time curve
  flowdeppressurecomponents.push_back(Teuchos::rcp(new IntConditionComponent("curve", true, true)));


  Teuchos::RCP<ConditionDefinition> lineflowdeppressure = Teuchos::rcp(new ConditionDefinition(
      "DESIGN LINE FLOW-DEPENDENT PRESSURE CONDITIONS", "LineFlowDepPressure",
      "LineFlowDepPressure", DRT::Condition::LineFlowDepPressure, true, DRT::Condition::Line));

  Teuchos::RCP<ConditionDefinition> surfflowdeppressure =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE FLOW-DEPENDENT PRESSURE CONDITIONS",
          "SurfaceFlowDepPressure", "SurfaceFlowDepPressure",
          DRT::Condition::SurfaceFlowDepPressure, true, DRT::Condition::Surface));

  // we attach all the components of this condition to this weak line DBC
  for (unsigned i = 0; i < flowdeppressurecomponents.size(); ++i)
  {
    lineflowdeppressure->AddComponent(flowdeppressurecomponents[i]);
    surfflowdeppressure->AddComponent(flowdeppressurecomponents[i]);
  }

  // and append it to the list of all conditions
  condlist.push_back(lineflowdeppressure);
  condlist.push_back(surfflowdeppressure);

  /*--------------------------------------------------------------------*/
  // Slip Supplemental Curved Boundary conditions

  std::vector<Teuchos::RCP<ConditionComponent>> slipsuppcomponents;

  slipsuppcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("USEUPDATEDNODEPOS")));
  slipsuppcomponents.push_back(
      Teuchos::rcp(new IntVectorConditionComponent("useupdatednodepos", 1)));

  Teuchos::RCP<ConditionDefinition> lineslipsupp = Teuchos::rcp(new ConditionDefinition(
      "DESIGN LINE SLIP SUPPLEMENTAL CURVED BOUNDARY CONDITIONS", "LineSlipSupp", "LineSlipSupp",
      DRT::Condition::LineSlipSupp, true, DRT::Condition::Line));

  Teuchos::RCP<ConditionDefinition> surfslipsupp = Teuchos::rcp(new ConditionDefinition(
      "DESIGN SURFACE SLIP SUPPLEMENTAL CURVED BOUNDARY CONDITIONS", "SurfaceSlipSupp",
      "SurfaceSlipSupp", DRT::Condition::SurfaceSlipSupp, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < slipsuppcomponents.size(); ++i)
  {
    lineslipsupp->AddComponent(slipsuppcomponents[i]);
    surfslipsupp->AddComponent(slipsuppcomponents[i]);
  }

  condlist.push_back(lineslipsupp);
  condlist.push_back(surfslipsupp);

  /*--------------------------------------------------------------------*/
  // Navier-slip boundary conditions

  std::vector<Teuchos::RCP<ConditionComponent>> navierslipcomponents;

  navierslipcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("SLIPCOEFFICIENT")));
  navierslipcomponents.push_back(Teuchos::rcp(new RealConditionComponent("slipcoefficient")));

  Teuchos::RCP<ConditionDefinition> linenavierslip = Teuchos::rcp(
      new ConditionDefinition("DESIGN LINE NAVIER-SLIP BOUNDARY CONDITIONS", "LineNavierSlip",
          "LineNavierSlip", DRT::Condition::LineNavierSlip, true, DRT::Condition::Line));

  Teuchos::RCP<ConditionDefinition> surfnavierslip = Teuchos::rcp(
      new ConditionDefinition("DESIGN SURF NAVIER-SLIP BOUNDARY CONDITIONS", "SurfNavierSlip",
          "SurfNavierSlip", DRT::Condition::SurfNavierSlip, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < navierslipcomponents.size(); ++i)
  {
    linenavierslip->AddComponent(navierslipcomponents[i]);
    surfnavierslip->AddComponent(navierslipcomponents[i]);
  }

  condlist.push_back(linenavierslip);
  condlist.push_back(surfnavierslip);

  /*--------------------------------------------------------------------*/
  // consistent outflow bcs for conservative element formulations

  Teuchos::RCP<ConditionDefinition> surfconsistentoutflowconsistency =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE CONSERVATIVE OUTFLOW CONSISTENCY",
          "SurfaceConservativeOutflowConsistency", "SurfaceConservativeOutflowConsistency",
          DRT::Condition::SurfaceConservativeOutflowConsistency, true, DRT::Condition::Surface));

  condlist.push_back(surfconsistentoutflowconsistency);

  /*--------------------------------------------------------------------*/
  // Neumann inflow for FLUID

  Teuchos::RCP<ConditionDefinition> linefluidneumanninflow = Teuchos::rcp(new ConditionDefinition(
      "FLUID NEUMANN INFLOW LINE CONDITIONS", "FluidNeumannInflow", "Line Fluid Neumann Inflow",
      DRT::Condition::FluidNeumannInflow, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surffluidneumanninflow = Teuchos::rcp(new ConditionDefinition(
      "FLUID NEUMANN INFLOW SURF CONDITIONS", "FluidNeumannInflow", "Surface Fluid Neumann Inflow",
      DRT::Condition::FluidNeumannInflow, true, DRT::Condition::Surface));

  condlist.push_back(linefluidneumanninflow);
  condlist.push_back(surffluidneumanninflow);

  /*--------------------------------------------------------------------*/
  // mixed/hybrid Dirichlet conditions

  std::vector<Teuchos::RCP<ConditionComponent>> mixhybDirichletcomponents;

  // we provide a vector of 3 values for velocities
  mixhybDirichletcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val", 3)));

  // and optional spatial functions
  mixhybDirichletcomponents.push_back(
      Teuchos::rcp(new IntVectorConditionComponent("funct", 3, false, false, true)));

  // characteristic velocity
  mixhybDirichletcomponents.push_back(Teuchos::rcp(new RealConditionComponent("u_C")));

  // the penalty parameter could be computed dynamically (using Spaldings
  // law of the wall) or using a fixed value (1)
  mixhybDirichletcomponents.push_back(
      Teuchos::rcp(new StringConditionComponent("Definition of penalty parameter", "constant",
          Teuchos::tuple<std::string>("constant", "Spalding"),
          Teuchos::tuple<std::string>("constant", "Spalding"))));

  // scaling factor for penalty parameter tauB
  mixhybDirichletcomponents.push_back(Teuchos::rcp(new RealConditionComponent("hB_divided_by")));

  // if Spaldings law is used, this defines the way how the traction at y is computed from utau
  mixhybDirichletcomponents.push_back(Teuchos::rcp(new StringConditionComponent("utau_computation",
      "at_wall", Teuchos::tuple<std::string>("at_wall", "viscous_tangent"),
      Teuchos::tuple<std::string>("at_wall", "viscous_tangent"))));


  Teuchos::RCP<ConditionDefinition> linemixhybDirichlet = Teuchos::rcp(new ConditionDefinition(
      "DESIGN LINE MIXED/HYBRID DIRICHLET CONDITIONS", "LineMixHybDirichlet", "LineMixHybDirichlet",
      DRT::Condition::LineMixHybDirichlet, true, DRT::Condition::Line));


  Teuchos::RCP<ConditionDefinition> surfmixhybDirichlet =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE MIXED/HYBRID DIRICHLET CONDITIONS",
          "SurfaceMixHybDirichlet", "SurfaceMixHybDirichlet",
          DRT::Condition::SurfaceMixHybDirichlet, true, DRT::Condition::Surface));

  // we attach all the components of this condition to this condition
  for (unsigned i = 0; i < mixhybDirichletcomponents.size(); ++i)
  {
    linemixhybDirichlet->AddComponent(mixhybDirichletcomponents[i]);
    surfmixhybDirichlet->AddComponent(mixhybDirichletcomponents[i]);
  }

  // and append it to the list of all conditions
  condlist.push_back(linemixhybDirichlet);
  condlist.push_back(surfmixhybDirichlet);

  /*--------------------------------------------------------------------*/
  // surface tension

  Teuchos::RCP<ConditionDefinition> surftension = Teuchos::rcp(new ConditionDefinition(
      "SURFACE TENSION CONDITIONS", "SurfaceStress", "Surface Stress (ideal water)",
      DRT::Condition::SurfaceTension, true, DRT::Condition::Surface));

  surftension->AddComponent(Teuchos::rcp(new IntConditionComponent("curve", true, true)));
  AddNamedReal(surftension, "gamma");

  condlist.push_back(surftension);

  /*--------------------------------------------------------------------*/
  // FREESURF

  std::vector<Teuchos::RCP<ConditionComponent>> freesurfcomponents;

  freesurfcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("FIELD")));
  freesurfcomponents.push_back(Teuchos::rcp(new StringConditionComponent("field", "fluid",
      Teuchos::tuple<std::string>("fluid", "ale"), Teuchos::tuple<std::string>("fluid", "ale"))));

  freesurfcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("COUPLING")));
  freesurfcomponents.push_back(Teuchos::rcp(new StringConditionComponent("coupling", "lagrange",
      Teuchos::tuple<std::string>("lagrange", "heightfunction", "sphereHeightFunction",
          "meantangentialvelocity", "meantangentialvelocityscaled"),
      Teuchos::tuple<std::string>("lagrange", "heightfunction", "sphereHeightFunction",
          "meantangentialvelocity", "meantangentialvelocityscaled"),
      true)));

  freesurfcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("VAL")));
  freesurfcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("val", 1)));

  freesurfcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("NODENORMALFUNCT")));
  freesurfcomponents.push_back(Teuchos::rcp(new IntConditionComponent("nodenormalfunct")));

  Teuchos::RCP<ConditionDefinition> linefreesurf = Teuchos::rcp(
      new ConditionDefinition("DESIGN FLUID FREE SURFACE LINE CONDITIONS", "FREESURFCoupling",
          "FREESURF Coupling", DRT::Condition::FREESURFCoupling, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surffreesurf = Teuchos::rcp(
      new ConditionDefinition("DESIGN FLUID FREE SURFACE SURF CONDITIONS", "FREESURFCoupling",
          "FREESURF Coupling", DRT::Condition::FREESURFCoupling, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < freesurfcomponents.size(); ++i)
  {
    linefreesurf->AddComponent(freesurfcomponents[i]);
    surffreesurf->AddComponent(freesurfcomponents[i]);
  }

  condlist.push_back(linefreesurf);
  condlist.push_back(surffreesurf);

  /*--------------------------------------------------------------------*/
  // fluid stress

  std::vector<Teuchos::RCP<ConditionComponent>> fluidstresscomponents;

  Teuchos::RCP<ConditionDefinition> linefluidstress =
      Teuchos::rcp(new ConditionDefinition("DESIGN FLUID STRESS CALC LINE CONDITIONS",
          "FluidStressCalc", "Line Fluid Stress Calculation", DRT::Condition::FluidStressCalc, true,
          DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surffluidstress =
      Teuchos::rcp(new ConditionDefinition("DESIGN FLUID STRESS CALC SURF CONDITIONS",
          "FluidStressCalc", "Surf Fluid Stress Calculation", DRT::Condition::FluidStressCalc, true,
          DRT::Condition::Surface));

  for (unsigned i = 0; i < fluidstresscomponents.size(); ++i)
  {
    linefluidstress->AddComponent(fluidstresscomponents[i]);
    surffluidstress->AddComponent(fluidstresscomponents[i]);
  }

  condlist.push_back(linefluidstress);
  condlist.push_back(surffluidstress);

  /*--------------------------------------------------------------------*/
  // lift & drag

  std::vector<Teuchos::RCP<ConditionComponent>> liftdragcomponents;

  liftdragcomponents.push_back(Teuchos::rcp(new IntConditionComponent("label")));
  liftdragcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("CENTER")));
  liftdragcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("centerCoord", 3)));
  // optional
  liftdragcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("AXIS", true)));
  liftdragcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("axis", 3, true)));

  Teuchos::RCP<ConditionDefinition> lineliftdrag =
      Teuchos::rcp(new ConditionDefinition("DESIGN FLUID LINE LIFT&DRAG", "LIFTDRAG",
          "Line LIFTDRAG", DRT::Condition::LineLIFTDRAG, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfliftdrag =
      Teuchos::rcp(new ConditionDefinition("DESIGN FLUID SURF LIFT&DRAG", "LIFTDRAG",
          "Surface LIFTDRAG", DRT::Condition::SurfLIFTDRAG, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < liftdragcomponents.size(); ++i)
  {
    lineliftdrag->AddComponent(liftdragcomponents[i]);
    surfliftdrag->AddComponent(liftdragcomponents[i]);
  }

  condlist.push_back(lineliftdrag);
  condlist.push_back(surfliftdrag);


  /*--------------------------------------------------------------------*/
  // flow rate through line

  std::vector<Teuchos::RCP<ConditionComponent>> lineflowratecomponents;
  lineflowratecomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  Teuchos::RCP<ConditionDefinition> lineflowrate =
      Teuchos::rcp(new ConditionDefinition("DESIGN FLOW RATE LINE CONDITIONS", "LineFlowRate",
          "Line Flow Rate", DRT::Condition::FlowRateThroughLine_2D, true, DRT::Condition::Line));

  for (unsigned i = 0; i < lineflowratecomponents.size(); ++i)
  {
    lineflowrate->AddComponent(lineflowratecomponents[i]);
  }
  condlist.push_back(lineflowrate);

  /*--------------------------------------------------------------------*/
  // flow rate through surface

  std::vector<Teuchos::RCP<ConditionComponent>> flowratecomponents;
  flowratecomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  Teuchos::RCP<ConditionDefinition> surfflowrate = Teuchos::rcp(new ConditionDefinition(
      "DESIGN FLOW RATE SURF CONDITIONS", "SurfFlowRate", "Surface Flow Rate",
      DRT::Condition::FlowRateThroughSurface_3D, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < flowratecomponents.size(); ++i)
  {
    surfflowrate->AddComponent(flowratecomponents[i]);
  }
  condlist.push_back(surfflowrate);

  /*--------------------------------------------------------------------*/
  // impuls rate through surface

  std::vector<Teuchos::RCP<ConditionComponent>> impulsratecomponents;
  impulsratecomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  Teuchos::RCP<ConditionDefinition> surfimpulsrate = Teuchos::rcp(new ConditionDefinition(
      "DESIGN IMPULS RATE SURF CONDITIONS", "SurfImpulsRate", "Surface Impuls Rate",
      DRT::Condition::ImpulsRateThroughSurface_3D, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < impulsratecomponents.size(); ++i)
  {
    surfimpulsrate->AddComponent(impulsratecomponents[i]);
  }
  condlist.push_back(surfimpulsrate);


  /*--------------------------------------------------------------------*/
  // Volumetric surface flow profile condition
  Teuchos::RCP<ConditionDefinition> volumetric_surface_flow_cond =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURF VOLUMETRIC FLOW CONDITIONS",
          "VolumetricSurfaceFlowCond", "volumetric surface flow condition",
          DRT::Condition::VolumetricSurfaceFlowCond, true, DRT::Condition::Surface));

  std::vector<Teuchos::RCP<ConditionComponent>> inflownormalcomponents;

  volumetric_surface_flow_cond->AddComponent(
      Teuchos::rcp(new IntConditionComponent("ConditionID")));

  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new StringConditionComponent(
      "ConditionType", "WOMERSLEY", Teuchos::tuple<std::string>("WOMERSLEY", "POLYNOMIAL"),
      Teuchos::tuple<std::string>("WOMERSLEY", "POLYNOMIAL"), true)));

  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new StringConditionComponent("prebiased",
      "NOTPREBIASED", Teuchos::tuple<std::string>("NOTPREBIASED", "PREBIASED", "FORCED"),
      Teuchos::tuple<std::string>("NOTPREBIASED", "PREBIASED", "FORCED"), true)));


  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new StringConditionComponent("FlowType",
      "InFlow", Teuchos::tuple<std::string>("InFlow", "OutFlow"),
      Teuchos::tuple<std::string>("InFlow", "OutFlow"), true)));

  volumetric_surface_flow_cond->AddComponent(
      Teuchos::rcp(new StringConditionComponent("CorrectionFlag", "WithOutCorrection",
          Teuchos::tuple<std::string>("WithOutCorrection", "WithCorrection"),
          Teuchos::tuple<std::string>("WithOutCorrection", "WithCorrection"), true)));
  AddNamedReal(volumetric_surface_flow_cond, "Period");
  AddNamedInt(volumetric_surface_flow_cond, "Order");
  AddNamedInt(volumetric_surface_flow_cond, "Harmonics");
  AddNamedReal(volumetric_surface_flow_cond, "Val");
  AddNamedInt(volumetric_surface_flow_cond, "Funct");

  volumetric_surface_flow_cond->AddComponent(
      Teuchos::rcp(new StringConditionComponent("NORMAL", "SelfEvaluateNormal",
          Teuchos::tuple<std::string>("SelfEvaluateNormal", "UsePrescribedNormal"),
          Teuchos::tuple<std::string>("SelfEvaluateNormal", "UsePrescribedNormal"), true)));


  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("n1")));
  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("n2")));
  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("n3")));


  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new StringConditionComponent(
      "CenterOfMass", "SelfEvaluateCenterOfMass",
      Teuchos::tuple<std::string>("SelfEvaluateCenterOfMass", "UsePrescribedCenterOfMass"),
      Teuchos::tuple<std::string>("SelfEvaluateCenterOfMass", "UsePrescribedCenterOfMass"), true)));

  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("c1")));
  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("c2")));
  volumetric_surface_flow_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("c3")));

  condlist.push_back(volumetric_surface_flow_cond);



  /*--------------------------------------------------------------------*/
  // Volumetric flow border nodes condition

  Teuchos::RCP<ConditionDefinition> volumetric_border_nodes_cond =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE VOLUMETRIC FLOW BORDER NODES",
          "VolumetricFlowBorderNodesCond", "volumetric flow border nodes condition",
          DRT::Condition::VolumetricFlowBorderNodes, true, DRT::Condition::Line));

  volumetric_border_nodes_cond->AddComponent(
      Teuchos::rcp(new IntConditionComponent("ConditionID")));


  condlist.push_back(volumetric_border_nodes_cond);

  /*--------------------------------------------------------------------*/
  // Volumetric surface total traction corrector
  Teuchos::RCP<ConditionDefinition> total_traction_correction_cond =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURF TOTAL TRACTION CORRECTION CONDITIONS",
          "TotalTractionCorrectionCond", "total traction correction condition",
          DRT::Condition::TotalTractionCorrectionCond, true, DRT::Condition::Surface));


  total_traction_correction_cond->AddComponent(
      Teuchos::rcp(new IntConditionComponent("ConditionID")));

  total_traction_correction_cond->AddComponent(Teuchos::rcp(new StringConditionComponent(
      "ConditionType", "POLYNOMIAL", Teuchos::tuple<std::string>("POLYNOMIAL", "WOMERSLEY"),
      Teuchos::tuple<std::string>("POLYNOMIAL", "WOMERSLEY"), true)));

  total_traction_correction_cond->AddComponent(
      Teuchos::rcp(new StringConditionComponent("prebiased", "NOTPREBIASED",
          Teuchos::tuple<std::string>("NOTPREBIASED", "PREBIASED", "FORCED"),
          Teuchos::tuple<std::string>("NOTPREBIASED", "PREBIASED", "FORCED"), true)));

  total_traction_correction_cond->AddComponent(Teuchos::rcp(new StringConditionComponent("FlowType",
      "InFlow", Teuchos::tuple<std::string>("InFlow", "OutFlow"),
      Teuchos::tuple<std::string>("InFlow", "OutFlow"), true)));

  total_traction_correction_cond->AddComponent(
      Teuchos::rcp(new StringConditionComponent("CorrectionFlag", "WithOutCorrection",
          Teuchos::tuple<std::string>("WithOutCorrection", "WithCorrection"),
          Teuchos::tuple<std::string>("WithOutCorrection", "WithCorrection"), true)));
  AddNamedReal(total_traction_correction_cond, "Period");
  AddNamedInt(total_traction_correction_cond, "Order");
  AddNamedInt(total_traction_correction_cond, "Harmonics");
  AddNamedReal(total_traction_correction_cond, "Val");
  AddNamedInt(total_traction_correction_cond, "Funct");

  total_traction_correction_cond->AddComponent(
      Teuchos::rcp(new StringConditionComponent("NORMAL", "SelfEvaluateNormal",
          Teuchos::tuple<std::string>("SelfEvaluateNormal", "UsePrescribedNormal"),
          Teuchos::tuple<std::string>("SelfEvaluateNormal", "UsePrescribedNormal"), true)));

  total_traction_correction_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("n1")));
  total_traction_correction_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("n2")));
  total_traction_correction_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("n3")));

  total_traction_correction_cond->AddComponent(Teuchos::rcp(new StringConditionComponent(
      "CenterOfMass", "SelfEvaluateCenterOfMass",
      Teuchos::tuple<std::string>("SelfEvaluateCenterOfMass", "UsePrescribedCenterOfMass"),
      Teuchos::tuple<std::string>("SelfEvaluateCenterOfMass", "UsePrescribedCenterOfMass"), true)));

  total_traction_correction_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("c1")));
  total_traction_correction_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("c2")));
  total_traction_correction_cond->AddComponent(Teuchos::rcp(new RealConditionComponent("c3")));

  condlist.push_back(total_traction_correction_cond);

  /*--------------------------------------------------------------------*/
  // Volumetric flow traction correction border nodes condition

  Teuchos::RCP<ConditionDefinition> traction_corrector_border_nodes_cond =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE TOTAL TRACTION CORRECTION BORDER NODES",
          "TotalTractionCorrectionBorderNodesCond",
          "total traction correction border nodes condition",
          DRT::Condition::TotalTractionCorrectionBorderNodes, true, DRT::Condition::Line));

  traction_corrector_border_nodes_cond->AddComponent(
      Teuchos::rcp(new IntConditionComponent("ConditionID")));


  condlist.push_back(traction_corrector_border_nodes_cond);



  /*--------------------------------------------------------------------*/
  // no penetration for darcy flow in porous media

  Teuchos::RCP<ConditionDefinition> nopenetration_surf = Teuchos::rcp(
      new ConditionDefinition("DESIGN SURFACE NORMAL NO PENETRATION CONDITION", "NoPenetration",
          "No Penetration", DRT::Condition::NoPenetration, true, DRT::Condition::Surface));

  condlist.push_back(nopenetration_surf);

  /*--------------------------------------------------------------------*/
  // no penetration for darcy flow in porous media

  Teuchos::RCP<ConditionDefinition> nopenetration_line = Teuchos::rcp(
      new ConditionDefinition("DESIGN LINE NORMAL NO PENETRATION CONDITION", "NoPenetration",
          "No Penetration", DRT::Condition::NoPenetration, true, DRT::Condition::Line));

  condlist.push_back(nopenetration_line);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of coupling terms in porous media

  Teuchos::RCP<ConditionDefinition> porocoupling_vol =
      Teuchos::rcp(new ConditionDefinition("DESIGN VOLUME POROCOUPLING CONDITION", "PoroCoupling",
          "Poro Coupling", DRT::Condition::PoroCoupling, true, DRT::Condition::Volume));

  condlist.push_back(porocoupling_vol);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of coupling terms in porous media

  Teuchos::RCP<ConditionDefinition> porocoupling_surf =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE POROCOUPLING CONDITION", "PoroCoupling",
          "Poro Coupling", DRT::Condition::PoroCoupling, true, DRT::Condition::Surface));

  condlist.push_back(porocoupling_surf);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in porous media problems

  Teuchos::RCP<ConditionDefinition> poropartint_surf =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURFACE PORO PARTIAL INTEGRATION", "PoroPartInt",
          "Poro Partial Integration", DRT::Condition::PoroPartInt, true, DRT::Condition::Surface));

  condlist.push_back(poropartint_surf);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in porous media problems

  Teuchos::RCP<ConditionDefinition> poropartint_line =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE PORO PARTIAL INTEGRATION", "PoroPartInt",
          "Poro Partial Integration", DRT::Condition::PoroPartInt, true, DRT::Condition::Line));

  condlist.push_back(poropartint_line);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in porous media problems

  Teuchos::RCP<ConditionDefinition> poropresint_surf = Teuchos::rcp(
      new ConditionDefinition("DESIGN SURFACE PORO PRESSURE INTEGRATION", "PoroPresInt",
          "Poro Pressure Integration", DRT::Condition::PoroPresInt, true, DRT::Condition::Surface));

  condlist.push_back(poropresint_surf);

  /*--------------------------------------------------------------------*/
  // condition for evaluation of boundary terms in porous media problems

  Teuchos::RCP<ConditionDefinition> poropresint_line =
      Teuchos::rcp(new ConditionDefinition("DESIGN LINE PORO PRESSURE INTEGRATION", "PoroPresInt",
          "Poro Pressure Integration", DRT::Condition::PoroPresInt, true, DRT::Condition::Line));

  condlist.push_back(poropresint_line);

  /*--------------------------------------------------------------------*/
  // Fluctuating Hydrodynamics Statistics on a surface

  std::vector<Teuchos::RCP<ConditionComponent>> flucthydrostatsurfcomponents;
  flucthydrostatsurfcomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  flucthydrostatsurfcomponents.push_back(
      Teuchos::rcp(new StringConditionComponent("evaluation type", "nodalbased",
          Teuchos::tuple<std::string>("elebased", "nodalbased", "ele_and_nodalbased"),
          Teuchos::tuple<std::string>("elebased", "nodalbased", "ele_and_nodalbased"))));


  Teuchos::RCP<ConditionDefinition> fluctHydro_statisticsSurf =
      Teuchos::rcp(new ConditionDefinition("DESIGN FLUCTHYDRO STATISTICS SURF CONDITIONS",
          "FluctHydroStatisticsSurf", "FluctHydro_StatisticsSurf",
          DRT::Condition::FluctHydro_StatisticsSurf, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < flucthydrostatsurfcomponents.size(); ++i)
    fluctHydro_statisticsSurf->AddComponent(flucthydrostatsurfcomponents[i]);

  condlist.push_back(fluctHydro_statisticsSurf);

  /*--------------------------------------------------------------------*/
  // Fluctuating Hydrodynamics Statistics on a line

  std::vector<Teuchos::RCP<ConditionComponent>> flucthydrostatlinecomponents;
  flucthydrostatlinecomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  flucthydrostatlinecomponents.push_back(
      Teuchos::rcp(new StringConditionComponent("evaluation type", "nodalbased",
          Teuchos::tuple<std::string>("elebased", "nodalbased", "ele_and_nodalbased"),
          Teuchos::tuple<std::string>("elebased", "nodalbased", "ele_and_nodalbased"))));

  Teuchos::RCP<ConditionDefinition> fluctHydro_statisticsLine =
      Teuchos::rcp(new ConditionDefinition("DESIGN FLUCTHYDRO STATISTICS LINE CONDITIONS",
          "FluctHydroStatisticsLine", "FluctHydro_StatisticsLine",
          DRT::Condition::FluctHydro_StatisticsLine, true, DRT::Condition::Line));

  for (unsigned i = 0; i < flucthydrostatlinecomponents.size(); ++i)
    fluctHydro_statisticsLine->AddComponent(flucthydrostatlinecomponents[i]);

  condlist.push_back(fluctHydro_statisticsLine);
}
