/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for scatra

\level 1


*/
/*----------------------------------------------------------------------*/
#include "baci_inpar_scatra.H"

#include "baci_inpar_bio.H"
#include "baci_inpar_fluid.H"
#include "baci_inpar_s2i.H"
#include "baci_inpar_validparameters.H"
#include "baci_lib_conditiondefinition.H"
#include "baci_linalg_equilibrate.H"
#include "baci_linalg_sparseoperator.H"

void INPAR::SCATRA::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& scatradyn = list->sublist(
      "SCALAR TRANSPORT DYNAMIC", false, "control parameters for scalar transport problems\n");

  setStringToIntegralParameter<int>("SOLVERTYPE", "linear_full", "type of scalar transport solver",
      tuple<std::string>("linear_full", "linear_incremental", "nonlinear",
          "nonlinear_multiscale_macrotomicro", "nonlinear_multiscale_macrotomicro_aitken",
          "nonlinear_multiscale_macrotomicro_aitken_dofsplit", "nonlinear_multiscale_microtomacro"),
      tuple<int>(solvertype_linear_full, solvertype_linear_incremental, solvertype_nonlinear,
          solvertype_nonlinear_multiscale_macrotomicro,
          solvertype_nonlinear_multiscale_macrotomicro_aitken,
          solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit,
          solvertype_nonlinear_multiscale_microtomacro),
      &scatradyn);

  setStringToIntegralParameter<int>("TIMEINTEGR", "One_Step_Theta", "Time Integration Scheme",
      tuple<std::string>("Stationary", "One_Step_Theta", "BDF2", "Gen_Alpha"),
      tuple<int>(timeint_stationary, timeint_one_step_theta, timeint_bdf2, timeint_gen_alpha),
      &scatradyn);

  DoubleParameter("MAXTIME", 1000.0, "Total simulation time", &scatradyn);
  IntParameter("NUMSTEP", 20, "Total number of time steps", &scatradyn);
  DoubleParameter("TIMESTEP", 0.1, "Time increment dt", &scatradyn);
  DoubleParameter("THETA", 0.5, "One-step-theta time integration factor", &scatradyn);
  DoubleParameter("ALPHA_M", 0.5, "Generalized-alpha time integration factor", &scatradyn);
  DoubleParameter("ALPHA_F", 0.5, "Generalized-alpha time integration factor", &scatradyn);
  DoubleParameter("GAMMA", 0.5, "Generalized-alpha time integration factor", &scatradyn);
  IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &scatradyn);
  IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &scatradyn);
  IntParameter("MATID", -1, "Material ID for automatic mesh generation", &scatradyn);

  setStringToIntegralParameter<int>("VELOCITYFIELD", "zero",
      "type of velocity field used for scalar transport problems",
      tuple<std::string>("zero", "function", "Navier_Stokes"),
      tuple<int>(velocity_zero, velocity_function, velocity_Navier_Stokes), &scatradyn);

  IntParameter("VELFUNCNO", -1, "function number for scalar transport velocity field", &scatradyn);

  {
    // a standard Teuchos::tuple can have at maximum 10 entries! We have to circumvent this here.
    Teuchos::Tuple<std::string, 13> name;
    Teuchos::Tuple<int, 13> label;
    name[0] = "zero_field";
    label[0] = initfield_zero_field;
    name[1] = "field_by_function";
    label[1] = initfield_field_by_function;
    name[2] = "field_by_condition";
    label[2] = initfield_field_by_condition;
    name[3] = "disturbed_field_by_function";
    label[3] = initfield_disturbed_field_by_function;
    name[4] = "1D_DISCONTPV";
    label[4] = initfield_discontprogvar_1D;
    name[5] = "FLAME_VORTEX_INTERACTION";
    label[5] = initfield_flame_vortex_interaction;
    name[6] = "RAYTAYMIXFRAC";
    label[6] = initfield_raytaymixfrac;
    name[7] = "L_shaped_domain";
    label[7] = initfield_Lshapeddomain;
    name[8] = "facing_flame_fronts";
    label[8] = initfield_facing_flame_fronts;
    name[9] = "oracles_flame";
    label[9] = initfield_oracles_flame;
    name[10] = "high_forced_hit";
    label[10] = initialfield_forced_hit_high_Sc;
    name[11] = "low_forced_hit";
    label[11] = initialfield_forced_hit_low_Sc;
    name[12] = "algebraic_field_dependence";
    label[12] = initialfield_algebraic_field_dependence;

    setStringToIntegralParameter<int>("INITIALFIELD", "zero_field",
        "Initial Field for scalar transport problem", name, label, &scatradyn);
  }

  IntParameter("INITFUNCNO", -1, "function number for scalar transport initial field", &scatradyn);

  BoolParameter("SPHERICALCOORDS", "No", "use of spherical coordinates", &scatradyn);

  setStringToIntegralParameter<int>("CALCERROR", "No",
      "compute error compared to analytical solution",
      tuple<std::string>("No", "Kwok_Wu", "ConcentricCylinders", "Electroneutrality",
          "error_by_function", "error_by_condition", "SphereDiffusion", "AnalyticSeries"),
      tuple<int>(calcerror_no, calcerror_Kwok_Wu, calcerror_cylinder, calcerror_electroneutrality,
          calcerror_byfunction, calcerror_bycondition, calcerror_spherediffusion,
          calcerror_AnalyticSeries),
      &scatradyn);

  IntParameter(
      "CALCERRORNO", -1, "function number for scalar transport error computation", &scatradyn);

  setStringToIntegralParameter<int>("CALCFLUX_DOMAIN", "No",
      "output of diffusive/total flux vectors inside domain",
      tuple<std::string>("No", "total", "diffusive"),
      tuple<int>(flux_none, flux_total, flux_diffusive), &scatradyn);

  BoolParameter("CALCFLUX_DOMAIN_LUMPED", "Yes",
      "perform approximate domain flux calculation involving matrix lumping", &scatradyn);

  setStringToIntegralParameter<int>("CALCFLUX_BOUNDARY", "No",
      "output of convective/diffusive/total flux vectors on boundary",
      tuple<std::string>("No", "total", "diffusive", "convective"),
      tuple<int>(flux_none, flux_total, flux_diffusive, flux_convective), &scatradyn);

  BoolParameter("CALCFLUX_BOUNDARY_LUMPED", "Yes",
      "perform approximate boundary flux calculation involving matrix lumping", &scatradyn);

  setNumericStringParameter("WRITEFLUX_IDS", "-1",
      "Write diffusive/total flux vector fields for these scalar fields only (starting with 1)",
      &scatradyn);

  setStringToIntegralParameter<int>("OUTPUTSCALARS", "none",
      "Output of total and mean values for transported scalars",
      tuple<std::string>("none", "entire_domain", "by_condition", "entire_domain_and_by_condition"),
      tuple<int>(outputscalars_none, outputscalars_entiredomain, outputscalars_condition,
          outputscalars_entiredomain_condition),
      &scatradyn);
  BoolParameter("OUTPUTSCALARSMEANGRAD", "No", "Output of mean gradient of scalars", &scatradyn);
  BoolParameter("OUTINTEGRREAC", "No", "Output of integral reaction values", &scatradyn);
  BoolParameter("OUTPUT_GMSH", "No", "Do you want to write Gmsh postprocessing files?", &scatradyn);

  BoolParameter("MATLAB_STATE_OUTPUT", "No",
      "Do you want to write the state solution to Matlab file?", &scatradyn);

  setStringToIntegralParameter<int>("CONVFORM", "convective", "form of convective term",
      tuple<std::string>("convective", "conservative"),
      tuple<int>(convform_convective, convform_conservative), &scatradyn);

  BoolParameter(
      "NEUMANNINFLOW", "no", "Flag to (de)activate potential Neumann inflow term(s)", &scatradyn);

  BoolParameter("CONV_HEAT_TRANS", "no",
      "Flag to (de)activate potential convective heat transfer boundary conditions", &scatradyn);

  BoolParameter(
      "SKIPINITDER", "no", "Flag to skip computation of initial time derivative", &scatradyn);

  setStringToIntegralParameter<int>("FSSUGRDIFF", "No", "fine-scale subgrid diffusivity",
      tuple<std::string>("No", "artificial", "Smagorinsky_all", "Smagorinsky_small"),
      tuple<int>(fssugrdiff_no, fssugrdiff_artificial, fssugrdiff_smagorinsky_all,
          fssugrdiff_smagorinsky_small),
      &scatradyn);

  // flag for output of performance statistics associated with nonlinear solver into *.csv file
  BoolParameter("ELECTROMAGNETICDIFFUSION", "No",
      "flag to activate electromagnetic diffusion problems", &scatradyn);

  // Current density source function for EMD problems
  IntParameter("EMDSOURCE", -1, "Current density source", &scatradyn);

  setStringToIntegralParameter<int>("MESHTYING", "no", "Flag to (de)activate mesh tying algorithm",
      tuple<std::string>("no", "Condensed_Smat", "Condensed_Bmat",
          "Condensed_Bmat_merged"),  // use the condensed_bmat_merged strategy
      tuple<int>(INPAR::FLUID::no_meshtying, INPAR::FLUID::condensed_smat,
          INPAR::FLUID::condensed_bmat,
          INPAR::FLUID::condensed_bmat_merged),  // use the condensed_bmat_merged strategy
      &scatradyn);

  // Type of coupling strategy between the two fields
  setStringToIntegralParameter<int>("FIELDCOUPLING", "matching",
      "Type of coupling strategy between fields", tuple<std::string>("matching", "volmortar"),
      tuple<int>(coupling_match, coupling_volmortar), &scatradyn);

  // linear solver id used for scalar transport/elch problems
  IntParameter(
      "LINEAR_SOLVER", -1, "number of linear solver used for scalar transport/elch...", &scatradyn);
  // linear solver id used for l2 projection problems (e.g. gradient projections)
  IntParameter("L2_PROJ_LINEAR_SOLVER", -1,
      "number of linear solver used for l2-projection sub-problems", &scatradyn);
  // IntParameter("SIMPLER_SOLVER",-1,"number of linear solver used for ELCH (solved with
  // SIMPLER)...",&scatradyn);

  // flag for equilibration of global system of equations
  setStringToIntegralParameter<CORE::LINALG::EquilibrationMethod>("EQUILIBRATION", "none",
      "flag for equilibration of global system of equations",
      tuple<std::string>("none", "rows_full", "rows_maindiag", "columns_full", "columns_maindiag",
          "rowsandcolumns_full", "rowsandcolumns_maindiag"),
      tuple<CORE::LINALG::EquilibrationMethod>(CORE::LINALG::EquilibrationMethod::none,
          CORE::LINALG::EquilibrationMethod::rows_full,
          CORE::LINALG::EquilibrationMethod::rows_maindiag,
          CORE::LINALG::EquilibrationMethod::columns_full,
          CORE::LINALG::EquilibrationMethod::columns_maindiag,
          CORE::LINALG::EquilibrationMethod::rowsandcolumns_full,
          CORE::LINALG::EquilibrationMethod::rowsandcolumns_maindiag),
      &scatradyn);

  // type of global system matrix in global system of equations
  setStringToIntegralParameter<CORE::LINALG::MatrixType>("MATRIXTYPE", "sparse",
      "type of global system matrix in global system of equations",
      tuple<std::string>("sparse", "block_condition", "block_condition_dof"),
      tuple<CORE::LINALG::MatrixType>(CORE::LINALG::MatrixType::sparse,
          CORE::LINALG::MatrixType::block_condition, CORE::LINALG::MatrixType::block_condition_dof),
      &scatradyn);

  // flag for natural convection effects
  BoolParameter("NATURAL_CONVECTION", "No", "Include natural convection effects", &scatradyn);

  // parameters for finite difference check
  setStringToIntegralParameter<int>("FDCHECK", "none",
      "flag for finite difference check: none, local, or global",
      tuple<std::string>("none",
          "global",           // perform finite difference check on time integrator level
          "global_extended",  // perform finite difference check on time integrator level for
                              // extended system matrix (e.g., involving Lagrange multipliers or
                              // interface layer thicknesses)
          "local"             // perform finite difference check on element level
          ),
      tuple<int>(fdcheck_none, fdcheck_global, fdcheck_global_extended, fdcheck_local), &scatradyn);
  DoubleParameter("FDCHECKEPS", 1.e-6,
      "dof perturbation magnitude for finite difference check (1.e-6 seems to work very well, "
      "whereas smaller values don't)",
      &scatradyn);
  DoubleParameter(
      "FDCHECKTOL", 1.e-6, "relative tolerance for finite difference check", &scatradyn);

  // parameter for optional computation of domain and boundary integrals, i.e., of surface areas and
  // volumes associated with specified nodesets
  setStringToIntegralParameter<int>("COMPUTEINTEGRALS", "none",
      "flag for optional computation of domain integrals",
      tuple<std::string>("none", "initial", "repeated"),
      tuple<int>(computeintegrals_none, computeintegrals_initial, computeintegrals_repeated),
      &scatradyn);

  // parameter for using p-adpativity and semi-implicit evaluation of the reaction term (at the
  // moment only used for HDG and cardiac monodomain problems)
  BoolParameter("PADAPTIVITY", "no", "Flag to (de)activate p-adativity", &scatradyn);
  DoubleParameter("PADAPTERRORTOL", 1e-6,
      "The error tolerance to calculate the variation of the elemental degree", &scatradyn);
  DoubleParameter("PADAPTERRORBASE", 1.66,
      "The error tolerance base to calculate the variation of the elemental degree", &scatradyn);
  IntParameter("PADAPTDEGREEMAX", 4, "The max. degree of the shape functions", &scatradyn);
  BoolParameter("SEMIIMPLICIT", "no",
      "Flag to (de)activate semi-implicit calculation of the reaction term", &scatradyn);

  // flag for output of performance statistics associated with linear solver into *.csv file
  BoolParameter("OUTPUTLINSOLVERSTATS", "No",
      "flag for output of performance statistics associated with linear solver into csv file",
      &scatradyn);

  // flag for output of performance statistics associated with nonlinear solver into *.csv file
  BoolParameter("OUTPUTNONLINSOLVERSTATS", "No",
      "flag for output of performance statistics associated with nonlinear solver into csv file",
      &scatradyn);

  // flag for point-based null space calculation
  BoolParameter(
      "NULLSPACE_POINTBASED", "No", "flag for point-based null space calculation", &scatradyn);

  // flag for adaptive time stepping
  BoolParameter("ADAPTIVE_TIMESTEPPING", "No", "flag for adaptive time stepping", &scatradyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scatra_nonlin = scatradyn.sublist(
      "NONLINEAR", false, "control parameters for solving nonlinear SCATRA problems\n");

  IntParameter("ITEMAX", 10, "max. number of nonlin. iterations", &scatra_nonlin);
  DoubleParameter("CONVTOL", 1e-6, "Tolerance for convergence check", &scatra_nonlin);
  IntParameter("ITEMAX_OUTER", 10,
      "Maximum number of outer iterations in partitioned coupling schemes (natural convection, "
      "multi-scale simulations etc.)",
      &scatra_nonlin);
  DoubleParameter("CONVTOL_OUTER", 1e-6,
      "Convergence check tolerance for outer loop in partitioned coupling schemes (natural "
      "convection, multi-scale simulations etc.)",
      &scatra_nonlin);
  BoolParameter("EXPLPREDICT", "no",
      "do an explicit predictor step before starting nonlinear iteration", &scatra_nonlin);
  DoubleParameter("ABSTOLRES", 1e-14,
      "Absolute tolerance for deciding if residual of nonlinear problem is already zero",
      &scatra_nonlin);

  // convergence criteria adaptivity
  BoolParameter("ADAPTCONV", "No",
      "Switch on adaptive control of linear solver tolerance for nonlinear solution",
      &scatra_nonlin);
  DoubleParameter("ADAPTCONV_BETTER", 0.1,
      "The linear solver shall be this much better than the current nonlinear residual in the "
      "nonlinear convergence limit",
      &scatra_nonlin);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scatradyn_stab = scatradyn.sublist("STABILIZATION", false,
      "control parameters for the stabilization of scalar transport problems");

  // this parameter governs type of stabilization
  setStringToIntegralParameter<int>("STABTYPE", "SUPG", "type of stabilization (if any)",
      tuple<std::string>("no_stabilization", "SUPG", "GLS", "USFEM", "centered", "upwind"),
      tuple<std::string>(
          "Do not use any stabilization -> only reasonable for low-Peclet-number flows", "Use SUPG",
          "Use GLS", "Use USFEM", "Use centered scheme", "Use upwinded scheme"),
      tuple<int>(stabtype_no_stabilization, stabtype_SUPG, stabtype_GLS, stabtype_USFEM,
          stabtype_hdg_centered, stabtype_hdg_upwind),
      &scatradyn_stab);

  // this parameter governs whether subgrid-scale velocity is included
  BoolParameter(
      "SUGRVEL", "no", "potential incorporation of subgrid-scale velocity", &scatradyn_stab);

  // this parameter governs whether all-scale subgrid diffusivity is included
  BoolParameter("ASSUGRDIFF", "no",
      "potential incorporation of all-scale subgrid diffusivity (a.k.a. discontinuity-capturing) "
      "term",
      &scatradyn_stab);

  // this parameter selects the tau definition applied
  setStringToIntegralParameter<int>("DEFINITION_TAU", "Franca_Valentin", "Definition of tau",
      tuple<std::string>("Taylor_Hughes_Zarins", "Taylor_Hughes_Zarins_wo_dt", "Franca_Valentin",
          "Franca_Valentin_wo_dt", "Shakib_Hughes_Codina", "Shakib_Hughes_Codina_wo_dt", "Codina",
          "Codina_wo_dt", "Franca_Madureira_Valentin", "Franca_Madureira_Valentin_wo_dt",
          "Exact_1D", "Zero", "Numerical_Value"),
      tuple<int>(tau_taylor_hughes_zarins, tau_taylor_hughes_zarins_wo_dt, tau_franca_valentin,
          tau_franca_valentin_wo_dt, tau_shakib_hughes_codina, tau_shakib_hughes_codina_wo_dt,
          tau_codina, tau_codina_wo_dt, tau_franca_madureira_valentin,
          tau_franca_madureira_valentin_wo_dt, tau_exact_1d, tau_zero, tau_numerical_value),
      &scatradyn_stab);

  // this parameter selects the characteristic element length for tau for all
  // stabilization parameter definitions requiring such a length
  setStringToIntegralParameter<int>("CHARELELENGTH", "streamlength",
      "Characteristic element length for tau",
      tuple<std::string>("streamlength", "volume_equivalent_diameter", "root_of_volume"),
      tuple<int>(streamlength, volume_equivalent_diameter, root_of_volume), &scatradyn_stab);

  // this parameter selects the all-scale subgrid-diffusivity definition applied
  setStringToIntegralParameter<int>("DEFINITION_ASSGD", "artificial_linear",
      "Definition of (all-scale) subgrid diffusivity",
      tuple<std::string>("artificial_linear", "artificial_linear_reinit",
          "Hughes_etal_86_nonlinear", "Tezduyar_Park_86_nonlinear",
          "Tezduyar_Park_86_nonlinear_wo_phizero", "doCarmo_Galeao_91_nonlinear",
          "Almeida_Silva_97_nonlinear", "YZbeta_nonlinear", "Codina_nonlinear"),
      tuple<std::string>("classical linear artificial subgrid-diffusivity",
          "simple linear artificial subgrid-diffusivity const*h",
          "nonlinear isotropic according to Hughes et al. (1986)",
          "nonlinear isotropic according to Tezduyar and Park (1986)",
          "nonlinear isotropic according to Tezduyar and Park (1986) without user parameter "
          "phi_zero",
          "nonlinear isotropic according to doCarmo and Galeao (1991)",
          "nonlinear isotropic according to Almeida and Silva (1997)", "nonlinear YZ beta model",
          "nonlinear isotropic according to Codina"),
      tuple<int>(assgd_artificial, assgd_lin_reinit, assgd_hughes, assgd_tezduyar,
          assgd_tezduyar_wo_phizero, assgd_docarmo, assgd_almeida, assgd_yzbeta, assgd_codina),
      &scatradyn_stab);

  // this parameter selects the location where tau is evaluated
  setStringToIntegralParameter<int>("EVALUATION_TAU", "element_center",
      "Location where tau is evaluated", tuple<std::string>("element_center", "integration_point"),
      tuple<std::string>("evaluate tau at element center", "evaluate tau at integration point"),
      tuple<int>(evaltau_element_center, evaltau_integration_point), &scatradyn_stab);

  // this parameter selects the location where the material law is evaluated
  // (does not fit here very well, but parameter transfer is easier)
  setStringToIntegralParameter<int>("EVALUATION_MAT", "element_center",
      "Location where material law is evaluated",
      tuple<std::string>("element_center", "integration_point"),
      tuple<std::string>(
          "evaluate material law at element center", "evaluate material law at integration point"),
      tuple<int>(evalmat_element_center, evalmat_integration_point), &scatradyn_stab);

  // this parameter selects methods for improving consistency of stabilization terms
  setStringToIntegralParameter<int>("CONSISTENCY", "no",
      "improvement of consistency for stabilization",
      tuple<std::string>("no", "L2_projection_lumped"),
      tuple<std::string>("inconsistent", "L2 projection with lumped mass matrix"),
      tuple<int>(consistency_no, consistency_l2_projection_lumped), &scatradyn_stab);

  // this parameter defines the numerical value, if stabilization with numerical values is used
  DoubleParameter("TAU_VALUE", 0.0, "Numerical value for tau for stabilization", &scatradyn_stab);

  // ----------------------------------------------------------------------
  // artery mesh tying
  Teuchos::ParameterList& scatradyn_art =
      scatradyn.sublist("ARTERY COUPLING", false, "Parameters for artery mesh tying");

  setStringToIntegralParameter<int>("ARTERY_COUPLING_METHOD", "None",
      "Coupling method for artery coupling.",
      tuple<std::string>("None", "Nodal", "GPTS", "MP", "NTP"),
      tuple<std::string>("none", "Nodal Coupling", "Gauss-Point-To-Segment Approach",
          "Mortar Penalty Approach", "1D node-to-point in 2D/3D Approach"),
      tuple<int>(INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::none,  // none
          INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::nodal,        // Nodal Coupling
          INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::gpts,  // Gauss-Point-To-Segment
                                                                          // Approach
          INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp,  // Mortar Penalty Approach//
          INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::ntp  // 1D node-to-point in 2D/3D

          ),
      &scatradyn_art);

  // penalty parameter
  DoubleParameter("PENALTY", 1000.0, "Penalty parameter for line-based coupling", &scatradyn_art);

  // coupled artery dofs for mesh tying
  setNumericStringParameter(
      "COUPLEDDOFS_ARTSCATRA", "-1.0", "coupled artery dofs for mesh tying", &scatradyn_art);

  // coupled porofluid dofs for mesh tying
  setNumericStringParameter(
      "COUPLEDDOFS_SCATRA", "-1.0", "coupled porofluid dofs for mesh tying", &scatradyn_art);

  // functions for coupling (arteryscatra part)
  setNumericStringParameter(
      "REACFUNCT_ART", "-1", "functions for coupling (arteryscatra part)", &scatradyn_art);

  // scale for coupling (arteryscatra part)
  setNumericStringParameter(
      "SCALEREAC_ART", "0", "scale for coupling (arteryscatra part)", &scatradyn_art);

  // functions for coupling (scatra part)
  setNumericStringParameter(
      "REACFUNCT_CONT", "-1", "functions for coupling (scatra part)", &scatradyn_art);

  // scale for coupling (scatra part)
  setNumericStringParameter(
      "SCALEREAC_CONT", "0", "scale for coupling (scatra part)", &scatradyn_art);
}



void INPAR::SCATRA::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // Boundary flux evaluation condition for scalar transport
  Teuchos::RCP<ConditionDefinition> linebndryfluxeval =
      Teuchos::rcp(new ConditionDefinition("SCATRA FLUX CALC LINE CONDITIONS", "ScaTraFluxCalc",
          "Scalar Transport Boundary Flux Calculation", DRT::Condition::ScaTraFluxCalc, true,
          DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfbndryfluxeval =
      Teuchos::rcp(new ConditionDefinition("SCATRA FLUX CALC SURF CONDITIONS", "ScaTraFluxCalc",
          "Scalar Transport Boundary Flux Calculation", DRT::Condition::ScaTraFluxCalc, true,
          DRT::Condition::Surface));
  condlist.emplace_back(linebndryfluxeval);
  condlist.emplace_back(surfbndryfluxeval);

  /*--------------------------------------------------------------------*/
  // conditions for calculation of total and mean values of transported scalars
  Teuchos::RCP<ConditionDefinition> totalandmeanscalarline =
      Teuchos::rcp(new ConditionDefinition("DESIGN TOTAL AND MEAN SCALAR LINE CONDITIONS",
          "TotalAndMeanScalar", "calculation of total and mean values of transported scalars",
          DRT::Condition::TotalAndMeanScalar, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> totalandmeanscalarsurf =
      Teuchos::rcp(new ConditionDefinition("DESIGN TOTAL AND MEAN SCALAR SURF CONDITIONS",
          "TotalAndMeanScalar", "calculation of total and mean values of transported scalars",
          DRT::Condition::TotalAndMeanScalar, true, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> totalandmeanscalarvol =
      Teuchos::rcp(new ConditionDefinition("DESIGN TOTAL AND MEAN SCALAR VOL CONDITIONS",
          "TotalAndMeanScalar", "calculation of total and mean values of transported scalars",
          DRT::Condition::TotalAndMeanScalar, true, DRT::Condition::Volume));

  // equip condition definitions with input file line components
  std::vector<Teuchos::RCP<::INPUT::LineComponent>> totalandmeanscalarcomponents;

  {
    totalandmeanscalarcomponents.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("ID")));
    totalandmeanscalarcomponents.emplace_back(
        Teuchos::rcp(new ::INPUT::IntComponent("ConditionID")));
  }

  // insert input file line components into condition definitions
  for (auto& totalandmeanscalarcomponent : totalandmeanscalarcomponents)
  {
    totalandmeanscalarline->AddComponent(totalandmeanscalarcomponent);
    totalandmeanscalarsurf->AddComponent(totalandmeanscalarcomponent);
    totalandmeanscalarvol->AddComponent(totalandmeanscalarcomponent);
  }

  // insert condition definitions into global list of valid condition definitions
  condlist.emplace_back(totalandmeanscalarline);
  condlist.emplace_back(totalandmeanscalarsurf);
  condlist.emplace_back(totalandmeanscalarvol);

  /*--------------------------------------------------------------------*/
  // conditions for calculation of relative error with reference to analytical solution
  Teuchos::RCP<ConditionDefinition> relerrorline =
      Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA RELATIVE ERROR LINE CONDITIONS",
          "ScatraRelError", "calculation of relative error with reference to analytical solution",
          DRT::Condition::ScatraRelError, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> relerrorsurf =
      Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA RELATIVE ERROR SURF CONDITIONS",
          "ScatraRelError", "calculation of relative error with reference to analytical solution",
          DRT::Condition::ScatraRelError, true, DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> relerrorvol =
      Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA RELATIVE ERROR VOL CONDITIONS",
          "ScatraRelError", "calculation of relative error with reference to analytical solution",
          DRT::Condition::ScatraRelError, true, DRT::Condition::Volume));

  // equip condition definitions with input file line components
  std::vector<Teuchos::RCP<::INPUT::LineComponent>> relerrorcomponents;

  {
    relerrorcomponents.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("ID")));
    relerrorcomponents.emplace_back(Teuchos::rcp(new ::INPUT::IntComponent("ConditionID")));
    relerrorcomponents.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("Function")));
    relerrorcomponents.emplace_back(Teuchos::rcp(new ::INPUT::IntComponent("FunctionID")));
  }

  // insert input file line components into condition definitions
  for (auto& relerrorcomponent : relerrorcomponents)
  {
    relerrorline->AddComponent(relerrorcomponent);
    relerrorsurf->AddComponent(relerrorcomponent);
    relerrorvol->AddComponent(relerrorcomponent);
  }

  // insert condition definitions into global list of valid condition definitions
  condlist.emplace_back(relerrorline);
  condlist.emplace_back(relerrorsurf);
  condlist.emplace_back(relerrorvol);

  /*--------------------------------------------------------------------*/
  // Coupling of different scalar transport fields

  std::vector<Teuchos::RCP<::INPUT::LineComponent>> scatracoupcomponents;

  scatracoupcomponents.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("NUMSCAL")));
  scatracoupcomponents.emplace_back(Teuchos::rcp(new ::INPUT::IntComponent("numscal")));
  scatracoupcomponents.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("ONOFF")));
  scatracoupcomponents.emplace_back(
      Teuchos::rcp(new ::INPUT::IntVectorComponent("onoff", ::INPUT::LengthFromInt("numscal"))));

  scatracoupcomponents.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("COUPID")));
  scatracoupcomponents.emplace_back(Teuchos::rcp(new ::INPUT::IntComponent("coupling id")));
  scatracoupcomponents.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("PERMCOEF")));
  scatracoupcomponents.emplace_back(
      Teuchos::rcp(new ::INPUT::RealComponent("permeability coefficient")));
  scatracoupcomponents.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("CONDUCT")));
  scatracoupcomponents.emplace_back(
      Teuchos::rcp(new ::INPUT::RealComponent("hydraulic conductivity")));
  scatracoupcomponents.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("FILTR")));
  scatracoupcomponents.emplace_back(
      Teuchos::rcp(new ::INPUT::RealComponent("filtration coefficient")));
  scatracoupcomponents.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("WSSONOFF")));
  scatracoupcomponents.emplace_back(Teuchos::rcp(new ::INPUT::IntComponent("wss onoff")));
  scatracoupcomponents.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("WSSCOEFFS")));
  scatracoupcomponents.emplace_back(
      Teuchos::rcp(new ::INPUT::RealVectorComponent("wss coeffs", 2)));


  Teuchos::RCP<ConditionDefinition> surfscatracoup = Teuchos::rcp(
      new ConditionDefinition("DESIGN SCATRA COUPLING SURF CONDITIONS", "ScaTraCoupling",
          "ScaTra Coupling", DRT::Condition::ScaTraCoupling, true, DRT::Condition::Surface));

  for (auto& scatracoupcomponent : scatracoupcomponents)
    surfscatracoup->AddComponent(scatracoupcomponent);

  condlist.emplace_back(surfscatracoup);

  /*--------------------------------------------------------------------*/
  // Robin boundary condition for scalar transport problems
  // line
  Teuchos::RCP<ConditionDefinition> scatrarobinline =
      Teuchos::rcp(new ConditionDefinition("DESIGN TRANSPORT ROBIN LINE CONDITIONS",
          "TransportRobin", "Scalar Transport Robin Boundary Condition",
          DRT::Condition::TransportRobin, true, DRT::Condition::Line));
  // surface
  Teuchos::RCP<ConditionDefinition> scatrarobinsurf =
      Teuchos::rcp(new ConditionDefinition("DESIGN TRANSPORT ROBIN SURF CONDITIONS",
          "TransportRobin", "Scalar Transport Robin Boundary Condition",
          DRT::Condition::TransportRobin, true, DRT::Condition::Surface));

  std::vector<Teuchos::RCP<::INPUT::LineComponent>> scatrarobincomponents;


  scatrarobincomponents.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("NUMSCAL")));
  scatrarobincomponents.emplace_back(Teuchos::rcp(new ::INPUT::IntComponent("numscal")));
  scatrarobincomponents.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("ONOFF")));
  scatrarobincomponents.emplace_back(
      Teuchos::rcp(new ::INPUT::IntVectorComponent("onoff", ::INPUT::LengthFromInt("numscal"))));

  scatrarobincomponents.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("PREFACTOR")));
  scatrarobincomponents.emplace_back(Teuchos::rcp(new ::INPUT::RealComponent("prefactor")));
  scatrarobincomponents.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("REFVALUE")));
  scatrarobincomponents.emplace_back(Teuchos::rcp(new ::INPUT::RealComponent("refvalue")));

  for (auto& scatrarobincomponent : scatrarobincomponents)
  {
    scatrarobinline->AddComponent(scatrarobincomponent);
    scatrarobinsurf->AddComponent(scatrarobincomponent);
  }

  condlist.emplace_back(scatrarobinline);
  condlist.emplace_back(scatrarobinsurf);

  /*--------------------------------------------------------------------*/
  // Neumann inflow for SCATRA

  Teuchos::RCP<ConditionDefinition> linetransportneumanninflow =
      Teuchos::rcp(new ConditionDefinition("TRANSPORT NEUMANN INFLOW LINE CONDITIONS",
          "TransportNeumannInflow", "Line Transport Neumann Inflow",
          DRT::Condition::TransportNeumannInflow, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surftransportneumanninflow =
      Teuchos::rcp(new ConditionDefinition("TRANSPORT NEUMANN INFLOW SURF CONDITIONS",
          "TransportNeumannInflow", "Surface Transport Neumann Inflow",
          DRT::Condition::TransportNeumannInflow, true, DRT::Condition::Surface));

  condlist.emplace_back(linetransportneumanninflow);
  condlist.emplace_back(surftransportneumanninflow);

  /*--------------------------------------------------------------------*/
  // Scatra convective heat transfer (Newton's law of heat transfer)

  std::vector<Teuchos::RCP<::INPUT::LineComponent>> transportthermoconvectcomponents;

  // decide here if approximation is sufficient
  // --> Tempn (old temperature T_n)
  // or if the exact solution is needed
  // --> Tempnp (current temperature solution T_n+1) with linearisation
  transportthermoconvectcomponents.emplace_back(Teuchos::rcp(new ::INPUT::SelectionComponent(
      "temperature state", "Tempnp", Teuchos::tuple<std::string>("Tempnp", "Tempn"),
      Teuchos::tuple<std::string>("Tempnp", "Tempn"))));
  // heat transfer coefficient h
  transportthermoconvectcomponents.emplace_back(
      Teuchos::rcp(new ::INPUT::SeparatorComponent("coeff")));
  transportthermoconvectcomponents.emplace_back(Teuchos::rcp(new ::INPUT::RealComponent("coeff")));
  // surrounding (fluid) temperature T_oo
  transportthermoconvectcomponents.emplace_back(
      Teuchos::rcp(new ::INPUT::SeparatorComponent("surtemp")));
  transportthermoconvectcomponents.emplace_back(
      Teuchos::rcp(new ::INPUT::RealComponent("surtemp")));
  // time curve to increase the surrounding (fluid) temperature T_oo in time
  transportthermoconvectcomponents.emplace_back(
      Teuchos::rcp(new ::INPUT::SeparatorComponent("surtempfunct")));
  transportthermoconvectcomponents.emplace_back(
      Teuchos::rcp(new ::INPUT::IntComponent("surtempfunct", {0, true, true})));
  // time curve to increase the complete boundary condition, i.e., the heat flux
  transportthermoconvectcomponents.emplace_back(
      Teuchos::rcp(new ::INPUT::SeparatorComponent("funct")));
  transportthermoconvectcomponents.emplace_back(
      Teuchos::rcp(new ::INPUT::IntComponent("funct", {0, true, true})));

  Teuchos::RCP<ConditionDefinition> linetransportthermoconvect =
      Teuchos::rcp(new ConditionDefinition("TRANSPORT THERMO CONVECTION LINE CONDITIONS",
          "TransportThermoConvections", "Line Transport Thermo Convections",
          DRT::Condition::TransportThermoConvections, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surftransportthermoconvect =
      Teuchos::rcp(new ConditionDefinition("TRANSPORT THERMO CONVECTION SURF CONDITIONS",
          "TransportThermoConvections", "Surface Transport Thermo Convections",
          DRT::Condition::TransportThermoConvections, true, DRT::Condition::Surface));

  for (auto& transportthermoconvectcomponent : transportthermoconvectcomponents)
  {
    linetransportthermoconvect->AddComponent(transportthermoconvectcomponent);
    surftransportthermoconvect->AddComponent(transportthermoconvectcomponent);
  }

  condlist.emplace_back(linetransportthermoconvect);
  condlist.emplace_back(surftransportthermoconvect);

  /*--------------------------------------------------------------------*/
  // conditions for calculation of calculation of heterogeneous reactions
  Teuchos::RCP<ConditionDefinition> scatraheteroreactionmasterline = Teuchos::rcp(
      new ConditionDefinition("DESIGN SCATRA HETEROGENEOUS REACTION LINE CONDITIONS / MASTER",
          "ScatraHeteroReactionMaster", "calculation of heterogeneous reactions",
          DRT::Condition::ScatraHeteroReactionCondMaster, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> scatraheteroreactionmastersurf = Teuchos::rcp(
      new ConditionDefinition("DESIGN SCATRA HETEROGENEOUS REACTION SURF CONDITIONS / MASTER",
          "ScatraHeteroReactionMaster", "calculation of heterogeneous reactions",
          DRT::Condition::ScatraHeteroReactionCondMaster, true, DRT::Condition::Surface));

  // insert condition definitions into global list of valid condition definitions
  condlist.emplace_back(scatraheteroreactionmasterline);
  condlist.emplace_back(scatraheteroreactionmastersurf);

  /*--------------------------------------------------------------------*/
  // conditions for calculation of calculation of heterogeneous reactions
  Teuchos::RCP<ConditionDefinition> scatraheteroreactionslaveline = Teuchos::rcp(
      new ConditionDefinition("DESIGN SCATRA HETEROGENEOUS REACTION LINE CONDITIONS / SLAVE",
          "ScatraHeteroReactionSlave", "calculation of heterogeneous reactions",
          DRT::Condition::ScatraHeteroReactionCondSlave, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> scatraheteroreactionslavesurf = Teuchos::rcp(
      new ConditionDefinition("DESIGN SCATRA HETEROGENEOUS REACTION SURF CONDITIONS / SLAVE",
          "ScatraHeteroReactionSlave", "calculation of heterogeneous reactions",
          DRT::Condition::ScatraHeteroReactionCondSlave, true, DRT::Condition::Surface));

  // insert condition definitions into global list of valid condition definitions
  condlist.emplace_back(scatraheteroreactionslaveline);
  condlist.emplace_back(scatraheteroreactionslavesurf);


  /*--------------------------------------------------------------------*/
  // scatra domain partitioning for block preconditioning of global system matrix
  // please note: this is currently only used in combination with scatra-scatra interface coupling
  // however the complete scatra matrix is subdivided into blocks which is not related to the
  // interface coupling at all
  {
    // partitioning of 2D domain into 2D subdomains
    auto scatrasurfpartitioning =
        Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA SURF CONDITIONS / PARTITIONING",
            "ScatraPartitioning", "Domain partitioning of scatra field",
            DRT::Condition::ScatraPartitioning, false, DRT::Condition::Surface));

    // partitioning of 3D domain into 3D subdomains
    auto scatravolpartitioning =
        Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA VOL CONDITIONS / PARTITIONING",
            "ScatraPartitioning", "Domain partitioning of scatra field",
            DRT::Condition::ScatraPartitioning, false, DRT::Condition::Volume));

    // insert condition definitions into global list of valid condition definitions
    condlist.emplace_back(scatrasurfpartitioning);
    condlist.emplace_back(scatravolpartitioning);
  }
}
