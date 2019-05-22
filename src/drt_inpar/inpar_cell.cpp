/*----------------------------------------------------------------------*/
/*!
\brief input parameters and conditions for cell migration

\maintainer Jonas Eichinger

\level 2
*/
/*----------------------------------------------------------------------*/


#include "inpar.H"
#include "inpar_cell.H"
#include "inpar_fluid.H"
#include "inpar_scatra.H"
#include "inpar_structure.H"

#include "drt_validparameters.H"
#include "../drt_lib/drt_conditiondefinition.H"


void INPAR::CELL::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using namespace INPAR::SCATRA;
  using namespace INPAR::STR;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;


  Teuchos::ParameterList& celldyn =
      list->sublist("CELL DYNAMIC", false, "Cell Migration Simulation");

  Teuchos::Tuple<std::string, 3> coupname;
  Teuchos::Tuple<int, 3> couplabel;

  coupname[0] = "basic_sequ_stagg";
  couplabel[0] = cell_basic_sequ_stagg;
  coupname[1] = "iter_stagg_fixed_rel_param";
  couplabel[1] = cell_iter_stagg_fixed_rel_param;
  coupname[2] = "iter_stagg_AITKEN_rel_param";
  couplabel[2] = cell_iter_stagg_AITKEN_rel_param;


  /* GENERAL PARAMETERS */

  setStringToIntegralParameter<int>("PSEUDO2D", "no",
      "True if a thin quasi 2-dimensional problem is modeled", tuple<std::string>("no", "yes"),
      tuple<int>(0, 1), &celldyn);

  setStringToIntegralParameter<int>("ARTIFICIAL_ECM_METHOD", "none",
      "How to treat artificial ECM solid phase", tuple<std::string>("none", "selective_assembly"),
      tuple<int>(0, 1), &celldyn);

  setStringToIntegralParameter<int>("ARTIFICIAL_FLUID_MOD", "simplified",
      "Full, simplified, or no modification of artificial flow",
      tuple<std::string>("none", "simplified", "full"),
      tuple<int>(mod_none, mod_simplified, mod_full), &celldyn);

  setStringToIntegralParameter<int>("SIMTYPE", "pureFSI", "Simulation Type",
      tuple<std::string>("pureFSI", "pureAdhesion", "pureCompression", "pureProtrusionFormation",
          "pureContraction", "pureEndoExocytosis", "Multiphysics"),
      tuple<int>(sim_type_pureFSI, sim_type_pureAdhesion, sim_type_pureConfinement,
          sim_type_pureProtrusionFormation, sim_type_pureContraction, sim_type_pureEndoExocytosis,
          sim_type_multiphysics),
      &celldyn);

  setStringToIntegralParameter<int>("MIGRATIONTYPE", "undefined",
      "Migration with or without ScaTra.",
      tuple<std::string>("undefined", "ameboid", "proteolytic"),
      tuple<int>(cell_migration_undefined, cell_migration_ameboid, cell_migration_proteolytic),
      &celldyn);

  setStringToIntegralParameter<int>("INITIALIZE_CELL", "no",
      "Use first time step as pre-simulation to initialize Cell", tuple<std::string>("no", "yes"),
      tuple<int>(0, 1), &celldyn);

  setStringToIntegralParameter<int>("ECM_INTERACTION", "no",
      "Switch on/off cell-ecm interaction model", tuple<std::string>("no", "yes"), tuple<int>(0, 1),
      &celldyn);

  setStringToIntegralParameter<int>("FLUID_INTERACTION", "no",
      "Switch on/off cell-fluid interaction model", tuple<std::string>("no", "yes"),
      tuple<int>(0, 1), &celldyn);

  setStringToIntegralParameter<int>("ADHESION_DYNAMICS", "no",
      "include adhesion dynamics in simulation", tuple<std::string>("no", "yes"), tuple<int>(0, 1),
      &celldyn);

  setStringToIntegralParameter<int>("SSI_CELL", "no",
      "build cell as ssi problem to enable intracellular and membrane transport",
      tuple<std::string>("no", "yes"), tuple<int>(0, 1), &celldyn);

  IntParameter("NUMSTEP", 200, "Total number of Timesteps", &celldyn);
  IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &celldyn);
  IntParameter("RESTARTEVRY", 9999, "Increment for writing restart", &celldyn);
  IntParameter(
      "INITIALIZATION_STEPS", 1, "Num of time steps for initialization (pre-simulation)", &celldyn);

  DoubleParameter("TIMESTEP", 0.1, "Time increment dt", &celldyn);
  DoubleParameter(
      "INITIAL_TIMESTEP", 0.1, "Time increment dt for first time step (pre-simulation)", &celldyn);
  DoubleParameter("MAXTIME", 1000.0, "Total simulation time", &celldyn);
  DoubleParameter(
      "kBT", 4.04530016e-6, "Thermal Energy (default [micrometer, mg, s] at 293K)", &celldyn);

  BoolParameter("DIFFTIMESTEPSIZE", "No", "use different step size for scatra and solid", &celldyn);


  /* GIVE DOF IDs OF DIFFERENT CHEMICAL SPECIES */

  Teuchos::ParameterList& scatradofdyn =
      celldyn.sublist("SCALAR TRANSPORT DOF IDS", false, "Control the dof ids of biomolecules");

  IntParameter(
      "ACTIN_MONOMER", -1, "ID of Actin monomer dof (start counting with 0)", &scatradofdyn);
  IntParameter("ROCK", -1, "ID of ROCK dof (start counting with 0)", &scatradofdyn);
  IntParameter(
      "RhoGEF", -1, "ID of Rho activating protein dof (start counting with 0)", &scatradofdyn);



  /* PROTEOLYSIS PARAMETERS */

  Teuchos::ParameterList& proteolysisdyn =
      celldyn.sublist("PROTEOLYSIS MODULE", false, "Control the models for proteolysis");

  setStringToIntegralParameter<int>("COUPALGO", "iter_stagg_fixed_rel_param",
      "Iteration Scheme over the fields", coupname, couplabel, &proteolysisdyn);

  setStringToIntegralParameter<int>("SEGREGATION", "undefined",
      "Segregation of chemical species by cell in volume or on surface",
      tuple<std::string>("undefined", "volume", "surface"),
      tuple<int>(segregation_undefined, segregation_volumetric, segregation_surface),
      &proteolysisdyn);

  setStringToIntegralParameter<int>("SEGREGATION_BY", "dirichlet",
      "Segregation via dirichlet or neumann condition",
      tuple<std::string>("undefined", "dirichlet", "neumann"),
      tuple<int>(segregation_by_undefined, segregation_by_dirichlet, segregation_by_neumann),
      &proteolysisdyn);

  setStringToIntegralParameter<int>("SEGREGATION_LAW", "undefined",
      "Segregation of chemical species by cell in volume or on surface",
      tuple<std::string>("undefined", "constant", "forcedependent"),
      tuple<int>(
          segregation_law_undefined, segregation_law_constant, segregation_law_forcedependent),
      &proteolysisdyn);

  DoubleParameter(
      "SEGREGATION_CONST", 0.0, "basic constant for segregation equation", &proteolysisdyn);



  /* CELL CONFINEMENT PARAMETERS */

  Teuchos::ParameterList& confinementdyn = celldyn.sublist("CONFINEMENT MODULE", false,
      "Control the models and algorithms for confinement of the cell to ecm void space");

  setStringToIntegralParameter<int>("COUPALGO", "iter_stagg_fixed_rel_param",
      "Iteration Scheme over the fields", coupname, couplabel, &confinementdyn);

  setStringToIntegralParameter<int>("COUPVARIABLE", "Displacement",
      "Coupling variable at the interface", tuple<std::string>("Displacement", "Force"),
      tuple<int>(0, 1), &confinementdyn);

  DoubleParameter("CONVTOL", 1e-6,
      "Tolerance for iteration over fields in case of partitioned scheme", &confinementdyn);
  DoubleParameter(
      "RELAX", 1.0, "fixed relaxation parameter for partitioned solver", &confinementdyn);
  DoubleParameter("MAXOMEGA", 0.0,
      "largest omega allowed for Aitken relaxation (0.0 means no constraint)", &confinementdyn);
  IntParameter("ITEMAX", 100, "Maximum number of iterations over fields", &confinementdyn);

  DoubleParameter(
      "PENALTY_INIT", 100.0, "Penalty parameter for the cell initialization", &confinementdyn);
  DoubleParameter("PENALTY_START", 100.0, "Penalty parameter at the beginning of each time step",
      &confinementdyn);
  DoubleParameter("ECM_FIBER_RADIUS", 0.3, "Average radius of ECM fibers", &confinementdyn);



  /* CELL ADHESION PARAMETERS */

  Teuchos::ParameterList& adhesiondyn = celldyn.sublist(
      "ADHESION MODULE", false, "Control the models and algorithms for cell-ecm adhesion");

  setStringToIntegralParameter<int>("COUPALGO", "iter_stagg_fixed_rel_param",
      "Iteration Scheme over the fields", coupname, couplabel, &adhesiondyn);


  setStringToIntegralParameter<int>("COUPMETHOD", "Dirichlet",
      "Coupling method at the adhesion interface", tuple<std::string>("Dirichlet", "Penalty"),
      tuple<int>(0, 1), &adhesiondyn);

  setStringToIntegralParameter<int>("COUPVARIABLE", "Force", "Coupling variable at the interface",
      tuple<std::string>("Displacement", "Force"), tuple<int>(0, 1), &adhesiondyn);


  DoubleParameter("CONVTOL", 1e-6,
      "Tolerance for iteration over fields in case of partitioned scheme", &adhesiondyn);
  DoubleParameter("RELAX", 1.0, "fixed relaxation parameter for partitioned solver", &adhesiondyn);
  DoubleParameter("MAXOMEGA", 0.0,
      "largest omega allowed for Aitken relaxation (0.0 means no constraint)", &adhesiondyn);
  IntParameter("ITEMAX", 100, "Maximum number of iterations over fields", &adhesiondyn);
  IntParameter("NUM_BOUNDSPECIES", -1, "Number of scalar fields which bind to ECM", &adhesiondyn);

  DoubleParameter("PENALTY", 1.0e6, "Penalty parameter for cell adhesion coupling", &adhesiondyn);
  DoubleParameter("RESET_CONC", 1.0e-10,
      "Minimum concentration defining adhesion site coupling (computational zero conc.)",
      &adhesiondyn);
  DoubleParameter(
      "ECM_FIBER_DIAMETER", -1.0, "Average ECM fiber diameter for Mikado model", &adhesiondyn);



  /* CELL PROTRUSION PARAMETERS */

  Teuchos::ParameterList& protdyn = celldyn.sublist(
      "PROTRUSION MODULE", false, "Control the models and algorithms for protrusion formation");

  setStringToIntegralParameter<int>("COUPALGO", "iter_stagg_fixed_rel_param",
      "Iteration Scheme over the fields", coupname, couplabel, &protdyn);

  setStringToIntegralParameter<int>("COUPVARIABLE", "Displacement",
      "Coupling variable at the interface", tuple<std::string>("Displacement", "Force"),
      tuple<int>(0, 1), &protdyn);

  setStringToIntegralParameter<int>("SSICOUPVARIABLE", "growth",
      "Check growth or usual ssi criterion", tuple<std::string>("undefined", "growth", "ssi"),
      tuple<int>(coup_growth_undefined, coup_growth_growth, coup_growth_ssi), &protdyn);

  DoubleParameter("CONVTOL", 1e-6,
      "Tolerance for iteration over fields in case of partitioned scheme", &protdyn);
  DoubleParameter("RELAX", 1.0, "fixed relaxation parameter for partitioned solver", &protdyn);
  DoubleParameter("MAXOMEGA", 0.0,
      "largest omega allowed for Aitken relaxation (0.0 means no constraint)", &protdyn);
  IntParameter("ITEMAX", 100, "Maximum number of iterations over fields", &protdyn);

  IntParameter("NUMDOF_ACTIN", -1, "Number of Scalar for Actin Monomer Concentration", &protdyn);
  IntParameter("NUMDOF_BARBEDENDS", -1, "Number of Scalar for Pointed End Concentration at surface",
      &protdyn);
  IntParameter("NUMDOF_BRANCHES", -1, "Number of Scalar for Branch Concentration", &protdyn);
  IntParameter("NUM_SURF_SCALARS", -1, "Total Number of Surface Scalars", &protdyn);
  IntParameter("NUM_FIL_ON_aBr", 3,
      "Number of filament barbed ends on reference length(2D)/area(3D) aBr", &protdyn);

  DoubleParameter("aBr", 37.1677531,
      "Minimal length (2D)/area (3D) between actin filament barbed ends\n"
      "above which branching can occur (default for simple 2D micro model).",
      &protdyn);
  DoubleParameter("RELAX_GROWTH", 1.0, "Fixed relaxation parameter for growth", &protdyn);
  DoubleParameter("K_ON", 10.0, "Actin polymerisation base rate", &protdyn);
  DoubleParameter("K_BR", 5.0, "Actin web branching base rate", &protdyn);
  DoubleParameter("ACTIN_MONOMER_SIZE", 2.7e-3, "Size of actin monomer", &protdyn);



  /* CELL-FLOW INTERACTION PARAMETERS */

  Teuchos::ParameterList& cfidyn = celldyn.sublist("FLOW INTERACTION MODULE", false,
      "Control the models and algorithms for cell-flow interaction");

  setStringToIntegralParameter<int>("COUPALGO", "iter_stagg_fixed_rel_param",
      "Iteration Scheme over the fields", coupname, couplabel, &cfidyn);


  setStringToIntegralParameter<int>("COUPVARIABLE", "Force", "Coupling variable at the interface",
      tuple<std::string>("Displacement", "Force"), tuple<int>(0, 1), &cfidyn);

  DoubleParameter("CONVTOL", 1e-6,
      "Tolerance for iteration over fields in case of partitioned scheme", &cfidyn);
  DoubleParameter("RELAX", 1.0, "fixed relaxation parameter for partitioned solver", &cfidyn);
  DoubleParameter("MAXOMEGA", 0.0,
      "largest omega allowed for Aitken relaxation (0.0 means no constraint)", &cfidyn);
  IntParameter("ITEMAX", 100, "Maximum number of iterations over fields", &cfidyn);



  /* ENDO-/EXOCYTOSIS PARAMETERS */


  Teuchos::ParameterList& endoexocytosisdyn =
      celldyn.sublist("ENDOEXOCYTOSIS MODULE", false, "Control the models for endo-/exocytosis");

  IntParameter("ENDOEXO_DELAY", 1,
      "Number of steps delay between internalization and begin of externalization",
      &endoexocytosisdyn);
  IntParameter(
      "EXO_DOF", 0, "Dof number of exocytosed species (starting with 0)", &endoexocytosisdyn);
  DoubleParameter("THETA", 1.0, "THETA use endo-/exocytosis", &endoexocytosisdyn);

  setStringToIntegralParameter<int>("EXODOMAIN", "surface", "exotysosis in surface of volume",
      tuple<std::string>("surface", "volume"), tuple<int>(exo_surface, exo_volume),
      &endoexocytosisdyn);


  /* CELL SCATRA PARAMS */

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& cellscatradyn = celldyn.sublist("SCALAR TRANSPORT", false,
      "control parameters for solving nonlinear SCATRA problems inside cell\n");

  setStringToIntegralParameter<int>("SOLVERTYPE", "nonlinear", "type of scalar transport solver",
      tuple<std::string>("linear_full", "linear_incremental", "nonlinear"),
      tuple<int>(solvertype_linear_full, solvertype_linear_incremental, solvertype_nonlinear),
      &cellscatradyn);

  setStringToIntegralParameter<int>("TIMEINTEGR", "One_Step_Theta", "Time Integration Scheme",
      tuple<std::string>("Stationary", "One_Step_Theta", "BDF2", "Gen_Alpha"),
      tuple<int>(timeint_stationary, timeint_one_step_theta, timeint_bdf2, timeint_gen_alpha),
      &cellscatradyn);

  DoubleParameter("THETA", 1.0, "One-step-theta time integration factor", &cellscatradyn);

  setStringToIntegralParameter<int>("VELOCITYFIELD", "Navier_Stokes",
      "type of velocity field used for scalar transport problems",
      tuple<std::string>("zero", "function", "Navier_Stokes"),
      tuple<int>(velocity_zero, velocity_function, velocity_Navier_Stokes), &cellscatradyn);

  IntParameter(
      "VELFUNCNO", -1, "function number for scalar transport velocity field", &cellscatradyn);

  {
    // a standard Teuchos::tuple can have at maximum 10 entries! We have to circumvent this here.
    Teuchos::Tuple<std::string, 12> name;
    Teuchos::Tuple<int, 12> label;
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

    setStringToIntegralParameter<int>("INITIALFIELD", "zero_field",
        "Initial Field for scalar transport problem", name, label, &cellscatradyn);
  }

  IntParameter(
      "INITFUNCNO", -1, "function number for scalar transport initial field", &cellscatradyn);

  BoolParameter("SPHERICALCOORDS", "No", "use of spherical coordinates", &cellscatradyn);

  setStringToIntegralParameter<int>("CALCERROR", "No",
      "compute error compared to analytical solution",
      tuple<std::string>("No", "Kwok_Wu", "ConcentricCylinders", "Electroneutrality",
          "error_by_function", "SphereDiffusion", "AnalyticSeries"),
      tuple<int>(calcerror_no, calcerror_Kwok_Wu, calcerror_cylinder, calcerror_electroneutrality,
          calcerror_byfunction, calcerror_spherediffusion, calcerror_AnalyticSeries),
      &cellscatradyn);

  IntParameter(
      "CALCERRORNO", -1, "function number for scalar transport error computation", &cellscatradyn);

  setStringToIntegralParameter<int>("CALCFLUX_DOMAIN", "No",
      "output of diffusive/total flux vectors inside domain",
      tuple<std::string>("No", "total", "diffusive"),
      tuple<int>(flux_none, flux_total, flux_diffusive), &cellscatradyn);

  setStringToIntegralParameter<int>("CALCFLUX_BOUNDARY", "No",
      "output of convective/diffusive/total flux vectors on boundary",
      tuple<std::string>("No", "total", "diffusive", "convective"),
      tuple<int>(flux_none, flux_total, flux_diffusive, flux_convective), &cellscatradyn);

  BoolParameter("CALCFLUX_DOMAIN_LUMPED", "Yes",
      "perform approximate domain flux calculation involving matrix lumping", &cellscatradyn);
  BoolParameter("CALCFLUX_BOUNDARY_LUMPED", "Yes",
      "perform approximate boundary flux calculation involving matrix lumping", &cellscatradyn);

  setNumericStringParameter("WRITEFLUX_IDS", "-1",
      "Write diffusive/total flux vector fields for these scalar fields only (starting with 1)",
      &cellscatradyn);

  setStringToIntegralParameter<int>("OUTPUTSCALARS", "none",
      "Output of total and mean values for transported scalars",
      tuple<std::string>("none", "entire_domain", "by_condition", "entire_domain_and_by_condition"),
      tuple<int>(outputscalars_none, outputscalars_entiredomain, outputscalars_condition,
          outputscalars_entiredomain_condition),
      &cellscatradyn);
  BoolParameter("OUTINTEGRREAC", "No", "Output of integral reaction values", &cellscatradyn);
  BoolParameter(
      "OUTPUT_GMSH", "No", "Do you want to write Gmsh postprocessing files?", &cellscatradyn);

  BoolParameter("MATLAB_STATE_OUTPUT", "No",
      "Do you want to write the state solution to Matlab file?", &cellscatradyn);

  setStringToIntegralParameter<int>("CONVFORM", "convective", "form of convective term",
      tuple<std::string>("convective", "conservative"),
      tuple<int>(convform_convective, convform_conservative), &cellscatradyn);

  BoolParameter("NEUMANNINFLOW", "no", "Flag to (de)activate potential Neumann inflow term(s)",
      &cellscatradyn);

  BoolParameter(
      "SKIPINITDER", "no", "Flag to skip computation of initial time derivative", &cellscatradyn);

  setStringToIntegralParameter<int>("MESHTYING", "Condensed_Smat",
      "Flag to (de)activate mesh tying algorithm",
      tuple<std::string>("no", "Condensed_Smat", "Condensed_Bmat",
          "Condensed_Bmat_merged"),  // use the condensed_bmat_merged strategy
      tuple<int>(INPAR::FLUID::no_meshtying, INPAR::FLUID::condensed_smat,
          INPAR::FLUID::condensed_bmat,
          INPAR::FLUID::condensed_bmat_merged),  // use the condensed_bmat_merged strategy
      &cellscatradyn);

  // linear solver id used for scalar transport/elch problems
  IntParameter("LINEAR_SOLVER", -1, "number of linear solver used for scalar transport/elch...",
      &cellscatradyn);
  // IntParameter("SIMPLER_SOLVER",-1,"number of linear solver used for ELCH (solved with
  // SIMPLER)...",&cellscatradyn);

  // flag for natural convection effects
  BoolParameter("NATURAL_CONVECTION", "No", "Include natural convection effects", &cellscatradyn);

  // parameters for finite difference check
  setStringToIntegralParameter<int>("FDCHECK", "none",
      "flag for finite difference check: none, local, or global",
      tuple<std::string>("none",
          "local",    // perform finite difference check on element level
          "global"),  // perform finite difference check on time integrator level
      tuple<int>(fdcheck_none, fdcheck_local, fdcheck_global), &cellscatradyn);
  DoubleParameter("FDCHECKEPS", 1.e-6,
      "dof perturbation magnitude for finite difference check (1.e-6 seems to work very well, "
      "whereas smaller values don't)",
      &cellscatradyn);
  DoubleParameter(
      "FDCHECKTOL", 1.e-6, "relative tolerance for finite difference check", &cellscatradyn);

  // parameter for optional computation of domain and boundary integrals, i.e., of surface areas and
  // volumes associated with specified nodesets
  setStringToIntegralParameter<int>("COMPUTEINTEGRALS", "none",
      "flag for optional computation of domain integrals",
      tuple<std::string>("none", "initial", "repeated"),
      tuple<int>(computeintegrals_none, computeintegrals_initial, computeintegrals_repeated),
      &cellscatradyn);
  // requested by base algorithm but not needed by cell migration
  setStringToIntegralParameter<int>("FSSUGRDIFF", "No", "fine-scale subgrid diffusivity",
      tuple<std::string>("No"), tuple<int>(fssugrdiff_no), &cellscatradyn);

  BoolParameter("CONV_HEAT_TRANS", "No",
      "Flag to (de)activate potential convective heat transfer boundary conditions",
      &cellscatradyn);

  // flag for output of performance statistics associated with linear solver into *.csv file
  BoolParameter("OUTPUTLINSOLVERSTATS", "No",
      "flag for output of performance statistics associated with linear solver into *.csv file",
      &cellscatradyn);

  // flag for output of performance statistics associated with nonlinear solver into *.csv file
  BoolParameter("OUTPUTNONLINSOLVERSTATS", "No",
      "flag for output of performance statistics associated with nonlinear solver into *.csv file",
      &cellscatradyn);

  // flag for point-based null space calculation
  BoolParameter(
      "NULLSPACE_POINTBASED", "No", "flag for point-based null space calculation", &cellscatradyn);

  // flag for adaptive time stepping
  BoolParameter("ADAPTIVE_TIMESTEPPING", "No", "flag for adaptive time stepping", &cellscatradyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& cellscatra_nonlin = cellscatradyn.sublist(
      "NONLINEAR", false, "control parameters for solving nonlinear CELL SCATRA problems\n");

  IntParameter("ITEMAX", 10, "max. number of nonlin. iterations", &cellscatra_nonlin);
  DoubleParameter("CONVTOL", 1e-6, "Tolerance for convergence check", &cellscatra_nonlin);
  IntParameter("ITEMAX_OUTER", 10,
      "Maximum number of outer iterations in partitioned coupling schemes (natural convection, "
      "multi-scale simulations etc.)",
      &cellscatra_nonlin);
  DoubleParameter("CONVTOL_OUTER", 1e-6,
      "Convergence check tolerance for outer loop in partitioned coupling schemes (natural "
      "convection, multi-scale simulations etc.)",
      &cellscatra_nonlin);
  BoolParameter("EXPLPREDICT", "no",
      "do an explicit predictor step before starting nonlinear iteration", &cellscatra_nonlin);
  DoubleParameter("ABSTOLRES", 1e-14,
      "Absolute tolerance for deciding if residual of nonlinear problem is already zero",
      &cellscatra_nonlin);

  // convergence criteria adaptivity
  BoolParameter("ADAPTCONV", "yes",
      "Switch on adaptive control of linear solver tolerance for nonlinear solution",
      &cellscatra_nonlin);
  DoubleParameter("ADAPTCONV_BETTER", 0.1,
      "The linear solver shall be this much better than the current nonlinear residual in the "
      "nonlinear convergence limit",
      &cellscatra_nonlin);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& cellscatradyn_stab = cellscatradyn.sublist("STABILIZATION", false,
      "control parameters for the stabilization of cell scalar transport problems");

  // this parameter governs type of stabilization
  setStringToIntegralParameter<int>("STABTYPE", "SUPG", "type of stabilization (if any)",
      tuple<std::string>("no_stabilization", "SUPG", "GLS", "USFEM"),
      tuple<std::string>(
          "Do not use any stabilization -> only reasonable for low-Peclet-number flows", "Use SUPG",
          "Use GLS", "Use USFEM"),
      tuple<int>(stabtype_no_stabilization, stabtype_SUPG, stabtype_GLS, stabtype_USFEM),
      &cellscatradyn_stab);

  // this parameter governs whether subgrid-scale velocity is included
  BoolParameter(
      "SUGRVEL", "no", "potential incorporation of subgrid-scale velocity", &cellscatradyn_stab);

  // this parameter governs whether all-scale subgrid diffusivity is included
  BoolParameter("ASSUGRDIFF", "no",
      "potential incorporation of all-scale subgrid diffusivity (a.k.a. discontinuity-capturing) "
      "term",
      &cellscatradyn_stab);

  // this parameter selects the tau definition applied
  setStringToIntegralParameter<int>("DEFINITION_TAU", "Franca_Valentin", "Definition of tau",
      tuple<std::string>("Taylor_Hughes_Zarins", "Taylor_Hughes_Zarins_wo_dt", "Franca_Valentin",
          "Franca_Valentin_wo_dt", "Shakib_Hughes_Codina", "Shakib_Hughes_Codina_wo_dt", "Codina",
          "Codina_wo_dt", "Franca_Madureira_Valentin", "Franca_Madureira_Valentin_wo_dt",
          "Exact_1D", "Zero"),
      tuple<int>(tau_taylor_hughes_zarins, tau_taylor_hughes_zarins_wo_dt, tau_franca_valentin,
          tau_franca_valentin_wo_dt, tau_shakib_hughes_codina, tau_shakib_hughes_codina_wo_dt,
          tau_codina, tau_codina_wo_dt, tau_franca_madureira_valentin,
          tau_franca_madureira_valentin_wo_dt, tau_exact_1d, tau_zero),
      &cellscatradyn_stab);

  // this parameter selects the characteristic element length for tau for all
  // stabilization parameter definitions requiring such a length
  setStringToIntegralParameter<int>("CHARELELENGTH", "streamlength",
      "Characteristic element length for tau",
      tuple<std::string>("streamlength", "volume_equivalent_diameter", "root_of_volume"),
      tuple<int>(streamlength, volume_equivalent_diameter, root_of_volume), &cellscatradyn_stab);

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
      &cellscatradyn_stab);

  // this parameter selects the location where tau is evaluated
  setStringToIntegralParameter<int>("EVALUATION_TAU", "integration_point",
      "Location where tau is evaluated", tuple<std::string>("element_center", "integration_point"),
      tuple<std::string>("evaluate tau at element center", "evaluate tau at integration point"),
      tuple<int>(evaltau_element_center, evaltau_integration_point), &cellscatradyn_stab);

  // this parameter selects the location where the material law is evaluated
  // (does not fit here very well, but parameter transfer is easier)
  setStringToIntegralParameter<int>("EVALUATION_MAT", "integration_point",
      "Location where material law is evaluated",
      tuple<std::string>("element_center", "integration_point"),
      tuple<std::string>(
          "evaluate material law at element center", "evaluate material law at integration point"),
      tuple<int>(evalmat_element_center, evalmat_integration_point), &cellscatradyn_stab);

  // this parameter selects methods for improving consistency of stabilization terms
  setStringToIntegralParameter<int>("CONSISTENCY", "no",
      "improvement of consistency for stabilization",
      tuple<std::string>("no", "L2_projection_lumped"),
      tuple<std::string>("inconsistent", "L2 projection with lumped mass matrix"),
      tuple<int>(consistency_no, consistency_l2_projection_lumped), &cellscatradyn_stab);



  /* CELL STRUCTURAL PARAMS */

  Teuchos::Array<std::string> yesnotuple =
      tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

  Teuchos::ParameterList& cellstructdyn =
      celldyn.sublist("STRUCTURAL DYNAMIC", false, "cell structure control parameters");

  /// copy already existing parameterlist instead of producing redundant code
  cellstructdyn = Teuchos::ParameterList(list->sublist("STRUCTURAL DYNAMIC", true));

}  // SetValidParameters


void INPAR::CELL::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // FOCAL ADHESION
  /*--------------------------------------------------------------------*/
  std::vector<Teuchos::RCP<ConditionComponent>> facomponents;

  facomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

  Teuchos::RCP<ConditionDefinition> linefa = Teuchos::rcp(
      new ConditionDefinition("DESIGN FOCAL ADHESION LINE CONDITIONS", "CellFocalAdhesion",
          "Cell Focal Adhesion", DRT::Condition::CellFocalAdhesion, true, DRT::Condition::Line));

  Teuchos::RCP<ConditionDefinition> surffa = Teuchos::rcp(
      new ConditionDefinition("DESIGN FOCAL ADHESION SURF CONDITIONS", "CellFocalAdhesion",
          "Cell Focal Adhesion", DRT::Condition::CellFocalAdhesion, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < facomponents.size(); ++i)
  {
    linefa->AddComponent(facomponents[i]);
    surffa->AddComponent(facomponents[i]);
  }

  condlist.push_back(linefa);
  condlist.push_back(surffa);

  /*--------------------------------------------------------------------*/
  // SURFACE FOR INTERNALIZATION AND EXTERNALIZATION
  /*--------------------------------------------------------------------*/
  /*--------------------------------------------------------------------*/
  // Boundary internalization condition for scalar transport
  std::vector<Teuchos::RCP<ConditionComponent>> internalizationcomponents;

  Teuchos::RCP<ConditionDefinition> lineinternalization =
      Teuchos::rcp(new ConditionDefinition("SCATRA CELL INTERNALIZATION LINE CONDITIONS",
          "ScaTraCellInt", "Scalar Transport Cell Internalization Calculation",
          DRT::Condition::ScaTraCellIntCalc, true, DRT::Condition::Line));

  Teuchos::RCP<ConditionDefinition> surfinternalization =
      Teuchos::rcp(new ConditionDefinition("SCATRA CELL INTERNALIZATION SURF CONDITIONS",
          "ScaTraCellInt", "Scalar Transport Cell Internalization Calculation",
          DRT::Condition::ScaTraCellIntCalc, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < internalizationcomponents.size(); ++i)
  {
    lineinternalization->AddComponent(internalizationcomponents[i]);
    surfinternalization->AddComponent(internalizationcomponents[i]);
  }


  condlist.push_back(lineinternalization);
  condlist.push_back(surfinternalization);


  /*--------------------------------------------------------------------*/
  // Boundary externalization condition for scalar transport
  std::vector<Teuchos::RCP<ConditionComponent>> externalizationcomponents;

  Teuchos::RCP<ConditionDefinition> lineexternalization =
      Teuchos::rcp(new ConditionDefinition("SCATRA CELL EXTERNALIZATION LINE CONDITIONS",
          "ScaTraCellExt", "Scalar Transport Cell Externalization Calculation",
          DRT::Condition::ScaTraCellExtCalc, true, DRT::Condition::Line));

  Teuchos::RCP<ConditionDefinition> surfexternalization =
      Teuchos::rcp(new ConditionDefinition("SCATRA CELL EXTERNALIZATION SURF CONDITIONS",
          "ScaTraCellExt", "Scalar Transport Cell Externalization Calculation",
          DRT::Condition::ScaTraCellExtCalc, true, DRT::Condition::Surface));

  Teuchos::RCP<ConditionDefinition> volexternalization =
      Teuchos::rcp(new ConditionDefinition("SCATRA CELL EXTERNALIZATION VOL CONDITIONS",
          "ScaTraCellExt", "Scalar Transport Cell Externalization Calculation",
          DRT::Condition::ScaTraCellExtCalc, true, DRT::Condition::Volume));

  externalizationcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ScalarID")));
  externalizationcomponents.push_back(Teuchos::rcp(new IntConditionComponent("ScalarID")));

  for (unsigned i = 0; i < externalizationcomponents.size(); ++i)
  {
    lineexternalization->AddComponent(externalizationcomponents[i]);
    surfexternalization->AddComponent(externalizationcomponents[i]);
    volexternalization->AddComponent(externalizationcomponents[i]);
  }

  condlist.push_back(lineexternalization);
  condlist.push_back(surfexternalization);
  condlist.push_back(volexternalization);
}
