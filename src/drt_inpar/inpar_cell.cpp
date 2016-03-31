/*----------------------------------------------------------------------*/
/*!
\file inpar_cell.cpp

<pre>
\maintainer Andreas Rauch
            rauch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 -15240
</pre>
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
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;


  Teuchos::ParameterList& celldyn = list->sublist("CELL DYNAMIC",false,
                                                  "Cell Migration Simulation");

  Teuchos::Tuple<std::string,3> coupname;
  Teuchos::Tuple<int,3> couplabel;

  coupname[ 0] = "basic_sequ_stagg";                              couplabel[ 0] = cell_basic_sequ_stagg;
  coupname[ 1] = "iter_stagg_fixed_rel_param";                    couplabel[ 1] = cell_iter_stagg_fixed_rel_param;
  coupname[ 2] = "iter_stagg_AITKEN_rel_param";                   couplabel[ 2] = cell_iter_stagg_AITKEN_rel_param;


  /* GENERAL PARAMETERS */

  setStringToIntegralParameter<int>(
                                    "PSEUDO2D","no",
                                    "True if a thin quasi 2-dimensional problem is modeled",
                                    tuple<std::string>(
                                                       "no",
                                                       "yes"),
                                    tuple<int>(
                                               0,
                                               1),
                                    &celldyn);

  setStringToIntegralParameter<int>(
                                    "SIMTYPE","pureFSI",
                                    "Simulation Type",
                                    tuple<std::string>(
                                                       "pureFSI",
                                                       "pureAdhesion",
                                                       "pureCompression",
                                                       "pureProtrusionFormation",
                                                       "pureContraction",
                                                       "Multiphysics"),
                                    tuple<int>(
                                               sim_type_pureFSI,
                                               sim_type_pureAdhesion,
                                               sim_type_pureConfinement,
                                               sim_type_pureProtrusionFormation,
                                               sim_type_pureContraction,
                                               sim_type_multiphysics),
                                    &celldyn);

  setStringToIntegralParameter<int>(
                                    "MIGRATIONTYPE","undefined",
                                    "Migration with or without ScaTra.",
                                    tuple<std::string>(
                                                       "undefined",
                                                       "ameboid",
                                                       "proteolytic"),
                                    tuple<int>(
                                               cell_migration_undefined,
                                               cell_migration_ameboid,
                                               cell_migration_proteolytic),
                                    &celldyn);

  setStringToIntegralParameter<int>(
                                    "INITIALIZE_CELL","no",
                                    "Use first time step as pre-simulation to initialize Cell",
                                    tuple<std::string>(
                                                       "no",
                                                       "yes"),
                                    tuple<int>(
                                               0,
                                               1),
                                    &celldyn);

  setStringToIntegralParameter<int>(
                                    "ECM_INTERACTION","no",
                                    "Switch on/off cell-ecm interaction model",
                                    tuple<std::string>(
                                                       "no",
                                                       "yes"),
                                    tuple<int>(
                                               0,
                                               1),
                                    &celldyn);

  setStringToIntegralParameter<int>(
                                    "FLUID_INTERACTION","yes",
                                    "Switch on/off cell-fluid interaction model",
                                    tuple<std::string>(
                                                       "no",
                                                       "yes"),
                                    tuple<int>(
                                               0,
                                               1),
                                    &celldyn);

  setStringToIntegralParameter<int>(
                                    "ADHESION_DYNAMICS","no",
                                    "include adhesion dynamics in simulation",
                                    tuple<std::string>(
                                                       "no",
                                                       "yes"),
                                    tuple<int>(
                                               0,
                                               1),
                                    &celldyn);

  setStringToIntegralParameter<int>(
                                    "SSI_CELL","no",
                                    "build cell as ssi problem to enable intracellular and membrane transport",
                                    tuple<std::string>(
                                                       "no",
                                                       "yes"),
                                    tuple<int>(
                                               0,
                                               1),
                                    &celldyn);

  IntParameter("NUMSTEP",200,"Total number of Timesteps",&celldyn);
  IntParameter("RESULTSEVRY",1,"Increment for writing solution",&celldyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&celldyn);
  IntParameter("INITIALIZATION_STEPS",1,"Num of time steps for initialization (pre-simulation)",&celldyn);

  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&celldyn);
  DoubleParameter("INITIAL_TIMESTEP",0.1,"Time increment dt for first time step (pre-simulation)",&celldyn);
  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&celldyn);




  /* PROTEOLYSIS PARAMETERS */


  Teuchos::ParameterList& proteolysisdyn = celldyn.sublist("PROTEOLYSIS MODULE",false,
                                                  "Control the models for proteolysis");

  setStringToIntegralParameter<int>("COUPALGO","iter_stagg_fixed_rel_param",
                                    "Iteration Scheme over the fields",
                                    coupname,
                                    couplabel,
                                    &proteolysisdyn);

  setStringToIntegralParameter<int>(
                                    "SEGREGATION","undefined",
                                    "Segregation of chemical species by cell in volume or on surface",
                                    tuple<std::string>(
                                                       "undefined",
                                                       "volume",
                                                       "surface"),
                                    tuple<int>(
                                               segregation_undefined,
                                               segregation_volumetric,
                                               segregation_surface),
                                    &proteolysisdyn);

  setStringToIntegralParameter<int>(
                                    "SEGREGATION_BY","dirichlet",
                                    "Segregation via dirichlet or neumann condition",
                                    tuple<std::string>(
                                                       "undefined",
                                                       "dirichlet",
                                                       "neumann"),
                                    tuple<int>(
                                               segregation_by_undefined,
                                               segregation_by_dirichlet,
                                               segregation_by_neumann),
                                    &proteolysisdyn);

  setStringToIntegralParameter<int>(
                                    "SEGREGATION_LAW","undefined",
                                    "Segregation of chemical species by cell in volume or on surface",
                                    tuple<std::string>(
                                                       "undefined",
                                                       "constant",
                                                       "forcedependent"),
                                    tuple<int>(
                                               segregation_law_undefined,
                                               segregation_law_constant,
                                               segregation_law_forcedependent),
                                    &proteolysisdyn);

  DoubleParameter("SEGREGATION_CONST",0.0,"basic constant for segregation equation",&proteolysisdyn);



  /* CELL CONFINEMENT PARAMETERS */

  Teuchos::ParameterList& confinementdyn = celldyn.sublist("CONFINEMENT MODULE",false,
                                                  "Control the models and algorithms for confinement of the cell to ecm void space");

  setStringToIntegralParameter<int>("COUPALGO","iter_stagg_fixed_rel_param",
                                    "Iteration Scheme over the fields",
                                    coupname,
                                    couplabel,
                                    &confinementdyn);

  DoubleParameter("PENALTY_INIT",100.0,"Penalty parameter for the cell initialization",&confinementdyn);
  DoubleParameter("PENALTY_START",100.0,"Penalty parameter at the beginning of each time step",&confinementdyn);
  DoubleParameter("ECM_FIBER_RADIUS",0.3,"Average radius of ECM fibers",&confinementdyn);



  /* CELL ADHESION PARAMETERS */

  Teuchos::ParameterList& adhesiondyn = celldyn.sublist("ADHESION MODULE",false,
                                                  "Control the models and algorithms for cell-ecm adhesion");

  setStringToIntegralParameter<int>("COUPALGO","iter_stagg_fixed_rel_param",
                                    "Iteration Scheme over the fields",
                                    coupname,
                                    couplabel,
                                    &adhesiondyn);



  /* CELL PROTRUSION PARAMETERS */

  Teuchos::ParameterList& protdyn = celldyn.sublist("PROTRUSION MODULE",false,
                                                  "Control the models and algorithms for protrusion formation");

  setStringToIntegralParameter<int>("COUPALGO","iter_stagg_fixed_rel_param",
                                    "Iteration Scheme over the fields",
                                    coupname,
                                    couplabel,
                                    &protdyn);



  /* CELL-FLOW INTERACTION PARAMETERS */

  Teuchos::ParameterList& cfidyn = celldyn.sublist("FLOW INTERACTION MODULE",false,
                                                  "Control the models and algorithms for cell-flow interaction");

  setStringToIntegralParameter<int>("COUPALGO","iter_stagg_fixed_rel_param",
                                    "Iteration Scheme over the fields",
                                    coupname,
                                    couplabel,
                                    &cfidyn);



  /* CELL SCATRA PARAMS */

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& cellscatradyn = celldyn.sublist(
      "SCALAR TRANSPORT",
      false,
      "control parameters for solving nonlinear SCATRA problems inside cell\n");

  setStringToIntegralParameter<int>("SOLVERTYPE","nonlinear",
                               "type of scalar transport solver",
                               tuple<std::string>(
                                 "linear_full",
                                 "linear_incremental",
                                 "nonlinear"
                                 ),
                               tuple<int>(
                                   solvertype_linear_full,
                                   solvertype_linear_incremental,
                                   solvertype_nonlinear),
                               &cellscatradyn);

  setStringToIntegralParameter<int>("TIMEINTEGR","One_Step_Theta",
                               "Time Integration Scheme",
                               tuple<std::string>(
                                 "Stationary",
                                 "One_Step_Theta",
                                 "BDF2",
                                 "Gen_Alpha"
                                 ),
                               tuple<int>(
                                   timeint_stationary,
                                   timeint_one_step_theta,
                                   timeint_bdf2,
                                   timeint_gen_alpha
                                 ),
                               &cellscatradyn);

  DoubleParameter("THETA",1.0,"One-step-theta time integration factor",&cellscatradyn);

  setStringToIntegralParameter<int>("VELOCITYFIELD","Navier_Stokes",
                               "type of velocity field used for scalar transport problems",
                               tuple<std::string>(
                                 "zero",
                                 "function",
                                 "function_and_curve",
                                 "Navier_Stokes"
                                 ),
                               tuple<int>(
                                   velocity_zero,
                                   velocity_function,
                                   velocity_function_and_curve,
                                   velocity_Navier_Stokes),
                               &cellscatradyn);

  IntParameter("VELFUNCNO",-1,"function number for scalar transport velocity field",&cellscatradyn);

  IntParameter("VELCURVENO",-1,"curve number for time-dependent scalar transport velocity field",&cellscatradyn);

  {
    // a standard Teuchos::tuple can have at maximum 10 entries! We have to circumvent this here.
    Teuchos::Tuple<std::string,12> name;
    Teuchos::Tuple<int,12> label;
    name[ 0] = "zero_field";                   label[ 0] = initfield_zero_field;
    name[ 1] = "field_by_function";            label[ 1] = initfield_field_by_function;
    name[ 2] = "field_by_condition";           label[ 2] = initfield_field_by_condition;
    name[ 3] = "disturbed_field_by_function";  label[ 3] = initfield_disturbed_field_by_function;
    name[ 4] = "1D_DISCONTPV";                 label[ 4] = initfield_discontprogvar_1D;
    name[ 5] = "FLAME_VORTEX_INTERACTION";     label[ 5] = initfield_flame_vortex_interaction;
    name[ 6] = "RAYTAYMIXFRAC";                label[ 6] = initfield_raytaymixfrac;
    name[ 7] = "L_shaped_domain";              label[ 7] = initfield_Lshapeddomain;
    name[ 8] = "facing_flame_fronts";          label[ 8] = initfield_facing_flame_fronts;
    name[ 9] = "oracles_flame";                label[ 9] = initfield_oracles_flame;
    name[10] = "high_forced_hit";              label[10] = initialfield_forced_hit_high_Sc;
    name[11] = "low_forced_hit";               label[11] = initialfield_forced_hit_low_Sc;

    setStringToIntegralParameter<int>(
        "INITIALFIELD",
        "zero_field",
        "Initial Field for scalar transport problem",
        name,
        label,
        &cellscatradyn);
  }

  IntParameter("INITFUNCNO",-1,"function number for scalar transport initial field",&cellscatradyn);

  setStringToIntegralParameter<int>("WRITEFLUX","No","output of diffusive/total flux vectors",
                               tuple<std::string>(
                                 "No",
                                 "totalflux_domain",
                                 "diffusiveflux_domain",
                                 "totalflux_boundary",
                                 "diffusiveflux_boundary",
                                 "convectiveflux_boundary"
                                 ),
                               tuple<int>(
                                   flux_no,
                                   flux_total_domain,
                                   flux_diffusive_domain,
                                   flux_total_boundary,
                                   flux_diffusive_boundary,
                                   flux_convective_boundary),
                               &cellscatradyn);


  setNumericStringParameter("WRITEFLUX_IDS","-1",
      "Write diffusive/total flux vector fields for these scalar fields only (starting with 1)",
      &cellscatradyn);

  BoolParameter("OUTPUTSCALARS","No","Output of total and mean values for transported scalars",&cellscatradyn);
  BoolParameter("OUTINTEGRREAC","No","Output of integral reaction values",&cellscatradyn);
  BoolParameter("OUTPUT_GMSH","No","Do you want to write Gmsh postprocessing files?",&cellscatradyn);

  BoolParameter("MATLAB_STATE_OUTPUT","No","Do you want to write the state solution to Matlab file?",&cellscatradyn);

  setStringToIntegralParameter<int>("CONVFORM","convective","form of convective term",
                               tuple<std::string>(
                                 "convective",
                                 "conservative"
                                 ),
                               tuple<int>(
                                 convform_convective,
                                 convform_conservative),
                               &cellscatradyn);

  BoolParameter("NEUMANNINFLOW",
      "no","Flag to (de)activate potential Neumann inflow term(s)",&cellscatradyn);

  BoolParameter("SKIPINITDER",
      "no","Flag to skip computation of initial time derivative",&cellscatradyn);

  setStringToIntegralParameter<int>("MESHTYING", "Condensed_Smat", "Flag to (de)activate mesh tying algorithm",
                                  tuple<std::string>(
                                      "no",
                                      "Condensed_Smat",
                                      "Condensed_Bmat",
                                      "Condensed_Bmat_merged"), //use the condensed_bmat_merged strategy
                                    tuple<int>(
                                        INPAR::FLUID::no_meshtying,
                                        INPAR::FLUID::condensed_smat,
                                        INPAR::FLUID::condensed_bmat,
                                        INPAR::FLUID::condensed_bmat_merged),   //use the condensed_bmat_merged strategy
                                    &cellscatradyn);

  // linear solver id used for scalar transport/elch problems
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for scalar transport/elch...",&cellscatradyn);
  //IntParameter("SIMPLER_SOLVER",-1,"number of linear solver used for ELCH (solved with SIMPLER)...",&cellscatradyn);

  // parameters for natural convection effects
  BoolParameter("NATURAL_CONVECTION","No","Include natural convection effects",&cellscatradyn);
  IntParameter("NATCONVITEMAX",10,"Maximum number of outer iterations for natural convection",&cellscatradyn);
  DoubleParameter("NATCONVCONVTOL",1e-6,"Convergence check tolerance for outer loop for natural convection",&cellscatradyn);

  // parameters for finite difference check
  setStringToIntegralParameter<int>("FDCHECK", "none", "flag for finite difference check: none, local, or global",
                                    tuple<std::string>(
                                      "none",
                                      "local",    // perform finite difference check on element level
                                      "global"),  // perform finite difference check on time integrator level
                                    tuple<int>(
                                        fdcheck_none,
                                        fdcheck_local,
                                        fdcheck_global),
                                    &cellscatradyn);
  DoubleParameter("FDCHECKEPS",1.e-6,"dof perturbation magnitude for finite difference check (1.e-6 seems to work very well, whereas smaller values don't)",&cellscatradyn);
  DoubleParameter("FDCHECKTOL",1.e-6,"relative tolerance for finite difference check",&cellscatradyn);

  // parameter for optional computation of domain and boundary integrals, i.e., of surface areas and volumes associated with specified nodesets
  setStringToIntegralParameter<int>(
      "COMPUTEINTEGRALS",
      "none",
      "flag for optional computation of domain integrals",
      tuple<std::string>(
          "none",
          "initial",
          "repeated"
          ),
      tuple<int>(
          computeintegrals_none,
          computeintegrals_initial,
          computeintegrals_repeated
          ),
      &cellscatradyn
      );
  // requested by base algorithm but not needed by cell migration
  setStringToIntegralParameter<int>("FSSUGRDIFF",
                               "No",
                               "fine-scale subgrid diffusivity",
                               tuple<std::string>(
                                 "No"),
                               tuple<int>(
                                   fssugrdiff_no),
                               &cellscatradyn);

  BoolParameter("CONV_HEAT_TRANS",
      "No","Flag to (de)activate potential convective heat transfer boundary conditions",&cellscatradyn);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& cellscatra_nonlin = cellscatradyn.sublist(
      "NONLINEAR",
      false,
      "control parameters for solving nonlinear CELL SCATRA problems\n");

  IntParameter("ITEMAX",10,"max. number of nonlin. iterations",&cellscatra_nonlin);
  DoubleParameter("CONVTOL",1e-6,"Tolerance for convergence check",&cellscatra_nonlin);
  BoolParameter("EXPLPREDICT","no","do an explicit predictor step before starting nonlinear iteration",&cellscatra_nonlin);
  DoubleParameter("ABSTOLRES",1e-14,"Absolute tolerance for deciding if residual of nonlinear problem is already zero",&cellscatra_nonlin);

  // convergence criteria adaptivity
  BoolParameter("ADAPTCONV","yes","Switch on adaptive control of linear solver tolerance for nonlinear solution",&cellscatra_nonlin);
  DoubleParameter("ADAPTCONV_BETTER",0.1,"The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",&cellscatra_nonlin);

/*----------------------------------------------------------------------*/
  Teuchos::ParameterList& cellscatradyn_stab = cellscatradyn.sublist("STABILIZATION",false,"control parameters for the stabilization of cell scalar transport problems");

  // this parameter governs type of stabilization
  setStringToIntegralParameter<int>("STABTYPE",
                                    "SUPG",
                                    "type of stabilization (if any)",
                               tuple<std::string>(
                                 "no_stabilization",
                                 "SUPG",
                                 "GLS",
                                 "USFEM"),
                               tuple<std::string>(
                                 "Do not use any stabilization -> only reasonable for low-Peclet-number flows",
                                 "Use SUPG",
                                 "Use GLS",
                                 "Use USFEM")  ,
                               tuple<int>(
                                   stabtype_no_stabilization,
                                   stabtype_SUPG,
                                   stabtype_GLS,
                                   stabtype_USFEM),
                               &cellscatradyn_stab);

  // this parameter governs whether subgrid-scale velocity is included
  BoolParameter("SUGRVEL","no","potential incorporation of subgrid-scale velocity",&cellscatradyn_stab);

  // this parameter governs whether all-scale subgrid diffusivity is included
  BoolParameter("ASSUGRDIFF","no",
      "potential incorporation of all-scale subgrid diffusivity (a.k.a. discontinuity-capturing) term",&cellscatradyn_stab);

  // this parameter selects the tau definition applied
  setStringToIntegralParameter<int>("DEFINITION_TAU",
                               "Franca_Valentin",
                               "Definition of tau",
                               tuple<std::string>(
                                 "Taylor_Hughes_Zarins",
                                 "Taylor_Hughes_Zarins_wo_dt",
                                 "Franca_Valentin",
                                 "Franca_Valentin_wo_dt",
                                 "Shakib_Hughes_Codina",
                                 "Shakib_Hughes_Codina_wo_dt",
                                 "Codina",
                                 "Codina_wo_dt",
                                 "Franca_Madureira_Valentin",
                                 "Franca_Madureira_Valentin_wo_dt",
                                 "Exact_1D",
                                 "Zero"),
                                tuple<int>(
                                    tau_taylor_hughes_zarins,
                                    tau_taylor_hughes_zarins_wo_dt,
                                    tau_franca_valentin,
                                    tau_franca_valentin_wo_dt,
                                    tau_shakib_hughes_codina,
                                    tau_shakib_hughes_codina_wo_dt,
                                    tau_codina,
                                    tau_codina_wo_dt,
                                    tau_franca_madureira_valentin,
                                    tau_franca_madureira_valentin_wo_dt,
                                    tau_exact_1d,
                                    tau_zero),
                               &cellscatradyn_stab);

  // this parameter selects the characteristic element length for tau for all
  // stabilization parameter definitions requiring such a length
  setStringToIntegralParameter<int>("CHARELELENGTH",
                               "streamlength",
                               "Characteristic element length for tau",
                               tuple<std::string>(
                                 "streamlength",
                                 "volume_equivalent_diameter",
                                 "root_of_volume"),
                               tuple<int>(
                                   streamlength,
                                   volume_equivalent_diameter,
                                   root_of_volume),
                               &cellscatradyn_stab);

  // this parameter selects the all-scale subgrid-diffusivity definition applied
  setStringToIntegralParameter<int>("DEFINITION_ASSGD",
                               "artificial_linear",
                               "Definition of (all-scale) subgrid diffusivity",
                               tuple<std::string>(
                                 "artificial_linear",
                                 "artificial_linear_reinit",
                                 "Hughes_etal_86_nonlinear",
                                 "Tezduyar_Park_86_nonlinear",
                                 "Tezduyar_Park_86_nonlinear_wo_phizero",
                                 "doCarmo_Galeao_91_nonlinear",
                                 "Almeida_Silva_97_nonlinear",
                                 "YZbeta_nonlinear",
                                 "Codina_nonlinear"),
                               tuple<std::string>(
                                 "classical linear artificial subgrid-diffusivity",
                                 "simple linear artificial subgrid-diffusivity const*h",
                                 "nonlinear isotropic according to Hughes et al. (1986)",
                                 "nonlinear isotropic according to Tezduyar and Park (1986)",
                                 "nonlinear isotropic according to Tezduyar and Park (1986) without user parameter phi_zero",
                                 "nonlinear isotropic according to doCarmo and Galeao (1991)",
                                 "nonlinear isotropic according to Almeida and Silva (1997)",
                                 "nonlinear YZ beta model",
                                 "nonlinear isotropic according to Codina")  ,
                                tuple<int>(
                                    assgd_artificial,
                                    assgd_lin_reinit,
                                    assgd_hughes,
                                    assgd_tezduyar,
                                    assgd_tezduyar_wo_phizero,
                                    assgd_docarmo,
                                    assgd_almeida,
                                    assgd_yzbeta,
                                    assgd_codina),
                               &cellscatradyn_stab);

  // this parameter selects the location where tau is evaluated
  setStringToIntegralParameter<int>("EVALUATION_TAU",
                               "integration_point",
                               "Location where tau is evaluated",
                               tuple<std::string>(
                                 "element_center",
                                 "integration_point"),
                               tuple<std::string>(
                                 "evaluate tau at element center",
                                 "evaluate tau at integration point")  ,
                                tuple<int>(
                                  evaltau_element_center,
                                  evaltau_integration_point),
                               &cellscatradyn_stab);

  // this parameter selects the location where the material law is evaluated
  // (does not fit here very well, but parameter transfer is easier)
  setStringToIntegralParameter<int>("EVALUATION_MAT",
                               "integration_point",
                               "Location where material law is evaluated",
                               tuple<std::string>(
                                 "element_center",
                                 "integration_point"),
                               tuple<std::string>(
                                 "evaluate material law at element center",
                                 "evaluate material law at integration point"),
                               tuple<int>(
                                 evalmat_element_center,
                                 evalmat_integration_point),
                               &cellscatradyn_stab);

  // this parameter selects methods for improving consistency of stabilization terms
  setStringToIntegralParameter<int>("CONSISTENCY",
                               "no",
                               "improvement of consistency for stabilization",
                               tuple<std::string>(
                                 "no",
                                 "L2_projection_lumped"),
                               tuple<std::string>(
                                 "inconsistent",
                                 "L2 projection with lumped mass matrix")  ,
                                tuple<int>(
                                  consistency_no,
                                  consistency_l2_projection_lumped),
                               &cellscatradyn_stab);



  /* CELL STRUCTURAL PARAMS */

  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  Teuchos::ParameterList& cellstructdyn = celldyn.sublist("STRUCTURAL DYNAMIC",false,"cell structure control parameters");

  setStringToIntegralParameter<int>("INT_STRATEGY","Old",
                               "global type of the used integration strategy",
                               tuple<std::string>(
                                 "Old",
                                 "Standard",
                                 "LOCA"),
                               tuple<int>(
                                 int_old,
                                 int_standard,
                                 int_loca),
                               &cellstructdyn);

  setStringToIntegralParameter<int>("DYNAMICTYP","OneStepTheta",
                               "type of the specific dynamic time integration scheme",
                               tuple<std::string>(
                                 "Statics",
                                 "GenAlpha",
                                 "OneStepTheta",
                                 "GEMM",
                                 "ExplicitEuler",
                                 "CentrDiff",
                                 "AdamsBashforth2",
                                 "EulerMaruyama",
                                 "EulerImpStoch",
                                 "StatMech"),
                               tuple<int>(
                                 dyna_statics,
                                 dyna_genalpha,
                                 dyna_onesteptheta,
                                 dyna_gemm,
                                 dyna_expleuler,
                                 dyna_centrdiff,
                                 dyna_ab2,
                                 dyna_euma,
                                 dyna_euimsto,
                                 dyna_statmech),
                               &cellstructdyn);

  setStringToIntegralParameter<int>("PRESTRESS","none","prestressing takes values none mulf id",
                               tuple<std::string>("none","None","NONE",
                                                  "mulf","Mulf","MULF",
                                                  "id","Id","ID"),
                               tuple<int>(prestress_none,prestress_none,prestress_none,
                                                            prestress_mulf,prestress_mulf,prestress_mulf,
                                                            prestress_id,prestress_id,prestress_id),
                               &cellstructdyn);

  DoubleParameter("PRESTRESSTIME",0.0,"time to switch from pre to post stressing",&cellstructdyn);

  // Output type
  IntParameter("RESEVRYERGY",0,"write system energies every requested step",&cellstructdyn);

  // Damping
  setStringToIntegralParameter<int>("DAMPING","No",
                               "type of damping: (1) Rayleigh damping matrix and use it from M_DAMP x M + K_DAMP x K, (2) Material based and calculated in elements",
                               tuple<std::string>(
                                 "no",
                                 "No",
                                 "NO",
                                 "yes",
                                 "Yes",
                                 "YES",
                                 "Rayleigh",
                                 "Material",
                                 "BrownianMotion"),
                               tuple<int>(
                                 damp_none,
                                 damp_none,
                                 damp_none,
                                 damp_rayleigh,
                                 damp_rayleigh,
                                 damp_rayleigh,
                                 damp_rayleigh,
                                 damp_material,
                                 damp_brownianmotion),
                               &cellstructdyn);
  DoubleParameter("M_DAMP",-1.0,"",&cellstructdyn);
  DoubleParameter("K_DAMP",-1.0,"",&cellstructdyn);

  DoubleParameter("TOLDISP",1.0E-10,
                  "tolerance in the displacement norm for the newton iteration",
                  &cellstructdyn);
  setStringToIntegralParameter<int>("NORM_DISP","Abs","type of norm for displacement convergence check",
                               tuple<std::string>(
                                 "Abs",
                                 "Rel",
                                 "Mix"),
                               tuple<int>(
                                 convnorm_abs,
                                 convnorm_rel,
                                 convnorm_mix),
                               &cellstructdyn);

  DoubleParameter("TOLRES",1.0E-08,
                  "tolerance in the residual norm for the newton iteration",
                  &cellstructdyn);
  setStringToIntegralParameter<int>("NORM_RESF","Abs","type of norm for residual convergence check",
                               tuple<std::string>(
                                 "Abs",
                                 "Rel",
                                 "Mix"),
                               tuple<int>(
                                 convnorm_abs,
                                 convnorm_rel,
                                 convnorm_mix),
                               &cellstructdyn);

  DoubleParameter("TOLPRE",1.0E-08,
                  "tolerance in pressure norm for the newton iteration",
                  &cellstructdyn);
  setStringToIntegralParameter<int>("NORM_PRES","Abs","type of norm for pressure convergence check",
                               tuple<std::string>(
                                 "Abs"),
                               tuple<int>(
                                 convnorm_abs),
                               &cellstructdyn);

  DoubleParameter("TOLINCO",1.0E-08,
                  "tolerance in the incompressible residual norm for the newton iteration",
                  &cellstructdyn);
  setStringToIntegralParameter<int>("NORM_INCO","Abs","type of norm for incompressible residual convergence check",
                               tuple<std::string>(
                                 "Abs"),
                               tuple<int>(
                                 convnorm_abs),
                               &cellstructdyn);

  setStringToIntegralParameter<int>("NORMCOMBI_DISPPRES","And","binary operator to combine pressure and displacement values",
                               tuple<std::string>(
                                 "And",
                                 "Or"),
                               tuple<int>(
                                 bop_and,
                                 bop_or),
                               &cellstructdyn);

  setStringToIntegralParameter<int>("NORMCOMBI_RESFINCO","And","binary operator to combine force and incompressible residual",
                               tuple<std::string>(
                                 "And",
                                 "Or"),
                               tuple<int>(
                                 bop_and,
                                 bop_or),
                               &cellstructdyn);

  setStringToIntegralParameter<int>("NORMCOMBI_RESFDISP","And","binary operator to combine displacement and residual force values",
                               tuple<std::string>(
                                 "And",
                                 "Or"),
                               tuple<int>(
                                 bop_and,
                                 bop_or),
                               &cellstructdyn);

  setStringToIntegralParameter<int>("STC_SCALING","no",
      "Scaled director conditioning for thin shell structures",
      tuple<std::string>(
        "no",
        "No",
        "NO",
        "Symmetric",
        "Right"),
      tuple<int>(
        stc_none,
        stc_none,
        stc_none,
        stc_currsym,
        stc_curr),
      &cellstructdyn);

  IntParameter("STC_LAYER",1,
               "number of STC layers for multilayer case",
               &cellstructdyn);

  DoubleParameter("PTCDT",0.1,
                  "pseudo time step for pseudo transient continuation (PTC) stabilized Newton procedure",
                  &cellstructdyn);

  DoubleParameter("TOLCONSTR",1.0E-08,
                  "tolerance in the constr error norm for the newton iteration",
                  &cellstructdyn);

  IntParameter("MAXITER",50,
               "maximum number of iterations allowed for Newton-Raphson iteration before failure",
               &cellstructdyn);
  IntParameter("MINITER",0,
               "minimum number of iterations to be done within Newton-Raphson loop",
               &cellstructdyn);
  setStringToIntegralParameter<int>("ITERNORM","L2","type of norm to be applied to residuals",
                               tuple<std::string>(
                                 "L1",
                                 "L2",
                                 "Rms",
                                 "Inf"),
                               tuple<int>(
                                 norm_l1,
                                 norm_l2,
                                 norm_rms,
                                 norm_inf),
                               &cellstructdyn);

  setStringToIntegralParameter<int>("DIVERCONT","stop","What to do with time integration when Newton-Raphson iteration failed",
                                tuple<std::string>(
                                  "stop",
                                  "continue",
                                  "repeat_step",
                                  "halve_step",
                                  "adapt_step",
                                  "rand_adapt_step",
                                  "rand_adapt_step_ele_err",
                                  "repeat_simulation"),
                                tuple<int>(
                                  divcont_stop,
                                  divcont_continue,
                                  divcont_repeat_step,
                                  divcont_halve_step,
                                  divcont_adapt_step,
                                  divcont_rand_adapt_step,
                                  divcont_rand_adapt_step_ele_err,
                                  divcont_repeat_simulation),
                                &cellstructdyn);

  setStringToIntegralParameter<int>("NLNSOL","fullnewton","Nonlinear solution technique",
                               tuple<std::string>(
                                 "vague",
                                 "fullnewton",
                                 "modnewton",
                                 "lsnewton",
                                 "ptc",
                                 "newtonlinuzawa",
                                 "augmentedlagrange",
                                 "NoxNewtonLineSearch",
                                 "noxgeneral",
                                 "noxnln",
                                 "NLNSOL"),
                               tuple<int>(
                                 soltech_vague,
                                 soltech_newtonfull,
                                 soltech_newtonmod,
                                 soltech_newtonls,
                                 soltech_ptc,
                                 soltech_newtonuzawalin,
                                 soltech_newtonuzawanonlin,
                                 soltech_noxnewtonlinesearch,
                                 soltech_noxgeneral,
                                 soltech_nox_nln,
                                 soltech_nlnsol),
                               &cellstructdyn);

  IntParameter("LSMAXITER",30,
               "maximum number of line search steps",
               &cellstructdyn);
  DoubleParameter("ALPHA_LS",0.5,
                  "step reduction factor alpha in (Newton) line search scheme",
                  &cellstructdyn);
  DoubleParameter("SIGMA_LS",1.e-4,
                  "sufficient descent factor in (Newton) line search scheme",
                  &cellstructdyn);

  setStringToIntegralParameter<int>("MATERIALTANGENT","analytical","way of evaluating the constitutive matrix",
                               tuple<std::string>(
                                 "analytical",
                                 "finitedifferences"),
                               tuple<int>(
                                 0,1),
                               &cellstructdyn);

  // Currently not used, but structure will be kept if someone wants to reimplement
  // AN 2013_05
  setStringToIntegralParameter<int>("CONTROLTYPE","load","load, disp, arc1, arc2 control",
                               tuple<std::string>(
                                 "load",
                                 "Load",
                                 "disp",
                                 "Disp",
                                 "Displacement",
                                 "arc1",
                                 "Arc1",
                                 "arc2",
                                 "Arc2"),
                               tuple<int>(
                                 control_load,
                                 control_load,
                                 control_disp,
                                 control_disp,
                                 control_disp,
                                 control_arc1,
                                 control_arc1,
                                 control_arc2,
                                 control_arc2),
                               &cellstructdyn);
  // Currently not used, but structure will be kept if someone wants to reimplement
  // AN 2013_05
  setNumericStringParameter("CONTROLNODE","-1 -1 -1",
                            "for methods other than load control: [node(fortran numbering)] [dof(c-numbering)] [curve(fortran numbering)]",
                            &cellstructdyn);

  setStringToIntegralParameter<int>("LOADLIN","No",
                                    "Use linearization of external follower load in Newton",
                                    yesnotuple,yesnovalue,&cellstructdyn);

  setStringToIntegralParameter<int>("MASSLIN","No","Application of nonlinear inertia terms",
  tuple<std::string>("No","no",
                     "Standard", "standard",
                     "Rotations", "rotations"),

  tuple<int>(ml_none,ml_none,
             ml_standard,ml_standard,
             ml_rotations,ml_rotations),
             &cellstructdyn);


// Since predicor "none" would be misleading, the usage of no predictor is called vague.
  setStringToIntegralParameter<int>("PREDICT","ConstDis","Type of predictor",
                               tuple<std::string>(
                                 "Vague",
                                 "ConstDis",
                                 "ConstVel",
                                 "ConstAcc",
                                 "ConstDisVelAcc",
                                 "TangDis",
                                 "ConstDisPres",
                                 "ConstDisVelAccPres"),
                               tuple<int>(
                                 pred_vague,
                                 pred_constdis,
                                 pred_constvel,
                                 pred_constacc,
                                 pred_constdisvelacc,
                                 pred_tangdis,
                                 pred_constdispres,
                                 pred_constdisvelaccpres),
                               &cellstructdyn);

  // Uzawa iteration for constraint systems
  DoubleParameter("UZAWAPARAM",1.0,"Parameter for Uzawa algorithm dealing with lagrange multipliers",&cellstructdyn);
  DoubleParameter("UZAWATOL",1.0E-8,"Tolerance for iterative solve with Uzawa algorithm",&cellstructdyn);
  IntParameter("UZAWAMAXITER",50,"maximum number of iterations allowed for uzawa algorithm before failure going to next newton step",&cellstructdyn);
  setStringToIntegralParameter<int>("UZAWAALGO","direct","",
                                 tuple<std::string>(
                                   "uzawa",
                                   "simple",
                                   "direct"),
                                 tuple<int>(
                                   consolve_uzawa,
                                   consolve_simple,
                                   consolve_direct),
                                 &cellstructdyn);

  // convergence criteria adaptivity
  setStringToIntegralParameter<int>("ADAPTCONV","No",
                               "Switch on adaptive control of linear solver tolerance for nonlinear solution",
                               yesnotuple,yesnovalue,&cellstructdyn);
  DoubleParameter("ADAPTCONV_BETTER",0.1,"The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",&cellstructdyn);

  setStringToIntegralParameter<int>("LUMPMASS","No",
                               "Lump the mass matrix for explicit time integration",
                               yesnotuple,yesnovalue,&cellstructdyn);

  setStringToIntegralParameter<int>("MODIFIEDEXPLEULER","Yes",
                               "Use the modified explicit Euler time integration scheme",
                               yesnotuple,yesnovalue,&cellstructdyn);

  // linear solver id used for structural problems
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for structural problems",&cellstructdyn);

  // flag decides if young's modulus is temperature dependent, so far only available
  // for temperature-dependent St.Venant Kirchhoff material
  setStringToIntegralParameter<int>("YOUNG_IS_TEMP_DEPENDENT","No",
                               "Use temperature-dependent Young's modulus",
                               yesnotuple,yesnovalue,&cellstructdyn);

  // where the geometry comes from
  setStringToIntegralParameter<int>(
    "GEOMETRY","full",
    "How the geometry is specified",
    tuple<std::string>(
      "full",
      "box",
      "file"),
    tuple<int>(
      INPAR::geometry_full,
      INPAR::geometry_box,
      INPAR::geometry_file),
    &cellstructdyn);

  /*----------------------------------------------------------------------*/
  /* parameters for one-step-theta cell structural integrator */
  Teuchos::ParameterList& onesteptheta = cellstructdyn.sublist("ONESTEPTHETA",false,"");

  DoubleParameter("THETA",1.0,"One-step-theta factor in (0,1]",&onesteptheta);


  /*----------------------------------------------------------------------*/
  /* timeadaptivity requested by base algorithm                           */

  Teuchos::ParameterList& cellstructtap = cellstructdyn.sublist("TIMEADAPTIVITY",false,"");

  setStringToIntegralParameter<int>(
    "KIND","NONE","Method for time step size adapivity",
    tuple<std::string>(
      "NONE",
      "ZienkiewiczXie",
      "AdamsBashforth2",
      "ExplicitEuler",
      "CentralDifference"),
    tuple<int>(
      INPAR::STR::timada_kind_none,
      INPAR::STR::timada_kind_zienxie,
      INPAR::STR::timada_kind_ab2,
      INPAR::STR::timada_kind_expleuler,
      INPAR::STR::timada_kind_centraldiff),
    &cellstructtap);

  DoubleParameter("OUTSYSPERIOD", 0.0, "Write system vectors (displacements, velocities, etc) every given period of time", &cellstructtap);
  DoubleParameter("OUTSTRPERIOD", 0.0, "Write stress/strain every given period of time", &cellstructtap);
  DoubleParameter("OUTENEPERIOD", 0.0, "Write energy every given period of time", &cellstructtap);
  DoubleParameter("OUTRESTPERIOD", 0.0, "Write restart data every given period of time", &cellstructtap);
  IntParameter("OUTSIZEEVERY", 0, "Write step size every given time step", &cellstructtap);

  DoubleParameter("STEPSIZEMAX", 0.0, "Limit maximally permitted time step size (>0)", &cellstructtap);
  DoubleParameter("STEPSIZEMIN", 0.0, "Limit minimally allowed time step size (>0)", &cellstructtap);
  DoubleParameter("SIZERATIOMAX", 0.0, "Limit maximally permitted change of time step size compared to previous size, important for multi-step schemes (>0)", &cellstructtap);
  DoubleParameter("SIZERATIOMIN", 0.0, "Limit minimally permitted change of time step size compared to previous size, important for multi-step schemes (>0)", &cellstructtap);
  DoubleParameter("SIZERATIOSCALE", 0.9, "This is a safety factor to scale theoretical optimal step size, should be lower than 1 and must be larger than 0", &cellstructtap);

  setStringToIntegralParameter<int>(
    "LOCERRNORM", "Vague", "Vector norm to treat error vector with",
    tuple<std::string>(
      "Vague",
      "L1",
      "L2",
      "Rms",
      "Inf"),
    tuple<int>(
      INPAR::STR::norm_vague,
      INPAR::STR::norm_l1,
      INPAR::STR::norm_l2,
      INPAR::STR::norm_rms,
      INPAR::STR::norm_inf),
    &cellstructtap);

  DoubleParameter("LOCERRTOL", 0.0, "Target local error tolerance (>0)", &cellstructtap);
  IntParameter("ADAPTSTEPMAX", 0, "Limit maximally allowed step size reduction attempts (>0)", &cellstructtap);

} // SetValidParameters


void INPAR::CELL::SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;

    /*--------------------------------------------------------------------*/
      // FOCAL ADHESION

    std::vector<Teuchos::RCP<ConditionComponent> > facomponents;

    facomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

    Teuchos::RCP<ConditionDefinition> linefa =
      Teuchos::rcp(new ConditionDefinition("DESIGN FOCAL ADHESION LINE CONDITIONS",
                                           "CellFocalAdhesion",
                                           "Cell Focal Adhesion",
                                           DRT::Condition::CellFocalAdhesion,
                                           true,
                                           DRT::Condition::Line));
    Teuchos::RCP<ConditionDefinition> surffa =
      Teuchos::rcp(new ConditionDefinition("DESIGN FOCAL ADHESION SURF CONDITIONS",
                                           "CellFocalAdhesion",
                                           "Cell Focal Adhesion",
                                           DRT::Condition::CellFocalAdhesion,
                                           true,
                                           DRT::Condition::Surface));

    for (unsigned i=0; i<facomponents.size(); ++i)
    {
      linefa->AddComponent(facomponents[i]);
      surffa->AddComponent(facomponents[i]);
    }

    condlist.push_back(linefa);
    condlist.push_back(surffa);


}
