/*----------------------------------------------------------------------*/
/*!
\file inpar_scatra.cpp

\brief Input parameters for scatra

\level 1

\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
*/
/*----------------------------------------------------------------------*/
#include "inpar_scatra.H"

#include "drt_validparameters.H"
#include "inpar_fluid.H"
#include "inpar_s2i.H"
#include "inpar_thermo.H"
#include "inpar_bio.H"

#include "../drt_lib/drt_conditiondefinition.H"

void INPAR::SCATRA::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& scatradyn = list->sublist(
      "SCALAR TRANSPORT DYNAMIC",
      false,
      "control parameters for scalar transport problems\n");

  setStringToIntegralParameter<int>("SOLVERTYPE","linear_full",
                               "type of scalar transport solver",
                               tuple<std::string>(
                                 "linear_full",
                                 "linear_incremental",
                                 "nonlinear",
                                 "nonlinear_multiscale_macrotomicro",
                                 "nonlinear_multiscale_macrotomicro_aitken",
                                 "nonlinear_multiscale_macrotomicro_aitken_dofsplit",
                                 "nonlinear_multiscale_microtomacro"
                                 ),
                               tuple<int>(
                                   solvertype_linear_full,
                                   solvertype_linear_incremental,
                                   solvertype_nonlinear,
                                   solvertype_nonlinear_multiscale_macrotomicro,
                                   solvertype_nonlinear_multiscale_macrotomicro_aitken,
                                   solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit,
                                   solvertype_nonlinear_multiscale_microtomacro
                                   ),
                               &scatradyn);

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
                               &scatradyn);

  DoubleParameter("MAXTIME",1000.0,"Total simulation time",&scatradyn);
  IntParameter("NUMSTEP",20,"Total number of time steps",&scatradyn);
  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&scatradyn);
  DoubleParameter("THETA",0.5,"One-step-theta time integration factor",&scatradyn);
  DoubleParameter("ALPHA_M",0.5,"Generalized-alpha time integration factor",&scatradyn);
  DoubleParameter("ALPHA_F",0.5,"Generalized-alpha time integration factor",&scatradyn);
  DoubleParameter("GAMMA",0.5,"Generalized-alpha time integration factor",&scatradyn);
  IntParameter("RESULTSEVRY",1,"Increment for writing solution",&scatradyn);
  IntParameter("RESTARTEVRY",1,"Increment for writing restart",&scatradyn);
  IntParameter("MATID",-1,"Material ID for automatic mesh generation",&scatradyn);

  setStringToIntegralParameter<int>("VELOCITYFIELD","zero",
                               "type of velocity field used for scalar transport problems",
                               tuple<std::string>(
                                 "zero",
                                 "function",
                                 "Navier_Stokes"
                                 ),
                               tuple<int>(
                                   velocity_zero,
                                   velocity_function,
                                   velocity_Navier_Stokes),
                               &scatradyn);

  IntParameter("VELFUNCNO",-1,"function number for scalar transport velocity field",&scatradyn);

  {
    // a standard Teuchos::tuple can have at maximum 10 entries! We have to circumvent this here.
    Teuchos::Tuple<std::string,13> name;
    Teuchos::Tuple<int,13> label;
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
    name[12] = "algebraic_field_dependence";   label[12] = initialfield_algebraic_field_dependence;

    setStringToIntegralParameter<int>(
        "INITIALFIELD",
        "zero_field",
        "Initial Field for scalar transport problem",
        name,
        label,
        &scatradyn);
  }

  IntParameter("INITFUNCNO",-1,"function number for scalar transport initial field",&scatradyn);

  BoolParameter("SPHERICALCOORDS","No","use of spherical coordinates",&scatradyn);

  setStringToIntegralParameter<int>("CALCERROR","No",
                               "compute error compared to analytical solution",
                               tuple<std::string>(
                                 "No",
                                 "Kwok_Wu",
                                 "ConcentricCylinders",
                                 "Electroneutrality",
                                 "error_by_function",
                                 "error_by_condition",
                                 "SphereDiffusion",
                                 "AnalyticSeries"
                                 ),
                               tuple<int>(
                                   calcerror_no,
                                   calcerror_Kwok_Wu,
                                   calcerror_cylinder,
                                   calcerror_electroneutrality,
                                   calcerror_byfunction,
                                   calcerror_bycondition,
                                   calcerror_spherediffusion,
                                   calcerror_AnalyticSeries
                                   ),
                               &scatradyn);

  IntParameter("CALCERRORNO",-1,"function number for scalar transport error computation",&scatradyn);

  setStringToIntegralParameter<int>(
      "CALCFLUX_DOMAIN",
      "No",
      "output of diffusive/total flux vectors inside domain",
      tuple<std::string>(
          "No",
          "total",
          "diffusive"
          ),
      tuple<int>(
          flux_none,
          flux_total,
          flux_diffusive
          ),
      &scatradyn
      );

  BoolParameter("CALCFLUX_DOMAIN_LUMPED","Yes","perform approximate domain flux calculation involving matrix lumping",&scatradyn);

  setStringToIntegralParameter<int>(
      "CALCFLUX_BOUNDARY",
      "No",
      "output of convective/diffusive/total flux vectors on boundary",
      tuple<std::string>(
          "No",
          "total",
          "diffusive",
          "convective"
          ),
      tuple<int>(
          flux_none,
          flux_total,
          flux_diffusive,
          flux_convective
          ),
      &scatradyn
      );

  BoolParameter("CALCFLUX_BOUNDARY_LUMPED","Yes","perform approximate boundary flux calculation involving matrix lumping",&scatradyn);

  setNumericStringParameter("WRITEFLUX_IDS","-1",
      "Write diffusive/total flux vector fields for these scalar fields only (starting with 1)",
      &scatradyn);

  setStringToIntegralParameter<int>("OUTPUTSCALARS","none","Output of total and mean values for transported scalars",
                               tuple<std::string>(
                                   "none",
                                   "entire_domain",
                                   "by_condition",
                                   "entire_domain_and_by_condition"
                                 ),
                               tuple<int>(
                                   outputscalars_none,
                                   outputscalars_entiredomain,
                                   outputscalars_condition,
                                   outputscalars_entiredomain_condition),
                               &scatradyn);
  BoolParameter("OUTINTEGRREAC","No","Output of integral reaction values",&scatradyn);
  BoolParameter("OUTPUT_GMSH","No","Do you want to write Gmsh postprocessing files?",&scatradyn);

  BoolParameter("MATLAB_STATE_OUTPUT","No","Do you want to write the state solution to Matlab file?",&scatradyn);

  setStringToIntegralParameter<int>("CONVFORM","convective","form of convective term",
                               tuple<std::string>(
                                 "convective",
                                 "conservative"
                                 ),
                               tuple<int>(
                                 convform_convective,
                                 convform_conservative),
                               &scatradyn);

  BoolParameter("NEUMANNINFLOW",
      "no","Flag to (de)activate potential Neumann inflow term(s)",&scatradyn);

  BoolParameter("CONV_HEAT_TRANS",
      "no","Flag to (de)activate potential convective heat transfer boundary conditions",&scatradyn);

  BoolParameter("SKIPINITDER",
      "no","Flag to skip computation of initial time derivative",&scatradyn);

  setStringToIntegralParameter<int>("FSSUGRDIFF",
                               "No",
                               "fine-scale subgrid diffusivity",
                               tuple<std::string>(
                                 "No",
                                 "artificial",
                                 "Smagorinsky_all",
                                 "Smagorinsky_small"
                                 ),
                               tuple<int>(
                                   fssugrdiff_no,
                                   fssugrdiff_artificial,
                                   fssugrdiff_smagorinsky_all,
                                   fssugrdiff_smagorinsky_small),
                               &scatradyn);

  setStringToIntegralParameter<int>("MESHTYING", "no", "Flag to (de)activate mesh tying algorithm",
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
                                    &scatradyn);

  // Type of coupling strategy between the two fields
  setStringToIntegralParameter<int>(
                              "FIELDCOUPLING","matching",
                              "Type of coupling strategy between fields",
                              tuple<std::string>(
                                "matching",
                                "volmortar"
                                ),
                              tuple<int>(
                                  coupling_match,
                                  coupling_volmortar
                                ),
                              &scatradyn);

  // linear solver id used for scalar transport/elch problems
  IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for scalar transport/elch...",&scatradyn);
  // linear solver id used for l2 projection problems (e.g. gradient projections)
  IntParameter("L2_PROJ_LINEAR_SOLVER",-1,"number of linear solver used for l2-projection sub-problems",&scatradyn);
  //IntParameter("SIMPLER_SOLVER",-1,"number of linear solver used for ELCH (solved with SIMPLER)...",&scatradyn);

  // flag for natural convection effects
  BoolParameter("NATURAL_CONVECTION","No","Include natural convection effects",&scatradyn);

  // parameters for finite difference check
  setStringToIntegralParameter<int>("FDCHECK", "none", "flag for finite difference check: none, local, or global",
                                    tuple<std::string>(
                                        "none",
                                        "global",            // perform finite difference check on time integrator level
                                        "global_extended",   // perform finite difference check on time integrator level for extended system matrix (e.g., involving Lagrange multipliers or interface layer thicknesses)
                                        "local"              // perform finite difference check on element level
                                        ),
                                    tuple<int>(
                                        fdcheck_none,
                                        fdcheck_global,
                                        fdcheck_global_extended,
                                        fdcheck_local
                                        ),
                                    &scatradyn);
  DoubleParameter("FDCHECKEPS",1.e-6,"dof perturbation magnitude for finite difference check (1.e-6 seems to work very well, whereas smaller values don't)",&scatradyn);
  DoubleParameter("FDCHECKTOL",1.e-6,"relative tolerance for finite difference check",&scatradyn);

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
      &scatradyn
      );

  // parameter for using p-adpativity and semi-implicit evaluation of the reaction term (at the moment only used for HDG and cardiac monodomain problems)
  BoolParameter("PADAPTIVITY", "no","Flag to (de)activate p-adativity",&scatradyn);
  DoubleParameter("PADAPTERRORTOL",1e-6,"The error tolerance to calculate the variation of the elemental degree",&scatradyn);
  DoubleParameter("PADAPTERRORBASE",1.66,"The error tolerance base to calculate the variation of the elemental degree",&scatradyn);
  IntParameter("PADAPTDEGREEMAX",4,"The max. degree of the shape functions",&scatradyn);
  BoolParameter("SEMIIMPLICIT", "no","Flag to (de)activate semi-implicit calculation of the reaction term",&scatradyn);

  // flag for output of performance statistics associated with linear solver into *.csv file
  BoolParameter("OUTPUTLINSOLVERSTATS","No","flag for output of performance statistics associated with linear solver into *.csv file",&scatradyn);

  // flag for output of performance statistics associated with nonlinear solver into *.csv file
  BoolParameter("OUTPUTNONLINSOLVERSTATS","No","flag for output of performance statistics associated with nonlinear solver into *.csv file",&scatradyn);

  // flag for point-based null space calculation
  BoolParameter("NULLSPACE_POINTBASED","No","flag for point-based null space calculation",&scatradyn);

  // flag for adaptive time stepping
  BoolParameter("ADAPTIVE_TIMESTEPPING","No","flag for adaptive time stepping",&scatradyn);
  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scatra_variational = scatradyn.sublist(
      "VARIATIONAL",
      false,
      "control parameters for solving problems under a VARIATIONAL setting\n");

  BoolParameter("SEMIMPLICITFUNCTIONAL","no","Flag to evaluate concentration implicitly or explicitly on the functional",&scatra_variational);
  BoolParameter("BLOCKPRECOND","NO","Switch to block-preconditioned family of solvers, only works with block preconditioners like CheapSIMPLE!",&scatra_variational);
  BoolParameter("ANALYTIC2PARAVIEW","NO","Flag to send the analytic solution to Paraview, if provided!",&scatra_variational);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scatra_nonlin = scatradyn.sublist(
      "NONLINEAR",
      false,
      "control parameters for solving nonlinear SCATRA problems\n");

  IntParameter("ITEMAX",10,"max. number of nonlin. iterations",&scatra_nonlin);
  DoubleParameter("CONVTOL",1e-6,"Tolerance for convergence check",&scatra_nonlin);
  IntParameter("ITEMAX_OUTER",10,"Maximum number of outer iterations in partitioned coupling schemes (natural convection, multi-scale simulations etc.)",&scatra_nonlin);
  DoubleParameter("CONVTOL_OUTER",1e-6,"Convergence check tolerance for outer loop in partitioned coupling schemes (natural convection, multi-scale simulations etc.)",&scatra_nonlin);
  BoolParameter("EXPLPREDICT","no","do an explicit predictor step before starting nonlinear iteration",&scatra_nonlin);
  DoubleParameter("ABSTOLRES",1e-14,"Absolute tolerance for deciding if residual of nonlinear problem is already zero",&scatra_nonlin);

  // convergence criteria adaptivity
  BoolParameter("ADAPTCONV","yes","Switch on adaptive control of linear solver tolerance for nonlinear solution",&scatra_nonlin);
  DoubleParameter("ADAPTCONV_BETTER",0.1,"The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",&scatra_nonlin);

/*----------------------------------------------------------------------*/
  Teuchos::ParameterList& scatradyn_stab = scatradyn.sublist("STABILIZATION",false,"control parameters for the stabilization of scalar transport problems");

  // this parameter governs type of stabilization
  setStringToIntegralParameter<int>("STABTYPE",
                                    "SUPG",
                                    "type of stabilization (if any)",
                               tuple<std::string>(
                                 "no_stabilization",
                                 "SUPG",
                                 "GLS",
                                 "USFEM",
                                 "centered",
                                 "upwind"),
                               tuple<std::string>(
                                 "Do not use any stabilization -> only reasonable for low-Peclet-number flows",
                                 "Use SUPG",
                                 "Use GLS",
                                 "Use USFEM",
                                 "Use centered scheme",
                                 "Use upwinded scheme")  ,
                               tuple<int>(
                                   stabtype_no_stabilization,
                                   stabtype_SUPG,
                                   stabtype_GLS,
                                   stabtype_USFEM,
                                   stabtype_hdg_centered,
                                   stabtype_hdg_upwind),
                               &scatradyn_stab);

  // this parameter governs whether subgrid-scale velocity is included
  BoolParameter("SUGRVEL","no","potential incorporation of subgrid-scale velocity",&scatradyn_stab);

  // this parameter governs whether all-scale subgrid diffusivity is included
  BoolParameter("ASSUGRDIFF","no",
      "potential incorporation of all-scale subgrid diffusivity (a.k.a. discontinuity-capturing) term",&scatradyn_stab);

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
                                 "Zero",
                                 "Numerical_Value"),
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
                                    tau_zero,
                                    tau_numerical_value),
                               &scatradyn_stab);

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
                               &scatradyn_stab);

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
                               &scatradyn_stab);

  // this parameter selects the location where tau is evaluated
  setStringToIntegralParameter<int>("EVALUATION_TAU",
                               "element_center",
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
                               &scatradyn_stab);

  // this parameter selects the location where the material law is evaluated
  // (does not fit here very well, but parameter transfer is easier)
  setStringToIntegralParameter<int>("EVALUATION_MAT",
                               "element_center",
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
                               &scatradyn_stab);

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
                               &scatradyn_stab);

   // this parameter defines the numerical value, if stabilization with numerical values is used
   DoubleParameter("TAU_VALUE",0.0,"Numerical value for tau for stabilization",&scatradyn_stab);

   // ----------------------------------------------------------------------
   // artery mesh tying
   Teuchos::ParameterList& scatradyn_art = scatradyn.sublist("ARTERY COUPLING",false,
     "Parameters for artery mesh tying"
     );

   setStringToIntegralParameter<int>("ARTERY_COUPLING_METHOD",
                                "None",
                                "Coupling method for artery coupling.",
                                tuple<std::string>(
                                  "None",
                                  "Nodal",
                                   "GPTS",
                                  "MP"),
                                tuple<std::string>(
                                  "none",
                                  "Nodal Coupling",
                                  "Gauss-Point-To-Segment Approach",
                                  "Mortar Penalty Approach"),
                                tuple<int>(
                                  INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::none,       // none
                                  INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::nodal,      // Nodal Coupling
                                  INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::gpts,       // Gauss-Point-To-Segment Approach
                                  INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp          // Mortar Penalty Approach
                                ),
                                &scatradyn_art);

   // penalty parameter
   DoubleParameter("PENALTY", 1000.0, "Penalty parameter for line-based coupling", &scatradyn_art);

   // coupled artery dofs for mesh tying
   setNumericStringParameter("COUPLEDDOFS_ARTSCATRA","-1.0",
                             "coupled artery dofs for mesh tying",
                             &scatradyn_art);

   // coupled porofluid dofs for mesh tying
   setNumericStringParameter("COUPLEDDOFS_SCATRA","-1.0",
                             "coupled porofluid dofs for mesh tying",
                             &scatradyn_art);

   // functions for coupling (arteryscatra part)
   setNumericStringParameter("REACFUNCT_ART","-1",
                             "functions for coupling (arteryscatra part)",
                             &scatradyn_art);

   // scale for coupling (arteryscatra part)
   setNumericStringParameter("SCALEREAC_ART","0",
                             "scale for coupling (arteryscatra part)",
                             &scatradyn_art);

   // functions for coupling (scatra part)
   setNumericStringParameter("REACFUNCT_CONT","-1",
                             "functions for coupling (scatra part)",
                             &scatradyn_art);

   // scale for coupling (scatra part)
   setNumericStringParameter("SCALEREAC_CONT","0",
                             "scale for coupling (scatra part)",
                             &scatradyn_art);

   // Flag if artery elements are evaluated in reference or current configuration
   BoolParameter("EVALUATE_IN_REF_CONFIG","yes","Flag if artery elements are evaluated in reference or current configuration",&scatradyn_art);
}



void INPAR::SCATRA::SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // Boundary flux evaluation condition for scalar transport
  Teuchos::RCP<ConditionDefinition> linebndryfluxeval =
    Teuchos::rcp(new ConditionDefinition("SCATRA FLUX CALC LINE CONDITIONS",
                                         "ScaTraFluxCalc",
                                         "Scalar Transport Boundary Flux Calculation",
                                         DRT::Condition::ScaTraFluxCalc,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surfbndryfluxeval =
    Teuchos::rcp(new ConditionDefinition("SCATRA FLUX CALC SURF CONDITIONS",
                                         "ScaTraFluxCalc",
                                         "Scalar Transport Boundary Flux Calculation",
                                         DRT::Condition::ScaTraFluxCalc,
                                         true,
                                         DRT::Condition::Surface));
  condlist.push_back(linebndryfluxeval);
  condlist.push_back(surfbndryfluxeval);

  /*--------------------------------------------------------------------*/
  // conditions for calculation of total and mean values of transported scalars
  Teuchos::RCP<ConditionDefinition> totalandmeanscalarline =
      Teuchos::rcp(new ConditionDefinition("DESIGN TOTAL AND MEAN SCALAR LINE CONDITIONS",
                                           "TotalAndMeanScalar",
                                           "calculation of total and mean values of transported scalars",
                                           DRT::Condition::TotalAndMeanScalar,
                                           true,
                                           DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> totalandmeanscalarsurf =
      Teuchos::rcp(new ConditionDefinition("DESIGN TOTAL AND MEAN SCALAR SURF CONDITIONS",
                                           "TotalAndMeanScalar",
                                           "calculation of total and mean values of transported scalars",
                                           DRT::Condition::TotalAndMeanScalar,
                                           true,
                                           DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> totalandmeanscalarvol =
      Teuchos::rcp(new ConditionDefinition("DESIGN TOTAL AND MEAN SCALAR VOL CONDITIONS",
                                           "TotalAndMeanScalar",
                                           "calculation of total and mean values of transported scalars",
                                           DRT::Condition::TotalAndMeanScalar,
                                           true,
                                           DRT::Condition::Volume));

  // equip condition definitions with input file line components
  std::vector<Teuchos::RCP<ConditionComponent> > totalandmeanscalarcomponents;

  {
    totalandmeanscalarcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ID")));
    totalandmeanscalarcomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));
  }

  // insert input file line components into condition definitions
  for (unsigned i=0; i<totalandmeanscalarcomponents.size(); ++i)
  {
    totalandmeanscalarline->AddComponent(totalandmeanscalarcomponents[i]);
    totalandmeanscalarsurf->AddComponent(totalandmeanscalarcomponents[i]);
    totalandmeanscalarvol->AddComponent(totalandmeanscalarcomponents[i]);
  }

  // insert condition definitions into global list of valid condition definitions
  condlist.push_back(totalandmeanscalarline);
  condlist.push_back(totalandmeanscalarsurf);
  condlist.push_back(totalandmeanscalarvol);

  /*--------------------------------------------------------------------*/
  // conditions for calculation of relative error with reference to analytical solution
  Teuchos::RCP<ConditionDefinition> relerrorline =
      Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA RELATIVE ERROR LINE CONDITIONS",
                                           "ScatraRelError",
                                           "calculation of relative error with reference to analytical solution",
                                           DRT::Condition::ScatraRelError,
                                           true,
                                           DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> relerrorsurf =
      Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA RELATIVE ERROR SURF CONDITIONS",
                                           "ScatraRelError",
                                           "calculation of relative error with reference to analytical solution",
                                           DRT::Condition::ScatraRelError,
                                           true,
                                           DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> relerrorvol =
      Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA RELATIVE ERROR VOL CONDITIONS",
                                           "ScatraRelError",
                                           "calculation of relative error with reference to analytical solution",
                                           DRT::Condition::ScatraRelError,
                                           true,
                                           DRT::Condition::Volume));

  // equip condition definitions with input file line components
  std::vector<Teuchos::RCP<ConditionComponent> > relerrorcomponents;

  {
    relerrorcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("ID")));
    relerrorcomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));
    relerrorcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("Function")));
    relerrorcomponents.push_back(Teuchos::rcp(new IntConditionComponent("FunctionID")));
  }

  // insert input file line components into condition definitions
  for (unsigned i=0; i<relerrorcomponents.size(); ++i)
  {
    relerrorline->AddComponent(relerrorcomponents[i]);
    relerrorsurf->AddComponent(relerrorcomponents[i]);
    relerrorvol->AddComponent(relerrorcomponents[i]);
  }

  // insert condition definitions into global list of valid condition definitions
  condlist.push_back(relerrorline);
  condlist.push_back(relerrorsurf);
  condlist.push_back(relerrorvol);

  /*--------------------------------------------------------------------*/
  // Coupling of different scalar transport fields

  std::vector<Teuchos::RCP<ConditionComponent> > scatracoupcomponents;

  std::vector<Teuchos::RCP<SeparatorConditionComponent> > KKintsepveccompstoich;
  KKintsepveccompstoich.push_back(Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
  // definition int vectors
  std::vector<Teuchos::RCP<IntVectorConditionComponent> > KKintveccompstoich;
  KKintveccompstoich.push_back(Teuchos::rcp(new IntVectorConditionComponent("onoff",2)));
  // definition separator for real vectors: length of the real vector is zero -> nothing is read
  std::vector<Teuchos::RCP<SeparatorConditionComponent> > KKrealsepveccompstoich;
  // definition real vectors: length of the real vector is zero -> nothing is read
  std::vector<Teuchos::RCP<RealVectorConditionComponent> > KKrealveccompstoich;


  scatracoupcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("NUMSCAL")) );
  scatracoupcomponents.push_back(Teuchos::rcp(new IntRealBundle(
                                 "intreal bundle numscal",
                                Teuchos::rcp(new IntConditionComponent("numscal")),
                                KKintsepveccompstoich,
                                KKintveccompstoich,
                                KKrealsepveccompstoich,
                                KKrealveccompstoich)) );
  scatracoupcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("COUPID")));
  scatracoupcomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));
  scatracoupcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("PERMCOEF")));
  scatracoupcomponents.push_back(Teuchos::rcp(new RealConditionComponent("permeability coefficient")));
  scatracoupcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("CONDUCT")));
  scatracoupcomponents.push_back(Teuchos::rcp(new RealConditionComponent("hydraulic conductivity")));
  scatracoupcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("FILTR")));
  scatracoupcomponents.push_back(Teuchos::rcp(new RealConditionComponent("filtration coefficient")));
  scatracoupcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("WSSONOFF")));
  scatracoupcomponents.push_back(Teuchos::rcp(new IntConditionComponent("wss onoff")));
  scatracoupcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("WSSCOEFFS")));
  scatracoupcomponents.push_back(Teuchos::rcp(new RealVectorConditionComponent("wss coeffs",2)));


  Teuchos::RCP<ConditionDefinition> surfscatracoup =
    Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA COUPLING SURF CONDITIONS",
                                         "ScaTraCoupling",
                                         "ScaTra Coupling",
                                         DRT::Condition::ScaTraCoupling,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<scatracoupcomponents.size(); ++i)
  {
    surfscatracoup->AddComponent(scatracoupcomponents[i]);
  }

  condlist.push_back(surfscatracoup);

  /*--------------------------------------------------------------------*/
  // Robin boundary condition for scalar transport problems
  // line
  Teuchos::RCP<ConditionDefinition> scatrarobinline =
    Teuchos::rcp(new ConditionDefinition("DESIGN TRANSPORT ROBIN LINE CONDITIONS",
                                         "TransportRobin",
                                         "Scalar Transport Robin Boundary Condition",
                                         DRT::Condition::TransportRobin,
                                         true,
                                         DRT::Condition::Line));
  // surface
  Teuchos::RCP<ConditionDefinition> scatrarobinsurf =
    Teuchos::rcp(new ConditionDefinition("DESIGN TRANSPORT ROBIN SURF CONDITIONS",
                                         "TransportRobin",
                                         "Scalar Transport Robin Boundary Condition",
                                         DRT::Condition::TransportRobin,
                                         true,
                                         DRT::Condition::Surface));

  std::vector<Teuchos::RCP<ConditionComponent> > scatrarobincomponents;

  std::vector<Teuchos::RCP<SeparatorConditionComponent> > Robinintsepveccompstoich;
  Robinintsepveccompstoich.push_back(Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
  // definition int vectors
  std::vector<Teuchos::RCP<IntVectorConditionComponent> > Robinintveccompstoich;
  Robinintveccompstoich.push_back(Teuchos::rcp(new IntVectorConditionComponent("onoff",2)));
  // definition separator for real vectors: length of the real vector is zero -> nothing is read
  std::vector<Teuchos::RCP<SeparatorConditionComponent> > Robinrealsepveccompstoich;
  // definition real vectors: length of the real vector is zero -> nothing is read
  std::vector<Teuchos::RCP<RealVectorConditionComponent> > Robinrealveccompstoich;


  scatrarobincomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("NUMSCAL")) );
  scatrarobincomponents.push_back(Teuchos::rcp(new IntRealBundle(
                                 "intreal bundle numscal",
                                Teuchos::rcp(new IntConditionComponent("numscal")),
                                Robinintsepveccompstoich,
                                Robinintveccompstoich,
                                Robinrealsepveccompstoich,
                                Robinrealveccompstoich)) );

  scatrarobincomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("PREFACTOR")));
  scatrarobincomponents.push_back(Teuchos::rcp(new RealConditionComponent("prefactor")));
  scatrarobincomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("REFVALUE")));
  scatrarobincomponents.push_back(Teuchos::rcp(new RealConditionComponent("refvalue")));

  for(unsigned i=0; i<scatrarobincomponents.size(); ++i)
  {
    scatrarobinline->AddComponent(scatrarobincomponents[i]);
    scatrarobinsurf->AddComponent(scatrarobincomponents[i]);
  }

  condlist.push_back(scatrarobinline);
  condlist.push_back(scatrarobinsurf);

  /*--------------------------------------------------------------------*/
  // Neumann inflow for SCATRA

  Teuchos::RCP<ConditionDefinition> linetransportneumanninflow =
    Teuchos::rcp(new ConditionDefinition("TRANSPORT NEUMANN INFLOW LINE CONDITIONS",
                                         "TransportNeumannInflow",
                                         "Line Transport Neumann Inflow",
                                         DRT::Condition::TransportNeumannInflow,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surftransportneumanninflow =
    Teuchos::rcp(new ConditionDefinition("TRANSPORT NEUMANN INFLOW SURF CONDITIONS",
                                         "TransportNeumannInflow",
                                         "Surface Transport Neumann Inflow",
                                         DRT::Condition::TransportNeumannInflow,
                                         true,
                                         DRT::Condition::Surface));

  condlist.push_back(linetransportneumanninflow);
  condlist.push_back(surftransportneumanninflow);

  /*--------------------------------------------------------------------*/
  // Scatra convective heat transfer (Newton's law of heat transfer)

  std::vector<Teuchos::RCP<ConditionComponent> > transportthermoconvectcomponents;

  // decide here if approximation is sufficient
  // --> Tempn (old temperature T_n)
  // or if the exact solution is needed
  // --> Tempnp (current temperature solution T_n+1) with linearisation
  transportthermoconvectcomponents.push_back(
    Teuchos::rcp(
      new StringConditionComponent(
        "temperature state","Tempnp",
        Teuchos::tuple<std::string>("Tempnp","Tempn"),
        Teuchos::tuple<std::string>("Tempnp","Tempn"))));
  // heat transfer coefficient h
  transportthermoconvectcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("coeff")));
  transportthermoconvectcomponents.push_back(Teuchos::rcp(new RealConditionComponent("coeff")));
  // surrounding (fluid) temperature T_oo
  transportthermoconvectcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("surtemp")));
  transportthermoconvectcomponents.push_back(Teuchos::rcp(new RealConditionComponent("surtemp")));
  // time curve to increase the surrounding (fluid) temperature T_oo in time
  transportthermoconvectcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("surtempfunct")));
  transportthermoconvectcomponents.push_back(Teuchos::rcp(new IntConditionComponent("surtempfunct",true,true)));
  // time curve to increase the complete boundary condition, i.e., the heat flux
  transportthermoconvectcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("funct")));
  transportthermoconvectcomponents.push_back(Teuchos::rcp(new IntConditionComponent("funct",true,true)));

  Teuchos::RCP<ConditionDefinition> linetransportthermoconvect =
    Teuchos::rcp(new ConditionDefinition("TRANSPORT THERMO CONVECTION LINE CONDITIONS",
                                         "TransportThermoConvections",
                                         "Line Transport Thermo Convections",
                                         DRT::Condition::TransportThermoConvections,
                                         true,
                                         DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surftransportthermoconvect =
    Teuchos::rcp(new ConditionDefinition("TRANSPORT THERMO CONVECTION SURF CONDITIONS",
                                         "TransportThermoConvections",
                                         "Surface Transport Thermo Convections",
                                         DRT::Condition::TransportThermoConvections,
                                         true,
                                         DRT::Condition::Surface));

  for (unsigned i=0; i<transportthermoconvectcomponents.size(); ++i)
  {
    linetransportthermoconvect->AddComponent(transportthermoconvectcomponents[i]);
    surftransportthermoconvect->AddComponent(transportthermoconvectcomponents[i]);
  }

  condlist.push_back(linetransportthermoconvect);
  condlist.push_back(surftransportthermoconvect);

  /*--------------------------------------------------------------------*/
  // Macro-micro coupling condition for micro scale in multi-scale scalar transport problems
  {
    // condition definition
    Teuchos::RCP<ConditionDefinition> multiscalecouplingpoint =
        Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA MULTI-SCALE COUPLING POINT CONDITIONS",
                                             "ScatraMultiScaleCoupling",
                                             "Scalar transport multi-scale coupling condition",
                                             DRT::Condition::ScatraMultiScaleCoupling,
                                             false,
                                             DRT::Condition::Point));

    // equip condition definition with input file line components
    std::vector<Teuchos::RCP<ConditionComponent> > multiscalecouplingcomponents;

    // kinetic models for macro-micro coupling
    std::vector<Teuchos::RCP<CondCompBundle> > kineticmodels;
    {
      {
        // constant permeability
        std::vector<Teuchos::RCP<ConditionComponent> > constperm;
        constperm.push_back(Teuchos::rcp(new SeparatorConditionComponent("numscal")));                // total number of existing scalars
        std::vector<Teuchos::RCP<SeparatorConditionComponent> > intsepcomp;                           // empty vector --> no separators for integer vectors needed
        std::vector<Teuchos::RCP<IntVectorConditionComponent> > intvectcomp;                          // empty vector --> no integer vectors needed
        std::vector<Teuchos::RCP<SeparatorConditionComponent> > realsepcomp;
        realsepcomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("permeabilities")));       // string separator in front of real permeability vector in input file line
        std::vector<Teuchos::RCP<RealVectorConditionComponent> > realvectcomp;
        realvectcomp.push_back(Teuchos::rcp(new RealVectorConditionComponent("permeabilities",0)));   // real vector of constant permeabilities
        constperm.push_back(Teuchos::rcp(new IntRealBundle(
            "permeabilities",
            Teuchos::rcp(new IntConditionComponent("numscal")),
            intsepcomp,
            intvectcomp,
            realsepcomp,
            realvectcomp
        )));

        kineticmodels.push_back(Teuchos::rcp(new CondCompBundle("ConstantPermeability",constperm,INPAR::S2I::kinetics_constperm)));
      }

      {
        // Butler-Volmer
        std::vector<Teuchos::RCP<ConditionComponent> > butlervolmer;
        butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("numscal")));            // total number of existing scalars
        std::vector<Teuchos::RCP<SeparatorConditionComponent> > intsepcomp;
        intsepcomp.push_back(Teuchos::rcp(new SeparatorConditionComponent("stoichiometries")));
        std::vector<Teuchos::RCP<IntVectorConditionComponent> > intvectcomp;                         // string separator in front of integer stoichiometry vector in input file line
        intvectcomp.push_back(Teuchos::rcp(new IntVectorConditionComponent("stoichiometries",0)));   // integer vector of stoichiometric coefficients
        std::vector<Teuchos::RCP<SeparatorConditionComponent> > realsepcomp;                         // empty vector --> no separators for real vectors needed
        std::vector<Teuchos::RCP<RealVectorConditionComponent> > realvectcomp;                       // empty vector --> no real vectors needed
        butlervolmer.push_back(Teuchos::rcp(new IntRealBundle(
            "stoichiometries",
            Teuchos::rcp(new IntConditionComponent("numscal")),
            intsepcomp,
            intvectcomp,
            realsepcomp,
            realvectcomp
            )));
        butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("e-")));
        butlervolmer.push_back(Teuchos::rcp(new IntConditionComponent("e-")));
        butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("k_r")));
        butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("k_r")));
        butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_a")));
        butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("alpha_a")));
        butlervolmer.push_back(Teuchos::rcp(new SeparatorConditionComponent("alpha_c")));
        butlervolmer.push_back(Teuchos::rcp(new RealConditionComponent("alpha_c")));

        kineticmodels.push_back(Teuchos::rcp(new CondCompBundle("Butler-Volmer",butlervolmer,INPAR::S2I::kinetics_butlervolmer)));
      }
    } // kinetic models for macro-micro coupling

    // insert kinetic models into vector with condition components
    multiscalecouplingcomponents.push_back(Teuchos::rcp(new SeparatorConditionComponent("KineticModel")));
    multiscalecouplingcomponents.push_back(Teuchos::rcp(new CondCompBundleSelector(
        "kinetic models for macro-micro coupling",
        Teuchos::rcp(new StringConditionComponent(
            "kinetic model",
            "ConstantPermeability",
            Teuchos::tuple<std::string>("ConstantPermeability","Butler-Volmer"),
            Teuchos::tuple<int>(INPAR::S2I::kinetics_constperm,INPAR::S2I::kinetics_butlervolmer))),
        kineticmodels)));

    // insert input file line components into condition definitions
    for(unsigned i=0; i<multiscalecouplingcomponents.size(); ++i)
      multiscalecouplingpoint->AddComponent(multiscalecouplingcomponents[i]);

    // insert condition definitions into global list of valid condition definitions
    condlist.push_back(multiscalecouplingpoint);
  }

  /*--------------------------------------------------------------------*/
  // conditions for calculation of calculation of heterogeneous reactions
  Teuchos::RCP<ConditionDefinition> scatraheteroreactionmasterline =
      Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA HETEROGENEOUS REACTION LINE CONDITIONS / MASTER",
                                           "ScatraHeteroReactionMaster",
                                           "calculation of heterogeneous reactions",
                                           DRT::Condition::ScatraHeteroReactionCondMaster,
                                           true,
                                           DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> scatraheteroreactionmastersurf =
      Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA HETEROGENEOUS REACTION SURF CONDITIONS / MASTER",
                                           "ScatraHeteroReactionMaster",
                                           "calculation of heterogeneous reactions",
                                           DRT::Condition::ScatraHeteroReactionCondMaster,
                                           true,
                                           DRT::Condition::Surface));

  // insert condition definitions into global list of valid condition definitions
  condlist.push_back(scatraheteroreactionmasterline);
  condlist.push_back(scatraheteroreactionmastersurf);

  /*--------------------------------------------------------------------*/
  // conditions for calculation of calculation of heterogeneous reactions
  Teuchos::RCP<ConditionDefinition> scatraheteroreactionslaveline =
      Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA HETEROGENEOUS REACTION LINE CONDITIONS / SLAVE",
                                           "ScatraHeteroReactionSlave",
                                           "calculation of heterogeneous reactions",
                                           DRT::Condition::ScatraHeteroReactionCondSlave,
                                           true,
                                           DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> scatraheteroreactionslavesurf =
      Teuchos::rcp(new ConditionDefinition("DESIGN SCATRA HETEROGENEOUS REACTION SURF CONDITIONS / SLAVE",
                                           "ScatraHeteroReactionSlave",
                                           "calculation of heterogeneous reactions",
                                           DRT::Condition::ScatraHeteroReactionCondSlave,
                                           true,
                                           DRT::Condition::Surface));

  // insert condition definitions into global list of valid condition definitions
  condlist.push_back(scatraheteroreactionslaveline);
  condlist.push_back(scatraheteroreactionslavesurf);

  return;
}
