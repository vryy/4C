/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameters for structure


\level 1
*/

/*----------------------------------------------------------------------*/

#include "baci_inpar_structure.H"

#include "baci_inpar.H"
#include "baci_inpar_validparameters.H"
#include "baci_lib_conditiondefinition.H"

namespace INPAR
{
  namespace STR
  {
    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void SetValidTimeAdaptivityParameters(Teuchos::ParameterList& list)
    {
      using namespace DRT::INPUT;
      using Teuchos::setStringToIntegralParameter;
      using Teuchos::tuple;

      setStringToIntegralParameter<int>("KIND", "None", "Method for time step size adaptivity",
          tuple<std::string>(
              "None", "ZienkiewiczXie", "AdamsBashforth2", "ExplicitEuler", "CentralDifference"),
          tuple<int>(INPAR::STR::timada_kind_none, INPAR::STR::timada_kind_zienxie,
              INPAR::STR::timada_kind_ab2, INPAR::STR::timada_kind_expleuler,
              INPAR::STR::timada_kind_centraldiff),
          &list);

      DoubleParameter("OUTSYSPERIOD", 0.0,
          "Write system vectors (displacements, velocities, etc) every given period of time",
          &list);
      DoubleParameter("OUTSTRPERIOD", 0.0, "Write stress/strain every given period of time", &list);
      DoubleParameter("OUTENEPERIOD", 0.0, "Write energy every given period of time", &list);
      DoubleParameter("OUTRESTPERIOD", 0.0, "Write restart data every given period of time", &list);
      IntParameter("OUTSIZEEVERY", 0, "Write step size every given time step", &list);

      DoubleParameter("STEPSIZEMAX", 0.0, "Limit maximally permitted time step size (>0)", &list);
      DoubleParameter("STEPSIZEMIN", 0.0, "Limit minimally allowed time step size (>0)", &list);
      DoubleParameter("SIZERATIOMAX", 0.0,
          "Limit maximally permitted change of time step size compared to previous size, important "
          "for multi-step schemes (>0)",
          &list);
      DoubleParameter("SIZERATIOMIN", 0.0,
          "Limit minimally permitted change of time step size compared to previous size, important "
          "for multi-step schemes (>0)",
          &list);
      DoubleParameter("SIZERATIOSCALE", 0.9,
          "This is a safety factor to scale theoretical optimal step size, should be lower than 1 "
          "and must be larger than 0",
          &list);

      setStringToIntegralParameter<int>("LOCERRNORM", "Vague",
          "Vector norm to treat error vector with",
          tuple<std::string>("Vague", "L1", "L2", "Rms", "Inf"),
          tuple<int>(INPAR::STR::norm_vague, INPAR::STR::norm_l1, INPAR::STR::norm_l2,
              INPAR::STR::norm_rms, INPAR::STR::norm_inf),
          &list);

      DoubleParameter("LOCERRTOL", 0.0, "Target local error tolerance (>0)", &list);
      IntParameter(
          "ADAPTSTEPMAX", 0, "Limit maximally allowed step size reduction attempts (>0)", &list);
    }



    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
    {
      using namespace DRT::INPUT;
      using Teuchos::setStringToIntegralParameter;
      using Teuchos::tuple;

      Teuchos::ParameterList& sdyn = list->sublist("STRUCTURAL DYNAMIC", false, "");

      setStringToIntegralParameter<int>("INT_STRATEGY", "Old",
          "global type of the used integration strategy", tuple<std::string>("Old", "Standard"),
          tuple<int>(int_old, int_standard), &sdyn);

      setStringToIntegralParameter<int>("DYNAMICTYP", "GenAlpha",
          "type of the specific dynamic time integration scheme",
          tuple<std::string>("Statics", "GenAlpha", "GenAlphaLieGroup", "OneStepTheta", "GEMM",
              "ExplicitEuler", "CentrDiff", "AdamsBashforth2", "AdamsBashforth4", "EulerMaruyama",
              "EulerImpStoch"),
          tuple<int>(dyna_statics, dyna_genalpha, dyna_genalpha_liegroup, dyna_onesteptheta,
              dyna_gemm, dyna_expleuler, dyna_centrdiff, dyna_ab2, dyna_ab4, dyna_euma,
              dyna_euimsto),
          &sdyn);

      setStringToIntegralParameter<INPAR::STR::PreStress>("PRESTRESS", "none",
          "prestressing takes values none mulf material_iterative",
          tuple<std::string>("none", "None", "NONE", "mulf", "Mulf", "MULF", "Material_Iterative",
              "MATERIAL_ITERATIVE", "material_iterative"),
          tuple<INPAR::STR::PreStress>(INPAR::STR::PreStress::none, INPAR::STR::PreStress::none,
              INPAR::STR::PreStress::none, INPAR::STR::PreStress::mulf, INPAR::STR::PreStress::mulf,
              INPAR::STR::PreStress::mulf, INPAR::STR::PreStress::material_iterative,
              INPAR::STR::PreStress::material_iterative, INPAR::STR::PreStress::material_iterative),
          &sdyn);

      DoubleParameter("PRESTRESSTIME", 0.0, "time to switch from pre to post stressing", &sdyn);

      DoubleParameter("PRESTRESSTOLDISP", 1e-9,
          "tolerance in the displacement norm during prestressing", &sdyn);
      IntParameter(
          "PRESTRESSMINLOADSTEPS", 0, "Minimum number of load steps during prestressing", &sdyn);

      // Output type
      IntParameter(
          "RESULTSEVRY", 1, "save displacements and contact forces every RESULTSEVRY steps", &sdyn);
      IntParameter("RESEVRYERGY", 0, "write system energies every requested step", &sdyn);
      IntParameter("RESTARTEVRY", 1, "write restart possibility every RESTARTEVRY steps", &sdyn);
      BoolParameter("CALC_ACC_ON_RESTART", "No",
          "Compute the initial state for a restart dynamics analysis", &sdyn);
      IntParameter("OUTPUT_STEP_OFFSET", 0,
          "An offset added to the current step to shift the steps to be written.", &sdyn);
      // Time loop control
      DoubleParameter("TIMESTEP", 0.05, "time step size", &sdyn);
      IntParameter("NUMSTEP", 200, "maximum number of steps", &sdyn);
      DoubleParameter("TIMEINIT", 0.0, "initial time", &sdyn);
      DoubleParameter("MAXTIME", 5.0, "maximum time", &sdyn);
      // Damping
      setStringToIntegralParameter<int>("DAMPING", "No",
          "type of damping: (1) Rayleigh damping matrix and use it from M_DAMP x M + K_DAMP x K, "
          "(2) Material based and calculated in elements",
          tuple<std::string>("no", "No", "NO", "yes", "Yes", "YES", "Rayleigh", "Material"),
          tuple<int>(damp_none, damp_none, damp_none, damp_rayleigh, damp_rayleigh, damp_rayleigh,
              damp_rayleigh, damp_material),
          &sdyn);
      DoubleParameter("M_DAMP", -1.0, "", &sdyn);
      DoubleParameter("K_DAMP", -1.0, "", &sdyn);

      DoubleParameter(
          "TOLDISP", 1.0E-10, "tolerance in the displacement norm for the newton iteration", &sdyn);
      setStringToIntegralParameter<int>("NORM_DISP", "Abs",
          "type of norm for displacement convergence check",
          tuple<std::string>("Abs", "Rel", "Mix"),
          tuple<int>(convnorm_abs, convnorm_rel, convnorm_mix), &sdyn);

      DoubleParameter(
          "TOLRES", 1.0E-08, "tolerance in the residual norm for the newton iteration", &sdyn);
      setStringToIntegralParameter<int>("NORM_RESF", "Abs",
          "type of norm for residual convergence check", tuple<std::string>("Abs", "Rel", "Mix"),
          tuple<int>(convnorm_abs, convnorm_rel, convnorm_mix), &sdyn);

      DoubleParameter(
          "TOLPRE", 1.0E-08, "tolerance in pressure norm for the newton iteration", &sdyn);
      setStringToIntegralParameter<int>("NORM_PRES", "Abs",
          "type of norm for pressure convergence check", tuple<std::string>("Abs"),
          tuple<int>(convnorm_abs), &sdyn);

      DoubleParameter("TOLINCO", 1.0E-08,
          "tolerance in the incompressible residual norm for the newton iteration", &sdyn);
      setStringToIntegralParameter<int>("NORM_INCO", "Abs",
          "type of norm for incompressible residual convergence check", tuple<std::string>("Abs"),
          tuple<int>(convnorm_abs), &sdyn);

      setStringToIntegralParameter<int>("NORMCOMBI_DISPPRES", "And",
          "binary operator to combine pressure and displacement values",
          tuple<std::string>("And", "Or"), tuple<int>(bop_and, bop_or), &sdyn);

      setStringToIntegralParameter<int>("NORMCOMBI_RESFINCO", "And",
          "binary operator to combine force and incompressible residual",
          tuple<std::string>("And", "Or"), tuple<int>(bop_and, bop_or), &sdyn);

      setStringToIntegralParameter<int>("NORMCOMBI_RESFDISP", "And",
          "binary operator to combine displacement and residual force values",
          tuple<std::string>("And", "Or"), tuple<int>(bop_and, bop_or), &sdyn);

      setStringToIntegralParameter<int>("STC_SCALING", "no",
          "Scaled director conditioning for thin shell structures",
          tuple<std::string>("no", "No", "NO", "Symmetric", "Right"),
          tuple<int>(stc_none, stc_none, stc_none, stc_currsym, stc_curr), &sdyn);

      IntParameter("STC_LAYER", 1, "number of STC layers for multilayer case", &sdyn);

      DoubleParameter("PTCDT", 0.1,
          "pseudo time step for pseudo transient continuation (PTC) stabilized Newton procedure",
          &sdyn);

      DoubleParameter("TOLCONSTR", 1.0E-08,
          "tolerance in the constr error norm for the newton iteration", &sdyn);

      DoubleParameter("TOLCONSTRINCR", 1.0E-08,
          "tolerance in the constr lm incr norm for the newton iteration", &sdyn);

      IntParameter("MAXITER", 50,
          "maximum number of iterations allowed for Newton-Raphson iteration before failure",
          &sdyn);
      IntParameter("MINITER", 0,
          "minimum number of iterations to be done within Newton-Raphson loop", &sdyn);
      setStringToIntegralParameter<int>("ITERNORM", "L2", "type of norm to be applied to residuals",
          tuple<std::string>("L1", "L2", "Rms", "Inf"),
          tuple<int>(norm_l1, norm_l2, norm_rms, norm_inf), &sdyn);

      setStringToIntegralParameter<int>("DIVERCONT", "stop",
          "What to do with time integration when Newton-Raphson iteration failed",
          tuple<std::string>("stop", "continue", "repeat_step", "halve_step", "adapt_step",
              "rand_adapt_step", "rand_adapt_step_ele_err", "repeat_simulation",
              "adapt_penaltycontact", "adapt_3D0Dptc_ele_err"),
          tuple<int>(divcont_stop, divcont_continue, divcont_repeat_step, divcont_halve_step,
              divcont_adapt_step, divcont_rand_adapt_step, divcont_rand_adapt_step_ele_err,
              divcont_repeat_simulation, divcont_adapt_penaltycontact,
              divcont_adapt_3D0Dptc_ele_err),
          &sdyn);

      IntParameter("MAXDIVCONREFINEMENTLEVEL", 10,
          "number of times timestep is halved in case nonlinear solver diverges", &sdyn);

      setStringToIntegralParameter<int>("NLNSOL", "fullnewton", "Nonlinear solution technique",
          tuple<std::string>("vague", "fullnewton", "modnewton", "lsnewton", "ptc",
              "newtonlinuzawa", "augmentedlagrange", "NoxNewtonLineSearch", "noxgeneral", "noxnln",
              "singlestep"),
          tuple<int>(soltech_vague, soltech_newtonfull, soltech_newtonmod, soltech_newtonls,
              soltech_ptc, soltech_newtonuzawalin, soltech_newtonuzawanonlin,
              soltech_noxnewtonlinesearch, soltech_noxgeneral, soltech_nox_nln, soltech_singlestep),
          &sdyn);

      IntParameter("LSMAXITER", 30, "maximum number of line search steps", &sdyn);
      DoubleParameter(
          "ALPHA_LS", 0.5, "step reduction factor alpha in (Newton) line search scheme", &sdyn);
      DoubleParameter(
          "SIGMA_LS", 1.e-4, "sufficient descent factor in (Newton) line search scheme", &sdyn);

      setStringToIntegralParameter<int>("MATERIALTANGENT", "analytical",
          "way of evaluating the constitutive matrix",
          tuple<std::string>("analytical", "finitedifferences"), tuple<int>(0, 1), &sdyn);

      BoolParameter(
          "LOADLIN", "No", "Use linearization of external follower load in Newton", &sdyn);

      setStringToIntegralParameter<int>("MASSLIN", "No", "Application of nonlinear inertia terms",
          tuple<std::string>("No", "no", "Standard", "standard", "Rotations", "rotations"),

          tuple<int>(ml_none, ml_none, ml_standard, ml_standard, ml_rotations, ml_rotations),
          &sdyn);

      BoolParameter("NEGLECTINERTIA", "No", "Neglect inertia", &sdyn);

      // Since predictor "none" would be misleading, the usage of no predictor is called vague.
      setStringToIntegralParameter<int>("PREDICT", "ConstDis", "Type of predictor",
          tuple<std::string>("Vague", "ConstDis", "ConstVel", "ConstAcc", "ConstDisVelAcc",
              "TangDis", "TangDisConstFext", "ConstDisPres", "ConstDisVelAccPres"),
          tuple<int>(pred_vague, pred_constdis, pred_constvel, pred_constacc, pred_constdisvelacc,
              pred_tangdis, pred_tangdis_constfext, pred_constdispres, pred_constdisvelaccpres),
          &sdyn);

      // Uzawa iteration for constraint systems
      DoubleParameter("UZAWAPARAM", 1.0,
          "Parameter for Uzawa algorithm dealing with lagrange multipliers", &sdyn);
      DoubleParameter(
          "UZAWATOL", 1.0E-8, "Tolerance for iterative solve with Uzawa algorithm", &sdyn);
      IntParameter("UZAWAMAXITER", 50,
          "maximum number of iterations allowed for uzawa algorithm before failure going to next "
          "newton step",
          &sdyn);
      setStringToIntegralParameter<int>("UZAWAALGO", "direct", "",
          tuple<std::string>("uzawa", "simple", "direct"),
          tuple<int>(consolve_uzawa, consolve_simple, consolve_direct), &sdyn);

      // convergence criteria adaptivity
      BoolParameter("ADAPTCONV", "No",
          "Switch on adaptive control of linear solver tolerance for nonlinear solution", &sdyn);
      DoubleParameter("ADAPTCONV_BETTER", 0.1,
          "The linear solver shall be this much better than the current nonlinear residual in the "
          "nonlinear convergence limit",
          &sdyn);

      BoolParameter("LUMPMASS", "No", "Lump the mass matrix for explicit time integration", &sdyn);

      BoolParameter("MODIFIEDEXPLEULER", "Yes",
          "Use the modified explicit Euler time integration scheme", &sdyn);

      // linear solver id used for structural problems
      IntParameter(
          "LINEAR_SOLVER", -1, "number of linear solver used for structural problems", &sdyn);

      // where the geometry comes from
      setStringToIntegralParameter<int>("GEOMETRY", "full", "How the geometry is specified",
          tuple<std::string>("full", "box", "file"),
          tuple<int>(INPAR::geometry_full, INPAR::geometry_box, INPAR::geometry_file), &sdyn);

      setStringToIntegralParameter<int>("MIDTIME_ENERGY_TYPE", "vague",
          "Specify the mid-averaging type for the structural energy contributions",
          tuple<std::string>("vague", "imrLike", "trLike"),
          tuple<int>(midavg_vague, midavg_imrlike, midavg_trlike), &sdyn);

      // Initial displacement
      setStringToIntegralParameter<int>("INITIALDISP", "zero_displacement",
          "Initial displacement for structure problem",
          tuple<std::string>("zero_displacement", "displacement_by_function"),
          tuple<int>(initdisp_zero_disp, initdisp_disp_by_function), &sdyn);

      // Function to evaluate initial displacement
      IntParameter("STARTFUNCNO", -1, "Function for Initial displacement", &sdyn);

      // Flag to (de)activate error calculations
      setStringToIntegralParameter<int>("CALCERROR", "no",
          "Flag to (de)activate error calculations", tuple<std::string>("no", "byfunct"),
          tuple<int>(no_error_calculation, byfunct), &sdyn);

      // Function number to calculate the error
      IntParameter("CALCERRORFUNCNO", -1, "Function for Error Calculation", &sdyn);

      /*--------------------------------------------------------------------*/
      /* parameters for time step size adaptivity in structural dynamics */
      Teuchos::ParameterList& tap = sdyn.sublist("TIMEADAPTIVITY", false, "");
      SetValidTimeAdaptivityParameters(tap);

      /*----------------------------------------------------------------------*/
      /* parameters for generalised-alpha structural integrator */
      Teuchos::ParameterList& genalpha = sdyn.sublist("GENALPHA", false, "");

      setStringToIntegralParameter<int>("GENAVG", "TrLike", "mid-average type of internal forces",
          tuple<std::string>("Vague", "ImrLike", "TrLike"),
          tuple<int>(midavg_vague, midavg_imrlike, midavg_trlike), &genalpha);
      DoubleParameter("BETA", -1.0, "Generalised-alpha factor in (0,1/2]", &genalpha);
      DoubleParameter("GAMMA", -1.0, "Generalised-alpha factor in (0,1]", &genalpha);
      DoubleParameter("ALPHA_M", -1.0, "Generalised-alpha factor in [0,1)", &genalpha);
      DoubleParameter("ALPHA_F", -1.0, "Generalised-alpha factor in [0,1)", &genalpha);
      DoubleParameter("RHO_INF", 1.0,
          "Spectral radius for generalised-alpha time integration, valid range is [0,1]",
          &genalpha);

      /*----------------------------------------------------------------------*/
      /* parameters for one-step-theta structural integrator */
      Teuchos::ParameterList& onesteptheta = sdyn.sublist("ONESTEPTHETA", false, "");

      DoubleParameter("THETA", 0.5, "One-step-theta factor in (0,1]", &onesteptheta);


      /*----------------------------------------------------------------------*/
      /* parameters for generalised-energy-momentum structural integrator */
      Teuchos::ParameterList& gemm = sdyn.sublist("GEMM", false, "");

      DoubleParameter("BETA", 0.25, "Generalised-alpha factor in (0,0.5]", &gemm);
      DoubleParameter("GAMMA", 0.5, "Generalised-alpha factor in (0,1]", &gemm);
      DoubleParameter("ALPHA_M", 0.5, "Generalised-alpha factor in [0,1)", &gemm);
      DoubleParameter("ALPHA_F", 0.5, "Generalised-alpha factor in [0,1)", &gemm);
      DoubleParameter("XI", 0.0, "generalisation factor in [0,1)", &gemm);
    }



    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
    {
      using namespace DRT::INPUT;

      /*--------------------------------------------------------------------*/

      // structural Robin spring dashpot boundary condition (spring and dashpot in parallel) - mhv
      // 08/16

      auto robinspringdashpotsurf =
          Teuchos::rcp(new ConditionDefinition("DESIGN SURF ROBIN SPRING DASHPOT CONDITIONS",
              "RobinSpringDashpot", "Robin Spring Dashpot", DRT::Condition::RobinSpringDashpot,
              true, DRT::Condition::Surface));

      auto robinspringdashpotpoint = Teuchos::rcp(new ConditionDefinition(
          "DESIGN POINT ROBIN SPRING DASHPOT CONDITIONS", "RobinSpringDashpot",
          "Robin Spring Dashpot", DRT::Condition::RobinSpringDashpot, true, DRT::Condition::Point));

      std::vector<Teuchos::RCP<::INPUT::LineComponent>> robinspringdashpotcomp;

      robinspringdashpotcomp.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("NUMDOF")));
      robinspringdashpotcomp.emplace_back(Teuchos::rcp(new ::INPUT::IntComponent("numdof")));

      robinspringdashpotcomp.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("ONOFF")));
      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new ::INPUT::IntVectorComponent("onoff", 3)));

      robinspringdashpotcomp.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("STIFF")));
      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new ::INPUT::RealVectorComponent("stiff", 3)));

      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new ::INPUT::SeparatorComponent("TIMEFUNCTSTIFF")));
      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new ::INPUT::IntVectorComponent("funct_stiff", 3)));

      robinspringdashpotcomp.emplace_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("VISCO")));
      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new ::INPUT::RealVectorComponent("visco", 3)));

      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new ::INPUT::SeparatorComponent("TIMEFUNCTVISCO")));
      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new ::INPUT::IntVectorComponent("funct_visco", 3)));

      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new ::INPUT::SeparatorComponent("DISPLOFFSET")));
      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new ::INPUT::RealVectorComponent("disploffset", 3)));

      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new ::INPUT::SeparatorComponent("TIMEFUNCTDISPLOFFSET")));
      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new ::INPUT::IntVectorComponent("funct_disploffset", 3)));

      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new ::INPUT::SeparatorComponent("FUNCTNONLINSTIFF")));
      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new ::INPUT::IntVectorComponent("funct_nonlinstiff", 3)));

      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new ::INPUT::SeparatorComponent("DIRECTION")));
      robinspringdashpotcomp.emplace_back(Teuchos::rcp(new ::INPUT::SelectionComponent("direction",
          "xyz", Teuchos::tuple<std::string>("xyz", "refsurfnormal", "cursurfnormal"),
          Teuchos::tuple<std::string>("xyz", "refsurfnormal", "cursurfnormal"), false)));

      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new ::INPUT::SeparatorComponent("COUPLING", "", true)));
      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new ::INPUT::IntVectorComponent("coupling id", 1, {0, true, true})));

      for (const auto& comp : robinspringdashpotcomp)
      {
        robinspringdashpotsurf->AddComponent(comp);
        robinspringdashpotpoint->AddComponent(comp);
      }

      condlist.emplace_back(robinspringdashpotsurf);
      condlist.emplace_back(robinspringdashpotpoint);

      /*--------------------------------------------------------------------*/
      // surface coupling for spring dashpot DIRECTION cursurfnormal
      // pfaller Apr15

      Teuchos::RCP<ConditionDefinition> springdashpotcoupcond = Teuchos::rcp(
          new ConditionDefinition("DESIGN SURF ROBIN SPRING DASHPOT COUPLING CONDITIONS",
              "RobinSpringDashpotCoupling", "RobinSpring Dashpot Coupling",
              DRT::Condition::RobinSpringDashpotCoupling, true, DRT::Condition::Surface));

      springdashpotcoupcond->AddComponent(Teuchos::rcp(new ::INPUT::IntComponent("coupling id")));

      condlist.push_back(springdashpotcoupcond);


      /*--------------------------------------------------------------------*/
      // surfactant

      Teuchos::RCP<ConditionDefinition> surfactant = Teuchos::rcp(new ConditionDefinition(
          "SURFACTANT CONDITIONS", "SurfaceStress", "Surface Stress (surfactant)",
          DRT::Condition::Surfactant, true, DRT::Condition::Surface));

      surfactant->AddComponent(Teuchos::rcp(new ::INPUT::IntComponent("funct", {0, true, true})));
      ::INPUT::AddNamedReal(surfactant, "k1xCbulk");
      ::INPUT::AddNamedReal(surfactant, "k2");
      ::INPUT::AddNamedReal(surfactant, "m1");
      ::INPUT::AddNamedReal(surfactant, "m2");
      ::INPUT::AddNamedReal(surfactant, "gamma_0");
      ::INPUT::AddNamedReal(surfactant, "gamma_min");

      condlist.push_back(surfactant);
    }

  }  // end of namespace STR
}  // end of namespace INPAR
