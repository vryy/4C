/*----------------------------------------------------------------------*/
/*! \file

\brief Input parameters for structure


\level 1
*/

/*----------------------------------------------------------------------*/

#include "4C_inpar_structure.hpp"

#include "4C_discretization_condition_definition.hpp"
#include "4C_io_geometry_type.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Inpar
{
  namespace STR
  {
    /*----------------------------------------------------------------------*/
    /*----------------------------------------------------------------------*/
    void SetValidTimeAdaptivityParameters(Teuchos::ParameterList& list)
    {
      using namespace Input;
      using Teuchos::setStringToIntegralParameter;
      using Teuchos::tuple;

      setStringToIntegralParameter<int>("KIND", "None", "Method for time step size adaptivity",
          tuple<std::string>("None", "ZienkiewiczXie", "JointExplicit",  //
              "AdamsBashforth2", "ExplicitEuler", "CentralDifference"),  // backward compatibility
          tuple<int>(Inpar::STR::timada_kind_none, Inpar::STR::timada_kind_zienxie,
              Inpar::STR::timada_kind_joint_explicit,  //
              Inpar::STR::timada_kind_ab2, Inpar::STR::timada_kind_expleuler,
              Inpar::STR::timada_kind_centraldiff),  // backward compatibility
          &list);

      Core::UTILS::DoubleParameter("OUTSYSPERIOD", 0.0,
          "Write system vectors (displacements, velocities, etc) every given period of time",
          &list);
      Core::UTILS::DoubleParameter(
          "OUTSTRPERIOD", 0.0, "Write stress/strain every given period of time", &list);
      Core::UTILS::DoubleParameter(
          "OUTENEPERIOD", 0.0, "Write energy every given period of time", &list);
      Core::UTILS::DoubleParameter(
          "OUTRESTPERIOD", 0.0, "Write restart data every given period of time", &list);
      Core::UTILS::IntParameter("OUTSIZEEVERY", 0, "Write step size every given time step", &list);

      Core::UTILS::DoubleParameter(
          "STEPSIZEMAX", 0.0, "Limit maximally permitted time step size (>0)", &list);
      Core::UTILS::DoubleParameter(
          "STEPSIZEMIN", 0.0, "Limit minimally allowed time step size (>0)", &list);
      Core::UTILS::DoubleParameter("SIZERATIOMAX", 0.0,
          "Limit maximally permitted change of time step size compared to previous size, important "
          "for multi-step schemes (>0)",
          &list);
      Core::UTILS::DoubleParameter("SIZERATIOMIN", 0.0,
          "Limit minimally permitted change of time step size compared to previous size, important "
          "for multi-step schemes (>0)",
          &list);
      Core::UTILS::DoubleParameter("SIZERATIOSCALE", 0.9,
          "This is a safety factor to scale theoretical optimal step size, should be lower than 1 "
          "and must be larger than 0",
          &list);

      setStringToIntegralParameter<int>("LOCERRNORM", "Vague",
          "Vector norm to treat error vector with",
          tuple<std::string>("Vague", "L1", "L2", "Rms", "Inf"),
          tuple<int>(Inpar::STR::norm_vague, Inpar::STR::norm_l1, Inpar::STR::norm_l2,
              Inpar::STR::norm_rms, Inpar::STR::norm_inf),
          &list);

      Core::UTILS::DoubleParameter("LOCERRTOL", 0.0, "Target local error tolerance (>0)", &list);
      Core::UTILS::IntParameter(
          "ADAPTSTEPMAX", 0, "Limit maximally allowed step size reduction attempts (>0)", &list);

      /// valid parameters for JOINT EXPLICIT

      Teuchos::ParameterList& jep = list.sublist("JOINT EXPLICIT", false, "");

      Core::UTILS::IntParameter(
          "LINEAR_SOLVER", -1, "number of linear solver used for auxiliary integrator", &jep);

      setStringToIntegralParameter<int>("INT_STRATEGY", "Standard",
          "global type of the used integration strategy", tuple<std::string>("Standard"),
          tuple<int>(int_standard), &jep);

      setStringToIntegralParameter<int>("DYNAMICTYP", "CentrDiff",
          "type of the specific auxiliary dynamic time integration scheme",
          tuple<std::string>("ExplicitEuler", "CentrDiff", "AdamsBashforth2", "AdamsBashforth4"),
          tuple<int>(dyna_expleuler, dyna_centrdiff, dyna_ab2, dyna_ab4), &jep);

      Core::UTILS::BoolParameter(
          "LUMPMASS", "No", "Lump the mass matrix for explicit time integration", &jep);

      setStringToIntegralParameter<int>("DAMPING", "No",
          "type of damping: (1) Rayleigh damping matrix and use it from M_DAMP x M + K_DAMP x K, "
          "(2) Material based and calculated in elements",
          tuple<std::string>("no", "No", "NO", "yes", "Yes", "YES", "Rayleigh", "Material"),
          tuple<int>(damp_none, damp_none, damp_none, damp_rayleigh, damp_rayleigh, damp_rayleigh,
              damp_rayleigh, damp_material),
          &jep);

      Core::UTILS::DoubleParameter("M_DAMP", -1.0, "", &jep);
      Core::UTILS::DoubleParameter("K_DAMP", -1.0, "", &jep);
    }



    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
    {
      using namespace Input;
      using Teuchos::setStringToIntegralParameter;
      using Teuchos::tuple;

      Teuchos::ParameterList& sdyn = list->sublist("STRUCTURAL DYNAMIC", false, "");

      setStringToIntegralParameter<int>("INT_STRATEGY", "Old",
          "global type of the used integration strategy", tuple<std::string>("Old", "Standard"),
          tuple<int>(int_old, int_standard), &sdyn);

      Core::UTILS::BoolParameter(
          "TIME_ADAPTIVITY", "No", "Enable adaptive time integration", &sdyn);

      setStringToIntegralParameter<int>("DYNAMICTYP", "GenAlpha",
          "type of the specific dynamic time integration scheme",
          tuple<std::string>("Statics", "GenAlpha", "GenAlphaLieGroup", "OneStepTheta",
              "ExplicitEuler", "CentrDiff", "AdamsBashforth2", "AdamsBashforth4", "EulerMaruyama",
              "EulerImpStoch"),
          tuple<int>(dyna_statics, dyna_genalpha, dyna_genalpha_liegroup, dyna_onesteptheta,
              dyna_expleuler, dyna_centrdiff, dyna_ab2, dyna_ab4, dyna_euma, dyna_euimsto),
          &sdyn);

      setStringToIntegralParameter<Inpar::STR::PreStress>("PRESTRESS", "none",
          "prestressing takes values none mulf material_iterative",
          tuple<std::string>("none", "None", "NONE", "mulf", "Mulf", "MULF", "Material_Iterative",
              "MATERIAL_ITERATIVE", "material_iterative"),
          tuple<Inpar::STR::PreStress>(Inpar::STR::PreStress::none, Inpar::STR::PreStress::none,
              Inpar::STR::PreStress::none, Inpar::STR::PreStress::mulf, Inpar::STR::PreStress::mulf,
              Inpar::STR::PreStress::mulf, Inpar::STR::PreStress::material_iterative,
              Inpar::STR::PreStress::material_iterative, Inpar::STR::PreStress::material_iterative),
          &sdyn);

      Core::UTILS::DoubleParameter(
          "PRESTRESSTIME", 0.0, "time to switch from pre to post stressing", &sdyn);

      Core::UTILS::DoubleParameter("PRESTRESSTOLDISP", 1e-9,
          "tolerance in the displacement norm during prestressing", &sdyn);
      Core::UTILS::IntParameter(
          "PRESTRESSMINLOADSTEPS", 0, "Minimum number of load steps during prestressing", &sdyn);

      // Output type
      Core::UTILS::IntParameter(
          "RESULTSEVRY", 1, "save displacements and contact forces every RESULTSEVRY steps", &sdyn);
      Core::UTILS::IntParameter(
          "RESEVRYERGY", 0, "write system energies every requested step", &sdyn);
      Core::UTILS::IntParameter(
          "RESTARTEVRY", 1, "write restart possibility every RESTARTEVRY steps", &sdyn);
      Core::UTILS::BoolParameter("CALC_ACC_ON_RESTART", "No",
          "Compute the initial state for a restart dynamics analysis", &sdyn);
      Core::UTILS::IntParameter("OUTPUT_STEP_OFFSET", 0,
          "An offset added to the current step to shift the steps to be written.", &sdyn);
      // Time loop control
      Core::UTILS::DoubleParameter("TIMESTEP", 0.05, "time step size", &sdyn);
      Core::UTILS::IntParameter("NUMSTEP", 200, "maximum number of steps", &sdyn);
      Core::UTILS::DoubleParameter("TIMEINIT", 0.0, "initial time", &sdyn);
      Core::UTILS::DoubleParameter("MAXTIME", 5.0, "maximum time", &sdyn);
      // Damping
      setStringToIntegralParameter<int>("DAMPING", "No",
          "type of damping: (1) Rayleigh damping matrix and use it from M_DAMP x M + K_DAMP x K, "
          "(2) Material based and calculated in elements",
          tuple<std::string>("no", "No", "NO", "yes", "Yes", "YES", "Rayleigh", "Material"),
          tuple<int>(damp_none, damp_none, damp_none, damp_rayleigh, damp_rayleigh, damp_rayleigh,
              damp_rayleigh, damp_material),
          &sdyn);
      Core::UTILS::DoubleParameter("M_DAMP", -1.0, "", &sdyn);
      Core::UTILS::DoubleParameter("K_DAMP", -1.0, "", &sdyn);

      Core::UTILS::DoubleParameter(
          "TOLDISP", 1.0E-10, "tolerance in the displacement norm for the newton iteration", &sdyn);
      setStringToIntegralParameter<int>("NORM_DISP", "Abs",
          "type of norm for displacement convergence check",
          tuple<std::string>("Abs", "Rel", "Mix"),
          tuple<int>(convnorm_abs, convnorm_rel, convnorm_mix), &sdyn);

      Core::UTILS::DoubleParameter(
          "TOLRES", 1.0E-08, "tolerance in the residual norm for the newton iteration", &sdyn);
      setStringToIntegralParameter<int>("NORM_RESF", "Abs",
          "type of norm for residual convergence check", tuple<std::string>("Abs", "Rel", "Mix"),
          tuple<int>(convnorm_abs, convnorm_rel, convnorm_mix), &sdyn);

      Core::UTILS::DoubleParameter(
          "TOLPRE", 1.0E-08, "tolerance in pressure norm for the newton iteration", &sdyn);
      setStringToIntegralParameter<int>("NORM_PRES", "Abs",
          "type of norm for pressure convergence check", tuple<std::string>("Abs"),
          tuple<int>(convnorm_abs), &sdyn);

      Core::UTILS::DoubleParameter("TOLINCO", 1.0E-08,
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

      Core::UTILS::IntParameter("STC_LAYER", 1, "number of STC layers for multilayer case", &sdyn);

      Core::UTILS::DoubleParameter("PTCDT", 0.1,
          "pseudo time step for pseudo transient continuation (PTC) stabilized Newton procedure",
          &sdyn);

      Core::UTILS::DoubleParameter("TOLCONSTR", 1.0E-08,
          "tolerance in the constr error norm for the newton iteration", &sdyn);

      Core::UTILS::DoubleParameter("TOLCONSTRINCR", 1.0E-08,
          "tolerance in the constr lm incr norm for the newton iteration", &sdyn);

      Core::UTILS::IntParameter("MAXITER", 50,
          "maximum number of iterations allowed for Newton-Raphson iteration before failure",
          &sdyn);
      Core::UTILS::IntParameter("MINITER", 0,
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

      Core::UTILS::IntParameter("MAXDIVCONREFINEMENTLEVEL", 10,
          "number of times timestep is halved in case nonlinear solver diverges", &sdyn);

      setStringToIntegralParameter<int>("NLNSOL", "fullnewton", "Nonlinear solution technique",
          tuple<std::string>("vague", "fullnewton", "modnewton", "lsnewton", "ptc",
              "newtonlinuzawa", "augmentedlagrange", "NoxNewtonLineSearch", "noxgeneral", "noxnln",
              "singlestep"),
          tuple<int>(soltech_vague, soltech_newtonfull, soltech_newtonmod, soltech_newtonls,
              soltech_ptc, soltech_newtonuzawalin, soltech_newtonuzawanonlin,
              soltech_noxnewtonlinesearch, soltech_noxgeneral, soltech_nox_nln, soltech_singlestep),
          &sdyn);

      Core::UTILS::IntParameter("LSMAXITER", 30, "maximum number of line search steps", &sdyn);
      Core::UTILS::DoubleParameter(
          "ALPHA_LS", 0.5, "step reduction factor alpha in (Newton) line search scheme", &sdyn);
      Core::UTILS::DoubleParameter(
          "SIGMA_LS", 1.e-4, "sufficient descent factor in (Newton) line search scheme", &sdyn);

      setStringToIntegralParameter<int>("MATERIALTANGENT", "analytical",
          "way of evaluating the constitutive matrix",
          tuple<std::string>("analytical", "finitedifferences"), tuple<int>(0, 1), &sdyn);

      Core::UTILS::BoolParameter(
          "LOADLIN", "No", "Use linearization of external follower load in Newton", &sdyn);

      setStringToIntegralParameter<int>("MASSLIN", "No", "Application of nonlinear inertia terms",
          tuple<std::string>("No", "no", "Standard", "standard", "Rotations", "rotations"),

          tuple<int>(ml_none, ml_none, ml_standard, ml_standard, ml_rotations, ml_rotations),
          &sdyn);

      Core::UTILS::BoolParameter("NEGLECTINERTIA", "No", "Neglect inertia", &sdyn);

      // Since predictor "none" would be misleading, the usage of no predictor is called vague.
      setStringToIntegralParameter<int>("PREDICT", "ConstDis", "Type of predictor",
          tuple<std::string>("Vague", "ConstDis", "ConstVel", "ConstAcc", "ConstDisVelAcc",
              "TangDis", "TangDisConstFext", "ConstDisPres", "ConstDisVelAccPres"),
          tuple<int>(pred_vague, pred_constdis, pred_constvel, pred_constacc, pred_constdisvelacc,
              pred_tangdis, pred_tangdis_constfext, pred_constdispres, pred_constdisvelaccpres),
          &sdyn);

      // Uzawa iteration for constraint systems
      Core::UTILS::DoubleParameter("UZAWAPARAM", 1.0,
          "Parameter for Uzawa algorithm dealing with lagrange multipliers", &sdyn);
      Core::UTILS::DoubleParameter(
          "UZAWATOL", 1.0E-8, "Tolerance for iterative solve with Uzawa algorithm", &sdyn);
      Core::UTILS::IntParameter("UZAWAMAXITER", 50,
          "maximum number of iterations allowed for uzawa algorithm before failure going to next "
          "newton step",
          &sdyn);
      setStringToIntegralParameter<int>("UZAWAALGO", "direct", "",
          tuple<std::string>("uzawa", "simple", "direct"),
          tuple<int>(consolve_uzawa, consolve_simple, consolve_direct), &sdyn);

      // convergence criteria adaptivity
      Core::UTILS::BoolParameter("ADAPTCONV", "No",
          "Switch on adaptive control of linear solver tolerance for nonlinear solution", &sdyn);
      Core::UTILS::DoubleParameter("ADAPTCONV_BETTER", 0.1,
          "The linear solver shall be this much better than the current nonlinear residual in the "
          "nonlinear convergence limit",
          &sdyn);

      Core::UTILS::BoolParameter(
          "LUMPMASS", "No", "Lump the mass matrix for explicit time integration", &sdyn);

      Core::UTILS::BoolParameter("MODIFIEDEXPLEULER", "Yes",
          "Use the modified explicit Euler time integration scheme", &sdyn);

      // linear solver id used for structural problems
      Core::UTILS::IntParameter(
          "LINEAR_SOLVER", -1, "number of linear solver used for structural problems", &sdyn);

      // where the geometry comes from
      setStringToIntegralParameter<int>("GEOMETRY", "full", "How the geometry is specified",
          tuple<std::string>("full", "box", "file"),
          tuple<int>(Core::IO::geometry_full, Core::IO::geometry_box, Core::IO::geometry_file),
          &sdyn);

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
      Core::UTILS::IntParameter("STARTFUNCNO", -1, "Function for Initial displacement", &sdyn);

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
      Core::UTILS::DoubleParameter("BETA", -1.0, "Generalised-alpha factor in (0,1/2]", &genalpha);
      Core::UTILS::DoubleParameter("GAMMA", -1.0, "Generalised-alpha factor in (0,1]", &genalpha);
      Core::UTILS::DoubleParameter("ALPHA_M", -1.0, "Generalised-alpha factor in [0,1)", &genalpha);
      Core::UTILS::DoubleParameter("ALPHA_F", -1.0, "Generalised-alpha factor in [0,1)", &genalpha);
      Core::UTILS::DoubleParameter("RHO_INF", 1.0,
          "Spectral radius for generalised-alpha time integration, valid range is [0,1]",
          &genalpha);

      /*----------------------------------------------------------------------*/
      /* parameters for one-step-theta structural integrator */
      Teuchos::ParameterList& onesteptheta = sdyn.sublist("ONESTEPTHETA", false, "");

      Core::UTILS::DoubleParameter("THETA", 0.5, "One-step-theta factor in (0,1]", &onesteptheta);
    }



    void SetValidConditions(
        std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
    {
      using namespace Input;

      /*--------------------------------------------------------------------*/

      // structural Robin spring dashpot boundary condition (spring and dashpot in parallel) - mhv
      // 08/16

      auto robinspringdashpotsurf = Teuchos::rcp(
          new Core::Conditions::ConditionDefinition("DESIGN SURF ROBIN SPRING DASHPOT CONDITIONS",
              "RobinSpringDashpot", "Robin Spring Dashpot", Core::Conditions::RobinSpringDashpot,
              true, Core::Conditions::geometry_type_surface));

      auto robinspringdashpotpoint = Teuchos::rcp(
          new Core::Conditions::ConditionDefinition("DESIGN POINT ROBIN SPRING DASHPOT CONDITIONS",
              "RobinSpringDashpot", "Robin Spring Dashpot", Core::Conditions::RobinSpringDashpot,
              true, Core::Conditions::geometry_type_point));

      std::vector<Teuchos::RCP<Input::LineComponent>> robinspringdashpotcomp;

      robinspringdashpotcomp.emplace_back(Teuchos::rcp(new Input::SeparatorComponent("NUMDOF")));
      robinspringdashpotcomp.emplace_back(Teuchos::rcp(new Input::IntComponent("numdof")));

      robinspringdashpotcomp.emplace_back(Teuchos::rcp(new Input::SeparatorComponent("ONOFF")));
      robinspringdashpotcomp.emplace_back(Teuchos::rcp(new Input::IntVectorComponent("onoff", 3)));

      robinspringdashpotcomp.emplace_back(Teuchos::rcp(new Input::SeparatorComponent("STIFF")));
      robinspringdashpotcomp.emplace_back(Teuchos::rcp(new Input::RealVectorComponent("stiff", 3)));

      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new Input::SeparatorComponent("TIMEFUNCTSTIFF")));
      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new Input::IntVectorComponent("funct_stiff", 3)));

      robinspringdashpotcomp.emplace_back(Teuchos::rcp(new Input::SeparatorComponent("VISCO")));
      robinspringdashpotcomp.emplace_back(Teuchos::rcp(new Input::RealVectorComponent("visco", 3)));

      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new Input::SeparatorComponent("TIMEFUNCTVISCO")));
      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new Input::IntVectorComponent("funct_visco", 3)));

      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new Input::SeparatorComponent("DISPLOFFSET")));
      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new Input::RealVectorComponent("disploffset", 3)));

      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new Input::SeparatorComponent("TIMEFUNCTDISPLOFFSET")));
      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new Input::IntVectorComponent("funct_disploffset", 3)));

      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new Input::SeparatorComponent("FUNCTNONLINSTIFF")));
      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new Input::IntVectorComponent("funct_nonlinstiff", 3)));

      robinspringdashpotcomp.emplace_back(Teuchos::rcp(new Input::SeparatorComponent("DIRECTION")));
      robinspringdashpotcomp.emplace_back(Teuchos::rcp(new Input::SelectionComponent("direction",
          "xyz", Teuchos::tuple<std::string>("xyz", "refsurfnormal", "cursurfnormal"),
          Teuchos::tuple<std::string>("xyz", "refsurfnormal", "cursurfnormal"), false)));

      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new Input::SeparatorComponent("COUPLING", "", true)));
      robinspringdashpotcomp.emplace_back(
          Teuchos::rcp(new Input::IntComponent("coupling id", {0, true, true})));

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

      Teuchos::RCP<Core::Conditions::ConditionDefinition> springdashpotcoupcond =
          Teuchos::rcp(new Core::Conditions::ConditionDefinition(
              "DESIGN SURF ROBIN SPRING DASHPOT COUPLING CONDITIONS", "RobinSpringDashpotCoupling",
              "RobinSpring Dashpot Coupling", Core::Conditions::RobinSpringDashpotCoupling, true,
              Core::Conditions::geometry_type_surface));

      springdashpotcoupcond->AddComponent(Teuchos::rcp(new Input::IntComponent("coupling id")));

      condlist.push_back(springdashpotcoupcond);


      /*--------------------------------------------------------------------*/
      // surfactant

      Teuchos::RCP<Core::Conditions::ConditionDefinition> surfactant =
          Teuchos::rcp(new Core::Conditions::ConditionDefinition("SURFACTANT CONDITIONS",
              "SurfaceStress", "Surface Stress (surfactant)", Core::Conditions::Surfactant, true,
              Core::Conditions::geometry_type_surface));

      surfactant->AddComponent(Teuchos::rcp(new Input::IntComponent("funct", {0, true, true})));
      Input::AddNamedReal(surfactant, "k1xCbulk");
      Input::AddNamedReal(surfactant, "k2");
      Input::AddNamedReal(surfactant, "m1");
      Input::AddNamedReal(surfactant, "m2");
      Input::AddNamedReal(surfactant, "gamma_0");
      Input::AddNamedReal(surfactant, "gamma_min");

      condlist.push_back(surfactant);
    }

  }  // end of namespace STR
}  // end of namespace Inpar

FOUR_C_NAMESPACE_CLOSE
