/*----------------------------------------------------------------------*/
/*!
\file inpar_structure.cpp

\brief Input parameters for combustion

\maintainer Alexander Seitz

\level 1
*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_structure.H"
#include "inpar.H"
#include "../drt_lib/drt_conditiondefinition.H"


namespace INPAR
{
  namespace STR
  {

  /*----------------------------------------------------------------------*/
  /*----------------------------------------------------------------------*/
  void SetValidTimeAdaptivityParameters(Teuchos::ParameterList& list)
  {
    using namespace DRT::INPUT;
    using Teuchos::tuple;
    using Teuchos::setStringToIntegralParameter;

    setStringToIntegralParameter<int>(
      "KIND","None","Method for time step size adapivity",
      tuple<std::string>(
        "None",
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
      &list);

    DoubleParameter("OUTSYSPERIOD", 0.0, "Write system vectors (displacements, velocities, etc) every given period of time", &list);
    DoubleParameter("OUTSTRPERIOD", 0.0, "Write stress/strain every given period of time", &list);
    DoubleParameter("OUTENEPERIOD", 0.0, "Write energy every given period of time", &list);
    DoubleParameter("OUTRESTPERIOD", 0.0, "Write restart data every given period of time", &list);
    IntParameter("OUTSIZEEVERY", 0, "Write step size every given time step", &list);

    DoubleParameter("STEPSIZEMAX", 0.0, "Limit maximally permitted time step size (>0)", &list);
    DoubleParameter("STEPSIZEMIN", 0.0, "Limit minimally allowed time step size (>0)", &list);
    DoubleParameter("SIZERATIOMAX", 0.0, "Limit maximally permitted change of time step size compared to previous size, important for multi-step schemes (>0)", &list);
    DoubleParameter("SIZERATIOMIN", 0.0, "Limit minimally permitted change of time step size compared to previous size, important for multi-step schemes (>0)", &list);
    DoubleParameter("SIZERATIOSCALE", 0.9, "This is a safety factor to scale theoretical optimal step size, should be lower than 1 and must be larger than 0", &list);

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
      &list);

    DoubleParameter("LOCERRTOL", 0.0, "Target local error tolerance (>0)", &list);
    IntParameter("ADAPTSTEPMAX", 0, "Limit maximally allowed step size reduction attempts (>0)", &list);
  }



  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
  {
    using namespace DRT::INPUT;
    using Teuchos::tuple;
    using Teuchos::setStringToIntegralParameter;

    Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
    Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

    Teuchos::ParameterList& sdyn = list->sublist("STRUCTURAL DYNAMIC",false,"");

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
                                 &sdyn);

    setStringToIntegralParameter<int>("DYNAMICTYP","GenAlpha",
                                 "type of the specific dynamic time integration scheme",
                                 tuple<std::string>(
                                   "Statics",
                                   "GenAlpha",
                                   "GenAlphaLieGroup",
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
                                   dyna_genalpha_liegroup,
                                   dyna_onesteptheta,
                                   dyna_gemm,
                                   dyna_expleuler,
                                   dyna_centrdiff,
                                   dyna_ab2,
                                   dyna_euma,
                                   dyna_euimsto,
                                   dyna_statmech),
                                 &sdyn);

    setStringToIntegralParameter<int>("PRESTRESS","none","prestressing takes values none mulf id",
                                 tuple<std::string>("none","None","NONE",
                                                    "mulf","Mulf","MULF",
                                                    "id","Id","ID"),
                                 tuple<int>(prestress_none,prestress_none,prestress_none,
                                                              prestress_mulf,prestress_mulf,prestress_mulf,
                                                              prestress_id,prestress_id,prestress_id),
                                 &sdyn);

    DoubleParameter("PRESTRESSTIME",0.0,"time to switch from pre to post stressing",&sdyn);

    // Output type
    IntParameter("RESULTSEVRY",1,"save displacements and contact forces every RESULTSEVRY steps",&sdyn);
    IntParameter("RESEVRYERGY",0,"write system energies every requested step",&sdyn);
    IntParameter("RESTARTEVRY",1,"write restart possibility every RESTARTEVRY steps",&sdyn);
    // Time loop control
    DoubleParameter("TIMESTEP",0.05,"time step size",&sdyn);
    IntParameter("NUMSTEP",200,"maximum number of steps",&sdyn);
    DoubleParameter("TIMEINIT",0.0,"initial time",&sdyn);
    DoubleParameter("MAXTIME",5.0,"maximum time",&sdyn);
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
                                 &sdyn);
    DoubleParameter("M_DAMP",-1.0,"",&sdyn);
    DoubleParameter("K_DAMP",-1.0,"",&sdyn);

    DoubleParameter("TOLDISP",1.0E-10,
                    "tolerance in the displacement norm for the newton iteration",
                    &sdyn);
    setStringToIntegralParameter<int>("NORM_DISP","Abs","type of norm for displacement convergence check",
                                 tuple<std::string>(
                                   "Abs",
                                   "Rel",
                                   "Mix"),
                                 tuple<int>(
                                   convnorm_abs,
                                   convnorm_rel,
                                   convnorm_mix),
                                 &sdyn);

    DoubleParameter("TOLRES",1.0E-08,
                    "tolerance in the residual norm for the newton iteration",
                    &sdyn);
    setStringToIntegralParameter<int>("NORM_RESF","Abs","type of norm for residual convergence check",
                                 tuple<std::string>(
                                   "Abs",
                                   "Rel",
                                   "Mix"),
                                 tuple<int>(
                                   convnorm_abs,
                                   convnorm_rel,
                                   convnorm_mix),
                                 &sdyn);

    DoubleParameter("TOLPRE",1.0E-08,
                    "tolerance in pressure norm for the newton iteration",
                    &sdyn);
    setStringToIntegralParameter<int>("NORM_PRES","Abs","type of norm for pressure convergence check",
                                 tuple<std::string>(
                                   "Abs"),
                                 tuple<int>(
                                   convnorm_abs),
                                 &sdyn);

    DoubleParameter("TOLINCO",1.0E-08,
                    "tolerance in the incompressible residual norm for the newton iteration",
                    &sdyn);
    setStringToIntegralParameter<int>("NORM_INCO","Abs","type of norm for incompressible residual convergence check",
                                 tuple<std::string>(
                                   "Abs"),
                                 tuple<int>(
                                   convnorm_abs),
                                 &sdyn);

    setStringToIntegralParameter<int>("NORMCOMBI_DISPPRES","And","binary operator to combine pressure and displacement values",
                                 tuple<std::string>(
                                   "And",
                                   "Or"),
                                 tuple<int>(
                                   bop_and,
                                   bop_or),
                                 &sdyn);

    setStringToIntegralParameter<int>("NORMCOMBI_RESFINCO","And","binary operator to combine force and incompressible residual",
                                 tuple<std::string>(
                                   "And",
                                   "Or"),
                                 tuple<int>(
                                   bop_and,
                                   bop_or),
                                 &sdyn);

    setStringToIntegralParameter<int>("NORMCOMBI_RESFDISP","And","binary operator to combine displacement and residual force values",
                                 tuple<std::string>(
                                   "And",
                                   "Or"),
                                 tuple<int>(
                                   bop_and,
                                   bop_or),
                                 &sdyn);

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
        &sdyn);

    IntParameter("STC_LAYER",1,
                 "number of STC layers for multilayer case",
                 &sdyn);

    DoubleParameter("PTCDT",0.1,
                    "pseudo time step for pseudo transient continuation (PTC) stabilized Newton procedure",
                    &sdyn);

    DoubleParameter("TOLCONSTR",1.0E-08,
                    "tolerance in the constr error norm for the newton iteration",
                    &sdyn);

    DoubleParameter("TOLCONSTRINCR",1.0E-08,
                    "tolerance in the constr lm incr norm for the newton iteration",
                    &sdyn);

    IntParameter("MAXITER",50,
                 "maximum number of iterations allowed for Newton-Raphson iteration before failure",
                 &sdyn);
    IntParameter("MINITER",0,
                 "minimum number of iterations to be done within Newton-Raphson loop",
                 &sdyn);
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
                                 &sdyn);

    setStringToIntegralParameter<int>("DIVERCONT","stop","What to do with time integration when Newton-Raphson iteration failed",
                                  tuple<std::string>(
                                    "stop",
                                    "continue",
                                    "repeat_step",
                                    "halve_step",
                                    "adapt_step",
                                    "rand_adapt_step",
                                    "rand_adapt_step_ele_err",
                                    "repeat_simulation",
                                    "adapt_penaltycontact"),
                                  tuple<int>(
                                    divcont_stop,
                                    divcont_continue,
                                    divcont_repeat_step,
                                    divcont_halve_step,
                                    divcont_adapt_step,
                                    divcont_rand_adapt_step,
                                    divcont_rand_adapt_step_ele_err,
                                    divcont_repeat_simulation,
                                    divcont_adapt_penaltycontact),
                                  &sdyn);

    IntParameter("MAXDIVCONREFINEMENTLEVEL",10,
                 "number of times timestep is halved in case nonlinear solver diverges",
                 &sdyn);

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
                                 &sdyn);

    IntParameter("LSMAXITER",30,
                 "maximum number of line search steps",
                 &sdyn);
    DoubleParameter("ALPHA_LS",0.5,
                    "step reduction factor alpha in (Newton) line search scheme",
                    &sdyn);
    DoubleParameter("SIGMA_LS",1.e-4,
                    "sufficient descent factor in (Newton) line search scheme",
                    &sdyn);

    setStringToIntegralParameter<int>("MATERIALTANGENT","analytical","way of evaluating the constitutive matrix",
                                 tuple<std::string>(
                                   "analytical",
                                   "finitedifferences"),
                                 tuple<int>(
                                   0,1),
                                 &sdyn);

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
                                 &sdyn);
    // Currently not used, but structure will be kept if someone wants to reimplement
    // AN 2013_05
    setNumericStringParameter("CONTROLNODE","-1 -1 -1",
                              "for methods other than load control: [node(fortran numbering)] [dof(c-numbering)] [curve(fortran numbering)]",
                              &sdyn);

    setStringToIntegralParameter<int>("LOADLIN","No",
                                      "Use linearization of external follower load in Newton",
                                      yesnotuple,yesnovalue,&sdyn);

    setStringToIntegralParameter<int>("MASSLIN","No","Application of nonlinear inertia terms",
    tuple<std::string>("No","no",
                       "Standard", "standard",
                       "Rotations", "rotations"),

    tuple<int>(ml_none,ml_none,
               ml_standard,ml_standard,
               ml_rotations,ml_rotations),
               &sdyn);


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
                                 &sdyn);

    // Uzawa iteration for constraint systems
    DoubleParameter("UZAWAPARAM",1.0,"Parameter for Uzawa algorithm dealing with lagrange multipliers",&sdyn);
    DoubleParameter("UZAWATOL",1.0E-8,"Tolerance for iterative solve with Uzawa algorithm",&sdyn);
    IntParameter("UZAWAMAXITER",50,"maximum number of iterations allowed for uzawa algorithm before failure going to next newton step",&sdyn);
    setStringToIntegralParameter<int>("UZAWAALGO","direct","",
                                   tuple<std::string>(
                                     "uzawa",
                                     "simple",
                                     "direct"),
                                   tuple<int>(
                                     consolve_uzawa,
                                     consolve_simple,
                                     consolve_direct),
                                   &sdyn);

    // convergence criteria adaptivity
    setStringToIntegralParameter<int>("ADAPTCONV","No",
                                 "Switch on adaptive control of linear solver tolerance for nonlinear solution",
                                 yesnotuple,yesnovalue,&sdyn);
    DoubleParameter("ADAPTCONV_BETTER",0.1,"The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",&sdyn);

    setStringToIntegralParameter<int>("LUMPMASS","No",
                                 "Lump the mass matrix for explicit time integration",
                                 yesnotuple,yesnovalue,&sdyn);

    setStringToIntegralParameter<int>("MODIFIEDEXPLEULER","Yes",
                                 "Use the modified explicit Euler time integration scheme",
                                 yesnotuple,yesnovalue,&sdyn);

    // linear solver id used for structural problems
    IntParameter("LINEAR_SOLVER",-1,"number of linear solver used for structural problems",&sdyn);

    // flag decides if young's modulus is temperature dependent, so far only available
    // for temperature-dependent St.Venant Kirchhoff material
    setStringToIntegralParameter<int>("YOUNG_IS_TEMP_DEPENDENT","No",
                                 "Use temperature-dependent Young's modulus",
                                 yesnotuple,yesnovalue,&sdyn);

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
      &sdyn);

    /*--------------------------------------------------------------------*/
    /* parameters for time step size adaptivity in structural dynamics */
    Teuchos::ParameterList& tap = sdyn.sublist("TIMEADAPTIVITY",false,"");
    SetValidTimeAdaptivityParameters(tap);

    /*----------------------------------------------------------------------*/
    /* parameters for generalised-alpha structural integrator */
    Teuchos::ParameterList& genalpha = sdyn.sublist("GENALPHA",false,"");

    setStringToIntegralParameter<int>("GENAVG","TrLike",
                                 "mid-average type of internal forces",
                                 tuple<std::string>(
                                   "Vague",
                                   "ImrLike",
                                   "TrLike"),
                                 tuple<int>(
                                   midavg_vague,
                                   midavg_imrlike,
                                   midavg_trlike),
                                 &genalpha);
    DoubleParameter("BETA",0.25,"Generalised-alpha factor in (0,1/2]",&genalpha);
    DoubleParameter("GAMMA",0.5,"Generalised-alpha factor in (0,1]",&genalpha);
    DoubleParameter("ALPHA_M",0.5,"Generalised-alpha factor in [0,1)",&genalpha);
    DoubleParameter("ALPHA_F",0.5,"Generalised-alpha factor in [0,1)",&genalpha);
    DoubleParameter("RHO_INF",-1.0,"Generalised-alpha factor in [0,1]",&genalpha);

    /*----------------------------------------------------------------------*/
    /* parameters for one-step-theta structural integrator */
    Teuchos::ParameterList& onesteptheta = sdyn.sublist("ONESTEPTHETA",false,"");

    DoubleParameter("THETA",0.5,"One-step-theta factor in (0,1]",&onesteptheta);


    /*----------------------------------------------------------------------*/
    /* parameters for generalised-energy-momentum structural integrator */
    Teuchos::ParameterList& gemm = sdyn.sublist("GEMM",false,"");

    DoubleParameter("BETA",0.25,"Generalised-alpha factor in (0,0.5]",&gemm);
    DoubleParameter("GAMMA",0.5,"Generalised-alpha factor in (0,1]",&gemm);
    DoubleParameter("ALPHA_M",0.5,"Generalised-alpha factor in [0,1)",&gemm);
    DoubleParameter("ALPHA_F",0.5,"Generalised-alpha factor in [0,1)",&gemm);
    DoubleParameter("XI",0.0,"generalisation factor in [0,1)",&gemm);
  }



  void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
  {
    using namespace DRT::INPUT;

    /*--------------------------------------------------------------------*/

    // structural Robin spring dashpot boundary condition (spring and dashpot in parallel) - mhv 08/16

    Teuchos::RCP<ConditionDefinition> robinspringdashpotcond =
      Teuchos::rcp(new ConditionDefinition("DESIGN SURF ROBIN SPRING DASHPOT CONDITIONS",
                                           "RobinSpringDashpot",
                                           "Robin Spring Dashpot",
                                           DRT::Condition::RobinSpringDashpot,
                                           true,
                                           DRT::Condition::Surface));


    robinspringdashpotcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("NUMDOF")));
    robinspringdashpotcond->AddComponent(Teuchos::rcp(new IntConditionComponent     ("numdof")));

    robinspringdashpotcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("ONOFF")));
    robinspringdashpotcond->AddComponent(Teuchos::rcp(new IntVectorConditionComponent("onoff", 3)));

    robinspringdashpotcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("STIFF")));
    robinspringdashpotcond->AddComponent(Teuchos::rcp(new RealVectorConditionComponent("stiff", 3)));

    robinspringdashpotcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("VISCO")));
    robinspringdashpotcond->AddComponent(Teuchos::rcp(new RealVectorConditionComponent("visco", 3)));

    robinspringdashpotcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("DISPLOFFSET")));
    robinspringdashpotcond->AddComponent(Teuchos::rcp(new RealVectorConditionComponent("disploffset", 3)));

    robinspringdashpotcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("DIRECTION")));
    robinspringdashpotcond->AddComponent(Teuchos::rcp(new StringConditionComponent("direction", "xyz",
                                                                                         Teuchos::tuple<std::string>("xyz","refsurfnormal","cursurfnormal"),
                                                                                         Teuchos::tuple<std::string>("xyz","refsurfnormal","cursurfnormal"),
                                                                                         false)));

    robinspringdashpotcond->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("COUPLING", true)));
    robinspringdashpotcond->AddComponent(Teuchos::rcp(new IntVectorConditionComponent("coupling id", 1, true, true)));

    condlist.push_back(robinspringdashpotcond);


    /*--------------------------------------------------------------------*/
    // surface coupling for spring dashpot DIRECTION cursurfnormal
    // pfaller Apr15

    Teuchos::RCP<ConditionDefinition> springdashpotcoupcond =
          Teuchos::rcp(new ConditionDefinition("DESIGN SURF ROBIN SPRING DASHPOT COUPLING CONDITIONS",
              "RobinSpringDashpotCoupling",
              "RobinSpring Dashpot Coupling",
              DRT::Condition::RobinSpringDashpotCoupling,
              true,
              DRT::Condition::Surface));

      springdashpotcoupcond->AddComponent(Teuchos::rcp(new IntConditionComponent("coupling id")));

      condlist.push_back(springdashpotcoupcond);


    /*--------------------------------------------------------------------*/
    // surfactant

    Teuchos::RCP<ConditionDefinition> surfactant =
      Teuchos::rcp(new ConditionDefinition("SURFACTANT CONDITIONS",
                                           "SurfaceStress",
                                           "Surface Stress (surfactant)",
                                           DRT::Condition::Surfactant,
                                           true,
                                           DRT::Condition::Surface));

    surfactant->AddComponent(Teuchos::rcp(new IntConditionComponent("curve",true,true)));
    AddNamedReal(surfactant,"k1xCbulk");
    AddNamedReal(surfactant,"k2");
    AddNamedReal(surfactant,"m1");
    AddNamedReal(surfactant,"m2");
    AddNamedReal(surfactant,"gamma_0");
    AddNamedReal(surfactant,"gamma_min");

    condlist.push_back(surfactant);
  }

  } // end of namespace STR
} // end of namespace INPAR

