/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for thermo problems

\level 1


*/

/*----------------------------------------------------------------------*/



#include "4C_inpar_thermo.hpp"

#include "4C_io_geometry_type.hpp"
#include "4C_lib_conditiondefinition.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void INPAR::THR::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& tdyn = list->sublist("THERMAL DYNAMIC", false, "");

  setStringToIntegralParameter<int>("DYNAMICTYP", "OneStepTheta",
      "type of time integration control",
      tuple<std::string>("Statics", "OneStepTheta", "GenAlpha", "ExplicitEuler"),
      tuple<int>(dyna_statics, dyna_onesteptheta, dyna_genalpha, dyna_expleuler), &tdyn);

  // output type
  CORE::UTILS::IntParameter("RESULTSEVRY", 1,
      "save temperature and other global quantities every RESULTSEVRY steps", &tdyn);
  CORE::UTILS::IntParameter("RESEVRYERGY", 0, "write system energies every requested step", &tdyn);
  CORE::UTILS::IntParameter(
      "RESTARTEVRY", 1, "write restart possibility every RESTARTEVRY steps", &tdyn);

  setStringToIntegralParameter<int>("INITIALFIELD", "zero_field",
      "Initial Field for thermal problem",
      tuple<std::string>("zero_field", "field_by_function", "field_by_condition"),
      tuple<int>(initfield_zero_field, initfield_field_by_function, initfield_field_by_condition),
      &tdyn);

  CORE::UTILS::IntParameter("INITFUNCNO", -1, "function number for thermal initial field", &tdyn);

  // Time loop control
  CORE::UTILS::DoubleParameter("TIMESTEP", 0.05, "time step size", &tdyn);
  CORE::UTILS::IntParameter("NUMSTEP", 200, "maximum number of steps", &tdyn);
  CORE::UTILS::DoubleParameter("MAXTIME", 5.0, "maximum time", &tdyn);

  // Iterationparameters
  CORE::UTILS::DoubleParameter(
      "TOLTEMP", 1.0E-10, "tolerance in the temperature norm of the Newton iteration", &tdyn);

  setStringToIntegralParameter<int>("NORM_TEMP", "Abs",
      "type of norm for temperature convergence check", tuple<std::string>("Abs", "Rel", "Mix"),
      tuple<int>(convnorm_abs, convnorm_rel, convnorm_mix), &tdyn);

  CORE::UTILS::DoubleParameter(
      "TOLRES", 1.0E-08, "tolerance in the residual norm for the Newton iteration", &tdyn);

  setStringToIntegralParameter<int>("NORM_RESF", "Abs",
      "type of norm for residual convergence check", tuple<std::string>("Abs", "Rel", "Mix"),
      tuple<int>(convnorm_abs, convnorm_rel, convnorm_mix), &tdyn);

  setStringToIntegralParameter<int>("NORMCOMBI_RESFTEMP", "And",
      "binary operator to combine temperature and residual force values",
      tuple<std::string>("And", "Or"), tuple<int>(bop_and, bop_or), &tdyn);

  CORE::UTILS::IntParameter("MAXITER", 50,
      "maximum number of iterations allowed for Newton-Raphson iteration before failure", &tdyn);

  CORE::UTILS::IntParameter(
      "MINITER", 0, "minimum number of iterations to be done within Newton-Raphson loop", &tdyn);

  setStringToIntegralParameter<int>("ITERNORM", "L2", "type of norm to be applied to residuals",
      tuple<std::string>("L1", "L2", "Rms", "Inf"),
      tuple<int>(norm_l1, norm_l2, norm_rms, norm_inf), &tdyn);

  setStringToIntegralParameter<int>("DIVERCONT", "stop",
      "What to do with time integration when Newton-Raphson iteration failed",
      tuple<std::string>("stop", "continue", "halve_step", "repeat_step", "repeat_simulation"),
      tuple<int>(divcont_stop, divcont_continue, divcont_halve_step, divcont_repeat_step,
          divcont_repeat_simulation),
      &tdyn);

  CORE::UTILS::IntParameter("MAXDIVCONREFINEMENTLEVEL", 10,
      "number of times timestep is halved in case nonlinear solver diverges", &tdyn);

  setStringToIntegralParameter<int>("NLNSOL", "fullnewton", "Nonlinear solution technique",
      tuple<std::string>("vague", "fullnewton"), tuple<int>(soltech_vague, soltech_newtonfull),
      &tdyn);

  setStringToIntegralParameter<int>("PREDICT", "ConstTemp",
      "Predictor of iterative solution techniques",
      tuple<std::string>("Vague", "ConstTemp", "ConstTempRate", "TangTemp"),
      tuple<int>(pred_vague, pred_consttemp, pred_consttemprate, pred_tangtemp), &tdyn);

  // convergence criteria solver adaptivity
  CORE::UTILS::BoolParameter("ADAPTCONV", "No",
      "Switch on adaptive control of linear solver tolerance for nonlinear solution", &tdyn);
  CORE::UTILS::DoubleParameter("ADAPTCONV_BETTER", 0.1,
      "The linear solver shall be this much better than the current nonlinear residual in the "
      "nonlinear convergence limit",
      &tdyn);

  CORE::UTILS::BoolParameter(
      "LUMPCAPA", "No", "Lump the capacity matrix for explicit time integration", &tdyn);

  // number of linear solver used for thermal problems
  CORE::UTILS::IntParameter(
      "LINEAR_SOLVER", -1, "number of linear solver used for thermal problems", &tdyn);

  // where the geometry comes from
  setStringToIntegralParameter<int>("GEOMETRY", "full", "How the geometry is specified",
      tuple<std::string>("full", "box", "file"),
      tuple<int>(IO::geometry_full, IO::geometry_box, IO::geometry_file), &tdyn);

  setStringToIntegralParameter<int>("CALCERROR", "No",
      "compute error compared to analytical solution", tuple<std::string>("No", "byfunct"),
      tuple<int>(no_error_calculation, calcerror_byfunct), &tdyn);
  CORE::UTILS::IntParameter("CALCERRORFUNCNO", -1, "Function for Error Calculation", &tdyn);

  /*----------------------------------------------------------------------*/
  /* parameters for generalised-alpha thermal integrator */
  Teuchos::ParameterList& tgenalpha = tdyn.sublist("GENALPHA", false, "");

  setStringToIntegralParameter<int>("GENAVG", "TrLike", "mid-average type of internal forces",
      tuple<std::string>("Vague", "ImrLike", "TrLike"),
      tuple<int>(midavg_vague, midavg_imrlike, midavg_trlike), &tgenalpha);

  // default values correspond to midpoint-rule
  CORE::UTILS::DoubleParameter("GAMMA", 0.5, "Generalised-alpha factor in (0,1]", &tgenalpha);
  CORE::UTILS::DoubleParameter("ALPHA_M", 0.5, "Generalised-alpha factor in [0.5,1)", &tgenalpha);
  CORE::UTILS::DoubleParameter("ALPHA_F", 0.5, "Generalised-alpha factor in [0.5,1)", &tgenalpha);
  CORE::UTILS::DoubleParameter("RHO_INF", -1.0, "Generalised-alpha factor in [0,1]", &tgenalpha);

  /*----------------------------------------------------------------------*/
  /* parameters for one-step-theta thermal integrator */
  Teuchos::ParameterList& tonesteptheta = tdyn.sublist("ONESTEPTHETA", false, "");

  CORE::UTILS::DoubleParameter("THETA", 0.5, "One-step-theta factor in (0,1]", &tonesteptheta);
}



void INPAR::THR::SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist)
{
  using namespace INPUT;

  /*--------------------------------------------------------------------*/
  // Convective heat transfer (Newton's law of heat transfer)

  std::vector<Teuchos::RCP<INPUT::LineComponent>> thermoconvectcomponents;

  // decide here if approximation is sufficient
  // --> Tempn (old temperature T_n)
  // or if the exact solution is needed
  // --> Tempnp (current temperature solution T_n+1) with linearisation
  thermoconvectcomponents.push_back(Teuchos::rcp(new INPUT::SelectionComponent("temperature state",
      "Tempnp", Teuchos::tuple<std::string>("Tempnp", "Tempn"),
      Teuchos::tuple<std::string>("Tempnp", "Tempn"))));
  // heat transfer coefficient h
  thermoconvectcomponents.push_back(Teuchos::rcp(new INPUT::SeparatorComponent("coeff")));
  thermoconvectcomponents.push_back(Teuchos::rcp(new INPUT::RealComponent("coeff")));
  // surrounding (fluid) temperature T_oo
  thermoconvectcomponents.push_back(Teuchos::rcp(new INPUT::SeparatorComponent("surtemp")));
  thermoconvectcomponents.push_back(Teuchos::rcp(new INPUT::RealComponent("surtemp")));
  // time curve to increase the surrounding (fluid) temperature T_oo in time
  thermoconvectcomponents.push_back(Teuchos::rcp(new INPUT::SeparatorComponent("surtempfunct")));
  thermoconvectcomponents.push_back(
      Teuchos::rcp(new INPUT::IntComponent("surtempfunct", {0, true, true})));
  // time curve to increase the complete boundary condition, i.e., the heat flux
  thermoconvectcomponents.push_back(Teuchos::rcp(new INPUT::SeparatorComponent("funct")));
  thermoconvectcomponents.push_back(
      Teuchos::rcp(new INPUT::IntComponent("funct", {0, true, true})));

  Teuchos::RCP<ConditionDefinition> linethermoconvect = Teuchos::rcp(new ConditionDefinition(
      "DESIGN THERMO CONVECTION LINE CONDITIONS", "ThermoConvections", "Line Thermo Convections",
      CORE::Conditions::ThermoConvections, true, CORE::Conditions::geometry_type_line));
  Teuchos::RCP<ConditionDefinition> surfthermoconvect = Teuchos::rcp(new ConditionDefinition(
      "DESIGN THERMO CONVECTION SURF CONDITIONS", "ThermoConvections", "Surface Thermo Convections",
      CORE::Conditions::ThermoConvections, true, CORE::Conditions::geometry_type_surface));

  for (unsigned i = 0; i < thermoconvectcomponents.size(); ++i)
  {
    linethermoconvect->AddComponent(thermoconvectcomponents[i]);
    surfthermoconvect->AddComponent(thermoconvectcomponents[i]);
  }

  condlist.push_back(linethermoconvect);
  condlist.push_back(surfthermoconvect);

  /*--------------------------------------------------------------------*/
  // Robin boundary conditions for heat transfer
  // NOTE: this condition must be
  Teuchos::RCP<ConditionDefinition> thermorobinline = Teuchos::rcp(new ConditionDefinition(
      "DESIGN THERMO ROBIN LINE CONDITIONS", "ThermoRobin", "Thermo Robin boundary condition",
      CORE::Conditions::ThermoRobin, true, CORE::Conditions::geometry_type_line));
  Teuchos::RCP<ConditionDefinition> thermorobinsurf = Teuchos::rcp(new ConditionDefinition(
      "DESIGN THERMO ROBIN SURF CONDITIONS", "ThermoRobin", "Thermo Robin boundary condition",
      CORE::Conditions::ThermoRobin, true, CORE::Conditions::geometry_type_surface));

  std::vector<Teuchos::RCP<INPUT::LineComponent>> thermorobincomponents;

  thermorobincomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("NUMSCAL")));
  thermorobincomponents.emplace_back(Teuchos::rcp(new INPUT::IntComponent("numscal")));
  thermorobincomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("ONOFF")));
  thermorobincomponents.emplace_back(
      Teuchos::rcp(new INPUT::IntVectorComponent("onoff", INPUT::LengthFromInt("numscal"))));

  thermorobincomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("PREFACTOR")));
  thermorobincomponents.emplace_back(Teuchos::rcp(new INPUT::RealComponent("prefactor")));
  thermorobincomponents.emplace_back(Teuchos::rcp(new INPUT::SeparatorComponent("REFVALUE")));
  thermorobincomponents.emplace_back(Teuchos::rcp(new INPUT::RealComponent("refvalue")));

  for (const auto& thermorobincomponent : thermorobincomponents)
  {
    thermorobinline->AddComponent(thermorobincomponent);
    thermorobinsurf->AddComponent(thermorobincomponent);
  }

  condlist.push_back(thermorobinline);
  condlist.push_back(thermorobinsurf);
}

FOUR_C_NAMESPACE_CLOSE
