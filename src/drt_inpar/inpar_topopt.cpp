/*----------------------------------------------------------------------*/
/*!
\level 3

\brief Input parameters for topopt

\maintainer Martin Winklmaier

*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_topopt.H"



void INPAR::TOPOPT::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& topoptcontrol = list->sublist("TOPOLOGY OPTIMIZATION CONTROL", false,
      "control parameters for topology optimization problems");

  DoubleParameter("MAXTIME", 10.0, "Total simulation time", &topoptcontrol);
  IntParameter("NUMSTEP", 100, "Total number of timesteps", &topoptcontrol);
  DoubleParameter("TIMESTEP", 0.1, "Time increment dt", &topoptcontrol);
  IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &topoptcontrol);
  IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &topoptcontrol);

  setStringToIntegralParameter<int>("DENS_TYPE", "node_based",
      "type of optimization = density = porosity field",
      tuple<std::string>("node_based", "element_based"),
      tuple<int>(dens_node_based, dens_ele_based), &topoptcontrol);

  setStringToIntegralParameter<int>("GRADIENT_TYPE", "adjoints", "basic type of adjoint equations",
      tuple<std::string>("adjoints", "FD1", "FD2"),
      tuple<int>(gradientByAdjoints, gradientByFD1, gradientByFD2), &topoptcontrol);

  setStringToIntegralParameter<int>("RESTART_ACTION", "Finished_Optimization_Step",
      "Startint field of Restart",
      tuple<std::string>("Fluid_Time_Step", "Adjoint_Time_Step", "Evaluated_Gradient",
          "Finished_Optimization_Step"),
      tuple<int>(fluid, adjoint, gradient, opti_step), &topoptcontrol);

  setStringToIntegralParameter<int>("CONV_CHECK_TYPE", "Residuum",
      "Convergence check due to given fields",
      tuple<std::string>("Increment", "Residuum", "Increment_and_Residuum"),
      tuple<int>(inc, res, inc_and_res), &topoptcontrol);

  DoubleParameter(
      "RESTOL", 1e-5, "Convergence tolerance of the objective function", &topoptcontrol);
  DoubleParameter(
      "INCTOL", 1e-5, "Convergence tolerance of the optimized variable", &topoptcontrol);

  setStringToIntegralParameter<int>("OBJECTIVE_DISSIPATION", "No",
      "Cdissipation part of the objective function", tuple<std::string>("No", "Yes", "Physical"),
      tuple<int>(obj_diss_no, obj_diss_yes, obj_diss_physical), &topoptcontrol);
  DoubleParameter(
      "DISSIPATION_FAC", -1.0, "factor for the dissipation part of the objective", &topoptcontrol);

  BoolParameter("OBJECTIVE_PRESSURE_DROP", "No", "pressure drop part of the objective function",
      &topoptcontrol);
  DoubleParameter("PRESSURE_DROP_FAC", -1.0,
      "factor for the mean pressure drop part of the objective", &topoptcontrol);

  IntParameter("NUM_OUTPUT_STEPS", 1, "number of output steps saved, 0 = all", &topoptcontrol);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& topoptoptimizer = topoptcontrol.sublist("TOPOLOGY OPTIMIZER", false,
      "control parameters for the optimizer of a topology optimization problem");

  setStringToIntegralParameter<int>("GCMMA_SOLVER", "original",
      "type of solver for gcmma subproblem", tuple<std::string>("original", "gaussian"),
      tuple<int>(INPAR::TOPOPT::gcmma_solver_orig, INPAR::TOPOPT::gcmma_solver_gauss),
      &topoptcontrol);

  DoubleParameter(
      "THETA", 0.5, "theta for temporal integration of objective function", &topoptoptimizer);
  IntParameter("MAX_ITER", 100, "Maximal number of optimization steps", &topoptoptimizer);
  IntParameter("MAX_GRAD_ITER", 100, "Maximal number of optimization steps containing the gradient",
      &topoptoptimizer);
  IntParameter(
      "MAX_INNER_ITER", 20, "Maximal number of inner optimization steps", &topoptoptimizer);
  IntParameter("MAX_SUB_ITER", 200, "Maximal iteration number within subproblem", &topoptoptimizer);
  IntParameter("MAX_INNER_SUB_ITER", 50, "Maximal iteration number within inner subproblem routine",
      &topoptoptimizer);
  IntParameter("MATID", -1, "Material ID for automatic mesh generation", &topoptoptimizer);

  setStringToIntegralParameter<int>("INITIALFIELD", "zero_field",
      "Initial Field for density = topology optimization's optimization variable",
      tuple<std::string>("zero_field", "field_by_function", "channel0", "channel05", "channel1",
          "channelstep0", "channelstep05", "channelstep1"),
      tuple<int>(initdensfield_zero_field, initdensfield_field_by_function,
          initdensfield_channelflow0, initdensfield_channelflow05, initdensfield_channelflow1,
          initdensfield_channelstepflow0, initdensfield_channelstepflow05,
          initdensfield_channelstepflow1),
      &topoptoptimizer);

  setStringToIntegralParameter<int>("TESTCASE", "test_no", "test case for optimizer",
      tuple<std::string>("test_no", "test_snake_one_constr", "test_snake_multiple_constr",
          "test_workflow_without_fluiddata", "test_channel", "test_channel_with_step",
          "test_cornerflow", "test_lin_poro", "test_quad_poro", "test_cubic_poro"),
      tuple<int>(optitest_no, optitest_snake_one_constr, optitest_snake_multiple_constr,
          optitest_workflow_without_fluiddata, optitest_channel, optitest_channel_with_step,
          optitest_cornerflow, optitest_lin_poro, optitest_quad_poro, optitest_cub_poro),
      &topoptoptimizer);

  IntParameter("INITFUNCNO", -1,
      "function number for initial density field in topology optimization", &topoptoptimizer);
  DoubleParameter("VOLUME_BOUNDARY", 0.7, "maximal percentage of fluid volume in background domain",
      &topoptoptimizer);
  DoubleParameter("TOL_KKT", 1.0e-5, "tolerance of optimization problem (for KKT-conditions)",
      &topoptoptimizer);
  DoubleParameter("TOL_SUB", 1.0e-3, "tolerance of subproblem", &topoptoptimizer);
  BoolParameter("TOL_SUB_ADAPTIV", "yes", "adaptive tolerance of subproblem", &topoptoptimizer);
  DoubleParameter("TOL_SUB_QUOT_FAC", 1.0e-03,
      "quotient factor for adaption of subproblem tolerance", &topoptoptimizer);
  DoubleParameter("TOL_SUB_MIN", 1.0e-10, "minimal tolerance of subproblem", &topoptoptimizer);
  BoolParameter("update_smooth", "no", "update smoothing parameter of impermeability function",
      &topoptoptimizer);
  IntParameter(
      "update_smooth_every_iter", 50, "update smoothing parameter every numiter", &topoptoptimizer);
  DoubleParameter("update_smooth_fac", 10.0, "update smoothing parameter factor", &topoptoptimizer);

  // output parameter
  BoolParameter("GMSH_OUTPUT", "No", "Write Gmsh files", &topoptoptimizer);
  IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &topoptoptimizer);

  // parameter of optimization algorithm GCMMA
  DoubleParameter("X_DIFF_MIN", 1.0e-5,
      "minimal difference of upper and lower boundary of optimization variable", &topoptoptimizer);
  DoubleParameter("RHOMIN", 1.0e-6, "minimal parameter value", &topoptoptimizer);
  DoubleParameter("FACMIN", 1.0e-10, "minimal parameter value", &topoptoptimizer);
  DoubleParameter("a_init", 0.0, "initial value for solver parameter", &topoptoptimizer);
  DoubleParameter("c_init", 1000.0, "initial value for solver parameter", &topoptoptimizer);
  DoubleParameter("fac_stepsize", -1.01,
      "factor for adjusting step size in every optimization step", &topoptoptimizer);
  DoubleParameter("RHO_FAC1", 1.1, "factor for updating rho", &topoptoptimizer);
  DoubleParameter("RHO_FAC2", 10.0, "factor for updating rho", &topoptoptimizer);
  DoubleParameter("asymptotes_fac1", 10.0, "factor for updating asymptotes", &topoptoptimizer);
  DoubleParameter("asymptotes_fac2", 0.01, "factor for updating asymptotes", &topoptoptimizer);
  DoubleParameter("fac_x_boundaries", 0.1,
      "unsensible factor for computation of boundaries for optimization variable",
      &topoptoptimizer);
  DoubleParameter("fac_sub_reg", 0.001, "regularisation factor in subproblem", &topoptoptimizer);
  DoubleParameter("GAMMA_UP", 2.3, "controls maximum stepsize increase", &topoptoptimizer);
  DoubleParameter("GAMMA_DOWN", 0.7, "controls maximum stepsize decrease", &topoptoptimizer);
  DoubleParameter("inc2_tol", 1.0e-10,
      "small tolerance for asymptotes update checked against delta x^2", &topoptoptimizer);
  DoubleParameter("GAMMA_UP", 2.3, "controls maximum stepsize increase", &topoptoptimizer);
  DoubleParameter("GAMMA_DOWN", 0.7, "controls maximum stepsize decrease", &topoptoptimizer);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& topoptadjointfluiddyn = topoptcontrol.sublist("TOPOLOGY ADJOINT FLUID",
      false, "control parameters for the adjoint fluid of a topology optimization problem");

  setStringToIntegralParameter<int>("ADJOINT_TYPE", "discrete_adjoint",
      "basic type of adjoint equations", tuple<std::string>("discrete_adjoint", "cont_adjoint"),
      tuple<int>(discrete_adjoint, cont_adjoint), &topoptadjointfluiddyn);

  setStringToIntegralParameter<int>("INITIALFIELD", "zero_field",
      "Initial field for adjoint problem", tuple<std::string>("zero_field", "field_by_function"),
      tuple<int>(initadjointfield_zero_field, initadjointfield_field_by_function),
      &topoptadjointfluiddyn);

  setStringToIntegralParameter<int>("TESTCASE", "test_no", "test case for adjoint problem",
      tuple<std::string>("test_no", "test_primal", "test_stat_const_vel_lin_pres",
          "test_stat_lin_vel_quad_pres", "test_stat_quad_vel_lin_pres",
          "test_stat_all_terms_all_constants", "test_instat_varying_theta",
          "test_instat_all_terms_all_constants", "test_instat_primal_and_dual"),
      tuple<int>(adjointtest_no, adjointtest_primal, adjointtest_stat_const_vel_lin_pres,
          adjointtest_stat_lin_vel_quad_pres, adjointtest_stat_quad_vel_lin_pres,
          adjointtest_stat_all_terms_all_constants, adjointtest_instat_varying_theta,
          adjointtest_instat_all_terms_all_constants, adjointtest_instat_primal_and_dual),
      &topoptadjointfluiddyn);

  IntParameter("INITFUNCNO", -1, "Function for initial field", &topoptadjointfluiddyn);
}
