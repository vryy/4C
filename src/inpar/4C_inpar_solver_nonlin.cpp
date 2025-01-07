// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_solver_nonlin.hpp"

#include "4C_solver_nonlin_nox_enum_lists.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Inpar::NlnSol::set_valid_parameters(Teuchos::ParameterList& list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  /*----------------------------------------------------------------------*
   * parameters for NOX - non-linear solution
   *----------------------------------------------------------------------*/
  Teuchos::ParameterList& snox = list.sublist("STRUCT NOX", false, "");

  {
    std::vector<std::string> nonlinear_solver_valid_input = {"Line Search Based",
        "Pseudo Transient", "Trust Region Based", "Inexact Trust Region Based", "Tensor Based",
        "Single Step"};

    Core::Utils::string_parameter("Nonlinear Solver", "Line Search Based",
        "Choose a nonlinear solver method.", &snox, nonlinear_solver_valid_input);
  }

  // sub-list direction
  Teuchos::ParameterList& direction = snox.sublist("Direction", false, "");

  {
    std::vector<std::string> newton_method_valid_input = {
        "Newton", "Steepest Descent", "NonlinearCG", "Broyden", "User Defined"};
    Core::Utils::string_parameter("Method", "Newton",
        "Choose a direction method for the nonlinear solver.", &direction,
        newton_method_valid_input);

    std::vector<std::string> user_defined_method_valid_input = {"Newton", "Modified Newton"};
    Core::Utils::string_parameter("User Defined Method", "Modified Newton",
        "Choose a user-defined direction method.", &direction, user_defined_method_valid_input);
  }

  // sub-sub-list "Newton"
  Teuchos::ParameterList& newton = direction.sublist("Newton", false, "");

  {
    std::vector<std::string> forcing_term_valid_input = {"Constant", "Type 1", "Type 2"};
    Core::Utils::string_parameter(
        "Forcing Term Method", "Constant", "", &newton, forcing_term_valid_input);

    Core::Utils::double_parameter(
        "Forcing Term Initial Tolerance", 0.1, "initial linear solver tolerance", &newton);
    Core::Utils::double_parameter("Forcing Term Minimum Tolerance", 1.0e-6, "", &newton);
    Core::Utils::double_parameter("Forcing Term Maximum Tolerance", 0.01, "", &newton);
    Core::Utils::double_parameter("Forcing Term Alpha", 1.5, "used only by \"Type 2\"", &newton);
    Core::Utils::double_parameter("Forcing Term Gamma", 0.9, "used only by \"Type 2\"", &newton);
    Core::Utils::bool_parameter("Rescue Bad Newton Solve", "Yes",
        "If set to true, we will use "
        "the computed direction even if the linear solve does not achieve the tolerance "
        "specified by the forcing term",
        &newton);
  }

  // sub-sub-list "Steepest Descent"
  Teuchos::ParameterList& steepestdescent = direction.sublist("Steepest Descent", false, "");

  {
    std::vector<std::string> scaling_type_valid_input = {
        "2-Norm", "Quadratic Model Min", "F 2-Norm", "None"};
    Core::Utils::string_parameter(
        "Scaling Type", "None", "", &steepestdescent, scaling_type_valid_input);
  }

  // sub-sub-sub-list "Modified Newton"
  Teuchos::ParameterList& modnewton = newton.sublist("Modified", false, "");
  {
    Core::Utils::double_parameter("Initial Primal Diagonal Correction", 1.0e-4,
        "Initial correction factor for the diagonal of the primal block.", &modnewton);

    Core::Utils::double_parameter("Minimal Primal Diagonal Correction", 1.0e-20,
        "Minimal correction factor for the diagonal of the primal block.", &modnewton);

    Core::Utils::double_parameter("Maximal Primal Diagonal Correction", 1.0e+40,
        "Maximal correction factor for the diagonal of the primal block.", &modnewton);

    Core::Utils::double_parameter("Primal Reduction Factor", 1.0 / 3.0,
        "Reduction factor for the adaption of the primal diagonal correction.", &modnewton);

    Core::Utils::double_parameter("Primal Accretion Factor", 8.0,
        "Accretion factor for the adaption of the primal diagonal correction.", &modnewton);

    Core::Utils::double_parameter("Primal High Accretion Factor", 100.0,
        "High accretion factor for the adaption of the primal diagonal correction.", &modnewton);

    Core::Utils::bool_parameter("Catch Floating Point Exceptions", "No",
        "Set to true, if"
        "floating point exceptions during the linear solver call should be "
        "caught by the algorithm.",
        &modnewton);
  }

  // sub-list "Pseudo Transient"
  Teuchos::ParameterList& ptc = snox.sublist("Pseudo Transient", false, "");

  {
    Core::Utils::double_parameter("deltaInit", -1.0,
        "Initial time step size. If its negative, the initial time step is calculated "
        "automatically.",
        &ptc);
    Core::Utils::double_parameter("deltaMax", std::numeric_limits<double>::max(),
        "Maximum time step size. "
        "If the new step size is greater than this value, the transient terms will be eliminated "
        "from the Newton iteration resulting in a full Newton solve.",
        &ptc);
    Core::Utils::double_parameter("deltaMin", 1.0e-5, "Minimum step size.", &ptc);
    Core::Utils::int_parameter(
        "Max Number of PTC Iterations", std::numeric_limits<int>::max(), "", &ptc);
    Core::Utils::double_parameter("SER_alpha", 1.0, "Exponent of SET.", &ptc);
    Core::Utils::double_parameter("ScalingFactor", 1.0, "Scaling Factor for ptc matrix.", &ptc);

    std::vector<std::string> time_step_control_valid_input = {"SER",
        "Switched Evolution Relaxation", "TTE", "Temporal Truncation Error", "MRR",
        "Model Reduction Ratio"};
    Core::Utils::string_parameter(
        "Time Step Control", "SER", "", &ptc, time_step_control_valid_input);

    std::vector<std::string> tsc_norm_type_valid_input = {"Two Norm", "One Norm", "Max Norm"};
    Core::Utils::string_parameter("Norm Type for TSC", "Max Norm",
        "Norm Type for the time step control", &ptc, tsc_norm_type_valid_input);

    std::vector<std::string> scaling_op_valid_input = {
        "Identity", "CFL Diagonal", "Lumped Mass", "Element based"};
    Core::Utils::string_parameter("Scaling Type", "Identity",
        "Type of the scaling matrix for the PTC method.", &ptc, scaling_op_valid_input);

    std::vector<std::string> build_scale_op_valid_input = {"every iter", "every timestep"};
    Core::Utils::string_parameter("Build scaling operator", "every timestep",
        "Build scaling operator in every iteration or timestep", &ptc, build_scale_op_valid_input);
  }

  // sub-list "Line Search"
  Teuchos::ParameterList& linesearch = snox.sublist("Line Search", false, "");

  {
    std::vector<std::string> method_valid_input = {
        "Full Step", "Backtrack", "Polynomial", "More'-Thuente", "User Defined"};
    Core::Utils::string_parameter("Method", "Full Step", "", &linesearch, method_valid_input);


    Teuchos::Array<std::string> checktypes =
        Teuchos::tuple<std::string>("Complete", "Minimal", "None");
    Teuchos::setStringToIntegralParameter<::NOX::StatusTest::CheckType>(
        "Inner Status Test Check Type", "Minimal",
        "Specify the check type for the inner status tests.", checktypes,
        Teuchos::tuple<::NOX::StatusTest::CheckType>(
            ::NOX::StatusTest::Complete, ::NOX::StatusTest::Minimal, ::NOX::StatusTest::None),
        &linesearch);
  }

  // sub-sub-list "Full Step"
  Teuchos::ParameterList& fullstep = linesearch.sublist("Full Step", false, "");

  {
    Core::Utils::double_parameter("Full Step", 1.0, "length of a full step", &fullstep);
  }

  // sub-sub-list "Backtrack"
  Teuchos::ParameterList& backtrack = linesearch.sublist("Backtrack", false, "");

  {
    Core::Utils::double_parameter("Default Step", 1.0, "starting step length", &backtrack);
    Core::Utils::double_parameter(
        "Minimum Step", 1.0e-12, "minimum acceptable step length", &backtrack);
    Core::Utils::double_parameter("Recovery Step", 1.0,
        "step to take when the line search fails (defaults to value for \"Default Step\")",
        &backtrack);
    Core::Utils::int_parameter(
        "Max Iters", 50, "maximum number of iterations (i.e., RHS computations)", &backtrack);
    Core::Utils::double_parameter("Reduction Factor", 0.5,
        "A multiplier between zero and one that reduces the step size between line search "
        "iterations",
        &backtrack);
    Core::Utils::bool_parameter("Allow Exceptions", "No",
        "Set to true, if exceptions during the force evaluation and backtracking routine should be "
        "allowed.",
        &backtrack);
  }

  // sub-sub-list "Polynomial"
  Teuchos::ParameterList& polynomial = linesearch.sublist("Polynomial", false, "");

  {
    Core::Utils::double_parameter("Default Step", 1.0, "Starting step length", &polynomial);
    Core::Utils::int_parameter("Max Iters", 100,
        "Maximum number of line search iterations. "
        "The search fails if the number of iterations exceeds this value",
        &polynomial);
    Core::Utils::double_parameter("Minimum Step", 1.0e-12,
        "Minimum acceptable step length. The search fails if the computed \f$\\lambda_k\f$ "
        "is less than this value",
        &polynomial);

    std::vector<std::string> recovery_step_type_valid_input = {"Constant", "Last Computed Step"};
    Core::Utils::string_parameter("Recovery Step Type", "Constant",
        "Determines the step size to take when the line search fails", &polynomial,
        recovery_step_type_valid_input);

    Core::Utils::double_parameter("Recovery Step", 1.0,
        "The value of the step to take when the line search fails. Only used if the \"Recovery "
        "Step Type\" is set to \"Constant\"",
        &polynomial);

    std::vector<std::string> interpolation_type_valid_input = {"Quadratic", "Quadratic3", "Cubic"};
    Core::Utils::string_parameter("Interpolation Type", "Cubic",
        "Type of interpolation that should be used", &polynomial, interpolation_type_valid_input);

    Core::Utils::double_parameter("Min Bounds Factor", 0.1,
        "Choice for \f$\\gamma_{\\min}\f$, i.e., the factor that limits the minimum size "
        "of the new step based on the previous step",
        &polynomial);
    Core::Utils::double_parameter("Max Bounds Factor", 0.5,
        "Choice for \f$\\gamma_{\\max}\f$, i.e., the factor that limits the maximum size "
        "of the new step based on the previous step",
        &polynomial);

    std::vector<std::string> sufficient_decrease_condition_valid_input = {
        "Armijo-Goldstein", "Ared/Pred", "None"};
    Core::Utils::string_parameter("Sufficient Decrease Condition", "Armijo-Goldstein",
        "Choice to use for the sufficient decrease condition", &polynomial,
        sufficient_decrease_condition_valid_input);

    Core::Utils::double_parameter(
        "Alpha Factor", 1.0e-4, "Parameter choice for sufficient decrease condition", &polynomial);
    Core::Utils::bool_parameter("Force Interpolation", "No",
        "Set to true if at least one interpolation step should be used. The default is false which "
        "means that the line search will stop if the default step length satisfies the convergence "
        "criteria",
        &polynomial);
    Core::Utils::bool_parameter("Use Counters", "Yes",
        "Set to true if we should use counters and then output the result to the parameter list as "
        "described in Output Parameters",
        &polynomial);
    Core::Utils::int_parameter("Maximum Iteration for Increase", 0,
        "Maximum index of the nonlinear iteration for which we allow a relative increase",
        &polynomial);
    Core::Utils::double_parameter("Allowed Relative Increase", 100, "", &polynomial);
  }

  // sub-sub-list "More'-Thuente"
  Teuchos::ParameterList& morethuente = linesearch.sublist("More'-Thuente", false, "");

  {
    Core::Utils::double_parameter("Sufficient Decrease", 1.0e-4,
        "The ftol in the sufficient decrease condition", &morethuente);
    Core::Utils::double_parameter(
        "Curvature Condition", 0.9999, "The gtol in the curvature condition", &morethuente);
    Core::Utils::double_parameter("Interval Width", 1.0e-15,
        "The maximum width of the interval containing the minimum of the modified function",
        &morethuente);
    Core::Utils::double_parameter(
        "Maximum Step", 1.0e6, "maximum allowable step length", &morethuente);
    Core::Utils::double_parameter(
        "Minimum Step", 1.0e-12, "minimum allowable step length", &morethuente);
    Core::Utils::int_parameter("Max Iters", 20,
        "maximum number of right-hand-side and corresponding Jacobian evaluations", &morethuente);
    Core::Utils::double_parameter("Default Step", 1.0, "starting step length", &morethuente);

    std::vector<std::string> recovery_step_type_valid_input = {"Constant", "Last Computed Step"};
    Core::Utils::string_parameter("Recovery Step Type", "Constant",
        "Determines the step size to take when the line search fails", &morethuente,
        recovery_step_type_valid_input);

    Core::Utils::double_parameter("Recovery Step", 1.0,
        "The value of the step to take when the line search fails. Only used if the \"Recovery "
        "Step Type\" is set to \"Constant\"",
        &morethuente);

    std::vector<std::string> sufficient_decrease_condition_valid_input = {
        "Armijo-Goldstein", "Ared/Pred", "None"};
    Core::Utils::string_parameter("Sufficient Decrease Condition", "Armijo-Goldstein",
        "Choice to use for the sufficient decrease condition", &morethuente,
        sufficient_decrease_condition_valid_input);

    Core::Utils::bool_parameter("Optimize Slope Calculation", "No",
        "Boolean value. If set to true the value of \f$s^T J^T F\f$ is estimated using a "
        "directional derivative in a call to ::NOX::LineSearch::Common::computeSlopeWithOutJac. "
        "If false the slope computation is computed with the "
        "::NOX::LineSearch::Common::computeSlope method. "
        "Setting this to true eliminates having to compute the Jacobian at each inner iteration of "
        "the More'-Thuente line search",
        &morethuente);
  }

  // sub-list "Trust Region"
  Teuchos::ParameterList& trustregion = snox.sublist("Trust Region", false, "");

  {
    Core::Utils::double_parameter("Minimum Trust Region Radius", 1.0e-6,
        "Minimum allowable trust region radius", &trustregion);
    Core::Utils::double_parameter("Maximum Trust Region Radius", 1.0e+9,
        "Maximum allowable trust region radius", &trustregion);
    Core::Utils::double_parameter("Minimum Improvement Ratio", 1.0e-4,
        "Minimum improvement ratio to accept the step", &trustregion);
    Core::Utils::double_parameter("Contraction Trigger Ratio", 0.1,
        "If the improvement ratio is less than this value, then the trust region is contracted by "
        "the amount specified by the \"Contraction Factor\". Must be larger than \"Minimum "
        "Improvement Ratio\"",
        &trustregion);
    Core::Utils::double_parameter("Contraction Factor", 0.25, "", &trustregion);
    Core::Utils::double_parameter("Expansion Trigger Ratio", 0.75,
        "If the improvement ratio is greater than this value, then the trust region is contracted "
        "by the amount specified by the \"Expansion Factor\"",
        &trustregion);
    Core::Utils::double_parameter("Expansion Factor", 4.0, "", &trustregion);
    Core::Utils::double_parameter("Recovery Step", 1.0, "", &trustregion);
  }

  // sub-list "Printing"
  Teuchos::ParameterList& printing = snox.sublist("Printing", false, "");

  {
    Core::Utils::bool_parameter("Error", "No", "", &printing);
    Core::Utils::bool_parameter("Warning", "Yes", "", &printing);
    Core::Utils::bool_parameter("Outer Iteration", "Yes", "", &printing);
    Core::Utils::bool_parameter("Inner Iteration", "Yes", "", &printing);
    Core::Utils::bool_parameter("Parameters", "No", "", &printing);
    Core::Utils::bool_parameter("Details", "No", "", &printing);
    Core::Utils::bool_parameter("Outer Iteration StatusTest", "Yes", "", &printing);
    Core::Utils::bool_parameter("Linear Solver Details", "No", "", &printing);
    Core::Utils::bool_parameter("Test Details", "No", "", &printing);
    Core::Utils::bool_parameter("Debug", "No", "", &printing);
  }

  // sub-list "Status Test"
  Teuchos::ParameterList& statusTest = snox.sublist("Status Test", false, "");

  {
    Core::Utils::string_parameter("XML File", "none",
        "Filename of XML file with configuration"
        " of nox status test",
        &statusTest);
  }

  // sub-list "Solver Options"
  Teuchos::ParameterList& solverOptions = snox.sublist("Solver Options", false, "");

  {
    Teuchos::Array<std::string> meritFct = Teuchos::tuple<std::string>("Sum of Squares");
    Teuchos::setStringToIntegralParameter<NOX::Nln::MeritFunction::MeritFctName>("Merit Function",
        "Sum of Squares", "", meritFct,
        Teuchos::tuple<NOX::Nln::MeritFunction::MeritFctName>(
            NOX::Nln::MeritFunction::mrtfct_sum_of_squares),
        &solverOptions);

    std::vector<std::string> status_test_check_type_valid_input = {"Complete", "Minimal", "None"};
    Core::Utils::string_parameter("Status Test Check Type", "Complete", "", &solverOptions,
        status_test_check_type_valid_input);
  }

  // sub-sub-sub-list "Linear Solver"
  Teuchos::ParameterList& linearSolver = newton.sublist("Linear Solver", false, "");

  {
    // convergence criteria adaptivity
    Core::Utils::bool_parameter("Adaptive Control", "No",
        "Switch on adaptive control of linear solver tolerance for nonlinear solution",
        &linearSolver);
    Core::Utils::double_parameter("Adaptive Control Objective", 0.1,
        "The linear solver shall be this much better than the current nonlinear residual in the "
        "nonlinear convergence limit",
        &linearSolver);
    Core::Utils::bool_parameter(
        "Zero Initial Guess", "Yes", "Zero out the delta X vector if requested.", &linearSolver);
    Core::Utils::bool_parameter("Computing Scaling Manually", "No",
        "Allows the manually scaling of your linear system (not supported at the moment).",
        &linearSolver);
    Core::Utils::bool_parameter("Output Solver Details", "Yes",
        "Switch the linear solver output on and off.", &linearSolver);
  }
}

FOUR_C_NAMESPACE_CLOSE
