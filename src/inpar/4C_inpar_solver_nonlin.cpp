/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for nonlinear solvers


\level 1
*/

/*----------------------------------------------------------------------*/



#include "4C_inpar_solver_nonlin.hpp"

#include "4C_solver_nonlin_nox_enum_lists.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Inpar::NlnSol::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  /*----------------------------------------------------------------------*
   * parameters for NOX - non-linear solution
   *----------------------------------------------------------------------*/
  Teuchos::ParameterList& snox = list->sublist("STRUCT NOX", false, "");

  {
    Teuchos::Array<std::string> st =
        Teuchos::tuple<std::string>("Line Search Based", "Pseudo Transient", "Trust Region Based",
            "Inexact Trust Region Based", "Tensor Based", "Single Step");
    Teuchos::setStringToIntegralParameter<int>("Nonlinear Solver", "Line Search Based", "", st,
        Teuchos::tuple<int>(0, 1, 2, 3, 4, 5), &snox);
  }

  // sub-list direction
  Teuchos::ParameterList& direction = snox.sublist("Direction", false, "");

  {
    Teuchos::Array<std::string> st = Teuchos::tuple<std::string>(
        "Newton", "Steepest Descent", "NonlinearCG", "Broyden", "User Defined");
    Teuchos::setStringToIntegralParameter<int>("Method", "Newton",
        "Choose a direction method for the nonlinear solver.", st,
        Teuchos::tuple<int>(0, 1, 2, 3, 4), &direction);

    Teuchos::Array<std::string> user_methods =
        Teuchos::tuple<std::string>("Newton", "Modified Newton");
    Teuchos::setStringToIntegralParameter<int>("User Defined Method", "Modified Newton",
        "Choose a user-defined direction method.", user_methods, Teuchos::tuple<int>(0, 1),
        &direction);
  }

  // sub-sub-list "Newton"
  Teuchos::ParameterList& newton = direction.sublist("Newton", false, "");

  {
    Teuchos::Array<std::string> forcingtermmethod =
        Teuchos::tuple<std::string>("Constant", "Type 1", "Type 2");
    Teuchos::setStringToIntegralParameter<int>("Forcing Term Method", "Constant", "",
        forcingtermmethod, Teuchos::tuple<int>(0, 1, 2), &newton);
    Core::UTILS::DoubleParameter(
        "Forcing Term Initial Tolerance", 0.1, "initial linear solver tolerance", &newton);
    Core::UTILS::DoubleParameter("Forcing Term Minimum Tolerance", 1.0e-6, "", &newton);
    Core::UTILS::DoubleParameter("Forcing Term Maximum Tolerance", 0.01, "", &newton);
    Core::UTILS::DoubleParameter("Forcing Term Alpha", 1.5, "used only by \"Type 2\"", &newton);
    Core::UTILS::DoubleParameter("Forcing Term Gamma", 0.9, "used only by \"Type 2\"", &newton);
    Core::UTILS::BoolParameter("Rescue Bad Newton Solve", "Yes",
        "If set to true, we will use "
        "the computed direction even if the linear solve does not achieve the tolerance "
        "specified by the forcing term",
        &newton);
  }

  // sub-sub-list "Steepest Descent"
  Teuchos::ParameterList& steepestdescent = direction.sublist("Steepest Descent", false, "");

  {
    Teuchos::Array<std::string> scalingtype =
        Teuchos::tuple<std::string>("2-Norm", "Quadratic Model Min", "F 2-Norm", "None");

    Teuchos::setStringToIntegralParameter<int>(
        "Scaling Type", "None", "", scalingtype, Teuchos::tuple<int>(0, 1, 2, 3), &steepestdescent);
  }

  // sub-sub-sub-list "Modified Newton"
  Teuchos::ParameterList& modnewton = newton.sublist("Modified", false, "");
  {
    Core::UTILS::DoubleParameter("Initial Primal Diagonal Correction", 1.0e-4,
        "Initial correction factor for the diagonal of the primal block.", &modnewton);

    Core::UTILS::DoubleParameter("Minimal Primal Diagonal Correction", 1.0e-20,
        "Minimal correction factor for the diagonal of the primal block.", &modnewton);

    Core::UTILS::DoubleParameter("Maximal Primal Diagonal Correction", 1.0e+40,
        "Maximal correction factor for the diagonal of the primal block.", &modnewton);

    Core::UTILS::DoubleParameter("Primal Reduction Factor", 1.0 / 3.0,
        "Reduction factor for the adaption of the primal diagonal correction.", &modnewton);

    Core::UTILS::DoubleParameter("Primal Accretion Factor", 8.0,
        "Accretion factor for the adaption of the primal diagonal correction.", &modnewton);

    Core::UTILS::DoubleParameter("Primal High Accretion Factor", 100.0,
        "High accretion factor for the adaption of the primal diagonal correction.", &modnewton);

    Teuchos::Array<std::string> defaultsteptests =
        Teuchos::tuple<std::string>("none", "Volume Change Control");
    Teuchos::setStringToIntegralParameter<NOX::Nln::Direction::DefaultStepTest>("Default Step Test",
        "none", "Choose a proper default step test strategy.", defaultsteptests,
        Teuchos::tuple<NOX::Nln::Direction::DefaultStepTest>(
            NOX::Nln::Direction::DefaultStepTest::none,
            NOX::Nln::Direction::DefaultStepTest::volume_change_control),
        &modnewton);

    Core::UTILS::BoolParameter("Catch Floating Point Exceptions", "No",
        "Set to true, if"
        "floating point exceptions during the linear solver call should be "
        "caught by the algorithm.",
        &modnewton);
  }

  // sub-list "Pseudo Transient"
  Teuchos::ParameterList& ptc = snox.sublist("Pseudo Transient", false, "");

  {
    Core::UTILS::DoubleParameter("deltaInit", -1.0,
        "Initial time step size. If its negative, the initial time step is calculated "
        "automatically.",
        &ptc);
    Core::UTILS::DoubleParameter("deltaMax", std::numeric_limits<double>::max(),
        "Maximum time step size. "
        "If the new step size is greater than this value, the transient terms will be eliminated "
        "from the Newton iteration resulting in a full Newton solve.",
        &ptc);
    Core::UTILS::DoubleParameter("deltaMin", 1.0e-5, "Minimum step size.", &ptc);
    Core::UTILS::IntParameter(
        "Max Number of PTC Iterations", std::numeric_limits<int>::max(), "", &ptc);
    Core::UTILS::DoubleParameter("SER_alpha", 1.0, "Exponent of SER.", &ptc);
    Core::UTILS::DoubleParameter("ScalingFactor", 1.0, "Scaling Factor for ptc matrix.", &ptc);
    Teuchos::Array<std::string> time_step_control =
        Teuchos::tuple<std::string>("SER", "Switched Evolution Relaxation", "TTE",
            "Temporal Truncation Error", "MRR", "Model Reduction Ratio");
    Teuchos::setStringToIntegralParameter<int>("Time Step Control", "SER", "", time_step_control,
        Teuchos::tuple<int>(0, 0, 1, 1, 2, 2), &ptc);
    Teuchos::Array<std::string> tsc_norm_type =
        Teuchos::tuple<std::string>("Two Norm", "One Norm", "Max Norm");
    Teuchos::setStringToIntegralParameter<int>("Norm Type for TSC", "Max Norm",
        "Norm Type for the time step control", tsc_norm_type, Teuchos::tuple<int>(0, 1, 2), &ptc);
    Teuchos::Array<std::string> scaling_op =
        Teuchos::tuple<std::string>("Identity", "CFL Diagonal", "Lumped Mass", "Element based");
    Teuchos::setStringToIntegralParameter<int>("Scaling Type", "Identity",
        "Type of the scaling matrix for the PTC method.", scaling_op,
        Teuchos::tuple<int>(0, 1, 2, 3), &ptc);
    // Build scaling Operator every iteration or every timestep
    Teuchos::Array<std::string> build_scale_op =
        Teuchos::tuple<std::string>("every iter", "every timestep");
    Teuchos::setStringToIntegralParameter<int>("Build scaling operator", "every timestep",
        "Build scaling operator in every iteration or timestep", build_scale_op,
        Teuchos::tuple<int>(0, 1), &ptc);
  }

  // sub-list "Line Search"
  Teuchos::ParameterList& linesearch = snox.sublist("Line Search", false, "");

  {
    Teuchos::Array<std::string> method = Teuchos::tuple<std::string>(
        "Full Step", "Backtrack", "Polynomial", "More'-Thuente", "User Defined");
    Teuchos::setStringToIntegralParameter<int>(
        "Method", "Full Step", "", method, Teuchos::tuple<int>(0, 1, 2, 3, 4), &linesearch);

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
    Core::UTILS::DoubleParameter("Full Step", 1.0, "length of a full step", &fullstep);
  }

  // sub-sub-list "Backtrack"
  Teuchos::ParameterList& backtrack = linesearch.sublist("Backtrack", false, "");

  {
    Core::UTILS::DoubleParameter("Default Step", 1.0, "starting step length", &backtrack);
    Core::UTILS::DoubleParameter(
        "Minimum Step", 1.0e-12, "minimum acceptable step length", &backtrack);
    Core::UTILS::DoubleParameter("Recovery Step", 1.0,
        "step to take when the line search fails (defaults to value for \"Default Step\")",
        &backtrack);
    Core::UTILS::IntParameter(
        "Max Iters", 50, "maximum number of iterations (i.e., RHS computations)", &backtrack);
    Core::UTILS::DoubleParameter("Reduction Factor", 0.5,
        "A multiplier between zero and one that reduces the step size between line search "
        "iterations",
        &backtrack);
    Core::UTILS::BoolParameter("Allow Exceptions", "No",
        "Set to true, if exceptions during the force evaluation and backtracking routine should be "
        "allowed.",
        &backtrack);
  }

  // sub-sub-list "Polynomial"
  Teuchos::ParameterList& polynomial = linesearch.sublist("Polynomial", false, "");

  {
    Core::UTILS::DoubleParameter("Default Step", 1.0, "Starting step length", &polynomial);
    Core::UTILS::IntParameter("Max Iters", 100,
        "Maximum number of line search iterations. "
        "The search fails if the number of iterations exceeds this value",
        &polynomial);
    Core::UTILS::DoubleParameter("Minimum Step", 1.0e-12,
        "Minimum acceptable step length. The search fails if the computed \f$\\lambda_k\f$ "
        "is less than this value",
        &polynomial);
    Teuchos::Array<std::string> recoverysteptype =
        Teuchos::tuple<std::string>("Constant", "Last Computed Step");
    Teuchos::setStringToIntegralParameter<int>("Recovery Step Type", "Constant",
        "Determines the step size to take when the line search fails", recoverysteptype,
        Teuchos::tuple<int>(0, 1), &polynomial);
    Core::UTILS::DoubleParameter("Recovery Step", 1.0,
        "The value of the step to take when the line search fails. Only used if the \"Recovery "
        "Step Type\" is set to \"Constant\"",
        &polynomial);
    Teuchos::Array<std::string> interpolationtype =
        Teuchos::tuple<std::string>("Quadratic", "Quadratic3", "Cubic");
    Teuchos::setStringToIntegralParameter<int>("Interpolation Type", "Cubic",
        "Type of interpolation that should be used", interpolationtype,
        Teuchos::tuple<int>(0, 1, 2), &polynomial);
    Core::UTILS::DoubleParameter("Min Bounds Factor", 0.1,
        "Choice for \f$\\gamma_{\\min}\f$, i.e., the factor that limits the minimum size "
        "of the new step based on the previous step",
        &polynomial);
    Core::UTILS::DoubleParameter("Max Bounds Factor", 0.5,
        "Choice for \f$\\gamma_{\\max}\f$, i.e., the factor that limits the maximum size "
        "of the new step based on the previous step",
        &polynomial);
    Teuchos::Array<std::string> sufficientdecreasecondition =
        Teuchos::tuple<std::string>("Armijo-Goldstein", "Ared/Pred", "None");
    Teuchos::setStringToIntegralParameter<int>("Sufficient Decrease Condition", "Armijo-Goldstein",
        "Choice to use for the sufficient decrease condition", sufficientdecreasecondition,
        Teuchos::tuple<int>(0, 1, 2), &polynomial);
    Core::UTILS::DoubleParameter(
        "Alpha Factor", 1.0e-4, "Parameter choice for sufficient decrease condition", &polynomial);
    Core::UTILS::BoolParameter("Force Interpolation", "No",
        "Set to true if at least one interpolation step should be used. The default is false which "
        "means that the line search will stop if the default step length satisfies the convergence "
        "criteria",
        &polynomial);
    Core::UTILS::BoolParameter("Use Counters", "Yes",
        "Set to true if we should use counters and then output the result to the paramter list as "
        "described in Output Parameters",
        &polynomial);
    Core::UTILS::IntParameter("Maximum Iteration for Increase", 0,
        "Maximum index of the nonlinear iteration for which we allow a relative increase",
        &polynomial);
    Core::UTILS::DoubleParameter("Allowed Relative Increase", 100, "", &polynomial);
  }

  // sub-sub-list "More'-Thuente"
  Teuchos::ParameterList& morethuente = linesearch.sublist("More'-Thuente", false, "");

  {
    Core::UTILS::DoubleParameter("Sufficient Decrease", 1.0e-4,
        "The ftol in the sufficient decrease condition", &morethuente);
    Core::UTILS::DoubleParameter(
        "Curvature Condition", 0.9999, "The gtol in the curvature condition", &morethuente);
    Core::UTILS::DoubleParameter("Interval Width", 1.0e-15,
        "The maximum width of the interval containing the minimum of the modified function",
        &morethuente);
    Core::UTILS::DoubleParameter(
        "Maximum Step", 1.0e6, "maximum allowable step length", &morethuente);
    Core::UTILS::DoubleParameter(
        "Minimum Step", 1.0e-12, "minimum allowable step length", &morethuente);
    Core::UTILS::IntParameter("Max Iters", 20,
        "maximum number of right-hand-side and corresponding Jacobian evaluations", &morethuente);
    Core::UTILS::DoubleParameter("Default Step", 1.0, "starting step length", &morethuente);
    Teuchos::Array<std::string> recoverysteptype =
        Teuchos::tuple<std::string>("Constant", "Last Computed Step");
    Teuchos::setStringToIntegralParameter<int>("Recovery Step Type", "Constant",
        "Determines the step size to take when the line search fails", recoverysteptype,
        Teuchos::tuple<int>(0, 1), &morethuente);
    Core::UTILS::DoubleParameter("Recovery Step", 1.0,
        "The value of the step to take when the line search fails. Only used if the \"Recovery "
        "Step Type\" is set to \"Constant\"",
        &morethuente);
    Teuchos::Array<std::string> sufficientdecreasecondition =
        Teuchos::tuple<std::string>("Armijo-Goldstein", "Ared/Pred", "None");
    Teuchos::setStringToIntegralParameter<int>("Sufficient Decrease Condition", "Armijo-Goldstein",
        "Choice to use for the sufficient decrease condition", sufficientdecreasecondition,
        Teuchos::tuple<int>(0, 1, 2), &morethuente);
    Core::UTILS::BoolParameter("Optimize Slope Calculation", "No",
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
    Core::UTILS::DoubleParameter("Minimum Trust Region Radius", 1.0e-6,
        "Minimum allowable trust region radius", &trustregion);
    Core::UTILS::DoubleParameter("Maximum Trust Region Radius", 1.0e+9,
        "Maximum allowable trust region radius", &trustregion);
    Core::UTILS::DoubleParameter("Minimum Improvement Ratio", 1.0e-4,
        "Minimum improvement ratio to accept the step", &trustregion);
    Core::UTILS::DoubleParameter("Contraction Trigger Ratio", 0.1,
        "If the improvement ratio is less than this value, then the trust region is contracted by "
        "the amount specified by the \"Contraction Factor\". Must be larger than \"Minimum "
        "Improvement Ratio\"",
        &trustregion);
    Core::UTILS::DoubleParameter("Contraction Factor", 0.25, "", &trustregion);
    Core::UTILS::DoubleParameter("Expansion Trigger Ratio", 0.75,
        "If the improvement ratio is greater than this value, then the trust region is contracted "
        "by the amount specified by the \"Expansion Factor\"",
        &trustregion);
    Core::UTILS::DoubleParameter("Expansion Factor", 4.0, "", &trustregion);
    Core::UTILS::DoubleParameter("Recovery Step", 1.0, "", &trustregion);
  }

  // sub-list "Printing"
  Teuchos::ParameterList& printing = snox.sublist("Printing", false, "");

  {
    Core::UTILS::BoolParameter("Error", "No", "", &printing);
    Core::UTILS::BoolParameter("Warning", "Yes", "", &printing);
    Core::UTILS::BoolParameter("Outer Iteration", "Yes", "", &printing);
    Core::UTILS::BoolParameter("Inner Iteration", "Yes", "", &printing);
    Core::UTILS::BoolParameter("Parameters", "No", "", &printing);
    Core::UTILS::BoolParameter("Details", "No", "", &printing);
    Core::UTILS::BoolParameter("Outer Iteration StatusTest", "Yes", "", &printing);
    Core::UTILS::BoolParameter("Linear Solver Details", "No", "", &printing);
    Core::UTILS::BoolParameter("Test Details", "No", "", &printing);
    Core::UTILS::BoolParameter("Debug", "No", "", &printing);
  }

  // sub-list "Status Test"
  Teuchos::ParameterList& statusTest = snox.sublist("Status Test", false, "");

  {
    Core::UTILS::StringParameter("XML File", "none",
        "Filename of XML file with configuration"
        " of nox status test",
        &statusTest);
  }

  // sub-list "Solver Options"
  Teuchos::ParameterList& solverOptions = snox.sublist("Solver Options", false, "");

  {
    Teuchos::Array<std::string> meritFct =
        Teuchos::tuple<std::string>("Sum of Squares", "Lagrangian", "Lagrangian Active");
    Teuchos::setStringToIntegralParameter<int>("Merit Function", "Sum of Squares", "", meritFct,
        Teuchos::tuple<int>(NOX::Nln::MeritFunction::mrtfct_sum_of_squares,
            NOX::Nln::MeritFunction::mrtfct_lagrangian,
            NOX::Nln::MeritFunction::mrtfct_lagrangian_active),
        &solverOptions);

    Teuchos::Array<std::string> scTestType =
        Teuchos::tuple<std::string>("Complete", "Minimal", "None");
    Teuchos::setStringToIntegralParameter<int>("Status Test Check Type", "Complete", "", scTestType,
        Teuchos::tuple<int>(0, 1, 2), &solverOptions);
  }

  // sub-sub-sub-list "Linear Solver"
  Teuchos::ParameterList& linearSolver = newton.sublist("Linear Solver", false, "");

  {
    // convergence criteria adaptivity
    Core::UTILS::BoolParameter("Adaptive Control", "No",
        "Switch on adaptive control of linear solver tolerance for nonlinear solution",
        &linearSolver);
    Core::UTILS::DoubleParameter("Adaptive Control Objective", 0.1,
        "The linear solver shall be this much better than the current nonlinear residual in the "
        "nonlinear convergence limit",
        &linearSolver);
    Core::UTILS::BoolParameter(
        "Zero Initial Guess", "Yes", "Zero out the delta X vector if requested.", &linearSolver);
    Core::UTILS::BoolParameter("Computing Scaling Manually", "No",
        "Allows the manually scaling of your linear system (not supported at the moment).",
        &linearSolver);
    Core::UTILS::BoolParameter("Output Solver Details", "Yes",
        "Switch the linear solver output on and off.", &linearSolver);
  }
}

FOUR_C_NAMESPACE_CLOSE
