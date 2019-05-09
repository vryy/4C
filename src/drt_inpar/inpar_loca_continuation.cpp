/*----------------------------------------------------------------------*/
/*!
\brief

\maintainer

\level

 */
/*----------------------------------------------------------------------*/

/*
 * inpar_loca_continuation.cpp
 *
 *  Created on: Nov 20, 2015
 *      Author: hiermeier
 */

#include "inpar_loca_continuation.H"
#include "drt_validparameters.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void INPAR::LOCA::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  /*--------------------------------------------------------------------------*
   * parameters for the LOCA continuation method
   *--------------------------------------------------------------------------*/
  Teuchos::ParameterList& loca = list->sublist("LOCA", false, "");
  DRT::INPUT::SetPrintEqualSign(loca, true);

  // sub-list stepper
  Teuchos::ParameterList& stepper = loca.sublist("Stepper", false,
      "Sub-list used by LOCA::Stepper to set parameters relevant for continuation run.");
  DRT::INPUT::SetPrintEqualSign(stepper, true);

  {
    Teuchos::Array<std::string> cm_str = Teuchos::tuple<std::string>("Natural", "Arc Length");
    Teuchos::setStringToIntegralParameter<int>("Continuation Method", "Arc Length",
        "Type of continuation to use.", cm_str, Teuchos::tuple<int>(0, 1), &stepper);

    Teuchos::Array<std::string> cp_str = Teuchos::tuple<std::string>("Unknown");
    Teuchos::setStringToIntegralParameter<int>("Continuation Parameter", "Unknown",
        "Name of continuation parameter.", cp_str, Teuchos::tuple<int>(0), &stepper);

    DRT::INPUT::DoubleParameter("Initial Value", -1.0e12,
        "(must be supplied) Initial value of continuation parameter.", &stepper);
    DRT::INPUT::DoubleParameter("Max Value", -1.0e12,
        "(must be supplied) Maximum value of continuation parameter.", &stepper);
    DRT::INPUT::DoubleParameter("Min Value", -1.0e12,
        "(must be supplied) Minimum value of continuation parameter.", &stepper);
    DRT::INPUT::IntParameter(
        "Max Steps", 100, "Maximum number of continuation steps (including failed step)", &stepper);
    DRT::INPUT::IntParameter("Max Nonlinear Iterations", 15,
        "Maximum number of nonlinear iterations per continuation step", &stepper);
    DRT::INPUT::BoolParameter("Skip Parameter Derivative", "Yes",
        "For natural continuation skip the "
        "df/dp computation which is usually unnecessary.",
        &stepper);
    DRT::INPUT::BoolParameter("Enable Arc Length Scaling", "Yes",
        "Enable arc-length scaling to equilibrate solution "
        "and parameter components to arc-length equations (see "
        "LOCA::MultiContinuation::ArcLengthGroup).",
        &stepper);
    DRT::INPUT::DoubleParameter("Goal Arc Length Parameter Contribution", 0.5,
        "Goal for parameter contribution to arc-length equation.", &stepper);
    DRT::INPUT::DoubleParameter("Max Arc Length Parameter Contribution", 0.8,
        "Max for parameter contribution to arc-length equation, triggering rescaling.", &stepper);
    DRT::INPUT::DoubleParameter("Initial Scale Factor", 1.0,
        "Initial scale factor for parameter term of arc-length equation.", &stepper);
    DRT::INPUT::DoubleParameter("Min Scale Factor", 1.0e-3,
        "Minimum scale factor for scaling parameter term of arc-length equation", &stepper);
    DRT::INPUT::BoolParameter("Enable Tangent Factor Step Size Scaling", "No",
        "Enable step size scaling "
        "by cosine between two consective tangent vector v_0 and v_1 to continuation curve |v_0 * "
        "v_1|^{alpha} where "
        "{alpha} is the tangent factor exponent.",
        &stepper);
    DRT::INPUT::DoubleParameter("Min Tangent Factor", 0.1,
        "Minimum cosine between two consecutive tangent vectors, below which the "
        "continuation step is failed.",
        &stepper);
    DRT::INPUT::DoubleParameter("Tangent Factor Exponent", 1.0,
        "Exponent on the cosine between two consecutive tangent vectors, which "
        "then modifies the step size.",
        &stepper);

    SetValidBorderedSolverMethod(stepper);

    DRT::INPUT::BoolParameter("Compute Eigenvalues", "No",
        "Flag for requesting eigenvalue calculation after each continuation step.", &stepper);
  }

  // sub-sub-list eigensolver
  Teuchos::ParameterList& eigensolver = stepper.sublist("Eigensolver", false,
      "Sub-list used by LOCA::Eigensolver::Factory to determine eigensolver strategies.");
  DRT::INPUT::SetPrintEqualSign(eigensolver, true);
  {
    Teuchos::Array<std::string> method_str = Teuchos::tuple<std::string>("Default", "Anasazi");
    Teuchos::setStringToIntegralParameter<int>(
        "Method", "Default", "", method_str, Teuchos::tuple<int>(0, 1), &eigensolver);

    /* Anasazi parameters */
    Teuchos::Array<std::string> op_str =
        Teuchos::tuple<std::string>("Jacobian Inverse", "Shift-Invert", "Cayley");
    Teuchos::setStringToIntegralParameter<int>("Operator", "Jacobian Inverse",
        "\"Jacobian Inverse\" - Eigenvalues of \f$J^{-1}\f$ | "
        "\"Shift-Invert\" - Eigenvalues of \f$(J-\\sigma M)^{-1}M\f$ | "
        "\"Cayley\" - Eigenvalues of \f$(J-\\sigma M)^{-1}(J - \\mu M)\f$ ",
        op_str, Teuchos::tuple<int>(0, 1, 2), &eigensolver);

    /* Shift-Invert */
    DRT::INPUT::DoubleParameter("Shift", 0.0, "\f$\\sigma\f$ of \"Shift-Invert\".", &eigensolver);
    /* Cayley */
    DRT::INPUT::DoubleParameter("Cayley Pole", 0.0, "\f$\\sigma\f$ of \"Cayley\".", &eigensolver);
    DRT::INPUT::DoubleParameter("Cayley Zero", 0.0, "\f$\\mu\f$ of \"Cayley\".", &eigensolver);

    DRT::INPUT::IntParameter("Block Size", 1, "Block Size", &eigensolver);
    DRT::INPUT::IntParameter("Num Blocks", 30,
        "Maximum number of blocks (equals the maximum length of the Arnoldi factorization.",
        &eigensolver);
    DRT::INPUT::IntParameter(
        "Num Eigenvalues", 4, "Number of requested eigenvalues.", &eigensolver);
    DRT::INPUT::DoubleParameter(
        "Convergence Tolerance", 1.0e-7, "Tolerance for the converged eigenvalues.", &eigensolver);
    DRT::INPUT::IntParameter(
        "Step Size", 1, "Checks convergence every so many steps.", &eigensolver);
    DRT::INPUT::IntParameter("Maximum Restarts", 1, "Number of restarts allowed.", &eigensolver);
    DRT::INPUT::BoolParameter("Symmetric", "No", "Is the operator symmetric?", &eigensolver);
    DRT::INPUT::BoolParameter("Normalize Eigenvectors with Mass Matrix", "No",
        "Option to normalize \f$v M v = 1\f$ instead of \f$v v =1\f$.", &eigensolver);

    Teuchos::Array<std::string> so_str =
        Teuchos::tuple<std::string>("LM", "Largest Magnitude", "LR", "Largest Real Component", "LI",
            "Largest Imaginary Component", "SM", "Smallest Magnitude", "SR",
            "Smallest Real Component", "SI", "Smallest Imaginary Component", "CA", "Cayley");
    Teuchos::setStringToIntegralParameter<int>("Sorting Order", "LM",
        "\"LM\" - Largest magnitude | "
        "\"LR\" - Largest real component | "
        "\"LI\" - Largest Imaginary component | "
        "\"SM\" - Smallest magnitude | "
        "\"SR\" - Smallest real component | "
        "\"SI\" - Smallest imaginary component | "
        "\"CA\" - Largest real part of inverse Cayley transformation",
        so_str, Teuchos::tuple<int>(0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6), &eigensolver);
  }

  // sub-list bifurcation
  Teuchos::ParameterList& bifurcation = loca.sublist("Bifurcation", false,
      "Sub-list used by LOCA::Bifurcation::Factory to determine what type of bifurcation "
      "calculation, if any to use.");
  DRT::INPUT::SetPrintEqualSign(bifurcation, true);

  {
    Teuchos::Array<std::string> type_str =
        Teuchos::tuple<std::string>("None", "Turning Point", "Pitchfork", "Hopf", "User Defined");
    Teuchos::setStringToIntegralParameter<int>("Type", "None", "Bifurcation method to use.",
        type_str, Teuchos::tuple<int>(0, 1, 2, 3, 4), &bifurcation);

    Teuchos::Array<std::string> formulation_str =
        Teuchos::tuple<std::string>("Moore-Spence", "Minimally Augmented");
    Teuchos::setStringToIntegralParameter<int>("Formulation", "Moore-Spence",
        "Name of the bifurcation formulation.", formulation_str, Teuchos::tuple<int>(0, 1),
        &bifurcation);

    Teuchos::Array<std::string> param_str = Teuchos::tuple<std::string>("Unknown");
    Teuchos::setStringToIntegralParameter<int>("Bifurcation Parameter", "Unknown",
        "(Must be supplied if \"Type\" is not \"None\") Name of bifurcation parameter.", param_str,
        Teuchos::tuple<int>(0), &bifurcation);

    DRT::INPUT::BoolParameter("Perturb Initial Solution", "No",
        "Flag indicating whether to apply an initial perturbation to the "
        "initial guess for the solution vector before starting bifurcation algorithm.",
        &bifurcation);
    DRT::INPUT::DoubleParameter("Relative Perturbation Size", 1.0e-3,
        "Size of relative perturbation of initial "
        "guess for solution vector.",
        &bifurcation);

    Teuchos::Array<std::string> bif_sol_method_str =
        Teuchos::tuple<std::string>("Salinger Bordering", "Phipps Bordering");
    Teuchos::setStringToIntegralParameter<int>("Solver Method", "Salinger Bordering",
        "Solver method", bif_sol_method_str, Teuchos::tuple<int>(0, 1), &bifurcation);

    // ToDo add the remaining parameters ...
  }

  // sub-list predictor
  Teuchos::ParameterList& predictor = loca.sublist("Predictor", false,
      "Sub-list used by LOCA::MultiPredictor::Factory to determine what type of predictor to use "
      "for each "
      "continuation step.");
  DRT::INPUT::SetPrintEqualSign(predictor, true);

  INPAR::LOCA::SetValidPredictorParams(predictor);

  // sub-sub-list first step predictor
  Teuchos::ParameterList& fs_predictor = predictor.sublist("First Step Predictor", false,
      "Sub-Sub-list used by the secant predictor to determine which predictor to use for the first "
      "continuation step.");
  DRT::INPUT::SetPrintEqualSign(fs_predictor, true);

  INPAR::LOCA::SetValidPredictorParams(fs_predictor);

  // sub-sub-list last step predictor
  Teuchos::ParameterList& ls_predictor = predictor.sublist("Last Step Predictor", false,
      "Sub-Sub-list used for last step of arc-length continuation to hit target (max or min) value "
      "exactly (usually \"Constant\" or \"Random\").");
  DRT::INPUT::SetPrintEqualSign(ls_predictor, true);

  INPAR::LOCA::SetValidPredictorParams(ls_predictor);

  // sub-list step size
  Teuchos::ParameterList& stepsize = loca.sublist("Step Size", false,
      "Sub-list used by LOCA::StepSize::Factory to determine step size control strategies.");
  DRT::INPUT::SetPrintEqualSign(stepsize, true);

  {
    Teuchos::Array<std::string> method_str = Teuchos::tuple<std::string>("Constant", "Adaptive");
    Teuchos::setStringToIntegralParameter<int>("Method", "Adaptive",
        "Step size control strategy to use. "
        "\"Constant\" - Use a constant step size in general, reducing the step size after a "
        "failure and increasing step size back up to original value after subsequent successes "
        "(see LOCA::StepSize::Constant) | "
        "\"Adaptive\" - Use an adaptive step size control strategy that adjusts step size "
        "according to the number of Newton iterations per step (see LOCA::StepSize::Adaptive)",
        method_str, Teuchos::tuple<int>(0, 1), &stepsize);

    DRT::INPUT::DoubleParameter(
        "Initial Step Size", 1.0, "Initial parameter step size.", &stepsize);
    DRT::INPUT::DoubleParameter(
        "Min Step Size", 1.0e-12, "Minimum parameter step size.", &stepsize);
    DRT::INPUT::DoubleParameter(
        "Max Step Size", 1.0e+12, "Maximum parameter step size.", &stepsize);
    DRT::INPUT::DoubleParameter("Failed Step Reduction Factor", 0.5,
        "Factor by which step size is reduced after a failed step.", &stepsize);
    DRT::INPUT::DoubleParameter("Successful Step Increase Factor", 1.26,
        "Factor by which step size is increased after a successful step when the step size is "
        "smaller than the initial step size (Constant step size method only).",
        &stepsize);
    DRT::INPUT::DoubleParameter("Aggressiveness", 0.5,
        "Aggressiveness factor in adaptive step size adjustment.", &stepsize);
  }

  // sub-list constraints
  Teuchos::ParameterList& constraints = loca.sublist(
      "Constraints", false, "Sub-list used to provide additional constraint equations.");
  DRT::INPUT::SetPrintEqualSign(constraints, true);

  SetValidBorderedSolverMethod(constraints);

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void INPAR::LOCA::SetValidPredictorParams(Teuchos::ParameterList& predictor)
{
  Teuchos::Array<std::string> method_str =
      Teuchos::tuple<std::string>("Constant", "Secant", "Tangent", "Random");
  Teuchos::setStringToIntegralParameter<int>("Method", "Secant",
      "Predictor method to use for computing the initial guess for each continuation step.",
      method_str, Teuchos::tuple<int>(0, 1, 2, 3), &predictor);

  DRT::INPUT::DoubleParameter(
      "Epsilon", 1.0e-3, "Relative size of perturbation for random predictor.", &predictor);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void INPAR::LOCA::SetValidBorderedSolverMethod(Teuchos::ParameterList& list)
{
  Teuchos::Array<std::string> bs_str =
      Teuchos::tuple<std::string>("Bordering", "Nested", "Householder", "Augmented");
  Teuchos::setStringToIntegralParameter<int>("Bordered Solver Method", "Householder",
      "Method for solving bordered system of equations", bs_str, Teuchos::tuple<int>(0, 1, 2, 3),
      &list);
}
