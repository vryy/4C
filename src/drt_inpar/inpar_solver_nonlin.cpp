/*----------------------------------------------------------------------*/
/*!
\file inpar_solver_nonlin.cpp

\maintainer Matthias Mayr, Michael Hiermeier

\brief Input parameters for nonlinear solvers
*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_solver_nonlin.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void INPAR::NLNSOL::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& nlnsol = list->sublist("NONLINEAR SOLVER", false, "Configuration of nonlinear solver package");

  DRT::INPUT::StringParameter("XML_FILE", "none", "Filename of XML file with configuration of nonlinear solver", &nlnsol);

  /*----------------------------------------------------------------------*
   * parameters for NOX - non-linear solution
   *----------------------------------------------------------------------*/
  Teuchos::ParameterList& snox = list->sublist("STRUCT NOX",false,"");
  SetPrintEqualSign(snox, true);

  {
    Teuchos::Array<std::string> st = Teuchos::tuple<std::string>(
        "Line Search Based",
        "Pseudo Transient",
        "Trust Region Based",
        "Inexact Trust Region Based",
        "Tensor Based");
    Teuchos::setStringToIntegralParameter<int>(
        "Nonlinear Solver","Line Search Based","",
        st,Teuchos::tuple<int>( 0, 1, 2, 3, 4 ),
        &snox);
  }

  // sub-list direction
  Teuchos::ParameterList& direction = snox.sublist("Direction",false,"");
  SetPrintEqualSign(direction,true);

  {
    Teuchos::Array<std::string> st = Teuchos::tuple<std::string>(
        "Newton",
        "Steepest Descent",
        "NonlinearCG",
        "Broyden");
    Teuchos::setStringToIntegralParameter<int>(
        "Method","Newton","",
        st,Teuchos::tuple<int>( 0, 1, 2, 3 ),
        &direction);
  }

  // sub-sub-list "Newton"
  Teuchos::ParameterList& newton = direction.sublist("Newton",false,"");
  SetPrintEqualSign(newton,true);

  {
    Teuchos::Array<std::string> forcingtermmethod = Teuchos::tuple<std::string>(
        "Constant",
        "Type 1",
        "Type 2");
    Teuchos::setStringToIntegralParameter<int>(
        "Forcing Term Method","Constant","",
        forcingtermmethod,Teuchos::tuple<int>( 0, 1, 2 ),
        &newton);
    DoubleParameter("Forcing Term Initial Tolerance",0.1,"initial linear solver tolerance",&newton);
    DoubleParameter("Forcing Term Minimum Tolerance",1.0e-6,"",&newton);
    DoubleParameter("Forcing Term Maximum Tolerance",0.01,"",&newton);
    DoubleParameter("Forcing Term Alpha",1.5,"used only by \"Type 2\"",&newton);
    DoubleParameter("Forcing Term Gamma",0.9,"used only by \"Type 2\"",&newton);
    BoolParameter("Rescue Bad Newton Solver","Yes","If set to true, we will use "
        "the computed direction even if the linear solve does not achieve the tolerance "
        "specified by the forcing term",&newton);
  }

  // sub-sub-list "Steepest Descent"
  Teuchos::ParameterList& steepestdescent = direction.sublist("Steepest Descent",false,"");
  SetPrintEqualSign(steepestdescent,true);

  {
    Teuchos::Array<std::string> scalingtype = Teuchos::tuple<std::string>(
        "2-Norm",
        "Quadratic Model Min",
        "F 2-Norm",
        "None");
    Teuchos::setStringToIntegralParameter<int>(
        "Scaling Type","None","",
        scalingtype,Teuchos::tuple<int>( 0, 1, 2, 3 ),
        &steepestdescent);
  }

  // sub-list "Pseudo Transient"
  Teuchos::ParameterList& ptc=snox.sublist("Pseudo Transient",false,"");
  SetPrintEqualSign(ptc,true);

  {
    DoubleParameter("deltaInit",-1.0,"Initial time step size. If its negative, the initial time step is calculated automatically.",&ptc);
    DoubleParameter("deltaMax",std::numeric_limits<double>::max(),"Maximum time step size. "
        "If the new step size is greater than this value, the transient terms will be eliminated "
        "from the Newton iteration resulting in a full Newton solve.",&ptc);
    DoubleParameter("deltaMin",1.0e-5,"Minimum step size.",&ptc);
    IntParameter("Max Number of PTC Iterations",std::numeric_limits<int>::max(),"",&ptc);
    Teuchos::Array<std::string> time_step_control = Teuchos::tuple<std::string>(
        "SER",
        "Switched Evolution Relaxation",
        "TTE",
        "Temporal Truncation Error",
        "MRR",
        "Model Reduction Ratio");
    Teuchos::setStringToIntegralParameter<int>(
        "Time Step Control","SER","",
        time_step_control,Teuchos::tuple<int>( 0, 0, 1, 1, 2, 2 ),
        &ptc);
    Teuchos::Array<std::string> tsc_norm_type = Teuchos::tuple<std::string>(
        "Two Norm",
        "One Norm",
        "Max Norm");
    Teuchos::setStringToIntegralParameter<int>(
        "Norm Type for TSC","Max Norm","Norm Type for the time step control",
        tsc_norm_type,Teuchos::tuple<int>( 0, 1, 2 ),
        &ptc);
    Teuchos::Array<std::string> scaling_op = Teuchos::tuple<std::string>(
        "Identity",
        "CFL Diagonal",
        "Lumped Mass");
    Teuchos::setStringToIntegralParameter<int>(
        "Scaling Type","Identity","Type of the scaling matrix for the PTC method.",
        scaling_op,Teuchos::tuple<int>( 0, 1, 2 ),
        &ptc);
  }

  // sub-list "Line Search"
  Teuchos::ParameterList& linesearch =snox.sublist("Line Search",false,"");
  SetPrintEqualSign(linesearch,true);

  {
    Teuchos::Array<std::string> method = Teuchos::tuple<std::string>(
        "Full Step",
        "Backtrack" ,
        "Polynomial",
        "More'-Thuente",
        "User Defined");
    Teuchos::setStringToIntegralParameter<int>(
        "Method","Full Step","",
        method,Teuchos::tuple<int>( 0, 1, 2, 3, 4 ),
        &linesearch);
  }

  // sub-sub-list "Full Step"
  Teuchos::ParameterList& fullstep = linesearch.sublist("Full Step",false,"");
  SetPrintEqualSign(fullstep,true);

  {
    DoubleParameter("Full Step",1.0,"length of a full step",&fullstep);
  }

  // sub-sub-list "Backtrack"
  Teuchos::ParameterList& backtrack = linesearch.sublist("Backtrack",false,"");
  SetPrintEqualSign(backtrack,true);

  {
    DoubleParameter("Default Step",1.0,"starting step length",&backtrack);
    DoubleParameter("Minimum Step",1.0e-12,"minimum acceptable step length",&backtrack);
    DoubleParameter("Recovery Step",1.0,"step to take when the line search fails (defaults to value for \"Default Step\")",&backtrack);
    IntParameter("Max Iters",50,"maximum number of iterations (i.e., RHS computations)",&backtrack);
    DoubleParameter("Reduction Factor",0.5,"A multiplier between zero and one that reduces the step size between line search iterations",&backtrack);
    BoolParameter("Allow Exceptions","No","Set to true, if exceptions during the force evaluation and backtracking routine should be allowed.",&backtrack);
  }

  // sub-sub-list "Polynomial"
  Teuchos::ParameterList& polynomial = linesearch.sublist("Polynomial",false,"");
  SetPrintEqualSign(polynomial,true);

  {
    DoubleParameter("Default Step",1.0,"Starting step length",&polynomial);
    IntParameter("Max Iters",100,"Maximum number of line search iterations. The search fails if the number of iterations exceeds this value",&polynomial);
    DoubleParameter("Minimum Step",1.0e-12,"Minimum acceptable step length. The search fails if the computed $lambda_k$ is less than this value",&polynomial);
    Teuchos::Array<std::string> recoverysteptype = Teuchos::tuple<std::string>(
        "Constant",
        "Last Computed Step");
    Teuchos::setStringToIntegralParameter<int>(
        "Recovery Step Type","Constant","Determines the step size to take when the line search fails",
        recoverysteptype,Teuchos::tuple<int>( 0, 1 ),
        &polynomial);
    DoubleParameter("Recovery Step",1.0,"The value of the step to take when the line search fails. Only used if the \"Recovery Step Type\" is set to \"Constant\"",&polynomial);
    Teuchos::Array<std::string> interpolationtype = Teuchos::tuple<std::string>(
        "Quadratic",
        "Quadratic3",
        "Cubic");
    Teuchos::setStringToIntegralParameter<int>(
        "Interpolation Type","Cubic","Type of interpolation that should be used",
        interpolationtype,Teuchos::tuple<int>( 0, 1, 2 ),
        &polynomial);
    DoubleParameter("Min Bounds Factor",0.1,"Choice for $gamma_{min}$, i.e., the factor that limits the minimum size of the new step based on the previous step",&polynomial);
    DoubleParameter("Max Bounds Factor",0.5,"Choice for $gamma_{max}$, i.e., the factor that limits the maximum size of the new step based on the previous step",&polynomial);
    Teuchos::Array<std::string> sufficientdecreasecondition = Teuchos::tuple<std::string>(
        "Armijo-Goldstein",
        "Ared/Pred",
        "None");
    Teuchos::setStringToIntegralParameter<int>(
        "Sufficient Decrease Condition","Armijo-Goldstein","Choice to use for the sufficient decrease condition",
        sufficientdecreasecondition,Teuchos::tuple<int>( 0, 1, 2 ),
        &polynomial);
    DoubleParameter("Alpha Factor",1.0e-4,"Parameter choice for sufficient decrease condition",&polynomial);
    BoolParameter("Force Interpolation","No","Set to true if at least one interpolation step should be used. The default is false which means that the line search will stop if the default step length satisfies the convergence criteria",&polynomial);
    BoolParameter("Use Counters","Yes","Set to true if we should use counters and then output the result to the paramter list as described in Output Parameters",&polynomial);
    IntParameter("Maximum Iteration for Increase",0,"Maximum index of the nonlinear iteration for which we allow a relative increase",&polynomial);
    DoubleParameter("Allowed Relative Increase",100,"",&polynomial);
  }

  // sub-sub-list "More'-Thuente"
  Teuchos::ParameterList& morethuente = linesearch.sublist("More'-Thuente",false,"");
  SetPrintEqualSign(morethuente,true);

  {
    DoubleParameter("Sufficient Decrease",1.0e-4,"The ftol in the sufficient decrease condition",&morethuente);
    DoubleParameter("Curvature Condition",0.9999,"The gtol in the curvature condition",&morethuente);
    DoubleParameter("Interval Width",1.0e-15,"The maximum width of the interval containing the minimum of the modified function",&morethuente);
    DoubleParameter("Maximum Step",1.0e6,"maximum allowable step length",&morethuente);
    DoubleParameter("Minimum Step",1.0e-12,"minimum allowable step length",&morethuente);
    IntParameter("Max Iters",20,"maximum number of right-hand-side and corresponding Jacobian evaluations",&morethuente);
    DoubleParameter("Default Step",1.0,"starting step length",&morethuente);
    Teuchos::Array<std::string> recoverysteptype = Teuchos::tuple<std::string>(
        "Constant",
        "Last Computed Step");
    Teuchos::setStringToIntegralParameter<int>(
        "Recovery Step Type","Constant","Determines the step size to take when the line search fails",
        recoverysteptype,Teuchos::tuple<int>( 0, 1 ),
        &morethuente);
    DoubleParameter("Recovery Step",1.0,"The value of the step to take when the line search fails. Only used if the \"Recovery Step Type\" is set to \"Constant\"",&morethuente);
    Teuchos::Array<std::string> sufficientdecreasecondition = Teuchos::tuple<std::string>(
        "Armijo-Goldstein",
        "Ared/Pred",
        "None");
    Teuchos::setStringToIntegralParameter<int>(
        "Sufficient Decrease Condition","Armijo-Goldstein","Choice to use for the sufficient decrease condition",
        sufficientdecreasecondition,Teuchos::tuple<int>( 0, 1, 2 ),
        &morethuente);
    BoolParameter("Optimize Slope Calculation","No","Boolean value. If set to true the value of $s^T J^T F$ is estimated using a directional derivative in a call to NOX::LineSearch::Common::computeSlopeWithOutJac. If false the slope computation is computed with the NOX::LineSearch::Common::computeSlope method. Setting this to true eliminates having to compute the Jacobian at each inner iteration of the More'-Thuente line search",&morethuente);
  }

  // sub-list "Trust Region"
  Teuchos::ParameterList& trustregion = snox.sublist("Trust Region",false,"");
  SetPrintEqualSign(trustregion,true);

  {
    DoubleParameter("Minimum Trust Region Radius",1.0e-6,"Minimum allowable trust region radius",&trustregion);
    DoubleParameter("Maximum Trust Region Radius",1.0e+9,"Maximum allowable trust region radius",&trustregion);
    DoubleParameter("Minimum Improvement Ratio",1.0e-4,"Minimum improvement ratio to accept the step",&trustregion);
    DoubleParameter("Contraction Trigger Ratio",0.1,"If the improvement ratio is less than this value, then the trust region is contracted by the amount specified by the \"Contraction Factor\". Must be larger than \"Minimum Improvement Ratio\"",&trustregion);
    DoubleParameter("Contraction Factor",0.25,"",&trustregion);
    DoubleParameter("Expansion Trigger Ratio",0.75,"If the improvement ratio is greater than this value, then the trust region is contracted by the amount specified by the \"Expansion Factor\"",&trustregion);
    DoubleParameter("Expansion Factor",4.0,"",&trustregion);
    DoubleParameter("Recovery Step",1.0,"",&trustregion);
  }

  // sub-list "Printing"
  Teuchos::ParameterList& printing = snox.sublist("Printing",false,"");
  SetPrintEqualSign(printing,true);

  {
    BoolParameter("Error","No","",&printing);
    BoolParameter("Warning","Yes","",&printing);
    BoolParameter("Outer Iteration","Yes","",&printing);
    BoolParameter("Inner Iteration","Yes","",&printing);
    BoolParameter("Parameters","No","",&printing);
    BoolParameter("Details","No","",&printing);
    BoolParameter("Outer Iteration StatusTest","Yes","",&printing);
    BoolParameter("Linear Solver Details","No","",&printing);
    BoolParameter("Test Details","No","",&printing);
    /*  // for LOCA
    BoolParameter("Stepper Iteration","No","",&printing);
    BoolParameter("Stepper Details","No","",&printing);
    BoolParameter("Stepper Parameters","Yes","",&printing);
     */
    BoolParameter("Debug","No","",&printing);
  }

  // sub-list "Status Test"
  Teuchos::ParameterList& statusTest = snox.sublist("Status Test",false,"");
  SetPrintEqualSign(statusTest,true);

  {
    StringParameter("XML File", "none", "Filename of XML file with configuration of nox status test", &statusTest);
  }

  // sub-list "Solver Options"
  Teuchos::ParameterList& solverOptions = snox.sublist("Solver Options",false,"");
  SetPrintEqualSign(solverOptions,true);

  {
    Teuchos::Array<std::string> meritFct = Teuchos::tuple<std::string>(
        "Sum of Squares",
        "Lagrangian");
    Teuchos::setStringToIntegralParameter<int>(
        "Merit Function","Sum of Squares","",
        meritFct,Teuchos::tuple<int>( 0, 1),
        &solverOptions);

    Teuchos::Array<std::string> scTestType = Teuchos::tuple<std::string>(
        "Complete",
        "Minimal",
        "None");
    Teuchos::setStringToIntegralParameter<int>(
        "Status Test Check Type","Complete","",
        scTestType,Teuchos::tuple<int>( 0, 1, 2),
        &solverOptions);
  }

  // sub-sub-sub-list "Linear Solver"
  Teuchos::ParameterList& linearSolver = newton.sublist("Linear Solver",false,"");
  SetPrintEqualSign(linearSolver,true);

  {
    // convergence criteria adaptivity
    BoolParameter("Adaptive Control","No",
        "Switch on adaptive control of linear solver tolerance for nonlinear solution",
        &linearSolver);
    DoubleParameter("Adaptive Control Objective",0.1,"The linear solver shall be this much better than the current nonlinear residual in the nonlinear convergence limit",&linearSolver);
    BoolParameter("Zero Initial Guess","Yes","Zero out the delta X vector if requested.",&linearSolver);
    BoolParameter("Computing Scaling Manually","No","Allows the manually scaling of your linear system (not supported at the moment).",&linearSolver);
    BoolParameter("Output Solver Details","Yes","Switch the linear solver output on and off.",&linearSolver);
  }

}
