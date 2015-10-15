/*-----------------------------------------------------------*/
/*!
\file nox_nln_solver_ptc.cpp

\maintainer Michael Hiermeier

\date Oct 6, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_solver_ptc.H"     // class definition

#include "nox_nln_linesearch_factory.H"
#include "nox_nln_direction_factory.H"
#include "nox_nln_aux.H"

#include "nox_nln_linearsystem.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_dserror.H"

#include <NOX_Solver_SolverUtils.H>
#include <NOX_Epetra_Group.H>
#include <NOX_Direction_Generic.H>
#include <NOX_StatusTest_Generic.H>
#include <NOX_LineSearch_Generic.H>

#include <Epetra_Vector.h>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Solver::PseudoTransient::PseudoTransient(
    const Teuchos::RCP<NOX::Abstract::Group>& grp,
    const Teuchos::RCP<NOX::StatusTest::Generic>& outerTests,
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests,
    const Teuchos::RCP<Teuchos::ParameterList>& params)
    : NOX::NLN::Solver::LineSearchBased(grp,outerTests,innerTests,params),
      iTestPtr_(innerTests),
      xDot_(Teuchos::null)
{
  init();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::init()
{
  checkType = NOX::Solver::parseStatusTestCheckType(paramsPtr->sublist("Solver Options"));

  // We use our own factories at this point to generate the line search
  // and the direction object.
  LineSearchBased::init(iTestPtr_);

  const Teuchos::ParameterList& p_ptc = paramsPtr->sublist("Pseudo Transient");
  deltaInit_ = p_ptc.get<double>("deltaInit");
  delta_ = deltaInit_;
  invDelta_ = 1.0/delta_;
  deltaMax_ = p_ptc.get<double>("deltMax");
  deltaMin_ = p_ptc.get<double>("deltaMin");
  pseudoTime_ = 0.0;

//  use_transient_residual_ =
//      p_ptc.get<bool>("Use Transient Residual in Direction Computation");

  maxPseudoTransientIterations_ =
      p_ptc.get<int>("Maximum Number of Pseudo Transient Iterations");

  // get the time step control type
  const std::string& control_str = p_ptc.get<std::string>("Time Step Control");
  tscType_ = String2TSCType(control_str);
  const std::string& norm_str = p_ptc.get<std::string>("Norm Type for Time Step Control");
  normType_ = NOX::NLN::AUX::String2NormType(norm_str);

  // get the scaling operator type
  const std::string& scaleop_str = p_ptc.get<std::string>("Scaling Operator");
  scaleOpType_ = String2ScaleOpType(scaleop_str);

  // ---------------------------------------------------------------------------
  // Begin: Setup the pre/post operator to modify the jacobian
  // ---------------------------------------------------------------------------
  prePostLinSysPtr_ = Teuchos::rcp(new NOX::NLN::LinSystem::PrePostOp::PseudoTransient(*this));
  // The name of the direction method and the corresponding sublist has to match!
  const std::string& dir_str = paramsPtr->sublist("Direction").get<std::string>("Method");
  if (paramsPtr->sublist("Direction").isSublist(dir_str))
  {
    // set the new pre/post operator for the linear system in the parameter list
    Teuchos::ParameterList& p_linsolver = paramsPtr->sublist("Direction").
        sublist(dir_str).sublist("Linear Solver");
    p_linsolver.set<Teuchos::RCP<NOX::NLN::LinSystem::PrePostOp::Generic> >
        ("User Defined Pre/Post Operator",prePostLinSysPtr_);
    /* Now the last thing to do is, that we have to reset the pre/post-operator in
     * the nln linear system class with the ptc pre/post-operator.
     */
    // Get the linear system
    Teuchos::RCP<NOX::Epetra::Group> epetra_soln =
        Teuchos::rcp_dynamic_cast<NOX::Epetra::Group>(solnPtr);
    Teuchos::RCP<NOX::NLN::LinearSystem>linsys =
        Teuchos::rcp_dynamic_cast<NOX::NLN::LinearSystem>(epetra_soln->getLinearSystem());
    if (linsys.is_null())
      throwError("init()","The linear system cast failed!");
    /* Call reset on the linear system
     * This resets also the pre/post-operator to the
     * NOX::NLN::LinSystem::PrePostOp::PseudoTransient() pre/post operator.
     */
    linsys->reset();
  }
  else
  {
    std::ostringstream msg;
    msg << "The name of the \"Method\" in the \"Direction\" sublist has to match the name "
        << "of the corresponding sub-sublist! There is no \"Direction->"
        << dir_str <<"\" sub-sublist!" << std::endl;
    throwError("init()",msg.str());
  }
  // ---------------------------------------------------------------------------
  // End: Setup the pre/post operator to modify the jacobian
  // ---------------------------------------------------------------------------


  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::reset(
    const NOX::Abstract::Vector& initialGuess,
    const Teuchos::RCP<NOX::StatusTest::Generic>& outerTests,
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests)
{
  solnPtr->setX(initialGuess);
  if (not outerTests.is_null())
    testPtr = outerTests;
  if (not innerTests.is_null())
    iTestPtr_ = innerTests;

  init();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::reset(
    const NOX::Abstract::Vector& initialGuess,
    const Teuchos::RCP<NOX::StatusTest::Generic>& outerTests)
{
  solnPtr->setX(initialGuess);
  if (not outerTests.is_null())
    testPtr = outerTests;

  init();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::StatusTest::StatusType NOX::NLN::Solver::PseudoTransient::step()
{
  prePostOperator.runPreIterate(*this);

  // On the first step do some initializations
  if (nIter == 0)
  {
    // Compute F of initital guess
    NOX::Abstract::Group::ReturnType rtype = solnPtr->computeF();
    if (rtype != NOX::Abstract::Group::Ok)
    {
      utilsPtr->out() << "NOX::Solver::PseudoTransient::init - "
              << "Unable to compute F" << std::endl;
      throw "NOX Error";
    }

    // Test the initial guess
    status = testPtr->checkStatus(*this, checkType);
    if ((status == NOX::StatusTest::Converged) and
        (utilsPtr->isPrintType(NOX::Utils::Warning)))
    {
      utilsPtr->out() << "Warning: NOX::Solver::PseudoTransient::init() - "
              << "The solution passed into the solver (either "
              << "through constructor or reset method) "
              << "is already converged!  The solver will not "
              << "attempt to solve this system since status is "
              << "flagged as converged." << std::endl;
    }
    // Update the pseudo time step size once at the beginning
    updatePseudoTimeStep();

    printUpdate();
  }

  // First check status
  if (status != NOX::StatusTest::Unconverged)
  {
    prePostOperator.runPostIterate(*this);
    printUpdate();
    return status;
  }

  // Copy pointers into temporary references
  NOX::Abstract::Group& soln = *solnPtr;
  NOX::StatusTest::Generic& otest = *testPtr;

  // Only necessary for the optional CFL scaling option
  computePseudoVelocity();

  /* Compute the direction for the update vector at the current
   * solution. Steady-state F is already computed so the only thing
   * to compute is J. The necessary non-linear solver dependent changes
   * to the Jacobian are done by using the runPreJacobianInverse() routine of the
   * pre/post-operator object of the NOX::NLN::LinearSystem class and
   * its derived classes.
   *
   * See the NOX::NLN::LinSystem::PrePostOP::PseudoTransient class for
   * more information.
   */
  bool ok = directionPtr->compute(*dirPtr,soln,*this);

  if (not ok)
  {
    utilsPtr->out() << "NOX::NLN::Solver::PseudoTransient::step - unable to calculate direction"
        << std::endl;
    status = NOX::StatusTest::Failed;
    prePostOperator.runPostIterate(*this);
    printUpdate();
    return status;
  }

  // Update iteration count.
  nIter ++;

  // Copy current soln to the old soln
  // (Changes the owner of the linear system. New owner is the old_soln_ptr_)
  *oldSolnPtr = *solnPtr;
  /* Update the pseudo time step size before the line search routine starts.
   * Otherwise, the line search routine could contradict the idea of the update
   * routine, since the line-search modifies the step-length in such a way, that
   * the given merit function value decreases. In many cases this coincides with
   * a reduced norm of the right-hand-side F and this will lead to an increased
   * pseudo time step. Therefore the pseudo time step would be increased even if
   * the actual behavior of the last solve step was bad.
   */
  // full step to get the trial point for the pseudo time step update routine
  solnPtr->computeX(*oldSolnPtr,*dirPtr,1.0);
  NOX::Abstract::Group::ReturnType rtype = solnPtr->computeF();
  if (rtype != NOX::Abstract::Group::Ok)
    throwError("compute","Unable to compute F!");
  // Update the pseudo time step with unmodified step length of the last Newton
  // step.
  updatePseudoTimeStep();

  /* Additional line search (optional).
   * You have to specify an inner status test to use it, otherwise
   * it will do nothing, but updating x with the specified default
   * step length (same behavior as for the "Full Step" case).
   */
  // Do line search and compute new soln.
  ok = lineSearchPtr->compute(soln, stepSize, *dirPtr, *this);
  if (!ok)
  {
    if (stepSize == 0.0)
    {
      utilsPtr->out() << "NOX::NLN::Solver::PseudoTransient::step - line search failed" << std::endl;
      status = NOX::StatusTest::Failed;
      prePostOperator.runPostIterate(*this);
      printUpdate();
      return status;
    }
    else if (utilsPtr->isPrintType(NOX::Utils::Warning))
      utilsPtr->out() << "NOX::NLN::Solver::PseudoTransient::step - using recovery step for line search"
              << std::endl;
  }

  // Compute F for new current solution.
  rtype = soln.computeF();
  if (rtype != NOX::Abstract::Group::Ok)
  {
    utilsPtr->out() << "NOX::Solver::PseudoTransient::iterate - unable to compute F" << std::endl;
    status = NOX::StatusTest::Failed;
    prePostOperator.runPostIterate(*this);
    printUpdate();
    return status;
  }

  // Evaluate the current status.
  status = otest.checkStatus(*this,checkType);

  prePostOperator.runPostIterate(*this);

  printUpdate();

  return status;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::StatusTest::StatusType NOX::NLN::Solver::PseudoTransient::solve()
{
  prePostOperator.runPreSolve(*this);

  // Iterate until converged or failed
  while (status == NOX::StatusTest::Unconverged)
    step();

  Teuchos::ParameterList& outputParams = paramsPtr->sublist("Output");
  outputParams.set("Nonlinear Iterations", nIter);
  outputParams.set("2-Norm of Residual", solnPtr->getNormF());

  prePostOperator.runPostSolve(*this);

  return status;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::computePseudoVelocity()
{
  if (nIter >= maxPseudoTransientIterations_ or
      scaleOpType_==scale_op_identity)
    return;

  if (xDot_.is_null())
    xDot_ = solnPtr->getX().clone(ShapeCopy);

  // computing xDot_ using finite differences
  if (nIter==0)
    xDot_->init(0.0);
  else
    xDot_->update(invDelta_,solnPtr->getX(),-invDelta_,oldSolnPtr->getX(),0.0);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::updatePseudoTimeStep()
{
  if (nIter < maxPseudoTransientIterations_)
  {
    deltaOld_ = delta_;

    if (nIter==0)
    {
      delta_ = deltaInit_;
      return;
    }

    switch (tscType_)
    {
      case tsc_ser:
        if (normType_ == NOX::Abstract::Vector::TwoNorm)
          delta_ = deltaOld_ * oldSolnPtr->getNormF() /
              solnPtr->getNormF();
        else
          delta_ = deltaOld_ * oldSolnPtr->getF().norm(normType_) /
              solnPtr->getF().norm(normType_);
        break;
      case tsc_tte:
        throwError("UpdatePseudoTimeStep()","The \"Temporal Truncation Error\" method is not yet implemented!");
        break;
      default:
        throwError("UpdatePseudoTimeStep()","Unknown time step control type.");
        break;
    }

    invDelta_ = 1.0 / delta_;
    if (delta_ > deltaMax_)
      invDelta_ = 0.0;

    if (delta_ < deltaMin_)
    {
      delta_ = deltaMin_;
      invDelta_ = 1.0 / delta_;
    }

    pseudoTime_ += delta_;
  }
  else
  {
    delta_ = std::numeric_limits<double>::max();
    invDelta_ = 0.0;
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const double& NOX::NLN::Solver::PseudoTransient::getInversePseudoTimeStep()
    const
{
  return invDelta_;
}

const enum NOX::NLN::Solver::PseudoTransient::ScaleOpType&
    NOX::NLN::Solver::PseudoTransient::getScalingOperatorType() const
{
  return scaleOpType_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::printUpdate()
{
  double normSoln = 0;
  double normStep = 0;

  // Print the status test parameters at each iteration if requested
  if ((status == NOX::StatusTest::Unconverged) &&
      (utilsPtr->isPrintType(NOX::Utils::OuterIterationStatusTest)))
  {
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Status Test Results --\n";
    testPtr->print(utilsPtr->out());
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
  }

  // All processes participate in the computation of these norms...
  if (utilsPtr->isPrintType(NOX::Utils::OuterIteration))
  {
    normSoln = solnPtr->getNormF();
    normStep = (nIter > 0) ? dirPtr->norm() : 0;
  }

  // ...But only the print process actually prints the result.
  if (utilsPtr->isPrintType(NOX::Utils::OuterIteration))
  {
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Nonlinear Solver Step " << nIter << " -- \n";
    utilsPtr->out() << "||F|| = " << utilsPtr->sciformat(normSoln);
    utilsPtr->out() << "  step = " << utilsPtr->sciformat(stepSize);
    utilsPtr->out() << "  dx = " << utilsPtr->sciformat(normStep);
    utilsPtr->out() << "  dt_ptc = " << utilsPtr->sciformat(delta_);
    if (status == NOX::StatusTest::Converged)
      utilsPtr->out() << " (Converged!)";
    if (status == NOX::StatusTest::Failed)
      utilsPtr->out() << " (Failed!)";
    utilsPtr->out() << "\n" << NOX::Utils::fill(72) << "\n" << std::endl;
  }

  // Print the final parameter values of the status test
  if ((status != NOX::StatusTest::Unconverged) &&
      (utilsPtr->isPrintType(NOX::Utils::OuterIteration)))
  {
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Final Status Test Results --\n";
    testPtr->print(utilsPtr->out());
    utilsPtr->out() << NOX::Utils::fill(72) << "\n";
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::throwError(
    const std::string& functionName,
    const std::string& errorMsg) const
{
  if (utilsPtr->isPrintType(NOX::Utils::Error))
  {
    std::ostringstream msg;
    msg << "ERROR - NOX::NLN::Sovler::PseudoTransient::" << functionName
     << " - " << errorMsg << std::endl;
    dserror(msg.str());
  }
  throw "NOX Error";
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::LinSystem::PrePostOp::PseudoTransient::PseudoTransient(
    const NOX::NLN::Solver::PseudoTransient& ptcsolver)
    : ptcsolver_(ptcsolver)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinSystem::PrePostOp::PseudoTransient::
    runPreApplyJacobianInverse(
    NOX::Abstract::Vector& rhs,
    LINALG::SparseOperator& jac,
    const NOX::NLN::LinearSystem& linsys)
{
  // get the type of the jacobian
  const enum NOX::NLN::LinearSystem::OperatorType& jactype =
      linsys.getJacobianOperatorType();

  switch (jactype)
  {
    case NOX::NLN::LinearSystem::LinalgSparseMatrix:
    {
      // First cast the LINALG::SparseOperator and do an additional safety
      // check.
      LINALG::SparseMatrix* jacPtr = dynamic_cast<LINALG::SparseMatrix*>(&jac);
      if (jacPtr == NULL)
        dserror("Something strange happened: The jacobian has not the "
            "operator type defined in the linear system object!");

      modifyJacobian(*jacPtr);

      break;
    }
    default:
    {
      dserror("Unsupported jacobian operator type: %s",
          NOX::NLN::LinearSystem::OperatorType2String(jactype).c_str());
      break;
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinSystem::PrePostOp::PseudoTransient::
    runPostApplyJacobianInverse(
    NOX::Abstract::Vector& rhs,
    LINALG::SparseOperator& jac,
    const NOX::NLN::LinearSystem& linsys)
{
  /* ToDo At the moment it seems unnecessary to undo the changes to the
   * jacobian. Nevertheless, if it becomes necessary, this would be the
   * right place.
   */
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinSystem::PrePostOp::PseudoTransient::modifyJacobian(
    LINALG::SparseMatrix& jac)
{
  // get the inverse pseudo time step
  const double& deltaInv = ptcsolver_.getInversePseudoTimeStep();
  const enum NOX::NLN::Solver::PseudoTransient::ScaleOpType scaleoptype =
      ptcsolver_.getScalingOperatorType();

  switch (scaleoptype)
  {
    case NOX::NLN::Solver::PseudoTransient::scale_op_identity:
    {
      /* Build the scaling operator V and multiply it with the inverse
       * pseudo time step. Finally, we modify the jacobian.
       *
       *        $ (\delta^{-1} \boldsymbol{I} + \boldsymbol{J})$
       */
      Teuchos::RCP<Epetra_Vector> v = LINALG::CreateVector(jac.RowMap(),false);
      v->PutScalar(deltaInv);
      // get the diagonal terms of the jacobian
      Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(jac.RowMap(),false);
      jac.ExtractDiagonalCopy(*diag);
      diag->Update(1.0,*v,1.0);
      // Finally modify the jacobian
      jac.ReplaceDiagonalValues(*diag);
      // NOTE: Since we modify only the diagonal terms the Dirichlet
      // conditions can be neglected.
      break;
    }
    case NOX::NLN::Solver::PseudoTransient::scale_op_cfl_diagonal:
    {
      dserror("Not yet implemented!");
      break;
    }
    default:
    {
      dserror("Unknown/unsupported scaling operator type!");
      break;
    }
  }

  return;
}

