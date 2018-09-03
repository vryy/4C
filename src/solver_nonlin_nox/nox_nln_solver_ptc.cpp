/*-----------------------------------------------------------*/
/*!
\file nox_nln_solver_ptc.cpp

\brief %NOX::NLN implementation of a pseudo transient non-linear
       solver.

\maintainer Michael Hiermeier

\date Oct 6, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_solver_ptc.H"  // class definition
#include "nox_nln_linesearch_factory.H"
#include "nox_nln_direction_factory.H"
#include "nox_nln_aux.H"
#include "nox_nln_group.H"
#include "nox_nln_linearsystem.H"
#include "nox_nln_statustest_normf.H"
#include "nox_nln_group_prepostoperator.H"
#include "nox_nln_linearsystem_prepostoperator.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_dserror.H"

#include "../drt_geometry/intersection_math.H"

#include <NOX_Solver_SolverUtils.H>
#include <NOX_Direction_Generic.H>
#include <NOX_StatusTest_Generic.H>
#include <NOX_LineSearch_Generic.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_GlobalData.H>
#include <NOX_MeritFunction_Generic.H>

#include <Epetra_Vector.h>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Solver::PseudoTransient::PseudoTransient(const Teuchos::RCP<NOX::Abstract::Group>& grp,
    const Teuchos::RCP<NOX::StatusTest::Generic>& outerTests,
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests,
    const Teuchos::RCP<Teuchos::ParameterList>& params)
    : NOX::NLN::Solver::LineSearchBased(grp, outerTests, innerTests, params),
      iTestPtr_(innerTests),
      xDot_(Teuchos::null),
      usePseudoTransientResidual_(false),
      calcDeltaInit_(false),
      isScalingOperator_(false),
      modelReductionRatio_(0.0)
{
  init();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::init()
{
  checkType = NOX::Solver::parseStatusTestCheckType(paramsPtr->sublist("Solver Options"));

  /* We use our own factories at this point to generate the line search
   * and the direction object. */
  LineSearchBased::init(iTestPtr_);

  const Teuchos::ParameterList& p_ptc = paramsPtr->sublist("Pseudo Transient");
  deltaInit_ = p_ptc.get<double>("deltaInit");
  if (deltaInit_ == 0.0)
    throwError("init()", "The initial pseudo time step is not allowed to be equal to 0.0!");
  else if (deltaInit_ < 0.0)
    calcDeltaInit_ = true;
  delta_ = deltaInit_;
  invDelta_ = 1.0 / delta_;
  deltaMax_ = p_ptc.get<double>("deltaMax");
  deltaMin_ = p_ptc.get<double>("deltaMin");
  SER_alpha_ = p_ptc.get<double>("SER_alpha");
  scaleFactor_ = p_ptc.get<double>("ScalingFactor");
  pseudoTime_ = 0.0;

  maxPseudoTransientIterations_ = p_ptc.get<int>("Max Number of PTC Iterations");

  // get the time step control type
  const std::string& control_str = p_ptc.get<std::string>("Time Step Control");
  tscType_ = String2TSCType(control_str);
  const std::string& norm_str = p_ptc.get<std::string>("Norm Type for TSC");
  normType_ = NOX::NLN::AUX::String2NormType(norm_str);

  // get the scaling operator type
  const std::string& scaleop_str = p_ptc.get<std::string>("Scaling Type");
  scaleOpType_ = String2ScaleOpType(scaleop_str);

  const std::string& build_scaling_op = p_ptc.get<std::string>("Build scaling operator");
  build_scaling_op_ = String2BuildOpType(build_scaling_op);

  // create the scaling operator
  createScalingOperator();

  // create the linear system pre/post operator
  createLinSystemPrePostOperator();

  // create the group pre/post operator
  createGroupPrePostOperator();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::reset(const NOX::Abstract::Vector& initialGuess,
    const Teuchos::RCP<NOX::StatusTest::Generic>& outerTests,
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests)
{
  solnPtr->setX(initialGuess);
  if (not outerTests.is_null()) testPtr = outerTests;
  if (not innerTests.is_null()) iTestPtr_ = innerTests;

  init();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::reset(const NOX::Abstract::Vector& initialGuess,
    const Teuchos::RCP<NOX::StatusTest::Generic>& outerTests)
{
  solnPtr->setX(initialGuess);
  if (not outerTests.is_null()) testPtr = outerTests;

  init();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::createScalingOperator()
{
  if (IsScalingOperator()) return;

  switch (scaleOpType_)
  {
    case NOX::NLN::Solver::PseudoTransient::scale_op_identity:
    {
      // identity matrix
      Teuchos::RCP<const NOX::Epetra::Vector> epetraXPtr =
          Teuchos::rcp_dynamic_cast<const NOX::Epetra::Vector>(solnPtr->getXPtr());
      if (epetraXPtr.is_null()) dserror("Cast to NOX::Epetra::Vector failed!");
      scalingDiagOpPtr_ =
          Teuchos::rcp(new Epetra_Vector(epetraXPtr->getEpetraVector().Map(), false));

      scalingDiagOpPtr_->PutScalar(1.0);

      isScalingOperator_ = true;

      break;
    }
    case NOX::NLN::Solver::PseudoTransient::scale_op_lumped_mass:
    {
      // get the lumped mass matrix
      scalingDiagOpPtr_ = Teuchos::rcp(new Epetra_Vector(
          *(Teuchos::rcp_dynamic_cast<NOX::NLN::Group>(solnPtr)->GetLumpedMassMatrixPtr())));
      break;
    }
    // get element based scaling operator
    case NOX::NLN::Solver::PseudoTransient::scale_op_element_based:
    {
      // Get element based contributions and assemble into SparseMatrix
      scalingMatrixOpPtr_ =
          Teuchos::rcp_dynamic_cast<NOX::NLN::Group>(solnPtr)->GetContributionsFromElementLevel();

      if (build_scaling_op_ == NOX::NLN::Solver::PseudoTransient::build_op_everyiter)
        isScalingOperator_ = false;
      else
        isScalingOperator_ = true;

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
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::createLinSystemPrePostOperator()
{
  prePostLinSysPtr_ = Teuchos::rcp(new NOX::NLN::LinSystem::PrePostOp::PseudoTransient(
      scalingDiagOpPtr_, scalingMatrixOpPtr_, *this));

  // The name of the direction method and the corresponding sublist has to match!
  const std::string& dir_str = paramsPtr->sublist("Direction").get<std::string>("Method");
  if (paramsPtr->sublist("Direction").isSublist(dir_str))
  {
    // set the new pre/post operator for the linear system in the parameter list
    Teuchos::ParameterList& p_linsolver =
        paramsPtr->sublist("Direction").sublist(dir_str).sublist("Linear Solver");
    // get the current map. If there is no map, return a new empty one. (reference)
    NOX::NLN::LinSystem::PrePostOperator::Map& prePostLinSystemMap =
        NOX::NLN::LinSystem::PrePostOp::GetMutableMap(p_linsolver);
    // insert/replace the old pointer in the map
    prePostLinSystemMap[NOX::NLN::LinSystem::prepost_ptc] = prePostLinSysPtr_;
    /* Now the last thing to do is, that we have to reset the pre/post-operator in
     * the nln linear system class with the ptc pre/post-operator. */
    // Get the linear system
    Teuchos::RCP<NOX::NLN::Group> nlnSolnPtr = Teuchos::rcp_dynamic_cast<NOX::NLN::Group>(solnPtr);
    if (nlnSolnPtr.is_null())
      throwError("createLinSystemPrePostOperator", "The group cast failed!");
    /* Call the reset function of the nox nln linear system to reset the
     * corresponding pre/post operator map. */
    nlnSolnPtr->ResetLinSysPrePostOperator(p_linsolver, true);
  }
  else
  {
    std::ostringstream msg;
    msg << "The name of the \"Method\" in the \"Direction\" sublist has to match the name "
        << "of the corresponding sub-sublist! There is no \"Direction->" << dir_str
        << "\" sub-sublist!" << std::endl;
    throwError("createLinSystemPrePostOperator()", msg.str());
  }
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::createGroupPrePostOperator()
{
  // ---------------------------------------------------------------------------
  // set the new pre/post operator for the nox nln group in the parameter list
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList& p_grpOpt = paramsPtr->sublist("Group Options");
  // get the current map. If there is no map, return a new empty one. (reference)
  NOX::NLN::GROUP::PrePostOperator::Map& prePostGroupMap =
      NOX::NLN::GROUP::PrePostOp::GetMutableMap(p_grpOpt);
  // insert or replace the old pointer in the map
  prePostGroupPtr_ = Teuchos::rcp(new NOX::NLN::GROUP::PrePostOp::PseudoTransient(
      scalingDiagOpPtr_, scalingMatrixOpPtr_, *this));

  // set the modified map
  prePostGroupMap[NOX::NLN::GROUP::prepost_ptc] = prePostGroupPtr_;
  /* Now the last thing to do is, that we have to reset the pre/post-operator in
   * the nln group class with the ptc pre/post-operator. */
  // Get the nln group
  Teuchos::RCP<NOX::NLN::Group> nlnSoln = Teuchos::rcp_dynamic_cast<NOX::NLN::Group>(solnPtr);
  if (nlnSoln.is_null()) throwError("createGroupPrePostOperator", "The group cast failed!");
  /* Call the reset function of the nox nln group to reset the corresponding
   * pre/post operator map. */
  nlnSoln->ResetPrePostOperator(p_grpOpt, true);

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
    if ((status == NOX::StatusTest::Converged) and (utilsPtr->isPrintType(NOX::Utils::Warning)))
    {
      utilsPtr->out() << "Warning: NOX::Solver::PseudoTransient::init() - "
                      << "The solution passed into the solver (either "
                      << "through constructor or reset method) "
                      << "is already converged! The solver will not "
                      << "attempt to solve this system since status is "
                      << "flagged as converged." << std::endl;
    }
    printUpdate();
  }
  else
  {
    // eventually create new scaling operator each iteration step
    createScalingOperator();

    /* Adjust the pseudo time step, such it fits to the chosen line search
     * step length. */
    adjustPseudoTimeStep();
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

  // Update the pseudo time step size.
  updatePseudoTimeStep();

  // Only necessary for the optional CFL scaling option
  computePseudoVelocity();

  /* Compute the direction for the update vector at the current
   * solution. Steady-state F is already computed so the only thing
   * to compute is J. The necessary non-linear solver dependent changes
   * to the Jacobian are done by using the runPostComputeJacobian() routine of the
   * pre/post-operator object of the NOX::NLN::LinearSystem class and
   * its derived classes.
   *
   * See the NOX::NLN::LinSystem::PrePostOP::PseudoTransient class for
   * more information. */
  bool ok = directionPtr->compute(*dirPtr, soln, *this);

  // Show the linear solver residual || Ax-b ||_2
  //  double lin_residual = 0.0;
  //  soln.getNormLastLinearSolveResidual(lin_residual);
  //  utilsPtr->out() << "Linear Solver Residual (L2-norm): " << lin_residual << std::endl;

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
  nIter++;

  // Copy current soln to the old soln
  // (Changes the owner of the linear system. New owner is the oldSolnPtr_)
  *oldSolnPtr = *solnPtr;

  /* Additional line search (optional).
   * You have to specify an inner status test to use it, otherwise
   * it will do nothing, but updating x with the specified default
   * step length (same behavior as for the "Full Step" case). */
  // Use the pseudo transient residual during the line search procedure.
  usePseudoTransientResidual_ = true;
  // Do line search and compute new soln.
  ok = lineSearchPtr->compute(soln, stepSize, *dirPtr, *this);
  /* call the computeF routine again, to be sure that it has been evaluated!
   * (this becomes necessary for a Full Step method call!) */
  solnPtr->computeF();
  usePseudoTransientResidual_ = false;
  // evaluate the model reduction ratio if desired
  if (ok) evalModelReductionRatio();

  if (!ok)
  {
    if (stepSize == 0.0)
    {
      utilsPtr->out() << "NOX::NLN::Solver::PseudoTransient::step - line search failed"
                      << std::endl;
      status = NOX::StatusTest::Failed;
      prePostOperator.runPostIterate(*this);
      printUpdate();
      return status;
    }
    else if (utilsPtr->isPrintType(NOX::Utils::Warning))
      utilsPtr->out()
          << "NOX::NLN::Solver::PseudoTransient::step - using recovery step for line search"
          << std::endl;
  }

  // Compute F for new current solution.
  NOX::Abstract::Group::ReturnType rtype = soln.computeF();
  if (rtype != NOX::Abstract::Group::Ok)
  {
    utilsPtr->out() << "NOX::Solver::PseudoTransient::iterate - unable to compute F" << std::endl;
    status = NOX::StatusTest::Failed;
    prePostOperator.runPostIterate(*this);
    printUpdate();
    return status;
  }

  // Evaluate the current status.
  status = otest.checkStatus(*this, checkType);

  prePostOperator.runPostIterate(*this);

  printUpdate();

  return status;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::evalModelReductionRatio()
{
  if (tscType_ != tsc_mrr) return;

  double fref = globalDataPtr->getMeritFunction()->computef(*oldSolnPtr);
  double fnew = globalDataPtr->getMeritFunction()->computef(*solnPtr);
  Teuchos::RCP<NOX::Abstract::Vector> sdirPtr = dirPtr->clone(DeepCopy);
  sdirPtr->scale(stepSize);
  double modelnew = globalDataPtr->getMeritFunction()->computeQuadraticModel(*sdirPtr, *oldSolnPtr);

  modelReductionRatio_ = (fref - fnew) / (fref - modelnew);
  if (utilsPtr->isPrintType(NOX::Utils::Details))
  {
    utilsPtr->out() << "(" << fref << " - " << fnew << ") / (" << fref << " - " << modelnew << ")"
                    << std::endl;
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::StatusTest::StatusType NOX::NLN::Solver::PseudoTransient::solve()
{
  isPtcSolve_ = true;
  prePostOperator.runPreSolve(*this);
  // Iterate until converged or failed
  while (status == NOX::StatusTest::Unconverged) step();

  Teuchos::ParameterList& outputParams = paramsPtr->sublist("Output");
  outputParams.set("Nonlinear Iterations", nIter);
  outputParams.set("2-Norm of Residual", solnPtr->getNormF());

  prePostOperator.runPostSolve(*this);
  isPtcSolve_ = false;
  return status;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::computePseudoVelocity()
{
  if (nIter >= maxPseudoTransientIterations_ or scaleOpType_ == scale_op_identity) return;

  if (xDot_.is_null()) xDot_ = solnPtr->getX().clone(ShapeCopy);

  // computing xDot_ using finite differences
  if (nIter == 0)
    xDot_->init(0.0);
  else
    xDot_->update(invDelta_, solnPtr->getX(), -invDelta_, oldSolnPtr->getX(), 0.0);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::adjustPseudoTimeStep()
{
  /* Check the step-length. If the step-length is not equal 1.0,
   * we adjust the pseudo time step size by using a least squares
   * approximation. */
  if (stepSize == 1.0 or nIter >= maxPseudoTransientIterations_) return;

  if (utilsPtr->isPrintType(NOX::Utils::Details))
  {
    utilsPtr->out() << "*--- Adjust Pseudo Time Step ---*" << std::endl
                    << "| Previous dt_ptc:    " << std::setw(9) << std::setprecision(3)
                    << std::scientific << delta_ << " |" << std::endl;
  }

  Teuchos::RCP<NOX::NLN::Group> oldNlnSolnPtr =
      Teuchos::rcp_dynamic_cast<NOX::NLN::Group>(oldSolnPtr);
  if (oldNlnSolnPtr.is_null())
    throwError("updatePseudoTimeStep", "Dynamic cast to NOX::NLN::Group failed!");
  oldNlnSolnPtr->adjustPseudoTimeStep(delta_, stepSize, *dirPtr, *this);

  if (utilsPtr->isPrintType(NOX::Utils::Details))
  {
    utilsPtr->out() << "| Corrected dt_ptc:   " << std::setw(9) << std::setprecision(3)
                    << std::scientific << delta_ << " |" << std::endl
                    << "*-------------------------------*" << std::endl;
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::updatePseudoTimeStep()
{
  if (nIter < maxPseudoTransientIterations_)
  {
    deltaOld_ = delta_;

    if (nIter == 0)
    {
      if (calcDeltaInit_)
      {
        double normF = 0.0;
        if (scaleOpType_ == scale_op_identity || scaleOpType_ == scale_op_element_based)
        {
          normF = solnPtr->getF().norm(normType_);
        }
        else
        {
          Teuchos::RCP<NOX::Abstract::Vector> scaledRHS = solnPtr->getF().clone(DeepCopy);
          NOX::Epetra::Vector& epetraScaledRHS = dynamic_cast<NOX::Epetra::Vector&>(*scaledRHS);
          epetraScaledRHS.getEpetraVector().ReciprocalMultiply(
              1.0, *scalingDiagOpPtr_, epetraScaledRHS.getEpetraVector(), 0.0);
          normF = epetraScaledRHS.norm(normType_);
        }
        if (normF > GEO::TOL12)
          deltaInit_ = 1.0 / (normF * normF);
        else
          deltaInit_ = 1.0;
      }
      delta_ = deltaInit_;
    }
    else
    {
      // -----------------------------------------------------------------
      // desired time stepping control method
      // -----------------------------------------------------------------
      switch (tscType_)
      {
        // ---------------------------------------------------------------
        // switched evolution relaxation
        // ---------------------------------------------------------------
        case tsc_ser:
        {
          if (normType_ == NOX::Abstract::Vector::TwoNorm)
            delta_ =
                deltaOld_ * std::pow((oldSolnPtr->getNormF() / solnPtr->getNormF()), SER_alpha_);
          else
            delta_ = deltaOld_ *
                     std::pow(oldSolnPtr->getF().norm(normType_) / solnPtr->getF().norm(normType_),
                         SER_alpha_);
          break;
        }
        // ---------------------------------------------------------------
        // temporal truncation error
        // ---------------------------------------------------------------
        case tsc_tte:
        {
          throwError("UpdatePseudoTimeStep()",
              "The \"Temporal Truncation Error\" method is not yet implemented!");
          break;
        }
        // ---------------------------------------------------------------
        // model reduction ratio
        // ---------------------------------------------------------------
        case tsc_mrr:
        {
          double normFOld = oldSolnPtr->getF().norm(normType_);
          double normF = solnPtr->getF().norm(normType_);
          double ratioF = normFOld / normF;
          // calculated the corrected old inverse mu
          double muinv = deltaOld_ * normFOld;

          /* If the normF Status test is already converged and no additional
           * step length modification has been necessary, we will just increase
           * the pseudo time step.
           * This is meant to speed up the convergence close to the solution. */
          if (stepSize == 1.0 and
              GetStatus<NOX::NLN::StatusTest::NormF>() == NOX::StatusTest::Converged)
            muinv *= 4.0;
          /* if the model performed badly: reduce the pseudo time step by a
           * factor tau_red, where 0.25 <= tau_red <= 0.8.
           * (ToDo supposed to become an input parameter) */
          else if (modelReductionRatio_ < 0.2)
            muinv *= std::min(std::max(0.25, ratioF), 0.8);
          /* if the model performed well: increase the pseudo time step by a
           * factor tau_inc, where 1.25 <= tau_inc <= 4.0.
           * (ToDo supposed to become an input parameter) */
          else if (modelReductionRatio_ > 0.8)
            muinv *= std::max(std::min(4.0, ratioF), 1.25);

          // update the pseudo time step
          delta_ = muinv / normF;
          if (utilsPtr->isPrintType(NOX::Utils::Details))
          {
            utilsPtr->out() << "*--- Model Reduction Ratio -----------------------------*"
                            << std::endl
                            << "| MRR: " << std::setw(9) << std::setprecision(3) << std::scientific
                            << modelReductionRatio_ << std::setw(41) << "|" << std::endl
                            << "| dt_ptc(k-1) --> (    muinv *  normFinv) --> dt_ptc(k) |"
                            << std::endl
                            << "| " << std::setw(11) << std::setprecision(3) << std::scientific
                            << deltaOld_ << " --> (" << std::setw(9) << std::setprecision(3)
                            << std::scientific << muinv << " * " << std::setw(9)
                            << std::setprecision(3) << std::scientific << 1.0 / normF << ") --> "
                            << std::setw(9) << std::setprecision(3) << std::scientific << delta_
                            << " |" << std::endl
                            << "*-------------------------------------------------------*"
                            << std::endl;
          }
          break;
        }
        default:
        {
          throwError("UpdatePseudoTimeStep()", "Unknown time step control type.");
          break;
        }
      }
    }

    // invert delta_ and check possible thresholds
    if (delta_ > deltaMax_)
      invDelta_ = 0.0;
    else
      invDelta_ = 1.0 / delta_;

    if (delta_ < deltaMin_)
    {
      delta_ = deltaMin_;
      invDelta_ = 1.0 / delta_;
    }

    // increase the pseudoTime_ counter
    pseudoTime_ += delta_;
  }
  /* if the maximum PTC iteration number is reached, we switch PTC off and
   * use a standard line search based solution procedure. */
  else
  {
    delta_ = std::numeric_limits<double>::infinity();
    invDelta_ = 0.0;
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const double& NOX::NLN::Solver::PseudoTransient::getInversePseudoTimeStep() const
{
  return invDelta_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const double& NOX::NLN::Solver::PseudoTransient::getScalingFactor() const { return scaleFactor_; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const enum NOX::NLN::Solver::PseudoTransient::ScaleOpType&
NOX::NLN::Solver::PseudoTransient::getScalingOperatorType() const
{
  return scaleOpType_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Vector& NOX::NLN::Solver::PseudoTransient::getScalingDiagOperator() const
{
  if (scalingDiagOpPtr_.is_null())
    throwError("getScalingDiagOperator", "The scaling operator is not initialized!");

  return *scalingDiagOpPtr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::NLN::Solver::PseudoTransient::usePseudoTransientResidual() const
{
  return usePseudoTransientResidual_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::NLN::Solver::PseudoTransient::isPtcSolve() const { return isPtcSolve_; }

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
    utilsPtr->out() << NOX::Utils::fill(82) << "\n";
    utilsPtr->out() << "-- Status Test Results --\n";
    testPtr->print(utilsPtr->out());
    utilsPtr->out() << NOX::Utils::fill(82) << "\n";
  }

  // All processes participate in the computation of these norms...
  if (utilsPtr->isPrintType(NOX::Utils::OuterIteration))
  {
    normSoln = solnPtr->getNormF();
    normStep = (nIter > 0) ? dirPtr->norm() : 0;
  }

  // ...But only the print process actually prints the result.
  // ------ standard output ------------------------------------------
  if (utilsPtr->isPrintType(NOX::Utils::OuterIteration) and
      (utilsPtr->isPrintType(NOX::Utils::OuterIterationStatusTest) or
          utilsPtr->isPrintType(NOX::Utils::InnerIteration)))
  {
    utilsPtr->out() << "\n" << NOX::Utils::fill(82) << "\n";
    utilsPtr->out() << "-- Nonlinear Solver Step " << nIter << " -- \n";
    utilsPtr->out() << "||F|| = " << utilsPtr->sciformat(normSoln);
    utilsPtr->out() << "  step = " << utilsPtr->sciformat(stepSize);
    utilsPtr->out() << "  dx = " << utilsPtr->sciformat(normStep);
    if (delta_ == -1.0)
      utilsPtr->out() << "  dt_ptc = auto";
    else
      utilsPtr->out() << "  dt_ptc = " << utilsPtr->sciformat(delta_);
    if (status == NOX::StatusTest::Converged) utilsPtr->out() << " (Converged!)";
    if (status == NOX::StatusTest::Failed) utilsPtr->out() << " (Failed!)";
    utilsPtr->out() << "\n" << NOX::Utils::fill(82) << "\n" << std::endl;
  }
  // ------ short output ---------------------------------------------
  else if (utilsPtr->isPrintType(NOX::Utils::OuterIteration))
  {
    // print the head line
    if (nIter == 0)
    {
      utilsPtr->out() << std::setw(4) << "#It" << std::setw(13) << "||F||_2" << std::setw(13)
                      << "step" << std::setw(13) << "||dx||_2" << std::setw(13) << "dt_ptc\n";
      utilsPtr->out() << NOX::Utils::fill(60, '^') << "\n";
    }
    utilsPtr->out() << std::setw(4) << nIter;
    utilsPtr->out() << "  " << utilsPtr->sciformat(normSoln);
    utilsPtr->out() << "  " << utilsPtr->sciformat(stepSize);
    utilsPtr->out() << "  " << utilsPtr->sciformat(normStep);
    if (delta_ == -1.0)
      utilsPtr->out() << std::setw(13) << "auto";
    else
      utilsPtr->out() << "  " << utilsPtr->sciformat(delta_);
    if (status == NOX::StatusTest::Converged) utilsPtr->out() << " (Converged!)";
    if (status == NOX::StatusTest::Failed) utilsPtr->out() << " (Failed!)";
    utilsPtr->out() << std::endl;
  }

  // Print the final parameter values of the status test
  if ((status != NOX::StatusTest::Unconverged) &&
      (utilsPtr->isPrintType(NOX::Utils::OuterIterationStatusTest)))
  {
    utilsPtr->out() << NOX::Utils::fill(82) << "\n";
    utilsPtr->out() << "-- Final Status Test Results --\n";
    testPtr->print(utilsPtr->out());
    utilsPtr->out() << NOX::Utils::fill(82) << "\n";
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::Solver::PseudoTransient::throwError(
    const std::string& functionName, const std::string& errorMsg) const
{
  if (utilsPtr->isPrintType(NOX::Utils::Error))
  {
    std::ostringstream msg;
    msg << "ERROR - NOX::NLN::Sovler::PseudoTransient::" << functionName << " - " << errorMsg
        << std::endl;
    dserror(msg.str());
  }
  throw "NOX Error";
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::LinSystem::PrePostOp::PseudoTransient::PseudoTransient(
    Teuchos::RCP<Epetra_Vector>& scalingDiagOpPtr,
    Teuchos::RCP<LINALG::SparseMatrix>& scalingMatrixOpPtr,
    const NOX::NLN::Solver::PseudoTransient& ptcsolver)
    : ptcsolver_(ptcsolver),
      scalingDiagOpPtr_(scalingDiagOpPtr),
      scalingMatrixOpPtr_(scalingMatrixOpPtr)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinSystem::PrePostOp::PseudoTransient::runPostComputeJacobian(
    LINALG::SparseOperator& jac, const Epetra_Vector& x, const NOX::NLN::LinearSystem& linsys)
{
  if (not ptcsolver_.isPtcSolve()) return;

  // get the type of the jacobian
  const enum NOX::NLN::LinSystem::OperatorType& jactype = linsys.getJacobianOperatorType();

  switch (jactype)
  {
    case NOX::NLN::LinSystem::LinalgSparseMatrix:
    {
      // First cast the LINALG::SparseOperator and do an additional sanity
      // check.
      LINALG::SparseMatrix* jacPtr = dynamic_cast<LINALG::SparseMatrix*>(&jac);
      if (jacPtr == NULL)
        dserror(
            "Something strange happened: The jacobian has not the "
            "operator type defined in the linear system object!");

      modifyJacobian(*jacPtr);

      break;
    }
    default:
    {
      dserror("Unsupported jacobian operator type: %s",
          NOX::NLN::LinSystem::OperatorType2String(jactype).c_str());
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinSystem::PrePostOp::PseudoTransient::runPostComputeFandJacobian(Epetra_Vector& rhs,
    LINALG::SparseOperator& jac, const Epetra_Vector& x, const NOX::NLN::LinearSystem& linsys)
{
  runPostComputeJacobian(jac, x, linsys);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinSystem::PrePostOp::PseudoTransient::modifyJacobian(LINALG::SparseMatrix& jac)
{
  // get the inverse pseudo time step
  const double& deltaInv = ptcsolver_.getInversePseudoTimeStep();
  const enum NOX::NLN::Solver::PseudoTransient::ScaleOpType& scaleoptype =
      ptcsolver_.getScalingOperatorType();
  const double& scaleFactor = ptcsolver_.getScalingFactor();

  switch (scaleoptype)
  {
    case NOX::NLN::Solver::PseudoTransient::scale_op_identity:
    case NOX::NLN::Solver::PseudoTransient::scale_op_lumped_mass:
    {
      /* Build the scaling operator V and multiply it with the inverse
       * pseudo time step. Finally, we modify the jacobian.
       *
       *        (\delta^{-1} \boldsymbol{I} + \boldsymbol{J}) */
      Teuchos::RCP<Epetra_Vector> v = Teuchos::rcp(new Epetra_Vector(*scalingDiagOpPtr_));
      // Scale v with scaling factor
      v->Scale(deltaInv * scaleFactor);
      // get the diagonal terms of the jacobian
      Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(jac.RowMap(), false);
      jac.ExtractDiagonalCopy(*diag);
      diag->Update(1.0, *v, 1.0);
      // Finally modify the jacobian
      jac.ReplaceDiagonalValues(*diag);
      break;
    }
    case NOX::NLN::Solver::PseudoTransient::scale_op_element_based:
    {
      /* Build the scaling operator V and multiply it with the inverse
       * pseudo time step. Finally, we modify the jacobian.
       *
       *        (\delta^{-1} \boldsymbol{V} + \boldsymbol{J}) */

      scalingMatrixOpPtr_->Complete();
      jac.Add(*scalingMatrixOpPtr_, false, scaleFactor * deltaInv, 1.0);
      jac.Complete();

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

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::GROUP::PrePostOp::PseudoTransient::PseudoTransient(
    Teuchos::RCP<Epetra_Vector>& scalingDiagOpPtr,
    Teuchos::RCP<LINALG::SparseMatrix>& scalingMatrixOpPtr,
    const NOX::NLN::Solver::PseudoTransient& ptcsolver)
    : ptcsolver_(ptcsolver),
      scalingDiagOpPtr_(scalingDiagOpPtr),
      scalingMatrixOpPtr_(scalingMatrixOpPtr),
      isPseudoTransientResidual_(false)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::Vector>
NOX::NLN::GROUP::PrePostOp::PseudoTransient::evalPseudoTransientFUpdate(const NOX::NLN::Group& grp)
{
  // get the current trial point
  Teuchos::RCP<NOX::Epetra::Vector> xUpdate =
      Teuchos::rcp_dynamic_cast<NOX::Epetra::Vector>(grp.getX().clone(DeepCopy));

  // get the old x vector
  const NOX::Abstract::Group& oldGrp = ptcsolver_.getPreviousSolutionGroup();
  const NOX::Epetra::Vector xOld = dynamic_cast<const NOX::Epetra::Vector&>(oldGrp.getX());

  /* Calculate the difference between the old and the new solution vector.
   * This is equivalent to the search direction scaled with the
   * step size. */
  xUpdate->update(-1.0, xOld, 1.0);

  const enum NOX::NLN::Solver::PseudoTransient::ScaleOpType& scaleoptype =
      ptcsolver_.getScalingOperatorType();
  switch (scaleoptype)
  {
    case NOX::NLN::Solver::PseudoTransient::scale_op_identity:
    {
      NOX::Epetra::Vector v = NOX::Epetra::Vector(scalingDiagOpPtr_);
      v.scale(ptcsolver_.getInversePseudoTimeStep());
      xUpdate->scale(v);

      break;
    }
    case NOX::NLN::Solver::PseudoTransient::scale_op_element_based:
    {
      scalingMatrixOpPtr_->Multiply(false, xUpdate->getEpetraVector(), xUpdate->getEpetraVector());
      xUpdate->scale(ptcsolver_.getInversePseudoTimeStep());

      break;
    }
    default:
    {
      dserror("Unknown scaling operator type!");
      break;
    }
  }
  return xUpdate;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::GROUP::PrePostOp::PseudoTransient::runPostComputeF(
    Epetra_Vector& F, const NOX::NLN::Group& grp)
{
  if (not ptcsolver_.isPtcSolve()) return;

  const bool usePseudoTransientResidual = ptcsolver_.usePseudoTransientResidual();

  /* If we need no pseudo transient residual or if the transient residual
   * has already been added, we can skip this function. */
  if (not usePseudoTransientResidual or isPseudoTransientResidual_) return;

  Teuchos::RCP<NOX::Epetra::Vector> fUpdate = evalPseudoTransientFUpdate(grp);

  // add the transient part
  F.Update(1.0, fUpdate->getEpetraVector(), 1.0);

  // set the flag
  isPseudoTransientResidual_ = true;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::GROUP::PrePostOp::PseudoTransient::runPreComputeF(
    Epetra_Vector& F, const NOX::NLN::Group& grp)
{
  if (not ptcsolver_.isPtcSolve()) return;

  // If the current rhs has not been calculated, yet.
  if (!grp.isF())
  {
    isPseudoTransientResidual_ = false;
    return;
  }

  /* Recalculate the static residual, if the current right hand side
   * has already been modified, though we need the static residual. */
  if (ptcsolver_.usePseudoTransientResidual() or !isPseudoTransientResidual_) return;

  Teuchos::RCP<NOX::Epetra::Vector> fUpdate = evalPseudoTransientFUpdate(grp);

  // subtract the transient part
  F.Update(-1.0, fUpdate->getEpetraVector(), 1.0);

  // set flag
  isPseudoTransientResidual_ = false;

  return;
}
