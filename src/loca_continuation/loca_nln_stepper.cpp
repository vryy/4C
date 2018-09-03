/*-----------------------------------------------------------*/
/*!
\file loca_nln_stepper.cpp

\maintainer Michael Hiermeier

\date Nov 16, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "loca_nln_stepper.H"
#include "../drt_lib/drt_dserror.H"

#include "../solver_nonlin_nox/nox_nln_inner_statustest_generic.H"
#include "../solver_nonlin_nox/nox_nln_solver_factory.H"

#include <LOCA_MultiContinuation_AbstractStrategy.H>
#include <LOCA_MultiContinuation_AbstractGroup.H>
#include <LOCA_MultiContinuation_ExtendedGroup.H>
#include <LOCA_GlobalData.H>
#include <LOCA_Factory.H>
#include <LOCA_Eigensolver_AbstractStrategy.H>
#include <LOCA_SaveEigenData_AbstractStrategy.H>
#include <LOCA_ErrorCheck.H>
#include <LOCA_Parameter_SublistParser.H>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LOCA::NLN::Stepper::Stepper(const Teuchos::RCP<LOCA::GlobalData>& loca_gdata,
    const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& initial_guess,
    const Teuchos::RCP<LOCA::StatusTest::Abstract>& loca_test,
    const Teuchos::RCP<NOX::StatusTest::Generic>& nox_otest,
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& nox_itest,
    const Teuchos::RCP<NOX::NLN::GlobalData>& nln_gdata,
    const Teuchos::RCP<Teuchos::ParameterList>& p)
    : LOCA::Stepper(loca_gdata, initial_guess, loca_test, nox_otest, p),
      noxInnerStatusTestPtr_(nox_itest),
      noxNlnGlobalDataPtr_(nln_gdata)
{
  // see base class constructor for more information
  resetSolver();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool LOCA::NLN::Stepper::reset(const Teuchos::RCP<LOCA::GlobalData>& loca_gdata,
    const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& initial_guess,
    const Teuchos::RCP<NOX::StatusTest::Generic>& nox_test,
    const Teuchos::RCP<Teuchos::ParameterList>& p)
{
  dserror("Not supported, since this method is deprecated!");
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool LOCA::NLN::Stepper::reset(const Teuchos::RCP<LOCA::GlobalData>& loca_gdata,
    const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& initial_guess,
    const Teuchos::RCP<LOCA::StatusTest::Abstract>& loca_test,
    const Teuchos::RCP<NOX::StatusTest::Generic>& nox_otest,
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& nox_itest,
    const Teuchos::RCP<NOX::NLN::GlobalData>& nln_gdata,
    const Teuchos::RCP<Teuchos::ParameterList>& p)
{
  noxInnerStatusTestPtr_ = nox_itest;
  noxNlnGlobalDataPtr_ = nln_gdata;

  bool check = LOCA::Stepper::reset(loca_gdata, initial_guess, loca_test, nox_otest, p);
  if (!check) dserror("LOCA::Stepper::reset() failed!");

  // reset the NOX::Solver object with own provided factory
  check = resetSolver();

  return check;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool LOCA::NLN::Stepper::resetSolver()
{
  solverPtr = NOX::NLN::Solver::BuildSolver(
      curGroupPtr, noxStatusTestPtr, noxInnerStatusTestPtr_, noxNlnGlobalDataPtr_);

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LOCA::Abstract::Iterator::IteratorStatus LOCA::NLN::Stepper::start()
{
  NOX::StatusTest::StatusType solverStatus;
  std::string callingFunction = "LOCA::NLN::Stepper::start()";

  // Allow continuation group to preprocess the step
  curGroupPtr->preProcessContinuationStep(LOCA::Abstract::Iterator::Successful);

  printStartStep();

  // Perform solve of initial conditions
  solverStatus = solverPtr->solve();

  // Allow continuation group to postprocess the step
  if (solverStatus == NOX::StatusTest::Converged)
    curGroupPtr->postProcessContinuationStep(LOCA::Abstract::Iterator::Successful);
  else
    curGroupPtr->postProcessContinuationStep(LOCA::Abstract::Iterator::Unsuccessful);

  // Set up continuation groups
  const LOCA::MultiContinuation::ExtendedGroup& constSolnGrp =
      dynamic_cast<const LOCA::MultiContinuation::ExtendedGroup&>(solverPtr->getSolutionGroup());
  Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> underlyingGroup =
      Teuchos::rcp_const_cast<LOCA::MultiContinuation::AbstractGroup>(
          constSolnGrp.getUnderlyingGroup());

  // Create continuation strategy
  curGroupPtr = globalData->locaFactory->createContinuationStrategy(
      parsedParams, stepperList, underlyingGroup, predictor, conParamIDs);

  // Do printing (stepNumber==0 case) after continuation group set up
  if (solverStatus == NOX::StatusTest::Failed)
    printEndStep(LOCA::Abstract::Iterator::Unsuccessful);
  else
    printEndStep(LOCA::Abstract::Iterator::Successful);

  // Set the initial step size
  curGroupPtr->setStepSize(stepSize);

  prevGroupPtr =
      Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::AbstractStrategy>(curGroupPtr->clone());

  // If nonlinear solve failed, return (this must be done after continuation
  // groups are created so Stepper::getSolutionGroup() functions correctly.
  if (solverStatus != NOX::StatusTest::Converged) return LOCA::Abstract::Iterator::Failed;

  // Save initial solution
  curGroupPtr->printSolution();

  // Compute eigenvalues/eigenvectors if requested
  if (calcEigenvalues)
  {
    Teuchos::RCP<std::vector<double>> evals_r;
    Teuchos::RCP<std::vector<double>> evals_i;
    Teuchos::RCP<NOX::Abstract::MultiVector> evecs_r;
    Teuchos::RCP<NOX::Abstract::MultiVector> evecs_i;
    eigensolver->computeEigenvalues(
        *curGroupPtr->getBaseLevelUnderlyingGroup(), evals_r, evals_i, evecs_r, evecs_i);

    saveEigenData->save(evals_r, evals_i, evecs_r, evecs_i);
  }

  // Compute predictor direction
  NOX::Abstract::Group::ReturnType predictorStatus = curGroupPtr->computePredictor();
  globalData->locaErrorCheck->checkReturnType(predictorStatus, callingFunction);
  curPredictorPtr = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(
      curGroupPtr->getPredictorTangent()[0].clone(NOX::DeepCopy));
  prevPredictorPtr = Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::ExtendedVector>(
      curGroupPtr->getPredictorTangent()[0].clone(NOX::ShapeCopy));

  // Create new solver using new continuation groups and combo status test
  resetSolver();

  return LOCA::Abstract::Iterator::NotFinished;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LOCA::Abstract::Iterator::IteratorStatus LOCA::NLN::Stepper::finish(
    LOCA::Abstract::Iterator::IteratorStatus itStatus)
{
  std::string callingFunction = "LOCA::NLN::Stepper::finish()";

  /* We don't need to check if the last step was successful since finish
   * is never called if it wasn't.  We might want to change that if there is
   * some post processing we want to do even if the run failed. */

  // Copy last solution
  curGroupPtr->copy(solverPtr->getSolutionGroup());

  // Return if iteration failed (reached max number of steps)
  if (itStatus == LOCA::Abstract::Iterator::Failed) return itStatus;

  bool do_target = stepperList->get("Hit Continuation Bound", true);
  if (!do_target) return LOCA::Abstract::Iterator::Finished;

  // Do one additional step using natural continuation to hit target value
  double value = curGroupPtr->getContinuationParameter();

  if (fabs(value - targetValue) > 1.0e-15 * (1.0 + fabs(targetValue)))
  {
    isTargetStep = true;

    // Save previous successful step information
    prevGroupPtr->copy(*curGroupPtr);

    // Get bifurcation group if there is one, or solution group if not
    Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> underlyingGrp =
        curGroupPtr->getUnderlyingGroup();

    // Create predictor strategy
    Teuchos::RCP<Teuchos::ParameterList> lastStepPredictorParams =
        parsedParams->getSublist("Last Step Predictor");
    // change default method to constant to avoid infinite stack recursion
    lastStepPredictorParams->get("Method", "Constant");
    predictor =
        globalData->locaFactory->createPredictorStrategy(parsedParams, lastStepPredictorParams);

    // Make a copy of the parameter list, change continuation method to
    // natural
    Teuchos::RCP<Teuchos::ParameterList> lastStepperParams =
        Teuchos::rcp(new Teuchos::ParameterList(*stepperList));
    lastStepperParams->set("Continuation Method", "Natural");

    // Create continuation strategy
    curGroupPtr = globalData->locaFactory->createContinuationStrategy(
        parsedParams, lastStepperParams, underlyingGrp, predictor, conParamIDs);

    // Set step size
    stepSize = targetValue - value;
    curGroupPtr->setStepSize(stepSize);

    // Get predictor direction
    NOX::Abstract::Group::ReturnType predictorStatus = curGroupPtr->computePredictor();
    globalData->locaErrorCheck->checkReturnType(predictorStatus, callingFunction);
    *curPredictorPtr = curGroupPtr->getPredictorTangent()[0];

    // Set previous solution vector in current solution group
    curGroupPtr->setPrevX(curGroupPtr->getX());

    // Take step in predictor direction
    curGroupPtr->computeX(*curGroupPtr, *curPredictorPtr, stepSize);

    // Allow continuation group to preprocess the step
    curGroupPtr->preProcessContinuationStep(LOCA::Abstract::Iterator::Successful);

    printStartStep();

    // Create new solver
    resetSolver();

    // Solve step
    NOX::StatusTest::StatusType solverStatus = solverPtr->solve();

    // Allow continuation group to postprocess the step
    if (solverStatus == NOX::StatusTest::Converged)
      curGroupPtr->postProcessContinuationStep(LOCA::Abstract::Iterator::Successful);
    else
      curGroupPtr->postProcessContinuationStep(LOCA::Abstract::Iterator::Unsuccessful);

    // Get solution
    curGroupPtr->copy(solverPtr->getSolutionGroup());

    if (solverStatus != NOX::StatusTest::Converged)
    {
      printEndStep(LOCA::Abstract::Iterator::Unsuccessful);
      return LOCA::Abstract::Iterator::Failed;
    }

    printEndStep(LOCA::Abstract::Iterator::Successful);

    curGroupPtr->printSolution();
  }

  return LOCA::Abstract::Iterator::Finished;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LOCA::Abstract::Iterator::StepStatus LOCA::NLN::Stepper::preprocess(
    LOCA::Abstract::Iterator::StepStatus stepStatus)
{
  if (stepStatus == LOCA::Abstract::Iterator::Unsuccessful)
  {
    // Restore previous step information
    curGroupPtr->copy(*prevGroupPtr);
  }
  else
  {
    // Save previous successful step information
    prevGroupPtr->copy(*curGroupPtr);
  }

  // Compute step size
  stepStatus = computeStepSize(stepStatus, stepSize);

  // Set step size in current solution group
  curGroupPtr->setStepSize(stepSize);

  // Set previous solution vector in current solution group
  curGroupPtr->setPrevX(prevGroupPtr->getX());

  // Take step in predictor direction
  curGroupPtr->computeX(*prevGroupPtr, *curPredictorPtr, stepSize);

  // Allow continuation group to preprocess the step
  curGroupPtr->preProcessContinuationStep(stepStatus);

  // Reset solver to compute new solution
  resetSolver();

  return stepStatus;
}
