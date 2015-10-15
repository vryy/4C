/*-----------------------------------------------------------*/
/*!
\file nox_nln_linesearch_backtrack.cpp

\maintainer Michael Hiermeier

\date Jun 11, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_linesearch_backtrack.H" // class definition

#include "../drt_lib/drt_dserror.H"

#include <NOX_Utils.H>
#include <NOX_GlobalData.H>
#include <NOX_Solver_Generic.H>
#include <NOX_Epetra_Group.H>

#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::LineSearch::Backtrack::Backtrack(
    const Teuchos::RCP<NOX::GlobalData>& gd,
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic> innerTests,
          Teuchos::ParameterList& params)
    : stepPtr_(NULL),
      checkType_(NOX::StatusTest::Complete),
      innerTestsPtr_(innerTests)
{
  reset(gd, params);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::LineSearch::Backtrack::reset
   (const Teuchos::RCP<NOX::GlobalData>& gd,
    Teuchos::ParameterList& params)
{
  Teuchos::ParameterList& p = params.sublist("Backtrack");

  utils_ = gd->getUtils();
  meritFunctionPtr_ = gd->getMeritFunction();

  lsIters_ = 0;
  searchDirectionPtr_ = Teuchos::null;

  status_ = NOX::NLN::INNER::StatusTest::status_unevaluated;

  defaultStep_ = p.get("Default Step", 1.0);
  reductionFactor_ = p.get("Reduction Factor", 0.5);
  if ((reductionFactor_ <= 0.0)  || (reductionFactor_ >= 1.0))
  {
    std::ostringstream msg;
    msg << "Invalid choice \"" << reductionFactor_ << "\" for \"Reduction Factor\"!\n"
        << "Value must be greater than zero and less than 1.0.";
    throwError("reset",msg.str());
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::LineSearch::Backtrack::reset()
{
  lsIters_ = 0;
  searchDirectionPtr_ = Teuchos::null;

  status_ = NOX::NLN::INNER::StatusTest::status_unevaluated;

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::LineSearch::Backtrack::compute
    (NOX::Abstract::Group& grp, double& step,
     const NOX::Abstract::Vector& dir,
     const NOX::Solver::Generic& s)
{
  // -------------------------------------------------
  // (re)set important line search parameters
  // -------------------------------------------------
  reset();
  // get the old solution group
  const NOX::Abstract::Group& oldGrp = s.getPreviousSolutionGroup();
  // update the search direction pointer
  searchDirectionPtr_ = Teuchos::rcpFromRef(dir);
  // set the step pointer to the inserted step variable
  stepPtr_ = &step;
  // reset the step length
  step = defaultStep_;
  // initialize the inner status test
  // ----------------------------------------------------------------------
  // BE CAREFUL HERE:
  // During the copy operation in NOX::Solver::LineSearchBased::step()
  // the current grp loses the ownership of the sharedLinearSystem. If we
  // want to access the jacobian, we have to use the oldGrp
  // (target of the copy process), instead.                hiermeier 08/15
  // ----------------------------------------------------------------------
  // ToDo Check what is happening in the PTC case, since we call computeX once
  // on the soln grp. Maybe this changes the ownership and we have to reset
  // the ownership to the old grp before we can call the compute routine of the
  // backtracking linesearch.
  const NOX::Epetra::Group* epetraOldGrpPtr =
        dynamic_cast<const NOX::Epetra::Group*>(&oldGrp);
  if (not epetraOldGrpPtr->isJacobian())
    throwError("compute()","Ownership changed unexpectedly!");

  status_ = innerTestsPtr_->CheckStatus(*this,oldGrp,checkType_);
  if ((status_ == NOX::NLN::INNER::StatusTest::status_converged) and
      (utils_->isPrintType(NOX::Utils::Warning)))
  {
    utils_->out() << "Warning: NOX::NLN::LineSearch::Backtrack::compute() - "
        << "The solution passed into the line search method fulfills the inner\n"
        << "stopping criterion already during the setup process.\n"
        << "The line search method will not attempt to find a capable "
        << "step length." << std::endl;
  }

  // increase iteration counter after initialization
  ++lsIters_;

  // -------------------------------------------------
  // update the solution vector and get a trial point
  // -------------------------------------------------
  grp.computeX(oldGrp, dir, step);

  NOX::Abstract::Group::ReturnType rtype = grp.computeF();
  if (rtype != NOX::Abstract::Group::Ok)
    throwError("compute","Unable to compute F!");

  // -------------------------------------------------
  // print header if desired
  // -------------------------------------------------
  if (utils_->isPrintType(NOX::Utils::InnerIteration))
  {
    utils_->out() << "\n" << NOX::Utils::fill(72) << "\n"
        << "-- Backtrack Line Search -- \n";
  }

  status_ = innerTestsPtr_->CheckStatus(*this,grp,checkType_);
  PrintUpdate();
  // -------------------------------------------------
  // inner backtracking loop
  // -------------------------------------------------
  while (status_ == NOX::NLN::INNER::StatusTest::status_step_too_long)
  {
    // -------------------------------------------------
    // reduce step length
    // -------------------------------------------------
    step *= reductionFactor_;

    // -------------------------------------------------
    // update the solution vector and get a trial point
    // -------------------------------------------------
    grp.computeX(oldGrp, dir, step);

    rtype = grp.computeF();
    if (rtype != NOX::Abstract::Group::Ok)
      throwError("compute","Unable to compute F!");
    status_ = innerTestsPtr_->CheckStatus(*this,grp,checkType_);
    PrintUpdate();
  }
  // -------------------------------------------------
  // print footer if desired
  // -------------------------------------------------
  if (utils_->isPrintType(NOX::Utils::InnerIteration))
  {
    utils_->out() << NOX::Utils::fill(72) << "\n";
  }

  if (status_ == NOX::NLN::INNER::StatusTest::status_step_too_short)
    throwError("compute()","The current step is too short and no "
        "restoration phase is implemented!");
  else if (status_ == NOX::NLN::INNER::StatusTest::status_no_descent_direction)
    throwError("compute()","The given search direction is no descent direction!");

  return (status_ == NOX::NLN::INNER::StatusTest::status_converged ? true : false);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const int& NOX::NLN::LineSearch::Backtrack::GetNumIterations() const
{
  return lsIters_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const NOX::MeritFunction::Generic& NOX::NLN::LineSearch::Backtrack::GetMeritFunction() const
{
  if (meritFunctionPtr_.is_null())
    throwError("GetMeritFunction","The merit function pointer is not initialized!");

  return *meritFunctionPtr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const NOX::Abstract::Vector& NOX::NLN::LineSearch::Backtrack::GetSearchDirection() const
{
  if (searchDirectionPtr_.is_null())
    throwError("GetSearchDirection","The search direction ptr is not initialized!");

  return *searchDirectionPtr_;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const double& NOX::NLN::LineSearch::Backtrack::GetStepLength() const
{
  if (stepPtr_==NULL)
    throwError("GetStepLength","Step pointer is NULL!");

  return *stepPtr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LineSearch::Backtrack::PrintUpdate()
{
  // Print the status test parameters at each iteration if requested
  if ((status_ == NOX::NLN::INNER::StatusTest::status_step_too_long) and
      (utils_->isPrintType(NOX::Utils::InnerIteration)))
  {
    utils_->out() << NOX::Utils::fill(72,'-') << "\n";
    utils_->out() << "-- Inner Status Test Results --\n";
    innerTestsPtr_->Print(utils_->out());
    utils_->out() << NOX::Utils::fill(72,'-') << "\n";
  }
  // Print the final parameter values of the status test
  if ((status_ != NOX::NLN::INNER::StatusTest::status_step_too_long) and
      (utils_->isPrintType(NOX::Utils::InnerIteration)))
  {
    utils_->out() << NOX::Utils::fill(72,'-') << "\n";
    utils_->out() << "-- Final Inner Status Test Results --\n";
    innerTestsPtr_->Print(utils_->out());
    utils_->out() << NOX::Utils::fill(72,'-') << "\n";
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LineSearch::Backtrack::throwError(
    const std::string& functionName,
    const std::string& errorMsg) const
{
    std::ostringstream msg;
    msg << "ERROR - NOX::NLN::LineSearch::Backtrack::" << functionName
        << " - " << errorMsg << std::endl;
    dserror(msg.str());
}
