/*-----------------------------------------------------------*/
/*!

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_inner_statustest_armijo.H"  // class definition
#include "nox_nln_linesearch_generic.H"

#include <NOX_Utils.H>
#include <NOX_MeritFunction_Generic.H>
#include <NOX_Abstract_Group.H>
#include <NOX_Epetra_Vector.H>
#include <Epetra_Vector.h>

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::Armijo::Armijo(
    const double& c_1, const bool& isMonotone, const std::size_t& maxHistSize)
    : status_(status_unevaluated),
      c_1_(c_1),
      fref_(0.0),
      fcurr_(0.0),
      slope_(0.0),
      step_(1.0),
      isMonotone_(isMonotone),
      maxHistSize_(maxHistSize),
      histVector_(std::deque<double>(0))
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool NOX::NLN::INNER::StatusTest::Armijo::Setup(
    const NOX::NLN::LineSearch::Generic& linesearch, const NOX::Abstract::Group& grp)
{
  const NOX::MeritFunction::Generic& mrtFct = linesearch.GetMeritFunction();

  // get the reference merit function value
  fref_ = mrtFct.computef(grp);

  // get the slope once (doesn't change during the inner iteration)
  slope_ = mrtFct.computeSlope(linesearch.GetSearchDirection(), grp);

  // return false if the search direction is no descent direction
  if (slope_ >= 0.0) return false;

  // -------------------------------------------
  // Non-monotone setup
  // -------------------------------------------
  if (not isMonotone_)
  {
    // Compare the current merit function value with the last
    // added value of the history vector and augment the history
    // vector only if the two values differ in more than
    // machine precision.
    if (fref_ != histVector_.front())
    {
      histVector_.push_front(fref_);
      // remove the last element if the maximum size of the history vector
      // is reached
      if (histVector_.size() > maxHistSize_) histVector_.pop_back();
    }
    // get the maximal merit function value of the last accepted steps for
    // the Armijo check
    fref_ = *std::max_element(histVector_.begin(), histVector_.end());
  }

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::StatusType NOX::NLN::INNER::StatusTest::Armijo::CheckStatus(
    const NOX::NLN::INNER::StatusTest::Interface::Required& interface,
    const NOX::Solver::Generic& solver, const NOX::Abstract::Group& grp,
    NOX::StatusTest::CheckType checkType)
{
  // check if it is a line search object
  // Amrijo rule plays only a role as inner status test for line search solvers
  const NOX::NLN::LineSearch::Generic* linesearch =
      dynamic_cast<const NOX::NLN::LineSearch::Generic*>(&interface);
  if (linesearch == NULL)
  {
    std::ostringstream msg;
    msg << "Dynamic cast to NOX::NLN::LineSearch::Generic failed!\n\n"
        << "The Armijo rule status test supports only Line Search problems!";
    throwError("CheckStatus", msg.str());
  }

  // setup for the current line search loop
  if (interface.GetNumIterations() == 0)
  {
    // If the search direction is no descent direction,
    // this function detects it and returns the corresponding
    // status.
    if (Setup(*linesearch, grp))
      status_ = status_unevaluated;
    else
      status_ = status_no_descent_direction;
  }
  else if (checkType == NOX::StatusTest::None)
  {
    fcurr_ = 0.0;
    status_ = status_unevaluated;
  }
  // If the setup call detected a non-descent direction,
  // we will not check the Armijo rule, since the check will
  // fail anyway.
  else if (status_ != status_no_descent_direction)
  {
    const NOX::MeritFunction::Generic& mrtFct = linesearch->GetMeritFunction();

    fcurr_ = mrtFct.computef(grp);

    step_ = linesearch->GetStepLength();

    // check the Armijo rule
    status_ = (fcurr_ < fref_ + c_1_ * step_ * slope_) ? status_converged : status_step_too_long;
  }

  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::StatusType NOX::NLN::INNER::StatusTest::Armijo::GetStatus() const
{
  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& NOX::NLN::INNER::StatusTest::Armijo::Print(std::ostream& stream, int indent) const
{
  std::string indent_string;
  indent_string.assign(indent, ' ');

  stream << indent_string;
  stream << status_;
  if (isMonotone_)
    stream << "Monotone";
  else
    stream << "Non-monotone";
  stream << " ";
  stream << "Armijo-Rule: ";
  stream << NOX::Utils::sciformat(fcurr_, 3) << " < "
         << NOX::Utils::sciformat(fref_ + c_1_ * step_ * slope_, 3) << "\n";

  stream << indent_string;
  stream << std::setw(13) << " ";
  stream << "(step = " << NOX::Utils::sciformat(step_, 3);
  stream << ", slope = " << NOX::Utils::sciformat(slope_, 3);
  if (not isMonotone_) stream << ", history = " << maxHistSize_;
  stream << ")\n";

  return stream;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::Armijo::throwError(
    const std::string& functionName, const std::string& errorMsg) const
{
  std::ostringstream msg;
  msg << "ERROR - NOX::NLN::INNER::StatusTest::Armijo::" << functionName << " - " << errorMsg
      << std::endl;
  dserror(msg.str());
}
