/*-----------------------------------------------------------*/
/*!
\file nox_nln_inner_statustest_upperbound.cpp

\brief inner status test that restricts value of update vector

\maintainer Maximilian Grill

\date 07/2016

\level 3
*/
/*-----------------------------------------------------------*/

#include "nox_nln_inner_statustest_upperbound.H"
#include "nox_nln_linesearch_generic.H"
#include "nox_nln_group.H"

#include <NOX_Utils.H>
#include <NOX_Abstract_Vector.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Solver_Generic.H>
#include <Epetra_Vector.h>


#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::UpperBound::UpperBound(const double& upperboundval,
    const NOX::Abstract::Vector::NormType normtype, const NOX::NLN::StatusTest::QuantityType qtype)
    : status_(status_unevaluated),
      normtype_(normtype),
      qtype_(qtype),
      upperboundval_(upperboundval),
      reduction_fac_(1.0),
      stepmaxval_(0.0)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::INNER::StatusTest::UpperBound::GetSearchDirectionLength(
    const NOX::NLN::LineSearch::Generic& linesearch, const NOX::Solver::Generic& solver,
    const NOX::Abstract::Group& grp) const
{
  const NOX::NLN::Group& nln_grp = dynamic_cast<const NOX::NLN::Group&>(grp);

  return nln_grp.GetTrialUpdateNorm(linesearch.GetSearchDirection(), normtype_, qtype_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::StatusType NOX::NLN::INNER::StatusTest::UpperBound::CheckStatus(
    const NOX::NLN::INNER::StatusTest::Interface::Required& interface,
    const NOX::Solver::Generic& solver, const NOX::Abstract::Group& grp,
    NOX::StatusTest::CheckType checkType)
{
  /* check if it is a line search object: upper bound for Newton step size only
   * makes sense as inner status test for line search solvers */
  const NOX::NLN::LineSearch::Generic* linesearch =
      dynamic_cast<const NOX::NLN::LineSearch::Generic*>(&interface);
  if (linesearch == NULL)
  {
    std::ostringstream msg;
    msg << "Dynamic cast to NOX::NLN::LineSearch::Generic failed!\n\n"
        << "The UpperBound rule status test supports only Line Search problems!";
    throwError("CheckStatus", msg.str());
  }

  /* we reduce the step length according to the upper bound criterion in the first
   * line search (i.e. inner) iteration and do nothing in all following iterations */
  if (interface.GetNumIterations() == 0)
  {
    const double dir_length = GetSearchDirectionLength(*linesearch, solver, grp);
    double steplength = linesearch->GetStepLength();

    // compute specified norm
    stepmaxval_ = steplength * dir_length;

    // check the value for specified upper bound
    status_ = (stepmaxval_ < upperboundval_) ? status_converged : status_step_too_long;

    /* in opposite to classical line search, there is no need to find the new step
     * length iteratively.
     * instead, we can immediately reduce the step length appropriately and avoid
     * several unnecessary evaluations of rhs (computeF calls). */
    if (status_ == status_step_too_long)
    {
      NOX::NLN::LineSearch::Generic* linesearch_mutable =
          const_cast<NOX::NLN::LineSearch::Generic*>(linesearch);

      /* the following is equivalent to dividing the step successively by two until
       * criterion is met. note: upperboundval_!=0 is checked in
       * NOX::NLN::INNER::StatusTest::Factory::BuildUpperBoundTest */
      reduction_fac_ =
          std::pow(0.5, std::ceil(std::log(stepmaxval_ / upperboundval_) / std::log(2)));

      steplength *= reduction_fac_;
      linesearch_mutable->SetStepLength(steplength);

      // adapt the stepmaxval_ variable accordingly to get correct output from Print()
      stepmaxval_ *= reduction_fac_;

      status_ = status_converged;
    }
    // nothing to do. reset the the reduction factor to its default value
    else
      reduction_fac_ = 1.0;
  }
  else if (checkType == NOX::StatusTest::None)
  {
    /*!
      The test can (and should, if possible) be skipped if
      checkType is NOX::StatusType::None. If the test is skipped, then
      the status should be set to NOX::StatusTest::Unevaluated.
    */
    status_ = status_unevaluated;
  }

  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::StatusType NOX::NLN::INNER::StatusTest::UpperBound::GetStatus() const
{
  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& NOX::NLN::INNER::StatusTest::UpperBound::Print(std::ostream& stream, int indent) const
{
  std::string indent_string;
  indent_string.assign(indent, ' ');

  stream << indent_string;
  stream << status_;
  stream << " ";
  stream << "Upper Bound for Newton step size: ";
  stream << NOX::Utils::sciformat(stepmaxval_, 3) << " < "
         << NOX::Utils::sciformat(upperboundval_, 3) << "\n";
  stream << indent_string << NOX::Utils::fill(13, ' ')
         << " (reduction factor = " << NOX::Utils::sciformat(reduction_fac_, 3) << ")\n";

  return stream;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::UpperBound::throwError(
    const std::string& functionName, const std::string& errorMsg) const
{
  std::ostringstream msg;
  msg << "ERROR - NOX::NLN::INNER::StatusTest::UpperBound::" << functionName << " - " << errorMsg
      << std::endl;
  dserror(msg.str());
}
