/*-----------------------------------------------------------*/
/*! \file

\brief inner status test that restricts value of update vector



\level 3
*/
/*-----------------------------------------------------------*/

#include "4C_solver_nonlin_nox_inner_statustest_upperbound.hpp"

#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_solver_nonlin_nox_linesearch_generic.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Vector.h>
#include <NOX_Abstract_Vector.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Solver_Generic.H>
#include <NOX_Utils.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Inner::StatusTest::UpperBound::UpperBound(const double& upperboundval,
    const ::NOX::Abstract::Vector::NormType normtype,
    const NOX::Nln::StatusTest::QuantityType qtype)
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
double NOX::Nln::Inner::StatusTest::UpperBound::get_search_direction_length(
    const NOX::Nln::LineSearch::Generic& linesearch, const ::NOX::Solver::Generic& solver,
    const ::NOX::Abstract::Group& grp) const
{
  const NOX::Nln::Group& nln_grp = dynamic_cast<const NOX::Nln::Group&>(grp);

  return nln_grp.GetTrialUpdateNorm(linesearch.GetSearchDirection(), normtype_, qtype_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Inner::StatusTest::StatusType NOX::Nln::Inner::StatusTest::UpperBound::CheckStatus(
    const NOX::Nln::Inner::StatusTest::Interface::Required& interface,
    const ::NOX::Solver::Generic& solver, const ::NOX::Abstract::Group& grp,
    ::NOX::StatusTest::CheckType checkType)
{
  /* check if it is a line search object: upper bound for Newton step size only
   * makes sense as inner status test for line search solvers */
  const NOX::Nln::LineSearch::Generic* linesearch =
      dynamic_cast<const NOX::Nln::LineSearch::Generic*>(&interface);
  if (linesearch == nullptr)
  {
    std::ostringstream msg;
    msg << "Dynamic cast to NOX::Nln::LineSearch::Generic failed!\n\n"
        << "The UpperBound rule status test supports only Line Search problems!";
    throw_error("CheckStatus", msg.str());
  }

  /* we reduce the step length according to the upper bound criterion in the first
   * line search (i.e. inner) iteration and do nothing in all following iterations */
  if (interface.GetNumIterations() == 0)
  {
    const double dir_length = get_search_direction_length(*linesearch, solver, grp);
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
      NOX::Nln::LineSearch::Generic* linesearch_mutable =
          const_cast<NOX::Nln::LineSearch::Generic*>(linesearch);

      /* the following is equivalent to dividing the step successively by two until
       * criterion is met. note: upperboundval_!=0 is checked in
       * NOX::Nln::Inner::StatusTest::Factory::build_upper_bound_test */
      reduction_fac_ =
          std::pow(0.5, std::ceil(std::log(stepmaxval_ / upperboundval_) / std::log(2)));

      steplength *= reduction_fac_;
      linesearch_mutable->SetStepLength(steplength);

      // adapt the stepmaxval_ variable accordingly to get correct output from print()
      stepmaxval_ *= reduction_fac_;

      status_ = status_converged;
    }
    // nothing to do. reset the the reduction factor to its default value
    else
      reduction_fac_ = 1.0;
  }
  else if (checkType == ::NOX::StatusTest::None)
  {
    /*!
      The test can (and should, if possible) be skipped if
      checkType is NOX::StatusType::None. If the test is skipped, then
      the status should be set to ::NOX::StatusTest::Unevaluated.
    */
    status_ = status_unevaluated;
  }

  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Inner::StatusTest::StatusType NOX::Nln::Inner::StatusTest::UpperBound::GetStatus() const
{
  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& NOX::Nln::Inner::StatusTest::UpperBound::print(std::ostream& stream, int indent) const
{
  std::string indent_string;
  indent_string.assign(indent, ' ');

  stream << indent_string;
  stream << status_;
  stream << " ";
  stream << "Upper Bound for Newton step size: ";
  stream << ::NOX::Utils::sciformat(stepmaxval_, 3) << " < "
         << ::NOX::Utils::sciformat(upperboundval_, 3) << "\n";
  stream << indent_string << ::NOX::Utils::fill(13, ' ')
         << " (reduction factor = " << ::NOX::Utils::sciformat(reduction_fac_, 3) << ")\n";

  return stream;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Inner::StatusTest::UpperBound::throw_error(
    const std::string& functionName, const std::string& errorMsg) const
{
  std::ostringstream msg;
  msg << "ERROR - NOX::Nln::Inner::StatusTest::UpperBound::" << functionName << " - " << errorMsg
      << std::endl;
  FOUR_C_THROW(msg.str());
}

FOUR_C_NAMESPACE_CLOSE
