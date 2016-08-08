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

#include <NOX_Utils.H>
#include <NOX_Abstract_Group.H>
#include <NOX_Abstract_Vector.H>
#include <NOX_Epetra_Vector.H>
#include <Epetra_Vector.h>


#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::UpperBound::UpperBound(
    const double& upperboundval,
    const NOX::Abstract::Vector::NormType& normtype) :
    status_(status_unevaluated),
    upperboundval_(upperboundval),
    normtype_(normtype),
    stepmaxval_(0.0)
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::UpperBound::Setup(
    const NOX::NLN::LineSearch::Generic& linesearch,
    const NOX::Abstract::Group& grp)
{
  /* ToDo setup a DofMap that specifies the entries of the update vector
   * which should be checked for upper bound (e.g. only structural position DoFs) */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::StatusType
    NOX::NLN::INNER::StatusTest::UpperBound::CheckStatus(
    const NOX::NLN::INNER::StatusTest::Interface::Required& interface,
    const NOX::Abstract::Group& grp,
    NOX::StatusTest::CheckType checkType)
{
  /* check if it is a line search object: upper bound for Newton step size only
   * makes sense as inner status test for line search solvers */
  const NOX::NLN::LineSearch::Generic* linesearch =
      dynamic_cast<const NOX::NLN::LineSearch::Generic*>(&interface);
  if(linesearch == NULL)
  {
    std::ostringstream msg;
    msg << "Dynamic cast to NOX::NLN::LineSearch::Generic failed!\n\n"
        << "The UpperBound rule status test supports only Line Search problems!";
    throwError("CheckStatus",msg.str());
  }

  // setup for the current line search loop
  if (interface.GetNumIterations()==0)
  {
    Setup(*linesearch,grp);
    status_ = status_unevaluated;
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
  else
  {
    // compute specified norm
    stepmaxval_ = linesearch->GetSearchDirection().norm(normtype_)
        * linesearch->GetStepLength();

    // check the value for specified upper bound
    status_ = (stepmaxval_ < upperboundval_) ? status_converged : status_step_too_long;
  }

  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::StatusType
    NOX::NLN::INNER::StatusTest::UpperBound::GetStatus() const
{
  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& NOX::NLN::INNER::StatusTest::UpperBound::Print(
    std::ostream& stream,
    int indent) const
{
  std::string indent_string;
  indent_string.assign(indent,' ');

  stream << indent_string;
  stream << status_;
  stream << " ";
  stream << "Upper Bound for Newton step size: ";
  stream << NOX::Utils::sciformat(stepmaxval_,3)
         << " < "
         << NOX::Utils::sciformat(upperboundval_,3) << "\n";

  return stream;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::INNER::StatusTest::UpperBound::throwError(
    const std::string& functionName,
    const std::string& errorMsg) const
{
  std::ostringstream msg;
  msg << "ERROR - NOX::NLN::INNER::StatusTest::UpperBound::" << functionName
      << " - " << errorMsg << std::endl;
  dserror(msg.str());
}
