/*-----------------------------------------------------------*/
/*!
\file nox_nln_statustest_activeset.cpp

\maintainer Michael Hiermeier

\date Oct 30, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_statustest_activeset.H"
#include "nox_nln_constraint_group.H"

#include <NOX_Solver_Generic.H>

#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::StatusTest::ActiveSet::ActiveSet(
    const NOX::NLN::StatusTest::QuantityType& qtype)
    : qtype_(qtype),
      status_(NOX::StatusTest::Unevaluated),
      cyclesize_(0),
      activesetsize_(0)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::StatusTest::StatusType NOX::NLN::StatusTest::ActiveSet::checkStatus(
    const NOX::Solver::Generic& problem,
    NOX::StatusTest::CheckType checkType)
{
  if (checkType == NOX::StatusTest::None)
  {
    status_ = NOX::StatusTest::Unevaluated;
    activesetsize_ = 0;
    cyclesize_ = 0;
  }
  else
  {
    // get the abstract solution group from the non-linear solver
    const NOX::Abstract::Group& grp = problem.getSolutionGroup();

    // check if the right hand side was already updated
    if (grp.isF())
      status_ = NOX::StatusTest::Unevaluated;
    else
    {
      // try to cast the nox group
      const NOX::NLN::CONSTRAINT::Group* cnlngrp =
          dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
      if (cnlngrp == NULL)
        dserror("NOX::NLN::CONSTRAINT::Group cast failed");

      // do the actual active set check
      status_ = cnlngrp->GetActiveSetStatus(qtype_,activesetsize_,cyclesize_);
    }
  }

  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::StatusTest::StatusType NOX::NLN::StatusTest::ActiveSet::getStatus() const
{
  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream& NOX::NLN::StatusTest::ActiveSet::print(std::ostream& stream,
    int indent) const
{
  std::string indent_string;
  indent_string.assign(indent,' ');

  stream << indent_string;
  stream << status_;
  stream << QuantityType2String(qtype_) << "-";
  stream << "Active-Set-Size = " << activesetsize_;
  stream << std::endl;
  // optional output
  if (cyclesize_>0)
    stream << indent_string << std::setw(13) << "CYCLING OF THE ACTIVE SET | CYCLE-SIZE IS "
    << cyclesize_ << std::endl;

  return stream;
}
