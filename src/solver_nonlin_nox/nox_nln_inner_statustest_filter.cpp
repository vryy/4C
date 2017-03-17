/*----------------------------------------------------------------------------*/
/*!
\file nox_nln_inner_statustest_filter.cpp

\brief Inner status test class for constraint problems. Filter
       techniques are based on ideas from multi-objective optimization:

       - Control of the two distinct goals of minimization of the objective
         function and satisfaction of the constraints.

       - Unlike merit functions, filter methods keep these two goals separate

\maintainer Michael Hiermeier

\date Mar 6, 2017

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "nox_nln_inner_statustest_filter.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::INNER::StatusTest::Filter::Filter()
    : status_( status_unevaluated )
{
  /* intentionally left blank */
}



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::NLN::INNER::StatusTest::StatusType
NOX::NLN::INNER::StatusTest::Filter::CheckStatus(
    const Interface::Required &  interface,
    const NOX::Abstract::Group & grp,
    NOX::StatusTest::CheckType   checkType)
{
  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::NLN::INNER::StatusTest::StatusType
NOX::NLN::INNER::StatusTest::Filter::GetStatus() const
{
  return status_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::ostream & NOX::NLN::INNER::StatusTest::Filter::Print(
    std::ostream& stream, int indent ) const
{
  return stream;
}
