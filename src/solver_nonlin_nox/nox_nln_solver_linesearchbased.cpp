/*-----------------------------------------------------------*/
/*!
\file nox_nln_solver_linesearchbased.cpp

\maintainer Michael Hiermeier

\date Jun 11, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_solver_linesearchbased.H"   // class definition
#include "nox_nln_inner_statustest_generic.H"
#include "nox_nln_linesearch_factory.H"
#include "nox_nln_direction_factory.H"
#include "nox_nln_group.H"

#include <NOX_StatusTest_Generic.H>
#include <NOX_Solver_SolverUtils.H>
#include <NOX_Direction_Factory.H>
#include <NOX_Direction_Generic.H>
#include <NOX_LineSearch_Generic.H>
#include <NOX_Utils.H>
#include <NOX_Abstract_Group.H>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Solver::LineSearchBased::LineSearchBased
   (const Teuchos::RCP<NOX::Abstract::Group>& grp,
    const Teuchos::RCP<NOX::StatusTest::Generic>& outerTests,
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests,
    const Teuchos::RCP<Teuchos::ParameterList>& params) :
    NOX::Solver::LineSearchBased(grp,outerTests,params)
{
  // call derived init() after base init() was called.
  init(innerTests);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::LineSearchBased::init(
    const Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests)
{
  // NOTE: We use different factories at this point!
  lineSearchPtr = NOX::NLN::LineSearch::
      BuildLineSearch(globalDataPtr,innerTests,paramsPtr->sublist("Line Search"));

  directionPtr = NOX::NLN::Direction::
    BuildDirection(globalDataPtr, paramsPtr->sublist("Direction"));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::LineSearchBased::printUpdate()
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

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::StatusTest::Generic& NOX::NLN::Solver::LineSearchBased::GetOuterStatusTest() const
{
  if (testPtr.is_null())
  {
    utilsPtr->err()
    << "The \"Status Test\" pointer is not initialized!" << std::endl;
    throw "NOX Error";
  }

  return (*testPtr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::Utils& NOX::NLN::Solver::LineSearchBased::GetUtils() const
{
  return *utilsPtr;
}

