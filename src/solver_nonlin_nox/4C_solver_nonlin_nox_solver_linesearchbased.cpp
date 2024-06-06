/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_solver_nonlin_nox_solver_linesearchbased.hpp"  // class definition

#include "4C_solver_nonlin_nox_aux.hpp"
#include "4C_solver_nonlin_nox_direction_factory.hpp"
#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_solver_nonlin_nox_inner_statustest_generic.hpp"
#include "4C_solver_nonlin_nox_linesearch_factory.hpp"
#include "4C_solver_nonlin_nox_statustest_activeset.hpp"
#include "4C_solver_nonlin_nox_statustest_normf.hpp"
#include "4C_solver_nonlin_nox_statustest_normupdate.hpp"
#include "4C_solver_nonlin_nox_statustest_normwrms.hpp"

#include <NOX_Abstract_Group.H>
#include <NOX_Direction_Factory.H>
#include <NOX_Direction_Generic.H>
#include <NOX_LineSearch_Generic.H>
#include <NOX_Solver_SolverUtils.H>
#include <NOX_Utils.H>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::Solver::LineSearchBased::LineSearchBased(const Teuchos::RCP<::NOX::Abstract::Group>& grp,
    const Teuchos::RCP<::NOX::StatusTest::Generic>& outerTests,
    const Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>& innerTests,
    const Teuchos::RCP<Teuchos::ParameterList>& params)
    : ::NOX::Solver::LineSearchBased(grp, outerTests, params)
{
  // call derived init() after base init() was called.
  init(innerTests);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Solver::LineSearchBased::init(
    const Teuchos::RCP<NOX::Nln::Inner::StatusTest::Generic>& innerTests)
{
  stepSize = 0.0;
  nIter = 0;
  status = ::NOX::StatusTest::Unconverged;
  checkType = ::NOX::Solver::parseStatusTestCheckType(paramsPtr->sublist("Solver Options"));

  // NOTE: We use different factories at this point!
  lineSearchPtr = NOX::Nln::LineSearch::BuildLineSearch(
      globalDataPtr, testPtr, innerTests, paramsPtr->sublist("Line Search"));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::Solver::LineSearchBased::printUpdate()
{
  double normSoln = 0;
  double normStep = 0;

  // Print the status test parameters at each iteration if requested
  if ((status == ::NOX::StatusTest::Unconverged) and
      (utilsPtr->isPrintType(::NOX::Utils::OuterIterationStatusTest)))
  {
    utilsPtr->out() << ::NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Status Test Results --\n";
    testPtr->print(utilsPtr->out());
    utilsPtr->out() << ::NOX::Utils::fill(72) << "\n";
  }

  // All processes participate in the computation of these norms...
  if (utilsPtr->isPrintType(::NOX::Utils::OuterIteration))
  {
    normSoln = solnPtr->getNormF();
    normStep = (nIter > 0) ? dirPtr->norm() : 0;
  }

  // ...But only the print process actually prints the result.
  // ------ standard output ------------------------------------------
  if (utilsPtr->isPrintType(::NOX::Utils::OuterIteration) and
      utilsPtr->isPrintType(::NOX::Utils::OuterIterationStatusTest))
  {
    utilsPtr->out() << "\n" << ::NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Nonlinear Solver Step " << nIter << " -- \n";
    utilsPtr->out() << "||F|| = " << utilsPtr->sciformat(normSoln);
    utilsPtr->out() << "  step = " << utilsPtr->sciformat(stepSize);
    utilsPtr->out() << "  dx = " << utilsPtr->sciformat(normStep);
    if (status == ::NOX::StatusTest::Converged) utilsPtr->out() << " (Converged!)";
    if (status == ::NOX::StatusTest::Failed) utilsPtr->out() << " (Failed!)";
    utilsPtr->out() << "\n" << ::NOX::Utils::fill(72) << "\n" << std::endl;
  }
  // ------ short output ---------------------------------------------
  else if (utilsPtr->isPrintType(::NOX::Utils::OuterIteration))
  {
    // print the head line
    if (nIter == 0)
    {
      utilsPtr->out() << std::setw(4) << "#It" << std::setw(13) << "||F||_2" << std::setw(13)
                      << "step" << std::setw(13) << "||dx||_2\n";
      utilsPtr->out() << ::NOX::Utils::fill(50, '^') << "\n";
    }
    utilsPtr->out() << std::setw(4) << nIter;
    utilsPtr->out() << "  " << utilsPtr->sciformat(normSoln);
    utilsPtr->out() << "  " << utilsPtr->sciformat(stepSize);
    utilsPtr->out() << "  " << utilsPtr->sciformat(normStep);
    if (status == ::NOX::StatusTest::Converged) utilsPtr->out() << " (Converged!)";
    if (status == ::NOX::StatusTest::Failed) utilsPtr->out() << " (Failed!)";
    utilsPtr->out() << std::endl;
  }

  // Print the final parameter values of the status test
  if ((status != ::NOX::StatusTest::Unconverged) and
      (utilsPtr->isPrintType(::NOX::Utils::OuterIteration) or
          utilsPtr->isPrintType(::NOX::Utils::OuterIterationStatusTest)))
  {
    utilsPtr->out() << ::NOX::Utils::fill(72) << "\n";
    utilsPtr->out() << "-- Final Status Test Results --\n";
    testPtr->print(utilsPtr->out());
    utilsPtr->out() << ::NOX::Utils::fill(72) << "\n";
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const ::NOX::StatusTest::Generic& NOX::Nln::Solver::LineSearchBased::GetOuterStatusTest() const
{
  if (testPtr.is_null())
  {
    utilsPtr->err() << "The \"Status Test\" pointer is not initialized!" << std::endl;
    throw "NOX Error";
  }

  return (*testPtr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class T>
::NOX::StatusTest::Generic* NOX::Nln::Solver::LineSearchBased::GetOuterStatusTest() const
{
  if (testPtr.is_null())
  {
    utilsPtr->err() << "The \"Status Test\" pointer is not initialized!" << std::endl;
    throw "NOX Error";
  }

  return NOX::Nln::Aux::GetOuterStatusTest<T>(*testPtr);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class T>
::NOX::StatusTest::Generic* NOX::Nln::Solver::LineSearchBased::get_outer_status_test_with_quantity(
    const NOX::Nln::StatusTest::QuantityType qtype) const
{
  if (testPtr.is_null())
  {
    utilsPtr->err() << "The \"Status Test\" pointer is not initialized!" << std::endl;
    throw "NOX Error";
  }

  return NOX::Nln::Aux::get_outer_status_test_with_quantity<T>(*testPtr, qtype);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::StatusTest::StatusType NOX::Nln::Solver::LineSearchBased::getStatus() const
{
  return status;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <class T>
::NOX::StatusTest::StatusType NOX::Nln::Solver::LineSearchBased::GetStatus() const
{
  ::NOX::StatusTest::StatusType gstatus = ::NOX::StatusTest::Unevaluated;
  int status = NOX::Nln::Aux::GetOuterStatus<T>(*testPtr);
  if (status != -100) gstatus = static_cast<::NOX::StatusTest::StatusType>(status);

  return gstatus;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const ::NOX::Utils& NOX::Nln::Solver::LineSearchBased::GetUtils() const { return *utilsPtr; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Direction::Generic& NOX::Nln::Solver::LineSearchBased::GetDirection() const
{
  if (directionPtr.is_null())
    FOUR_C_THROW("nullptr ptr: The direction pointer is not yet initialized.");

  return *directionPtr;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template ::NOX::StatusTest::StatusType
NOX::Nln::Solver::LineSearchBased::GetStatus<NOX::Nln::StatusTest::NormF>() const;
template ::NOX::StatusTest::StatusType
NOX::Nln::Solver::LineSearchBased::GetStatus<NOX::Nln::StatusTest::NormUpdate>() const;
template ::NOX::StatusTest::StatusType
NOX::Nln::Solver::LineSearchBased::GetStatus<NOX::Nln::StatusTest::NormWRMS>() const;
template ::NOX::StatusTest::StatusType
NOX::Nln::Solver::LineSearchBased::GetStatus<NOX::Nln::StatusTest::ActiveSet>() const;

template ::NOX::StatusTest::Generic*
NOX::Nln::Solver::LineSearchBased::GetOuterStatusTest<NOX::Nln::StatusTest::ActiveSet>() const;
template ::NOX::StatusTest::Generic*
NOX::Nln::Solver::LineSearchBased::get_outer_status_test_with_quantity<NOX::Nln::StatusTest::NormF>(
    const NOX::Nln::StatusTest::QuantityType qtype) const;

FOUR_C_NAMESPACE_CLOSE
