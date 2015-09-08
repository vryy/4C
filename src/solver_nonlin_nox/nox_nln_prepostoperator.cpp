/*-----------------------------------------------------------*/
/*!
\file nox_nln_prepostoperator.cpp

\maintainer Michael Hiermeier

\date Jul 29, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_prepostoperator.H"    // class definition
#include "nox_nln_solver_linesearchbased.H"
#include "nox_nln_statustest_normf.H"
#include "nox_nln_aux.H"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::PrePostOperator::PrePostOperator()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::PrePostOperator::runPreIterate(const NOX::Solver::Generic& nlnSolver)
{
  const NOX::NLN::Solver::LineSearchBased* lsSolver =
      dynamic_cast<const NOX::NLN::Solver::LineSearchBased*>(&nlnSolver);

  if (lsSolver == 0)
    dserror("runPreItertate - non-linear solver cast failed!");

  // set the current number of nonlinear iterations
  // this is necessary for the linear solver in some cases (e.g. contact)
  const Teuchos::RCP<Teuchos::ParameterList>& params = lsSolver->GetMutableListPtr();
  Teuchos::ParameterList& linearSolverParams = params->sublist("Direction",true).
      sublist("Newton",true).sublist("Linear Solver",true);
  linearSolverParams.set<int>("Number of Nonlinear Iterations",lsSolver->getNumIterations());


  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::PrePostOperator::runPreSolve(const NOX::Solver::Generic& nlnSolver)
{
  const NOX::NLN::Solver::LineSearchBased* lsSolver =
      dynamic_cast<const NOX::NLN::Solver::LineSearchBased*>(&nlnSolver);

  if (lsSolver == 0)
    dserror("runPreSolve - non-linear solver cast failed!");

  // set the wanted tolerance for the linear solver
  const Teuchos::RCP<Teuchos::ParameterList>& params = lsSolver->GetMutableListPtr();
  Teuchos::ParameterList& linearSolverParams = params->sublist("Direction",true).
      sublist("Newton",true).sublist("Linear Solver",true);

  // Find and get the "specified tolerance" of the structural normF test in the statusTest object
  const NOX::StatusTest::Generic& statusTest = lsSolver->GetOuterStatusTest();
  double wanted = NOX::NLN::AUX::GetNormFClassVariable(
      statusTest,NOX::NLN::StatusTest::quantity_structure,"SpecifiedTolerance");

  if (wanted==-1.0)
  {
    lsSolver->GetUtils().err() << "\n*** WARNING ***\n"
                               << "There is no NOX::NLN::StatusTest::NormF test for the structural\n"
                               << "components. The \"Wanted Tolerance\" for the sublist\n"
                               << "\"Linear Solver\" was set to its default value 1.0e-6!\n";

    linearSolverParams.set<double>("Wanted Tolerance",1.0e-6);
  }
  else
    linearSolverParams.set<double>("Wanted Tolerance",wanted);

  return;
}

