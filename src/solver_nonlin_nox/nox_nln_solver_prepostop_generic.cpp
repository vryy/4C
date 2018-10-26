/*-----------------------------------------------------------*/
/*!
\file nox_nln_solver_prepostop_generic.cpp

\maintainer Michael Hiermeier

\date Jul 29, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_solver_prepostop_generic.H"  // class definition
#include "nox_nln_solver_linesearchbased.H"
#include "nox_nln_statustest_normf.H"
#include "nox_nln_aux.H"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Solver::PrePostOp::Generic::Generic()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PrePostOp::Generic::runPreIterate(const NOX::Solver::Generic& nlnSolver)
{
  // ToDo Use the getListPtr() routine of the NOX::Solver::Generic class
  // and do a Teuchos::rcp_const_cast() instead!
  const NOX::NLN::Solver::LineSearchBased* lsSolver =
      dynamic_cast<const NOX::NLN::Solver::LineSearchBased*>(&nlnSolver);

  if (lsSolver == 0) dserror("runPreItertate - non-linear solver cast failed!");

  // Set the current number of nonlinear iterations
  // this is necessary for the linear solver in some cases (e.g. contact)
  const Teuchos::RCP<Teuchos::ParameterList>& params = lsSolver->GetMutableListPtr();
  const std::string dir_method_str(NOX::NLN::AUX::GetDirectionMethodListName(*params));
  if (params->sublist("Direction").isSublist(dir_method_str))
    if (params->sublist("Direction").sublist(dir_method_str).isSublist("Linear Solver"))
    {
      Teuchos::ParameterList& linearSolverParams =
          params->sublist("Direction").sublist(dir_method_str).sublist("Linear Solver");
      linearSolverParams.set<int>("Number of Nonlinear Iterations", lsSolver->getNumIterations());
    }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PrePostOp::Generic::runPreSolve(const NOX::Solver::Generic& nlnSolver)
{
  // ToDo Use the getListPtr() routine of the NOX::Solver::Generic class
  // and do a Teuchos::rcp_const_cast() instead!
  const NOX::NLN::Solver::LineSearchBased* lsSolver =
      dynamic_cast<const NOX::NLN::Solver::LineSearchBased*>(&nlnSolver);

  if (lsSolver == 0) dserror("runPreSolve - non-linear solver cast failed!");

  // set the wanted tolerance for the linear solver
  const Teuchos::RCP<Teuchos::ParameterList>& params = lsSolver->GetMutableListPtr();
  const std::string dir_method_str(NOX::NLN::AUX::GetDirectionMethodListName(*params));
  if (params->sublist("Direction").isSublist(dir_method_str))
    if (params->sublist("Direction").sublist(dir_method_str).isSublist("Linear Solver"))
    {
      Teuchos::ParameterList& linearSolverParams =
          params->sublist("Direction").sublist(dir_method_str).sublist("Linear Solver");

      // Find and get the "specified tolerance" of the structural normF test in the statusTest
      // object
      const NOX::StatusTest::Generic& statusTest = lsSolver->GetOuterStatusTest();
      double wanted = NOX::NLN::AUX::GetNormFClassVariable(
          statusTest, NOX::NLN::StatusTest::quantity_structure, "SpecifiedTolerance");
      if (wanted == -1.0)
        wanted = NOX::NLN::AUX::GetNormFClassVariable(
            statusTest, NOX::NLN::StatusTest::quantity_levelset_reinit, "SpecifiedTolerance");

      if (wanted == -1.0)
      {
        if (lsSolver->GetUtils().isPrintType(NOX::Utils::Warning))
        {
          lsSolver->GetUtils().out()
              << "\n*** WARNING ***\n"
              << "There is no NOX::NLN::StatusTest::NormF test for the primal field \n"
              << "components. The \"Wanted Tolerance\" for the sublist\n"
              << "\"Linear Solver\" was set to its default value 1.0e-6!\n";
        }
        linearSolverParams.set<double>("Wanted Tolerance", 1.0e-6);
      }
      else
        linearSolverParams.set<double>("Wanted Tolerance", wanted);
    }
  return;
}
