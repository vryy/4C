/*-----------------------------------------------------------*/
/*! \file



\level 3

*/
/*-----------------------------------------------------------*/

#include "baci_solver_nonlin_nox_solver_prepostop_generic.hpp"  // class definition

#include "baci_solver_nonlin_nox_aux.hpp"
#include "baci_solver_nonlin_nox_solver_linesearchbased.hpp"
#include "baci_solver_nonlin_nox_statustest_normf.hpp"
#include "baci_utils_exceptions.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Solver::PrePostOp::Generic::Generic()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PrePostOp::Generic::runPreIterate(const ::NOX::Solver::Generic& nlnSolver)
{
  // ToDo Use the getListPtr() routine of the ::NOX::Solver::Generic class
  // and do a Teuchos::rcp_const_cast() instead!
  const NOX::NLN::Solver::LineSearchBased* lsSolver =
      dynamic_cast<const NOX::NLN::Solver::LineSearchBased*>(&nlnSolver);

  if (lsSolver != nullptr)
  {
    // Set the current number of nonlinear iterations
    // this is necessary for the linear solver in some cases (e.g. contact)
    const Teuchos::RCP<Teuchos::ParameterList>& params = lsSolver->GetListPtr();
    const std::string dir_method_str(NOX::NLN::AUX::GetDirectionMethodListName(*params));
    if (params->sublist("Direction").isSublist(dir_method_str))
    {
      if (params->sublist("Direction").sublist(dir_method_str).isSublist("Linear Solver"))
      {
        Teuchos::ParameterList& linearSolverParams =
            params->sublist("Direction").sublist(dir_method_str).sublist("Linear Solver");
        linearSolverParams.set<int>("Number of Nonlinear Iterations", lsSolver->getNumIterations());
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Solver::PrePostOp::Generic::runPreSolve(const ::NOX::Solver::Generic& nlnSolver)
{
  // ToDo Use the getListPtr() routine of the ::NOX::Solver::Generic class
  // and do a Teuchos::rcp_const_cast() instead!
  const NOX::NLN::Solver::LineSearchBased* lsSolver =
      dynamic_cast<const NOX::NLN::Solver::LineSearchBased*>(&nlnSolver);

  if (lsSolver != nullptr)
  {
    // set the wanted tolerance for the linear solver
    const Teuchos::RCP<Teuchos::ParameterList>& params = lsSolver->GetListPtr();
    const std::string dir_method_str(NOX::NLN::AUX::GetDirectionMethodListName(*params));
    if (params->sublist("Direction").isSublist(dir_method_str))
    {
      if (params->sublist("Direction").sublist(dir_method_str).isSublist("Linear Solver"))
      {
        Teuchos::ParameterList& linearSolverParams =
            params->sublist("Direction").sublist(dir_method_str).sublist("Linear Solver");

        // Find and get the "specified tolerance" of the structural normF test in the statusTest
        // object
        const ::NOX::StatusTest::Generic& statusTest = lsSolver->GetOuterStatusTest();
        double specified_tol = NOX::NLN::AUX::GetNormFClassVariable(
            statusTest, NOX::NLN::StatusTest::quantity_structure, "SpecifiedTolerance");
        if (specified_tol == -1.0)
          specified_tol = NOX::NLN::AUX::GetNormFClassVariable(
              statusTest, NOX::NLN::StatusTest::quantity_levelset_reinit, "SpecifiedTolerance");

        if (specified_tol == -1.0)
        {
          if (lsSolver->GetUtils().isPrintType(::NOX::Utils::Warning))
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
          linearSolverParams.set<double>("Wanted Tolerance", specified_tol);
      }
    }
  }
}

BACI_NAMESPACE_CLOSE
