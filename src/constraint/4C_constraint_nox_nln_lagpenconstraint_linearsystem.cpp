/*-----------------------------------------------------------*/
/*! \file

\brief Derived class which manages the special requirements to the linear
       solver for structural-constraint problems.


\date Jul 15, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_constraint_nox_nln_lagpenconstraint_linearsystem.hpp"

#include "4C_linalg_sparseoperator.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_solver_nonlin_nox_interface_jacobian.hpp"
#include "4C_solver_nonlin_nox_interface_required.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::LAGPENCONSTRAINT::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams, const SolverMap& solvers,
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
    const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& J,
    const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec,
    const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& M, const ::NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<::NOX::Epetra::Scaling> scalingObject)
    : NOX::NLN::LinearSystem(printParams, linearSolverParams, solvers, iReq, iJac, J, iPrec, M,
          cloneVector, scalingObject),
      i_constr_(iConstr),
      i_constr_prec_(iConstrPrec)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::LAGPENCONSTRAINT::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams, const SolverMap& solvers,
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
    const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& J,
    const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec,
    const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& M, const ::NOX::Epetra::Vector& cloneVector)
    : NOX::NLN::LinearSystem(
          printParams, linearSolverParams, solvers, iReq, iJac, J, iPrec, M, cloneVector),
      i_constr_(iConstr),
      i_constr_prec_(iConstrPrec)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::SolverParams NOX::NLN::LAGPENCONSTRAINT::LinearSystem::SetSolverOptions(
    Teuchos::ParameterList& p, Teuchos::RCP<CORE::LINALG::Solver>& solverPtr,
    const NOX::NLN::SolutionType& solverType)
{
  CORE::LINALG::SolverParams solver_params;

  bool isAdaptiveControl = p.get<bool>("Adaptive Control");
  double adaptiveControlObjective = p.get<double>("Adaptive Control Objective");

  if (isAdaptiveControl)
  {
    // dynamic cast of the required/rhs interface
    Teuchos::RCP<NOX::NLN::Interface::Required> iNlnReq =
        Teuchos::rcp_dynamic_cast<NOX::NLN::Interface::Required>(reqInterfacePtr_);
    if (iNlnReq.is_null()) throwError("setSolverOptions", "required interface cast failed");

    double worst = iNlnReq->CalcRefNormForce();
    // This value has to be specified in the PrePostOperator object of
    // the non-linear solver (i.e. runPreSolve())
    double wanted = p.get<double>("Wanted Tolerance");
    solver_params.nonlin_tolerance = wanted;
    solver_params.nonlin_residual = worst;
    solver_params.lin_tol_better = adaptiveControlObjective;
  }
  return solver_params;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::SolutionType NOX::NLN::LAGPENCONSTRAINT::LinearSystem::GetActiveLinSolver(
    const std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>& solvers,
    Teuchos::RCP<CORE::LINALG::Solver>& currSolver)
{
  // check input
  if (solvers.size() > 2)
    throwError("GetCurrentLinSolver",
        "There have to be exactly two CORE::LINALG::Solvers (structure + contact)!");
  // ---------------------------------------------------------------------
  // Solving a saddle point system
  // Lagrange multiplier constraints -> SaddlePoint
  // ---------------------------------------------------------------------
  NOX::NLN::CONSTRAINT::PrecInterfaceMap::const_iterator cit;
  bool issaddlepoint = false;
  for (cit = i_constr_prec_.begin(); cit != i_constr_prec_.end(); ++cit)
  {
    if (cit->second->IsSaddlePointSystem())
    {
      issaddlepoint = true;
      break;
    }
  }

  if (issaddlepoint)
  {
    currSolver = solvers.at(NOX::NLN::sol_lag_pen_constraint);
    return NOX::NLN::sol_lag_pen_constraint;
  }
  // ----------------------------------------------------------------------
  // Solving a standard displacement-based system
  // Penalty-enforced constraints
  // ----------------------------------------------------------------------
  if (!issaddlepoint)
  {
    currSolver = solvers.at(NOX::NLN::sol_structure);
    return NOX::NLN::sol_structure;
  }

  // default return
  currSolver = solvers.at(NOX::NLN::sol_lag_pen_constraint);
  return NOX::NLN::sol_lag_pen_constraint;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LAGPENCONSTRAINT::LinearSystem::throwError(
    const std::string& functionName, const std::string& errorMsg) const
{
  if (utils_.isPrintType(::NOX::Utils::Error))
  {
    utils_.out() << "NOX::LAGPENCONSTRAINT::LinearSystem::" << functionName << " - " << errorMsg
                 << std::endl;
  }
  throw "NOX Error";
}

FOUR_C_NAMESPACE_CLOSE
