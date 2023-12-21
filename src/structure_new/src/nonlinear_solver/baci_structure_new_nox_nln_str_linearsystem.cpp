/*-----------------------------------------------------------*/
/*! \file

\brief Derived class which manages the special requirements to the linear
       solver for structural problems.


\level 3

*/
/*-----------------------------------------------------------*/

#include "baci_structure_new_nox_nln_str_linearsystem.H"

#include "baci_linalg_sparseoperator.H"
#include "baci_linear_solver_method_linalg.H"
#include "baci_solver_nonlin_nox_interface_jacobian.H"
#include "baci_solver_nonlin_nox_interface_required.H"

#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::STR::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>& solvers,
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& J,
    const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& M, const ::NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<::NOX::Epetra::Scaling> scalingObject)
    : NOX::NLN::LinearSystem(printParams, linearSolverParams, solvers, iReq, iJac, J, iPrec, M,
          cloneVector, scalingObject)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::STR::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>& solvers,
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& J,
    const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& M, const ::NOX::Epetra::Vector& cloneVector)
    : NOX::NLN::LinearSystem(
          printParams, linearSolverParams, solvers, iReq, iJac, J, iPrec, M, cloneVector)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::STR::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>& solvers,
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& J, const ::NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<::NOX::Epetra::Scaling> scalingObject)
    : NOX::NLN::LinearSystem(
          printParams, linearSolverParams, solvers, iReq, iJac, J, cloneVector, scalingObject)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::STR::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>& solvers,
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& J, const ::NOX::Epetra::Vector& cloneVector)
    : NOX::NLN::LinearSystem(printParams, linearSolverParams, solvers, iReq, iJac, J, cloneVector)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::SolverParams NOX::NLN::STR::LinearSystem::SetSolverOptions(Teuchos::ParameterList& p,
    Teuchos::RCP<CORE::LINALG::Solver>& solverPtr, const NOX::NLN::SolutionType& solverType)
{
  CORE::LINALG::SolverParams solver_params;
  bool isAdaptiveControl = p.get<bool>("Adaptive Control");

  if (isAdaptiveControl)
  {
    // dynamic cast of the required/rhs interface
    Teuchos::RCP<NOX::NLN::Interface::Required> iNlnReq =
        Teuchos::rcp_dynamic_cast<NOX::NLN::Interface::Required>(reqInterfacePtr_, true);

    solver_params.nonlin_tolerance = p.get<double>("Wanted Tolerance");
    solver_params.nonlin_residual = iNlnReq->CalcRefNormForce();
    solver_params.lin_tol_better = p.get<double>("Adaptive Control Objective");
  }

  return solver_params;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::SolutionType NOX::NLN::STR::LinearSystem::GetActiveLinSolver(
    const std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>& solvers,
    Teuchos::RCP<CORE::LINALG::Solver>& currSolver)
{
  // check input
  if (solvers.size() > 1) dserror("There has to be exactly one CORE::LINALG::Solver (structure)!");

  currSolver = solvers.at(NOX::NLN::sol_structure);
  return NOX::NLN::sol_structure;
}

BACI_NAMESPACE_CLOSE
