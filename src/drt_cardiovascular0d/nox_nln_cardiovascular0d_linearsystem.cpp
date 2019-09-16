/*-----------------------------------------------------------*/
/*! \file

\brief Derived class which manages the special requirements to the linear
       solver for 0D cardiovascular-structural problems.

\maintainer Amadeus Gebauer

\date Jul 15, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_cardiovascular0d_linearsystem.H"

#include "../solver_nonlin_nox/nox_nln_interface_jacobian.H"
#include "../solver_nonlin_nox/nox_nln_interface_required.H"

#include <Teuchos_ParameterList.hpp>

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_sparseoperator.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::CARDIOVASCULAR0D::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& solvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<LINALG::SparseOperator>& J,
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
    const Teuchos::RCP<LINALG::SparseOperator>& M, const NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject)
    : NOX::NLN::LinearSystem(printParams, linearSolverParams, solvers, iReq, iJac, J, iPrec, M,
          cloneVector, scalingObject)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::CARDIOVASCULAR0D::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& solvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<LINALG::SparseOperator>& J,
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
    const Teuchos::RCP<LINALG::SparseOperator>& M, const NOX::Epetra::Vector& cloneVector)
    : NOX::NLN::LinearSystem(
          printParams, linearSolverParams, solvers, iReq, iJac, J, iPrec, M, cloneVector)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::CARDIOVASCULAR0D::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& solvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<LINALG::SparseOperator>& J, const NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject)
    : NOX::NLN::LinearSystem(
          printParams, linearSolverParams, solvers, iReq, iJac, J, cloneVector, scalingObject)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::CARDIOVASCULAR0D::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& solvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<LINALG::SparseOperator>& J, const NOX::Epetra::Vector& cloneVector)
    : NOX::NLN::LinearSystem(printParams, linearSolverParams, solvers, iReq, iJac, J, cloneVector)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::CARDIOVASCULAR0D::LinearSystem::SetSolverOptions(Teuchos::ParameterList& p,
    Teuchos::RCP<LINALG::Solver>& solverPtr, const NOX::NLN::SolutionType& solverType)
{
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
    solverPtr->AdaptTolerance(wanted, worst, adaptiveControlObjective);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::SolutionType NOX::NLN::CARDIOVASCULAR0D::LinearSystem::GetActiveLinSolver(
    const std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& solvers,
    Teuchos::RCP<LINALG::Solver>& currSolver)
{
  //  // check input
  //  if (solvers.size()>1)
  //    throwError("GetCurrentLinSolver","There has to be exactly one "
  //        "LINALG::Solver (structure)!");

  currSolver = solvers.at(NOX::NLN::sol_cardiovascular0d);
  return NOX::NLN::sol_cardiovascular0d;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::CARDIOVASCULAR0D::LinearSystem::throwError(
    const std::string& functionName, const std::string& errorMsg) const
{
  if (utils_.isPrintType(NOX::Utils::Error))
  {
    utils_.out() << "NOX::NLN::CARDIOVASCULAR0D::LinearSystem::" << functionName << " - "
                 << errorMsg << std::endl;
  }
  throw "NOX Error";
}
