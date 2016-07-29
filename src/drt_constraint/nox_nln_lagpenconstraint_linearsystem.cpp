/*-----------------------------------------------------------*/
/*!
\file nox_nln_lagpenconstraint_linearsystem.cpp

\brief Derived class which manages the special requirements to the linear
       solver for structural-constraint problems.

\maintainer Marc Hirschvogel

\date Jul 15, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_lagpenconstraint_linearsystem.H"

#include "../solver_nonlin_nox/nox_nln_interface_jacobian.H"
#include "../solver_nonlin_nox/nox_nln_interface_required.H"

#include <Teuchos_ParameterList.hpp>

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_sparseoperator.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::LAGPENCONSTRAINT::LinearSystem::LinearSystem(
    Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const SolverMap& solvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr,
    const Teuchos::RCP<LINALG::SparseOperator>& J,
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
    const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec,
    const Teuchos::RCP<LINALG::SparseOperator>& M,
    const NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject)
    : NOX::NLN::LinearSystem(printParams,linearSolverParams,solvers,
        iReq,iJac,J,iPrec,M,cloneVector,scalingObject),
      iConstr_(iConstr),
      iConstrPrec_(iConstrPrec)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::LAGPENCONSTRAINT::LinearSystem::LinearSystem(
    Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const SolverMap& solvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iConstr,
    const Teuchos::RCP<LINALG::SparseOperator>& J,
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
    const NOX::NLN::CONSTRAINT::PrecInterfaceMap& iConstrPrec,
    const Teuchos::RCP<LINALG::SparseOperator>& M,
    const NOX::Epetra::Vector& cloneVector)
    : NOX::NLN::LinearSystem(printParams,linearSolverParams,solvers,
        iReq,iJac,J,iPrec,M,cloneVector),
      iConstr_(iConstr),
      iConstrPrec_(iConstrPrec)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LAGPENCONSTRAINT::LinearSystem::SetSolverOptions(
    Teuchos::ParameterList& p,
    Teuchos::RCP<LINALG::Solver>& solverPtr,
    const NOX::NLN::SolutionType& solverType)
{

  bool isAdaptiveControl          = p.get<bool>("Adaptive Control");
  double adaptiveControlObjective = p.get<double>("Adaptive Control Objective");
  // This value is specified in the underlying time integrator
  // (i.e. RunPreNoxNlnSolve())
//  int step = p.get<int>("Current Time Step");
  // This value is specified in the PrePostOperator object of
  // the non-linear solver (i.e. runPreIterate())
//  int nlnIter = p.get<int>("Number of Nonlinear Iterations");

  if (isAdaptiveControl)
  {
    // dynamic cast of the required/rhs interface
    Teuchos::RCP<NOX::NLN::Interface::Required> iNlnReq =
        Teuchos::rcp_dynamic_cast<NOX::NLN::Interface::Required>(reqInterfacePtr_);
    if (iNlnReq.is_null())
      throwError("setSolverOptions","required interface cast failed");

    double worst  = iNlnReq->CalcRefNormForce();
    // This value has to be specified in the PrePostOperator object of
    // the non-linear solver (i.e. runPreSolve())
    double wanted = p.get<double>("Wanted Tolerance");
    solverPtr->AdaptTolerance(wanted,worst,adaptiveControlObjective);
  }

  // nothing more to do for a pure structural solver
  if (solverType == NOX::NLN::sol_structure)
    return;


  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::SolutionType NOX::NLN::LAGPENCONSTRAINT::LinearSystem::GetActiveLinSolver(
    const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& solvers,
    Teuchos::RCP<LINALG::Solver>& currSolver)
{

  // check input
  if (solvers.size()>2)
    throwError("GetCurrentLinSolver",
        "There have to be exactly two LINALG::Solvers (structure + contact)!");
  // ---------------------------------------------------------------------
  // Solving a saddle point system
  // (1) Standard / Dual Lagrange multipliers -> SaddlePoint
  // (2) Direct Augmented Lagrange strategy
  // ---------------------------------------------------------------------
  NOX::NLN::CONSTRAINT::PrecInterfaceMap::const_iterator cit;
  bool issaddlepoint = false;
  for (cit=iConstrPrec_.begin();cit!=iConstrPrec_.end();++cit)
  {
    if (cit->second->IsSaddlePointSystem())
    {
      issaddlepoint = true;
      break;
    }
  }
  // ---------------------------------------------------------------------
  // Solving a purely displacement based system
  // (1) Dual (not Standard) Lagrange multipliers -> Condensed
  // (2) Penalty and Uzawa Augmented Lagrange strategies
  // ---------------------------------------------------------------------
  bool iscondensed = false;
  for (cit=iConstrPrec_.begin();cit!=iConstrPrec_.end();++cit)
  {
    if (cit->second->IsCondensedSystem())
    {
      iscondensed = true;
      break;
    }
  }

  if (issaddlepoint or iscondensed)
  {
    currSolver =  solvers.at(NOX::NLN::sol_lag_pen_constraint);
    return NOX::NLN::sol_lag_pen_constraint;
  }
  // ----------------------------------------------------------------------
  // check if constraint contributions are present,
  // if not we make a standard solver call to speed things up
  // ----------------------------------------------------------------------
  if (!issaddlepoint and !iscondensed)
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
    const std::string& functionName,
    const std::string& errorMsg) const
{
  if (utils_.isPrintType(NOX::Utils::Error)) {
    utils_.out() << "NOX::LAGPENCONSTRAINT::LinearSystem::" << functionName
     << " - " << errorMsg << std::endl;
  }
  throw "NOX Error";
}
