/*-----------------------------------------------------------*/
/*!
\file nox_nln_contact_linearsystem.cpp

\maintainer Michael Hiermeier

\date Jul 14, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_contact_linearsystem.H"     // base class
#include "nox_nln_contact_interface_preconditioner.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_sparseoperator.H"

#include "../solver_nonlin_nox/nox_nln_interface_jacobian.H"
#include "../solver_nonlin_nox/nox_nln_interface_required.H"

#include "../drt_mortar/mortar_strategy_base.H"
#include "../drt_contact/contact_abstract_strategy.H"

#include "../drt_inpar/inpar_contact.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::CONTACT::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& solvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<LINALG::SparseOperator>& J,
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
    const Teuchos::RCP<LINALG::SparseOperator>& M,
    const NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject)
    : NOX::NLN::LinearSystem(printParams,linearSolverParams,solvers,iReq,iJac,J,iPrec,M,cloneVector,scalingObject),
      strat_(Teuchos::null),
      isMeshtying_(Teuchos::null),
      isContact_(Teuchos::null)
{
  Init();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::CONTACT::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& solvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<LINALG::SparseOperator>& J,
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
    const Teuchos::RCP<LINALG::SparseOperator>& M,
    const NOX::Epetra::Vector& cloneVector)
    : NOX::NLN::LinearSystem(printParams,linearSolverParams,solvers,iReq,iJac,J,iPrec,M,cloneVector),
      strat_(Teuchos::null),
      isMeshtying_(Teuchos::null),
      isContact_(Teuchos::null)
{
  Init();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::CONTACT::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& solvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<LINALG::SparseOperator>& J,
    const NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject)
    : NOX::NLN::LinearSystem(printParams,linearSolverParams,solvers,iReq,iJac,J,cloneVector,scalingObject),
      strat_(Teuchos::null),
      isMeshtying_(Teuchos::null),
      isContact_(Teuchos::null)
{
  Init();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::CONTACT::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& solvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<LINALG::SparseOperator>& J,
    const NOX::Epetra::Vector& cloneVector)
    : NOX::NLN::LinearSystem(printParams,linearSolverParams,solvers,iReq,iJac,J,cloneVector),
      strat_(Teuchos::null),
      isMeshtying_(Teuchos::null),
      isContact_(Teuchos::null)
{
  Init();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::CONTACT::LinearSystem::Init()
{
  if (precInterfacePtr_.is_null())
    throwError("LinearSystem","the preconditioner interface is not initialized!");

  Teuchos::RCP<NOX::NLN::CONTACT::Interface::Preconditioner> iCoPrec =
      Teuchos::rcp_dynamic_cast<NOX::NLN::CONTACT::Interface::Preconditioner>(precInterfacePtr_);

  if (iCoPrec.is_null())
    throwError("NOX::CONTACT::LinearSystem::LinearSystem",
        "NOX::CONTACT::Interface::Preconditioner was not initialized!");

  // initialize booleans
  isMeshtying_ = iCoPrec->HaveMeshtying();
  isContact_   = iCoPrec->HaveContact();

  // initialize strategy pointer
  strat_ = Teuchos::rcpFromRef(iCoPrec->GetStrategy());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::CONTACT::LinearSystem::SetSolverOptions(
    Teuchos::ParameterList& p,
    Teuchos::RCP<LINALG::Solver>& solverPtr,
    const NOX::NLN::SolutionType& solverType)
{
  bool isAdaptiveControl          = p.get<bool>("Adaptive Control");
  double adaptiveControlObjective = p.get<double>("Adaptive Control Objective");
  // This value is specified in the underlying time integrator
  // (i.e. RunPreNoxNlnSolve())
  int step = p.get<int>("Current Time Step");
  // This value is specified in the PrePostOperator object of
  // the non-linear solver (i.e. runPreIterate())
  int nlnIter = p.get<int>("Number of Nonlinear Iterations");

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

  // update information about active slave dofs
  // ---------------------------------------------------------------------
  // feed solver/preconditioner with additional information about the
  // contact/meshtying problem
  // ---------------------------------------------------------------------
  {
    //TODO: maps for merged meshtying and contact problem !!!

    // feed Aztec or Belos based solvers with contact information
    if (solverPtr->Params().isSublist("Aztec Parameters")
        or solverPtr->Params().isSublist("Belos Parameters"))
    {
      Teuchos::ParameterList& mueluParams = solverPtr->Params().sublist("Aztec Parameters");
      Teuchos::RCP<Epetra_Map> masterDofMap;
      Teuchos::RCP<Epetra_Map> slaveDofMap;
      Teuchos::RCP<Epetra_Map> innerDofMap;
      Teuchos::RCP<Epetra_Map> activeDofMap;
      GetStrategy().CollectMapsForPreconditioner(masterDofMap, slaveDofMap, innerDofMap, activeDofMap);
      Teuchos::ParameterList & linSystemProps = mueluParams.sublist("Linear System properties");
      linSystemProps.set<Teuchos::RCP<Epetra_Map> >("contact masterDofMap",masterDofMap);
      linSystemProps.set<Teuchos::RCP<Epetra_Map> >("contact slaveDofMap",slaveDofMap);
      linSystemProps.set<Teuchos::RCP<Epetra_Map> >("contact innerDofMap",innerDofMap);
      linSystemProps.set<Teuchos::RCP<Epetra_Map> >("contact activeDofMap",activeDofMap);
      Teuchos::RCP< ::CONTACT::CoAbstractStrategy> costrat = Teuchos::rcp_dynamic_cast< ::CONTACT::CoAbstractStrategy>(strat_);
      // contact or contact/meshtying
      if (costrat != Teuchos::null) linSystemProps.set<std::string>("ProblemType", "contact");
      // only meshtying
      else                         linSystemProps.set<std::string>("ProblemType", "meshtying");
      linSystemProps.set<int>("time step",step);
      // increase counter by one (historical reasons)
      linSystemProps.set<int>("iter",nlnIter+1);
    }
  } // end: feed solver with contact/meshtying information

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::SolutionType NOX::NLN::CONTACT::LinearSystem::GetCurrentLinSolver(
    const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& solvers,
    Teuchos::RCP<LINALG::Solver>& currSolver)
{
  // check input
  if (solvers.size()>2)
    throwError("GetCurrentLinSolver","There have to be exactly two LINALG::Solvers (structure + contact)!");

  INPAR::CONTACT::SolvingStrategy soltype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(strat_->Params(),"STRATEGY");
  INPAR::CONTACT::SystemType      systype =
      DRT::INPUT::IntegralValue<INPAR::CONTACT::SystemType>(strat_->Params(),"SYSTEM");

  // ---------------------------------------------------------------------
  // Solving a saddle point system
  // (1) Standard / Dual Lagrange multipliers -> SaddlePoint
  // (2) Direct Augmented Lagrange strategy
  // ---------------------------------------------------------------------
  if ((soltype==INPAR::CONTACT::solution_lagmult || soltype==INPAR::CONTACT::solution_augmented)
      && systype!=INPAR::CONTACT::system_condensed)
  {
    // check if contact contributions are present,
    // if not we make a standard solver call to speed things up
    if (!GetStrategy().IsInContact() && !GetStrategy().WasInContact() && !GetStrategy().WasInContactLastTimeStep())
    {
      currSolver = solvers.at(NOX::NLN::sol_structure);
      return NOX::NLN::sol_structure;
    }
    else
    {
      currSolver =  solvers.at(NOX::NLN::sol_contact);
      return NOX::NLN::sol_contact;
    }
  }
  // ---------------------------------------------------------------------
  // Solving a purely displacement based system
  // (1) Dual (not Standard) Lagrange multipliers -> Condensed
  // (2) Penalty and Uzawa Augmented Lagrange strategies
  // ---------------------------------------------------------------------
  if(isMeshtying_)
  {
    currSolver = solvers.at(NOX::NLN::sol_contact);
    return NOX::NLN::sol_contact;
  }
  else if (isContact_)
    // check if contact contributions are present,
    // if not we make a standard solver call to speed things up
    if (!GetStrategy().IsInContact() and
        !GetStrategy().WasInContact() and
        !GetStrategy().WasInContactLastTimeStep())
    {
      currSolver = solvers.at(NOX::NLN::sol_structure);
      return NOX::NLN::sol_structure;
    }

  // default return
  currSolver = solvers.at(NOX::NLN::sol_contact);
  return NOX::NLN::sol_contact;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MORTAR::StrategyBase& NOX::NLN::CONTACT::LinearSystem::GetStrategy()
{
  return *strat_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::CONTACT::LinearSystem::throwError(
    const std::string& functionName,
    const std::string& errorMsg) const
{
  if (utils_.isPrintType(NOX::Utils::Error)) {
    utils_.out() << "NOX::CONTACT::LinearSystem::" << functionName
     << " - " << errorMsg << std::endl;
  }
  throw "NOX Error";
}
