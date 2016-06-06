/*----------------------------------------------------------------------*/
/*!
\file nox_nln_contact_linearsystem.cpp

\brief Derived class which manages the special requirements to the linear
       solver for contact problems.

\level 3

\maintainer Michael Hiermeier

\date Jul 14, 2015

*/
/*----------------------------------------------------------------------*/
#include "nox_nln_contact_linearsystem.H"                    // base class

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_blocksparsematrix.H"

#include "../solver_nonlin_nox/nox_nln_interface_jacobian.H"
#include "../solver_nonlin_nox/nox_nln_interface_required.H"

#include "../drt_mortar/mortar_strategy_base.H"
#include "../drt_contact/contact_abstract_strategy.H"

#include "../drt_inpar/inpar_contact.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::CONTACT::LinearSystem::LinearSystem(
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
NOX::NLN::CONTACT::LinearSystem::LinearSystem(
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
    // TODO: maps for merged meshtying and contact problem !!!
    // feed Aztec or Belos based solvers with contact information
    if (solverPtr->Params().isSublist("Aztec Parameters")
        or solverPtr->Params().isSublist("Belos Parameters"))
    {
      if (iConstrPrec_.size()>1)
        dserror("Currently only one constraint preconditioner interface can be handled! \n "
            "Needs to be extended!");

      Teuchos::ParameterList& mueluParams = solverPtr->Params().sublist("Aztec Parameters");
      // vector entries:
      // (0) masterDofMap
      // (1) slaveDofMap
      // (2) innerDofMap
      // (3) activeDofMap
      std::vector<Teuchos::RCP<Epetra_Map> > prec_maps(4,Teuchos::null);
      iConstrPrec_.begin()->second->FillMapsForPreconditioner(prec_maps);
      Teuchos::ParameterList & linSystemProps = mueluParams.sublist("Linear System properties");
      linSystemProps.set<Teuchos::RCP<Epetra_Map> >("contact masterDofMap",prec_maps[0]);
      linSystemProps.set<Teuchos::RCP<Epetra_Map> >("contact slaveDofMap",prec_maps[1]);
      linSystemProps.set<Teuchos::RCP<Epetra_Map> >("contact innerDofMap",prec_maps[2]);
      linSystemProps.set<Teuchos::RCP<Epetra_Map> >("contact activeDofMap",prec_maps[3]);
      // contact or contact/meshtying
      if (iConstrPrec_.begin()->first == NOX::NLN::sol_contact)
        linSystemProps.set<std::string>("ProblemType", "contact");
      // only meshtying
      else if (iConstrPrec_.begin()->first == NOX::NLN::sol_meshtying)
        linSystemProps.set<std::string>("ProblemType", "meshtying");
      else
        dserror("Currently we support only a pure meshtying OR a pure contact problem!");

      linSystemProps.set<int>("time step",step);
      // increase counter by one (historical reasons)
      linSystemProps.set<int>("iter",nlnIter+1);
    }
  } // end: feed solver with contact/meshtying information

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::SolutionType NOX::NLN::CONTACT::LinearSystem::GetActiveLinSolver(
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
    currSolver =  solvers.at(NOX::NLN::sol_contact);
    return NOX::NLN::sol_contact;
  }
  // ----------------------------------------------------------------------
  // check if contact contributions are present,
  // if not we make a standard solver call to speed things up
  // ----------------------------------------------------------------------
  if (!issaddlepoint and !iscondensed)
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
