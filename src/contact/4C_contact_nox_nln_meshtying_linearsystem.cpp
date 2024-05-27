/*----------------------------------------------------------------------*/
/*! \file
\brief Derived class which manages the special requirements to the linear
       solver for mesh tying problems.

\level 3


*/
/*----------------------------------------------------------------------*/
#include "4C_contact_nox_nln_meshtying_linearsystem.hpp"  // base class

#include "4C_contact_abstract_strategy.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_mortar_strategy_base.hpp"
#include "4C_solver_nonlin_nox_aux.hpp"
#include "4C_solver_nonlin_nox_interface_jacobian.hpp"
#include "4C_solver_nonlin_nox_interface_required.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::MESHTYING::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
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
NOX::NLN::MESHTYING::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
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
CORE::LINALG::SolverParams NOX::NLN::MESHTYING::LinearSystem::set_solver_options(
    Teuchos::ParameterList& p, Teuchos::RCP<CORE::LINALG::Solver>& solverPtr,
    const NOX::NLN::SolutionType& solverType)
{
  CORE::LINALG::SolverParams solver_params;

  bool isAdaptiveControl = p.get<bool>("Adaptive Control");
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
    if (iNlnReq.is_null()) throw_error("setSolverOptions", "required interface cast failed");

    double worst = iNlnReq->CalcRefNormForce();
    // This value has to be specified in the PrePostOperator object of
    // the non-linear solver (i.e. runPreSolve())
    double wanted = p.get<double>("Wanted Tolerance");
    solver_params.nonlin_tolerance = wanted;
    solver_params.nonlin_residual = worst;
    solver_params.lin_tol_better = adaptiveControlObjective;
  }

  // nothing more to do for a pure structural solver
  if (solverType == NOX::NLN::sol_structure) return solver_params;

  // update information about active slave dofs
  // ---------------------------------------------------------------------
  // feed solver/preconditioner with additional information about the
  // contact/meshtying problem
  // ---------------------------------------------------------------------
  {
    // TODO: maps for merged meshtying and contact problem !!!
    // feed Belos based solvers with contact information
    if (solverPtr->Params().isSublist("Belos Parameters"))
    {
      if (i_constr_prec_.size() > 1)
        FOUR_C_THROW(
            "Currently only one constraint preconditioner interface can be handled! \n "
            "Needs to be extended!");

      Teuchos::ParameterList& mueluParams = solverPtr->Params().sublist("Belos Parameters");

      // vector entries:
      // (0) masterDofMap
      // (1) slaveDofMap
      // (2) innerDofMap
      // (3) activeDofMap
      std::vector<Teuchos::RCP<Epetra_Map>> prec_maps(4, Teuchos::null);
      i_constr_prec_.begin()->second->fill_maps_for_preconditioner(prec_maps);
      mueluParams.set<Teuchos::RCP<Epetra_Map>>("contact masterDofMap", prec_maps[0]);
      mueluParams.set<Teuchos::RCP<Epetra_Map>>("contact slaveDofMap", prec_maps[1]);
      mueluParams.set<Teuchos::RCP<Epetra_Map>>("contact innerDofMap", prec_maps[2]);
      mueluParams.set<Teuchos::RCP<Epetra_Map>>("contact activeDofMap", prec_maps[3]);
      // contact or contact/meshtying
      if (i_constr_prec_.begin()->first == NOX::NLN::sol_contact)
        mueluParams.set<std::string>("GLOBAL::ProblemType", "contact");
      // only meshtying
      else if (i_constr_prec_.begin()->first == NOX::NLN::sol_meshtying)
        mueluParams.set<std::string>("GLOBAL::ProblemType", "meshtying");
      else
        FOUR_C_THROW("Currently we support only a pure meshtying OR a pure contact problem!");

      mueluParams.set<int>("time step", step);
      // increase counter by one (historical reasons)
      mueluParams.set<int>("iter", nlnIter + 1);
    }
  }  // end: feed solver with contact/meshtying information

  return solver_params;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::SolutionType NOX::NLN::MESHTYING::LinearSystem::get_active_lin_solver(
    const std::map<NOX::NLN::SolutionType, Teuchos::RCP<CORE::LINALG::Solver>>& solvers,
    Teuchos::RCP<CORE::LINALG::Solver>& currSolver)
{
  currSolver = solvers.at(NOX::NLN::sol_meshtying);
  return NOX::NLN::sol_meshtying;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::MESHTYING::LinearSystem::throw_error(
    const std::string& functionName, const std::string& errorMsg) const
{
  if (utils_.isPrintType(::NOX::Utils::Error))
  {
    utils_.out() << "NOX::CONTACT::LinearSystem::" << functionName << " - " << errorMsg
                 << std::endl;
  }
  throw "NOX Error";
}

FOUR_C_NAMESPACE_CLOSE
