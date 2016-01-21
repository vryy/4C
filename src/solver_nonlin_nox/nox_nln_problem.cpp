/*-----------------------------------------------------------*/
/*!
\file nox_nln_problem.cpp

\maintainer Michael Hiermeier

\date Jun 30, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_PROBLEM_CPP_
#define NOX_NLN_PROBLEM_CPP_

#include "nox_nln_problem.H"    // class definition
#include "nox_nln_globaldata.H"
#include "nox_nln_constraint_interface_required.H"
#include "nox_nln_interface_jacobian.H"
#include "nox_nln_linearsystem.H"
#include "nox_nln_linearsystem_factory.H"
#include "nox_nln_constraint_group.H"
#include "nox_nln_inner_statustest_factory.H"

#include <Teuchos_ParameterList.hpp>

#include <Epetra_Operator.h>
#include "../linalg/linalg_solver.H"

#include <NOX_Utils.H>
#include <NOX_Epetra_Scaling.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_StatusTest_Generic.H>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Problem::Problem(Teuchos::RCP<NOX::NLN::GlobalData>& noxNlnGlobalData,
    Teuchos::RCP<NOX::Epetra::Vector>& x,
    Teuchos::RCP<LINALG::SparseOperator> A) :
  noxNlnGlobalData_(noxNlnGlobalData)
{
  Initialize(x,A);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Problem::Initialize(Teuchos::RCP<NOX::Epetra::Vector>& x,
                                      Teuchos::RCP<LINALG::SparseOperator> A)
{
  // in the standard case, we use the input rhs and matrix
  // ToDo Check if CreateView is sufficient
  xVector_ = x;
  jac_     = A;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem> NOX::NLN::Problem::CreateLinearSystem() const
{
  const std::map<NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& linSolvers =
      noxNlnGlobalData_->GetLinSolvers();
  const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq =
      noxNlnGlobalData_->GetRequiredInterface();
  const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac =
      noxNlnGlobalData_->GetJacobianInterface();
  const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec =
      noxNlnGlobalData_->GetPreconditionerInterface();

  Teuchos::ParameterList& params = noxNlnGlobalData_->GetNlnParameterList();
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys;
  // printing parameters
  Teuchos::ParameterList& printParams = params.sublist("Printing",true);
  // linear solver parameters
  Teuchos::ParameterList& lsParams    = params.sublist("Direction",true).
      sublist("Newton",true).sublist("Linear Solver",true);

  // preconditioner
  // not used at the moment
  Teuchos::RCP<LINALG::SparseOperator> preconditioner = Teuchos::null;

  // scaling
  // not used at the moment
  Teuchos::RCP<NOX::Epetra::Scaling> s = Teuchos::null;

  // build the linear system --> factory call
  linSys = NOX::NLN::LinSystem::BuildLinearSystem
      (printParams,lsParams,linSolvers,iReq,iJac,jac_,iPrec,preconditioner,*xVector_,s);

  return linSys;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Abstract::Group> NOX::NLN::Problem::CreateGroup(
    const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys
    ) const
{
  Teuchos::RCP<NOX::Abstract::Group> noxgrp = Teuchos::null;

  Teuchos::ParameterList& params = noxNlnGlobalData_->GetNlnParameterList();
  const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq =
      noxNlnGlobalData_->GetRequiredInterface();
  if (noxNlnGlobalData_->GetIsConstrained())
  {
    const ConstraintMap& iconstr =
          noxNlnGlobalData_->GetConstraintInterfaces();
    noxgrp = Teuchos::rcp(new NOX::NLN::CONSTRAINT::Group(params.sublist("Printing"),
        params.sublist("Group Options"),iReq,*xVector_,linSys,iconstr));
  }
  else
  {
    noxgrp = Teuchos::rcp(new NOX::NLN::Group(params.sublist("Printing"),
        params.sublist("Group Options"),iReq,*xVector_,linSys));
  }

  return noxgrp;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Problem::CreateStatusTests(
    Teuchos::RCP<NOX::StatusTest::Generic>& outerTests,
    Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests) const
{
  Teuchos::ParameterList& p = noxNlnGlobalData_->GetNlnParameterList();

  // A "Outer Status Test" has to be supplied by the user
  Teuchos::ParameterList& oParams =
      p.sublist("Status Test",true).sublist("Outer Status Test",true);
  outerTests = NOX::NLN::StatusTest::BuildOuterStatusTests(oParams, noxNlnGlobalData_->GetNoxUtils());

  // A "Inner Status Test" is optional in some cases.
  // Check if there is a "Inner Status Test" sublist and if it is filled.
  if (p.sublist("Status Test",true).isSublist("Inner Status Test")
      and p.sublist("Status Test").sublist("Inner Status Test").numParams() != 0)
  {
    Teuchos::ParameterList& iParams =
        p.sublist("Status Test",true).sublist("Inner Status Test");
    innerTests = NOX::NLN::INNER::StatusTest::BuildInnerStatusTests(iParams, noxNlnGlobalData_->GetNoxUtils());
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Problem::CheckFinalStatus(
    const NOX::StatusTest::StatusType& finalStatus) const
{
  if (finalStatus != NOX::StatusTest::Converged)
  {
    throwError("CheckFinalStatus()",
        "The nonlinear solver did not converge!");
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Problem::throwError(
    const std::string& functionName,
    const std::string& errorMsg) const
{
  std::ostringstream msg;
  msg << "ERROR - NOX::NLN::NoxProblem::" << functionName
      << " - " << errorMsg << std::endl;
  dserror(msg.str());
}

#endif /* NOX_NLN_PROBLEM_CPP_ */
