/*-----------------------------------------------------------*/
/*!
\file nox_nln_problem.cpp

\brief This class manages some of the necessary factory calls
       if a %NOX::NLN solver is supposed to be used. Therefore a
       lean function call becomes possible.

\maintainer Michael Hiermeier

\date Jun 30, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_PROBLEM_CPP_
#define NOX_NLN_PROBLEM_CPP_

#include "nox_nln_problem.H"  // class definition
#include "nox_nln_globaldata.H"
#include "nox_nln_constraint_interface_required.H"
#include "nox_nln_interface_jacobian.H"
#include "nox_nln_linearsystem.H"
#include "nox_nln_linearsystem_factory.H"
#include "nox_nln_constraint_group.H"
#include "nox_nln_inner_statustest_factory.H"
#include "nox_nln_aux.H"

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
NOX::NLN::Problem::Problem(const Teuchos::RCP<NOX::NLN::GlobalData>& noxNlnGlobalData)
    : isinit_(false),
      isjac_(false),
      noxNlnGlobalData_(noxNlnGlobalData),
      xVector_(NULL),
      jac_(NULL),
      precMat_(Teuchos::null)
{
  /* intentionally left blank */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Problem::Problem(const Teuchos::RCP<NOX::NLN::GlobalData>& noxNlnGlobalData,
    const Teuchos::RCP<NOX::Epetra::Vector>& x, const Teuchos::RCP<LINALG::SparseOperator>& A)
    : isinit_(false),
      noxNlnGlobalData_(noxNlnGlobalData),
      xVector_(NULL),
      jac_(NULL),
      precMat_(Teuchos::null)
{
  Initialize(x, A);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Problem::Initialize(
    const Teuchos::RCP<NOX::Epetra::Vector>& x, const Teuchos::RCP<LINALG::SparseOperator>& A)
{
  // in the standard case, we use the input rhs and matrix
  // ToDo Check if CreateView is sufficient
  if (x.is_null()) dserror("You have to provide a state vector pointer unequal Teuchos::null!");

  xVector_ = &x;
  isjac_ = (not A.is_null());
  jac_ = &A;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem> NOX::NLN::Problem::CreateLinearSystem() const
{
  CheckInit();
  if (not IsJac())
    dserror(
        "You have to set a jacobian first, before you can create a "
        "linear system!");

  const NOX::NLN::LinSystem::LinearSystemType linsystype =
      NOX::NLN::AUX::GetLinearSystemType(noxNlnGlobalData_->GetLinSolvers());
  Teuchos::RCP<NOX::Epetra::Scaling> scalingObject = noxNlnGlobalData_->GetScalingObject();

  // build the linear system --> factory call
  return NOX::NLN::LinSystem::BuildLinearSystem(
      linsystype, *noxNlnGlobalData_, *jac_, *xVector_, precMat_, scalingObject);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Abstract::Group> NOX::NLN::Problem::CreateGroup(
    const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys) const
{
  CheckInit();
  Teuchos::RCP<NOX::Abstract::Group> noxgrp = Teuchos::null;

  Teuchos::ParameterList& params = noxNlnGlobalData_->GetNlnParameterList();
  const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq =
      noxNlnGlobalData_->GetRequiredInterface();
  if (noxNlnGlobalData_->GetIsConstrained())
  {
    const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iconstr =
        noxNlnGlobalData_->GetConstraintInterfaces();
    noxgrp = Teuchos::rcp(new NOX::NLN::CONSTRAINT::Group(params.sublist("Printing"),
        params.sublist("Group Options"), iReq, **xVector_, linSys, iconstr));
  }
  else
  {
    noxgrp = Teuchos::rcp(new NOX::NLN::Group(
        params.sublist("Printing"), params.sublist("Group Options"), iReq, **xVector_, linSys));
  }

  return noxgrp;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Problem::CreateOuterStatusTest(
    Teuchos::RCP<NOX::StatusTest::Generic>& outerTests) const
{
  Teuchos::ParameterList& p = noxNlnGlobalData_->GetNlnParameterList();

  // A "Outer Status Test" has to be supplied by the user
  Teuchos::ParameterList& oParams =
      p.sublist("Status Test", true).sublist("Outer Status Test", true);
  outerTests =
      NOX::NLN::StatusTest::BuildOuterStatusTests(oParams, noxNlnGlobalData_->GetNoxUtils());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Problem::CreateStatusTests(Teuchos::RCP<NOX::StatusTest::Generic>& outerTests,
    Teuchos::RCP<NOX::NLN::INNER::StatusTest::Generic>& innerTests) const
{
  CreateOuterStatusTest(outerTests);

  // A "Inner Status Test" is optional in some cases.
  // Check if there is a "Inner Status Test" sublist and if it is filled.
  Teuchos::ParameterList& p = noxNlnGlobalData_->GetNlnParameterList();
  if (p.sublist("Status Test", true).isSublist("Inner Status Test") and
      p.sublist("Status Test").sublist("Inner Status Test").numParams() != 0)
  {
    Teuchos::ParameterList& iParams = p.sublist("Status Test", true).sublist("Inner Status Test");
    innerTests = NOX::NLN::INNER::StatusTest::BuildInnerStatusTests(
        iParams, noxNlnGlobalData_->GetNoxUtils());
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Problem::CheckFinalStatus(const NOX::StatusTest::StatusType& finalStatus) const
{
  if (finalStatus != NOX::StatusTest::Converged)
  {
    dserror("The nonlinear solver did not converge!");
  }

  return;
}

#endif /* NOX_NLN_PROBLEM_CPP_ */
