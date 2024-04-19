/*-----------------------------------------------------------*/
/*! \file

\brief This class manages some of the necessary factory calls
       if a %NOX::NLN solver is supposed to be used. Therefore a
       lean function call becomes possible.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SOLVER_NONLIN_NOX_SOLVER_NONLIN_NOX_PROBLEM_CPP
#define SOLVER_NONLIN_NOX_SOLVER_NONLIN_NOX_PROBLEM_CPP

#include "baci_solver_nonlin_nox_problem.hpp"  // class definition

#include "baci_linear_solver_method_linalg.hpp"
#include "baci_solver_nonlin_nox_aux.hpp"
#include "baci_solver_nonlin_nox_constraint_group.hpp"
#include "baci_solver_nonlin_nox_constraint_interface_required.hpp"
#include "baci_solver_nonlin_nox_globaldata.hpp"
#include "baci_solver_nonlin_nox_inner_statustest_factory.hpp"
#include "baci_solver_nonlin_nox_interface_jacobian.hpp"
#include "baci_solver_nonlin_nox_linearsystem.hpp"
#include "baci_solver_nonlin_nox_linearsystem_factory.hpp"
#include "baci_solver_nonlin_nox_singlestep_group.hpp"

#include <Epetra_Operator.h>
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_Scaling.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_StatusTest_Generic.H>
#include <NOX_Utils.H>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Problem::Problem(const Teuchos::RCP<NOX::NLN::GlobalData>& noxNlnGlobalData)
    : isinit_(false),
      isjac_(false),
      noxNlnGlobalData_(noxNlnGlobalData),
      xVector_(nullptr),
      jac_(nullptr),
      precMat_(Teuchos::null)
{
  /* intentionally left blank */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Problem::Problem(const Teuchos::RCP<NOX::NLN::GlobalData>& noxNlnGlobalData,
    const Teuchos::RCP<::NOX::Epetra::Vector>& x,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& A)
    : isinit_(false),
      noxNlnGlobalData_(noxNlnGlobalData),
      xVector_(nullptr),
      jac_(nullptr),
      precMat_(Teuchos::null)
{
  Initialize(x, A);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::Problem::Initialize(const Teuchos::RCP<::NOX::Epetra::Vector>& x,
    const Teuchos::RCP<CORE::LINALG::SparseOperator>& A)
{
  // in the standard case, we use the input rhs and matrix
  // ToDo Check if CreateView is sufficient
  if (x.is_null())
    FOUR_C_THROW("You have to provide a state vector pointer unequal Teuchos::null!");

  xVector_ = &x;
  isjac_ = (not A.is_null());
  jac_ = &A;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Epetra::LinearSystem> NOX::NLN::Problem::CreateLinearSystem() const
{
  CheckInit();
  if (not IsJac())
    FOUR_C_THROW(
        "You have to set a jacobian first, before you can create a "
        "linear system!");

  const NOX::NLN::LinSystem::LinearSystemType linsystype =
      NOX::NLN::AUX::GetLinearSystemType(noxNlnGlobalData_->GetLinSolvers());
  Teuchos::RCP<::NOX::Epetra::Scaling> scalingObject = noxNlnGlobalData_->GetScalingObject();

  // build the linear system --> factory call
  return NOX::NLN::LinSystem::BuildLinearSystem(
      linsystype, *noxNlnGlobalData_, *jac_, *xVector_, precMat_, scalingObject);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Abstract::Group> NOX::NLN::Problem::CreateGroup(
    const Teuchos::RCP<::NOX::Epetra::LinearSystem>& linSys) const
{
  CheckInit();
  Teuchos::RCP<::NOX::Abstract::Group> noxgrp = Teuchos::null;

  Teuchos::ParameterList& params = noxNlnGlobalData_->GetNlnParameterList();
  const Teuchos::RCP<::NOX::Epetra::Interface::Required>& iReq =
      noxNlnGlobalData_->GetRequiredInterface();
  const std::string nlnSolver = params.get<std::string>("Nonlinear Solver", "");
  if (noxNlnGlobalData_->GetIsConstrained())
  {
    const NOX::NLN::CONSTRAINT::ReqInterfaceMap& iconstr =
        noxNlnGlobalData_->GetConstraintInterfaces();
    noxgrp = Teuchos::rcp(new NOX::NLN::CONSTRAINT::Group(params.sublist("Printing"),
        params.sublist("Group Options"), iReq, **xVector_, linSys, iconstr));
  }
  else if (nlnSolver.compare("Single Step") == 0)
  {
    std::cout << "Single Step Group is selected" << std::endl;
    noxgrp = Teuchos::rcp(new NOX::NLN::SINGLESTEP::Group(
        params.sublist("Printing"), params.sublist("Group Options"), iReq, **xVector_, linSys));
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
    Teuchos::RCP<::NOX::StatusTest::Generic>& outerTests) const
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
void NOX::NLN::Problem::CreateStatusTests(Teuchos::RCP<::NOX::StatusTest::Generic>& outerTests,
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
void NOX::NLN::Problem::CheckFinalStatus(const ::NOX::StatusTest::StatusType& finalStatus) const
{
  if (finalStatus != ::NOX::StatusTest::Converged)
  {
    FOUR_C_THROW("The nonlinear solver did not converge!");
  }

  return;
}

#endif

FOUR_C_NAMESPACE_CLOSE
