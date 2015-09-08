/*-----------------------------------------------------------*/
/*!
\file nox_nln_constraint_problem.cpp

\maintainer Michael Hiermeier

\date Jun 30, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_constraint_problem.H"
#include "nox_nln_constraint_group.H"

#include <Epetra_Vector.h>

#include <Teuchos_ParameterList.hpp>

#include <NOX_Epetra_Vector.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::CONSTRAINT::NoxProblem::NoxProblem(Teuchos::RCP<NOX::NLN::GlobalData>& nlnGlobalData,
    Teuchos::RCP<Epetra_Vector>& x,
    Teuchos::RCP<LINALG::SparseOperator> A11) :
    NOX::NLN::NoxProblem(nlnGlobalData,x,A11)
{
  Initialize(x,A11);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::CONSTRAINT::NoxProblem::Initialize(Teuchos::RCP<Epetra_Vector>& x1,
                                                  Teuchos::RCP<LINALG::SparseOperator> A11)
{
  // create a default saddle point system, if necessary
  if (not GetConstrPtr()->IsCondensed())
  {
    xVector_ = GetConstrPtr()->BuildSaddlePointXVector(*x1);
    jac_ = GetConstrPtr()->InitializeSaddlePointSystem(A11);
  }
  else
  {
    // ToDo Check if CreateView is sufficient
    xVector_ =  Teuchos::rcp(new NOX::Epetra::Vector(x1,NOX::Epetra::Vector::CreateCopy));
    jac_ = A11;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Abstract::Group> NOX::NLN::CONSTRAINT::NoxProblem::CreateNoxGroup(
    const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys) const
{
  Teuchos::ParameterList& params = nlnGlobalData_->GetNlnParameterList();
  const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq =
      nlnGlobalData_->GetRequiredInterface();

  return Teuchos::rcp(new NOX::NLN::CONSTRAINT::Group(params.sublist("Printing"),
    iReq,*xVector_, linSys,GetConstrPtr()));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required>
NOX::NLN::CONSTRAINT::NoxProblem::GetConstrPtr() const
{
  return nlnGlobalData_->GetConstraintInterface();
}
