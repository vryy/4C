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
  // Initialize is called in the base class.
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
    iReq,*xVector_,linSys,GetConstrInterfaces()));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::map<enum NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> >&
NOX::NLN::CONSTRAINT::NoxProblem::GetConstrInterfaces() const
{
  return nlnGlobalData_->GetConstraintInterfaces();
}
