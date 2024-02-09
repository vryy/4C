/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/

#include "baci_linear_solver_preconditioner_point.hpp"

#include "baci_utils_exceptions.hpp"

#include <Epetra_CrsMatrix.h>

BACI_NAMESPACE_OPEN

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
CORE::LINEAR_SOLVER::InfNormPreconditioner::InfNormPreconditioner(
    Teuchos::RCP<PreconditionerType> preconditioner)
    : preconditioner_(preconditioner)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void CORE::LINEAR_SOLVER::InfNormPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(matrix);
  if (A == nullptr) dserror("CrsMatrix expected");

  // do infnorm scaling
  rowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
  colsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
  A->InvRowSums(*rowsum_);
  A->InvColSums(*colsum_);

  SetupLinearProblem(matrix, x, b);

  lp_.LeftScale(*rowsum_);
  lp_.RightScale(*colsum_);

  preconditioner_->Setup(create, matrix, x, b);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void CORE::LINEAR_SOLVER::InfNormPreconditioner::Finish(
    Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  preconditioner_->Finish(matrix, x, b);

  // undo scaling
  Epetra_Vector invrowsum(rowsum_->Map(), false);
  invrowsum.Reciprocal(*rowsum_);
  rowsum_ = Teuchos::null;
  Epetra_Vector invcolsum(colsum_->Map(), false);
  invcolsum.Reciprocal(*colsum_);
  colsum_ = Teuchos::null;
  lp_.LeftScale(invrowsum);
  lp_.RightScale(invcolsum);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
CORE::LINEAR_SOLVER::SymDiagPreconditioner::SymDiagPreconditioner(
    Teuchos::RCP<PreconditionerType> preconditioner)
    : preconditioner_(preconditioner)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void CORE::LINEAR_SOLVER::SymDiagPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(matrix);
  if (A == nullptr) dserror("CrsMatrix expected");

  SetupLinearProblem(matrix, x, b);

  // do symmetric diagonal scaling
  Epetra_Vector invdiag(A->RowMap(), false);
  diag_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
  A->ExtractDiagonalCopy(*diag_);
  invdiag.Reciprocal(*diag_);

  lp_.LeftScale(invdiag);
  lp_.RightScale(invdiag);

  preconditioner_->Setup(create, matrix, x, b);
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void CORE::LINEAR_SOLVER::SymDiagPreconditioner::Finish(
    Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  preconditioner_->Finish(matrix, x, b);

  lp_.LeftScale(*diag_);
  lp_.RightScale(*diag_);
  diag_ = Teuchos::null;
}

BACI_NAMESPACE_CLOSE
