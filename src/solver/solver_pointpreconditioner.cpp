/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/

#include <Epetra_CrsMatrix.h>

#include "drt_dserror.H"

#include "solver_pointpreconditioner.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::NonePreconditioner::NonePreconditioner(FILE* outfile, Teuchos::ParameterList& list)
    : PreconditionerType(outfile)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::NonePreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  SetupLinearProblem(matrix, x, b);

  if (create)
  {
    prec_ = Teuchos::rcp(new NoneOperator(matrix));
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::InfNormPreconditioner::InfNormPreconditioner(
    Teuchos::RCP<PreconditionerType> preconditioner)
    : PreconditionerType(NULL), preconditioner_(preconditioner)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::InfNormPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(matrix);
  if (A == NULL) dserror("CrsMatrix expected");

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
void LINALG::SOLVER::InfNormPreconditioner::Finish(
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
LINALG::SOLVER::SymDiagPreconditioner::SymDiagPreconditioner(
    Teuchos::RCP<PreconditionerType> preconditioner)
    : PreconditionerType(NULL), preconditioner_(preconditioner)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::SymDiagPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(matrix);
  if (A == NULL) dserror("CrsMatrix expected");

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
void LINALG::SOLVER::SymDiagPreconditioner::Finish(
    Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  preconditioner_->Finish(matrix, x, b);

  lp_.LeftScale(*diag_);
  lp_.RightScale(*diag_);
  diag_ = Teuchos::null;
}
