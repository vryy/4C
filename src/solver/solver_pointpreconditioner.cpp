/*!----------------------------------------------------------------------

\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
Created on: Jul 4, 2011
*----------------------------------------------------------------------*/

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

#include "../drt_lib/drt_dserror.H"

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
LINALG::SOLVER::DWindPreconditioner::DWindPreconditioner(
    FILE* outfile, Teuchos::RCP<PreconditionerType> preconditioner, Teuchos::ParameterList& azlist)
    : PreconditionerType(outfile), preconditioner_(preconditioner), azlist_(azlist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::DWindPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(matrix);
  if (A == NULL) dserror("CrsMatrix expected");

  if (create)
  {
    double tau = azlist_.get<double>("downwinding tau", 1.0);
    int nv = azlist_.get<int>("downwinding nv", 1);
    int np = azlist_.get<int>("downwinding np", 0);
    Teuchos::RCP<Epetra_CrsMatrix> fool = Teuchos::rcp(A, false);
    dwind_ = Teuchos::rcp(
        new LINALG::DownwindMatrix(fool, nv, np, tau, azlist_.get<int>("AZ_output", 0)));
  }

  dwA_ = dwind_->Permute(A);
  dwx_ = dwind_->Permute(x);
  dwb_ = dwind_->Permute(b);

  // build the actual preconditioner based on the permuted matrix
  preconditioner_->Setup(create, &*dwA_, &*dwx_, &*dwb_);

  Epetra_LinearProblem& lp = preconditioner_->LinearProblem();

  SetupLinearProblem(lp.GetOperator(), lp.GetLHS(), lp.GetRHS());
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::DWindPreconditioner::Finish(
    Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  preconditioner_->Finish(&*dwA_, &*dwx_, &*dwb_);

  // undo reordering of lhs, don't care for rhs
  dwind_->InvPermute(&*dwx_, x);
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
