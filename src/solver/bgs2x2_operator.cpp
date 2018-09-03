/*!----------------------------------------------------------------------
\file bgs2x2_operator.cpp

<pre>
\brief Declaration
\level 1
\maintainer Lena Wiechert
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/

#include "ml_MultiLevelPreconditioner.h"

#include "../linalg/linalg_precond.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"

#include "bgs2x2_operator.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::BGS2x2_Operator::BGS2x2_Operator(Teuchos::RCP<Epetra_Operator> A,
    const Teuchos::ParameterList& list1, const Teuchos::ParameterList& list2, int global_iter,
    double global_omega, int block1_iter, double block1_omega, int block2_iter, double block2_omega,
    bool fliporder, FILE* outfile)
    : outfile_(outfile),
      list1_(list1),
      list2_(list2),
      global_iter_(global_iter),
      global_omega_(global_omega),
      block1_iter_(block1_iter),
      block1_omega_(block1_omega),
      block2_iter_(block2_iter),
      block2_omega_(block2_omega)
{
  if (!fliporder)
  {
    firstind_ = 0;
    secind_ = 1;
  }
  else
  {
    firstind_ = 1;
    secind_ = 0;

    // switch parameter lists according to fliporder
    list2_ = list1;
    list1_ = list2;
  }

  A_ = Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(A);
  if (A_ != Teuchos::null)
  {
    // Make a shallow copy of the block matrix as the preconditioners on the
    // blocks will be reused and the next assembly will replace the block
    // matrices.
    A_ = A_->Clone(View);
    mmex_ = A_->RangeExtractor();
  }
  else
  {
    dserror("BGS2x2: provided operator is not a BlockSparseMatrix!");
  }

  SetupBlockPreconditioners();

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::BGS2x2_Operator::SetupBlockPreconditioners()
{
  Teuchos::RCP<Teuchos::ParameterList> rcplist1 = Teuchos::rcp(&list1_, false);
  Teuchos::RCP<LINALG::Solver> s1 =
      Teuchos::rcp(new LINALG::Solver(rcplist1, A_->Comm(), outfile_));
  solver1_ = Teuchos::rcp(new LINALG::Preconditioner(s1));
  const LINALG::SparseMatrix& Op11 = A_->Matrix(firstind_, firstind_);
  solver1_->Setup(Op11.EpetraMatrix());

  Teuchos::RCP<Teuchos::ParameterList> rcplist2 = Teuchos::rcp(&list2_, false);
  Teuchos::RCP<LINALG::Solver> s2 =
      Teuchos::rcp(new LINALG::Solver(rcplist2, A_->Comm(), outfile_));
  solver2_ = Teuchos::rcp(new LINALG::Preconditioner(s2));
  const LINALG::SparseMatrix& Op22 = A_->Matrix(secind_, secind_);
  solver2_->Setup(Op22.EpetraMatrix());

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::BGS2x2_Operator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  Teuchos::RCP<Epetra_MultiVector> y1 = mmex_.ExtractVector(Y, firstind_);
  Teuchos::RCP<Epetra_MultiVector> y2 = mmex_.ExtractVector(Y, secind_);

  Teuchos::RCP<Epetra_MultiVector> z1 =
      Teuchos::rcp(new Epetra_MultiVector(y1->Map(), y1->NumVectors()));
  Teuchos::RCP<Epetra_MultiVector> z2 =
      Teuchos::rcp(new Epetra_MultiVector(y2->Map(), y2->NumVectors()));

  Teuchos::RCP<Epetra_MultiVector> tmpx1 =
      Teuchos::rcp(new Epetra_MultiVector(A_->DomainMap(firstind_), y1->NumVectors()));
  Teuchos::RCP<Epetra_MultiVector> tmpx2 =
      Teuchos::rcp(new Epetra_MultiVector(A_->DomainMap(secind_), y2->NumVectors()));

  const LINALG::SparseMatrix& Op11 = A_->Matrix(firstind_, firstind_);
  const LINALG::SparseMatrix& Op22 = A_->Matrix(secind_, secind_);
  const LINALG::SparseMatrix& Op12 = A_->Matrix(firstind_, secind_);
  const LINALG::SparseMatrix& Op21 = A_->Matrix(secind_, firstind_);

  // outer Richardson loop
  for (int run = 0; run < global_iter_; ++run)
  {
    Teuchos::RCP<Epetra_MultiVector> x1 = A_->DomainExtractor().ExtractVector(X, firstind_);
    Teuchos::RCP<Epetra_MultiVector> x2 = A_->DomainExtractor().ExtractVector(X, secind_);

    // ----------------------------------------------------------------
    // first block

    if (run > 0)
    {
      Op11.Multiply(false, *y1, *tmpx1);
      x1->Update(-1.0, *tmpx1, 1.0);
      Op12.Multiply(false, *y2, *tmpx1);
      x1->Update(-1.0, *tmpx1, 1.0);
    }

    solver1_->Solve(Op11.EpetraMatrix(), z1, x1, true);

    LocalBlockRichardson(solver1_, Op11, x1, z1, tmpx1, block1_iter_, block1_omega_);

    if (run > 0)
    {
      y1->Update(global_omega_, *z1, 1.0);
    }
    else
    {
      y1->Update(global_omega_, *z1, 0.0);
    }

    // ----------------------------------------------------------------
    // second block

    if (run > 0)
    {
      Op22.Multiply(false, *y2, *tmpx2);
      x2->Update(-1.0, *tmpx2, 1.0);
    }

    Op21.Multiply(false, *y1, *tmpx2);
    x2->Update(-1.0, *tmpx2, 1.0);

    solver2_->Solve(Op22.EpetraMatrix(), z2, x2, true);

    LocalBlockRichardson(solver2_, Op22, x2, z2, tmpx2, block2_iter_, block2_omega_);

    if (run > 0)
    {
      y2->Update(global_omega_, *z2, 1.0);
    }
    else
    {
      y2->Update(global_omega_, *z2, 0.0);
    }
  }

  mmex_.InsertVector(*y1, firstind_, Y);
  mmex_.InsertVector(*y2, secind_, Y);

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::BGS2x2_Operator::LocalBlockRichardson(Teuchos::RCP<LINALG::Preconditioner> solver,
    const LINALG::SparseMatrix& Op, Teuchos::RCP<Epetra_MultiVector> x,
    Teuchos::RCP<Epetra_MultiVector> y, Teuchos::RCP<Epetra_MultiVector> tmpx, int iter,
    double omega) const
{
  if (iter > 0)
  {
    y->Scale(omega);
    Teuchos::RCP<Epetra_MultiVector> tmpy =
        Teuchos::rcp(new Epetra_MultiVector(y->Map(), y->NumVectors()));

    for (int i = 0; i < iter; ++i)
    {
      Op.EpetraMatrix()->Multiply(false, *y, *tmpx);
      tmpx->Update(1.0, *x, -1.0);

      solver->Solve(Op.EpetraMatrix(), tmpy, tmpx, false);
      y->Update(omega, *tmpy, 1.0);
    }
  }
}
