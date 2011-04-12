/*!----------------------------------------------------------------------
\file bgs2x2_operator.cpp

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#include "bgs2x2_operator.H"
#include "linalg_solver.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::BGS2x2_Operator::BGS2x2_Operator(RCP<Epetra_Operator> A,
                                         const ParameterList& list1,
                                         const ParameterList& list2,
                                         int global_iter,
                                         double global_omega,
                                         int block1_iter,
                                         double block1_omega,
                                         int block2_iter,
                                         double block2_omega,
                                         bool fliporder,
                                         FILE* outfile)
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
  }

  A_ = rcp_dynamic_cast<BlockSparseMatrixBase>(A);
  if (A_!=null)
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

  Teuchos::RCP<ParameterList> rcplist1 = rcp(&list1_,false);
  Teuchos::RCP<LINALG::Solver> s1 = rcp(new LINALG::Solver(rcplist1,A_->Comm(),outfile_));
  solver1_ = Teuchos::rcp(new LINALG::Preconditioner(s1));

  Teuchos::RCP<ParameterList> rcplist2 = rcp(&list2_,false);
  Teuchos::RCP<LINALG::Solver> s2 = rcp(new LINALG::Solver(rcplist2,A_->Comm(),outfile_));
  solver2_ = Teuchos::rcp(new LINALG::Preconditioner(s2));

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::BGS2x2_Operator::SetupPreconditioner()
{
  const LINALG::SparseMatrix& Op11 = A_->Matrix(firstind_,firstind_);
  solver1_->Setup(Op11.EpetraMatrix());

  const LINALG::SparseMatrix& Op22 = A_->Matrix(secind_,secind_);
  solver2_->Setup(Op22.EpetraMatrix());

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::BGS2x2_Operator::ApplyInverse(const Epetra_MultiVector& X,
                                          Epetra_MultiVector& Y) const
{
  const Epetra_Vector &x = Teuchos::dyn_cast<const Epetra_Vector>(X);
  Epetra_Vector &y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  Teuchos::RCP<Epetra_Vector> y1 = mmex_.ExtractVector(y,firstind_);
  Teuchos::RCP<Epetra_Vector> y2 = mmex_.ExtractVector(y,secind_);

  Teuchos::RCP<Epetra_Vector> z1 = Teuchos::rcp(new Epetra_Vector(y1->Map()));
  Teuchos::RCP<Epetra_Vector> z2 = Teuchos::rcp(new Epetra_Vector(y2->Map()));

  Teuchos::RCP<Epetra_Vector> tmpx1 = Teuchos::rcp(new Epetra_Vector(A_->DomainMap(firstind_)));
  Teuchos::RCP<Epetra_Vector> tmpx2 = Teuchos::rcp(new Epetra_Vector(A_->DomainMap(secind_)));

  const LINALG::SparseMatrix& Op11 = A_->Matrix(firstind_,firstind_);
  const LINALG::SparseMatrix& Op22 = A_->Matrix(secind_,secind_);
  const LINALG::SparseMatrix& Op12 = A_->Matrix(firstind_,secind_);
  const LINALG::SparseMatrix& Op21 = A_->Matrix(secind_,firstind_);

  // outer Richardson loop
  for (int run=0; run<global_iter_; ++run)
  {
    Teuchos::RCP<Epetra_Vector> x1 = A_->DomainExtractor().ExtractVector(x,firstind_);
    Teuchos::RCP<Epetra_Vector> x2 = A_->DomainExtractor().ExtractVector(x,secind_);

    // ----------------------------------------------------------------
    // first block

    if (run>0)
    {
      Op11.Multiply(false,*y1,*tmpx1);
      x1->Update(-1.0,*tmpx1,1.0);
      Op12.Multiply(false,*y2,*tmpx1);
      x1->Update(-1.0,*tmpx1,1.0);
    }

    solver1_->Solve(Op11.EpetraMatrix(),z1,x1,true);

    LocalBlockRichardson(solver1_,Op11,x1,z1,tmpx1,block1_iter_,block1_omega_);

    if (run>0)
    {
      y1->Update(global_omega_,*z1,1.0);
    }
    else
    {
      y1->Update(global_omega_,*z1,0.0);
    }

    // ----------------------------------------------------------------
    // second block

    if (run>0)
    {
      Op22.Multiply(false,*y2,*tmpx2);
      x2->Update(-1.0,*tmpx2,1.0);
    }

    Op21.Multiply(false,*y1,*tmpx2);
    x2->Update(-1.0,*tmpx2,1.0);

    solver2_->Solve(Op22.EpetraMatrix(),z2,x2,true);

    LocalBlockRichardson(solver2_,Op22,x2,z2,tmpx2,block2_iter_,block2_omega_);

    if (run>0)
    {
      y2->Update(global_omega_,*z2,1.0);
    }
    else
    {
      y2->Update(global_omega_,*z2,0.0);
    }
  }

  mmex_.InsertVector(*y1,firstind_,y);
  mmex_.InsertVector(*y2,secind_,y);

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::BGS2x2_Operator::LocalBlockRichardson(Teuchos::RCP<LINALG::Preconditioner> solver,
                                                   const LINALG::SparseMatrix& Op,
                                                   Teuchos::RCP<Epetra_Vector> x,
                                                   Teuchos::RCP<Epetra_Vector> y,
                                                   Teuchos::RCP<Epetra_Vector> tmpx,
                                                   int iter,
                                                   double omega) const
{
  if (iter > 0)
  {
    y->Scale(omega);
    Teuchos::RCP<Epetra_Vector> tmpy = Teuchos::rcp(new Epetra_Vector(y->Map()));

    for (int i=0; i < iter; ++i)
    {
      Op.EpetraMatrix()->Multiply(false,*y,*tmpx);
      tmpx->Update(1.0,*x,-1.0);

      solver->Solve(Op.EpetraMatrix(),tmpy,tmpx,false);
      y->Update(omega,*tmpy,1.0);
    }
  }
}


#endif  // #ifdef CCADISCRET
