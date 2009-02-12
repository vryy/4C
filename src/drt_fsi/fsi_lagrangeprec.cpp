#ifdef CCADISCRET

#include "fsi_lagrangeprec.H"
#include <Epetra_Time.h>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::LagrangianBlockMatrix::LagrangianBlockMatrix(const LINALG::MultiMapExtractor& maps,
                                                  Teuchos::RCP<LINALG::Solver> structuresolver,
                                                  Teuchos::RCP<LINALG::Solver> fluidsolver,
                                                  Teuchos::RCP<LINALG::Solver> alesolver)
  : LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(maps,maps,81,false,true)
{
  structuresolver_ = Teuchos::rcp(new LINALG::Preconditioner(structuresolver));
  fluidsolver_ = Teuchos::rcp(new LINALG::Preconditioner(fluidsolver));
  alesolver_ = Teuchos::rcp(new LINALG::Preconditioner(alesolver));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int FSI::LagrangianBlockMatrix::ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
#ifdef BLOCKMATRIXMERGE
  MergeSolve(X,Y);
#else
  SAFLowerGS(X,Y);
#endif

  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::LagrangianBlockMatrix::SAFLowerGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  // Extract matrix blocks

  const LINALG::SparseMatrix& S    = Matrix(0,0);
  const LINALG::SparseMatrix& F    = Matrix(1,1);
  const LINALG::SparseMatrix& Fg   = Matrix(1,2);
  const LINALG::SparseMatrix& Aii  = Matrix(2,2);
  //const LINALG::SparseMatrix& Aig  = Matrix(2,1); <-- this is dropped momentarily

  const LINALG::SparseMatrix& CSF  = Matrix(3,0);
  const LINALG::SparseMatrix& CFS  = Matrix(3,1);

  //const LINALG::SparseMatrix& CSFT = Matrix(0,3); <-- this is the one we drop here
  const LINALG::SparseMatrix& CFST = Matrix(1,3);

  // Extract vector blocks

  // RHS
  const Epetra_Vector &x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  Teuchos::RCP<Epetra_Vector> sx = DomainExtractor().ExtractVector(x,0);
  Teuchos::RCP<Epetra_Vector> fx = DomainExtractor().ExtractVector(x,1);
  Teuchos::RCP<Epetra_Vector> ax = DomainExtractor().ExtractVector(x,2);
  Teuchos::RCP<Epetra_Vector> lx = DomainExtractor().ExtractVector(x,3);

  // initial guess
  Epetra_Vector &y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  Teuchos::RCP<Epetra_Vector> sy = RangeExtractor().ExtractVector(y,0);
  Teuchos::RCP<Epetra_Vector> fy = RangeExtractor().ExtractVector(y,1);
  Teuchos::RCP<Epetra_Vector> ay = RangeExtractor().ExtractVector(y,2);
  Teuchos::RCP<Epetra_Vector> ly = RangeExtractor().ExtractVector(y,3);
  Teuchos::RCP<Epetra_Vector> lytemp = RangeExtractor().ExtractVector(y,3);

  // block preconditioner

  structuresolver_->Solve(S.EpetraMatrix(),sy,sx,true);
  alesolver_->Solve(Aii.EpetraMatrix(),ay,ax,true);
  fluidsolver_->Solve(F.EpetraMatrix(),fy,fx,true);

  // build solution vector

  RangeExtractor().InsertVector(*sy,0,y);
  RangeExtractor().InsertVector(*fy,1,y);
  RangeExtractor().InsertVector(*ay,2,y);
  RangeExtractor().InsertVector(*ly,3,y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::LagrangianBlockMatrix::MergeSolve(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  // this is really evil :)

  Teuchos::RCP<LINALG::SparseMatrix> sparse = Merge();

  const Epetra_Vector &x = Teuchos::dyn_cast<const Epetra_Vector>(X);
  Epetra_Vector &y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  fluidsolver_->Solve(sparse->EpetraMatrix(),
                      Teuchos::rcp(&y,false),
                      Teuchos::rcp(new Epetra_Vector(x)),
                      true);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* FSI::LagrangianBlockMatrix::Label() const
{
  return "FSI::LagrangianBlockMatrix";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::LagrangianBlockMatrix::SetupPreconditioner()
{
  const LINALG::SparseMatrix& structInnerOp = Matrix(0,0);
  const LINALG::SparseMatrix& fluidInnerOp  = Matrix(1,1);
  const LINALG::SparseMatrix& aleInnerOp    = Matrix(2,2);

  structuresolver_->Setup(structInnerOp.EpetraMatrix());
  fluidsolver_    ->Setup(fluidInnerOp .EpetraMatrix());
  alesolver_      ->Setup(aleInnerOp   .EpetraMatrix());
}


#endif
