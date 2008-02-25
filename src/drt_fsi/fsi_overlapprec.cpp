#ifdef CCADISCRET

#include "fsi_overlapprec.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::OverlappingBlockMatrix::OverlappingBlockMatrix(const LINALG::MultiMapExtractor& maps,
                                                    Teuchos::RCP<LINALG::Solver> structuresolver,
                                                    Teuchos::RCP<LINALG::Solver> fluidsolver,
                                                    Teuchos::RCP<LINALG::Solver> alesolver)
  : LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(maps,maps,81,false,true),
    structuresolver_(structuresolver),
    fluidsolver_(fluidsolver),
    alesolver_(alesolver)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int FSI::OverlappingBlockMatrix::ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  const LINALG::SparseMatrix& structInnerOp = Matrix(0,0);
  const LINALG::SparseMatrix& fluidInnerOp  = Matrix(1,1);
  const LINALG::SparseMatrix& fluidBoundOp  = Matrix(1,0);
  const LINALG::SparseMatrix& aleInnerOp    = Matrix(2,2);
  const LINALG::SparseMatrix& aleBoundOp    = Matrix(2,0);

  // Extract vector blocks
  // RHS

  const Epetra_Vector &x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  Teuchos::RCP<Epetra_Vector> sx = DomainExtractor().ExtractVector(x,0);
  Teuchos::RCP<Epetra_Vector> fx = DomainExtractor().ExtractVector(x,1);
  Teuchos::RCP<Epetra_Vector> ax = DomainExtractor().ExtractVector(x,2);

  // initial guess

  Epetra_Vector &y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  Teuchos::RCP<Epetra_Vector> sy = RangeExtractor().ExtractVector(y,0);
  Teuchos::RCP<Epetra_Vector> fy = RangeExtractor().ExtractVector(y,1);
  Teuchos::RCP<Epetra_Vector> ay = RangeExtractor().ExtractVector(y,2);

  Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));
  Teuchos::RCP<Epetra_Vector> tmpax = Teuchos::rcp(new Epetra_Vector(DomainMap(2)));

  // Solve structure equations for sy with the rhs sx

  structuresolver_->Solve(structInnerOp.EpetraMatrix(),sy,sx,true);

  // Solve fluid equations for fy with the rhs fx - F(I,Gamma) sy

  fluidBoundOp.Multiply(false,*sy,*tmpfx);
  fx->Update(-1.0,*tmpfx,1.0);
  fluidsolver_->Solve(fluidInnerOp.EpetraMatrix(),fy,fx,true);

  // Solve ale equations for ay with the rhs ax - A(I,Gamma) sy

  aleBoundOp.Multiply(false,*sy,*tmpax);
  ax->Update(-1.0,*tmpax,1.0);
  alesolver_->Solve(aleInnerOp.EpetraMatrix(),ay,ax,true);

  // build solution vector

  RangeExtractor().InsertVector(*sy,0,y);
  RangeExtractor().InsertVector(*fy,1,y);
  RangeExtractor().InsertVector(*ay,2,y);

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* FSI::OverlappingBlockMatrix::Label() const
{
  return "FSI::OverlappingBlockMatrix";
}


#endif
