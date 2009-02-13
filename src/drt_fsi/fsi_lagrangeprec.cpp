#ifdef CCADISCRET

#include "fsi_lagrangeprec.H"
#include <Epetra_Time.h>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::LagrangianBlockMatrix::LagrangianBlockMatrix(const LINALG::MultiMapExtractor& maps,
                                                  ADAPTER::Structure& structure,
                                                  ADAPTER::Fluid& fluid,
                                                  ADAPTER::Ale& ale,
                                                  int symmetric,
                                                  double omega,
                                                  int iterations,
                                                  double somega,
                                                  int siterations,
                                                  double fomega,
                                                  int fiterations,
                                                  FILE* err)
  : BlockPreconditioningMatrix(maps,
                               structure,
                               fluid,
                               ale,
                               symmetric,
                               omega,
                               iterations,
                               somega,
                               siterations,
                               fomega,
                               fiterations,
                               err)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::LagrangianBlockMatrix::SGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
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
const char* FSI::LagrangianBlockMatrix::Label() const
{
  return "FSI::LagrangianBlockMatrix";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::LagrangianBlockMatrix::SetupPreconditioner()
{
#ifdef BLOCKMATRIXMERGE
  BlockPreconditioningMatrix::SetupPreconditioner();
#else
  const LINALG::SparseMatrix& structInnerOp = Matrix(0,0);
  const LINALG::SparseMatrix& fluidInnerOp  = Matrix(1,1);
  const LINALG::SparseMatrix& aleInnerOp    = Matrix(2,2);

  structuresolver_->Setup(structInnerOp.EpetraMatrix());
  fluidsolver_    ->Setup(fluidInnerOp .EpetraMatrix());
  if (constalesolver_==Teuchos::null)
    alesolver_    ->Setup(aleInnerOp   .EpetraMatrix());
#endif
}


#endif
