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
  const LINALG::SparseMatrix& Aig  = Matrix(2,1);

  const LINALG::SparseMatrix& CSF  = Matrix(3,0);
  const LINALG::SparseMatrix& CFS  = Matrix(3,1);

  const LINALG::SparseMatrix& CSFT = Matrix(0,3);
  const LINALG::SparseMatrix& CFST = Matrix(1,3);

  // Extract vector blocks

  // RHS
  const Epetra_Vector &x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  // initial guess
  Epetra_Vector &y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  // these are assumed to be zero
  Teuchos::RCP<Epetra_Vector> sy = RangeExtractor().ExtractVector(y,0);
  Teuchos::RCP<Epetra_Vector> fy = RangeExtractor().ExtractVector(y,1);
  Teuchos::RCP<Epetra_Vector> ay = RangeExtractor().ExtractVector(y,2);
  Teuchos::RCP<Epetra_Vector> ly = RangeExtractor().ExtractVector(y,3);

  Teuchos::RCP<Epetra_Vector> sz = Teuchos::rcp(new Epetra_Vector(sy->Map()));
  Teuchos::RCP<Epetra_Vector> fz = Teuchos::rcp(new Epetra_Vector(fy->Map()));
  Teuchos::RCP<Epetra_Vector> az = Teuchos::rcp(new Epetra_Vector(ay->Map()));

  Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));
  Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));
  Teuchos::RCP<Epetra_Vector> tmpax = Teuchos::rcp(new Epetra_Vector(DomainMap(2)));

  // block preconditioner

  if (symmetric_)
    dserror("symmetric Gauss-Seidel not implemented");

  // outer Richardson loop
  for (int run=0; run<iterations_; ++run)
  {
    Teuchos::RCP<Epetra_Vector> sx = DomainExtractor().ExtractVector(x,0);
    Teuchos::RCP<Epetra_Vector> fx = DomainExtractor().ExtractVector(x,1);
    Teuchos::RCP<Epetra_Vector> ax = DomainExtractor().ExtractVector(x,2);
    Teuchos::RCP<Epetra_Vector> lx = DomainExtractor().ExtractVector(x,3);

    // structure

    if (run>0)
    {
      S.Multiply(false,*sy,*tmpsx);
      sx->Update(-1.0,*tmpsx,1.0);
      CSFT.Multiply(false,*ly,*tmpsx);
      sx->Update(-1.0,*tmpsx,1.0);
    }

    structuresolver_->Solve(S.EpetraMatrix(),sz,sx,true);
    LocalBlockRichardson(structuresolver_,S,sx,sz,tmpsx,siterations_,somega_,err_,Comm());

    if (run>0)
    {
      sy->Update(omega_,*sz,1.0);
    }
    else
    {
      sy->Update(omega_,*sz,0.0);
    }

    // lagrange coupling

    // ale

    if (run>0)
    {
      Aii.Multiply(false,*ay,*tmpax);
      ax->Update(-1.0,*tmpax,1.0);
      Aig.Multiply(false,*fy,*tmpax);
      ax->Update(-1.0,*tmpax,1.0);
    }

    alesolver_->Solve(Aii.EpetraMatrix(),ay,ax,true);

    if (run>0)
    {
      ay->Update(omega_,*az,1.0);
    }
    else
    {
      ay->Update(omega_,*az,0.0);
    }

    // fluid

    if (run>0)
    {
      F.Multiply(false,*fy,*tmpfx);
      fx->Update(-1.0,*tmpfx,1.0);
    }

    Fg.Multiply(false,*ay,*tmpfx);
    fx->Update(-1.0,*tmpfx,1.0);
    CFST.Multiply(false,*ly,*tmpfx);
    fx->Update(-1.0,*tmpfx,1.0);

    fluidsolver_->Solve(F.EpetraMatrix(),fy,fx,true);
    LocalBlockRichardson(fluidsolver_,F,fx,fz,tmpfx,fiterations_,fomega_,err_,Comm());

    if (run>0)
    {
      fy->Update(omega_,*fz,1.0);
    }
    else
    {
      fy->Update(omega_,*fz,0.0);
    }
  }

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
