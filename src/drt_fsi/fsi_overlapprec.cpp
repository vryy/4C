#ifdef CCADISCRET

#include "fsi_overlapprec.H"
#include <Epetra_Time.h>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
static void LocalBlockRichardson(Teuchos::RCP<LINALG::Preconditioner> solver,
                                 const LINALG::SparseMatrix& innerOp,
                                 Teuchos::RCP<Epetra_Vector> x,
                                 Teuchos::RCP<Epetra_Vector> y,
                                 Teuchos::RCP<Epetra_Vector> tmpx,
                                 int iterations,
                                 double omega,
                                 FILE* err,
                                 const Epetra_Comm& comm)
{
  if (iterations > 0)
  {
    y->Scale(omega);
    Teuchos::RCP<Epetra_Vector> tmpy = Teuchos::rcp(new Epetra_Vector(y->Map()));
    if (err!=NULL)
      if (comm.MyPID()==0)
        fprintf(err,"    fluid richardson (%d,%f):",iterations,omega);
    for (int i=0; i<iterations; ++i)
    {
      innerOp.EpetraMatrix()->Multiply(false,*y,*tmpx);
      tmpx->Update(1.0,*x,-1.0);

      if (err!=NULL)
      {
        double n;
        tmpx->Norm2(&n);
        if (comm.MyPID()==0)
          fprintf(err," %e",n);
      }

      solver->Solve(innerOp.EpetraMatrix(),tmpy,tmpx,false);
      y->Update(omega,*tmpy,1.0);
    }
    if (err!=NULL)
      if (comm.MyPID()==0)
        fprintf(err,"\n");
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::OverlappingBlockMatrix::OverlappingBlockMatrix(const LINALG::MultiMapExtractor& maps,
                                                    ADAPTER::Structure& structure,
                                                    ADAPTER::Fluid& fluid,
                                                    ADAPTER::Ale& ale,
                                                    bool structuresplit,
                                                    int symmetric,
                                                    double omega,
                                                    int iterations,
                                                    double somega,
                                                    int siterations,
                                                    double fomega,
                                                    int fiterations,
                                                    FILE* err)
  : LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(maps,maps,81,false,true),
    structuresplit_(structuresplit),
    symmetric_(symmetric),
    omega_(omega),
    iterations_(iterations),
    somega_(somega),
    siterations_(siterations),
    fomega_(fomega),
    fiterations_(fiterations),
    err_(err)
{
  fluidsolver_ = Teuchos::rcp(new LINALG::Preconditioner(fluid.LinearSolver()));

#ifndef BLOCKMATRIXMERGE
  structuresolver_ = Teuchos::rcp(new LINALG::Preconditioner(structure.LinearSolver()));

  constalesolver_ = ale.ConstPreconditioner();
  if (constalesolver_==Teuchos::null)
    alesolver_ = Teuchos::rcp(new LINALG::Preconditioner(ale.LinearSolver()));
  else
    alesolver_ = constalesolver_;
#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int FSI::OverlappingBlockMatrix::ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  if (UseTranspose())
    dserror("no transpose preconditioning");

#ifdef BLOCKMATRIXMERGE
  MergeSolve(X, Y);
#else
  SGS(X, Y);
#endif

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrix::MergeSolve(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
#ifdef BLOCKMATRIXMERGE
  const Epetra_Vector &x = Teuchos::dyn_cast<const Epetra_Vector>(X);
  Epetra_Vector &y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  fluidsolver_->Solve(sparse_->EpetraMatrix(),
                      Teuchos::rcp(&y,false),
                      Teuchos::rcp(new Epetra_Vector(x)),
                      true);
#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrix::SAFLowerGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  const LINALG::SparseMatrix& structInnerOp = Matrix(0,0);
  const LINALG::SparseMatrix& fluidInnerOp  = Matrix(1,1);
  const LINALG::SparseMatrix& fluidMeshOp   = Matrix(1,2);
  const LINALG::SparseMatrix& fluidBoundOp  = Matrix(1,0);
  const LINALG::SparseMatrix& aleInnerOp    = Matrix(2,2);

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

  {
    // Solve structure equations for sy with the rhs sx
#if 0
    if (Comm().MyPID()==0)
      std::cout << "    structural stime: " << std::flush;

    Epetra_Time ts(Comm());
#endif

    structuresolver_->Solve(structInnerOp.EpetraMatrix(),sy,sx,true);

    // do Richardson iteration
    Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));
    LocalBlockRichardson(structuresolver_,structInnerOp,sx,sy,tmpsx,siterations_,somega_,err_,Comm());

#if 0
    if (Comm().MyPID()==0)
      std::cout << std::scientific << ts.ElapsedTime() << std::flush;
#endif
  }

  {
    // Solve ale equations for ay with the rhs ax - A(I,Gamma) sy

#if 0
    if (Comm().MyPID()==0)
      std::cout << "    ale stime: " << std::flush;

    Epetra_Time ta(Comm());
#endif

    Teuchos::RCP<Epetra_Vector> tmpax = Teuchos::rcp(new Epetra_Vector(DomainMap(2)));
    if (structuresplit_)
    {
//       const LINALG::SparseMatrix& aleBoundOp    = Matrix(2,1);
//       aleBoundOp.Multiply(false,*fy,*tmpax);
//       ax->Update(-1.0,*tmpax,1.0);
    }
    else
    {
      const LINALG::SparseMatrix& aleBoundOp    = Matrix(2,0);
      aleBoundOp.Multiply(false,*sy,*tmpax);
      ax->Update(-1.0,*tmpax,1.0);
    }
    alesolver_->Solve(aleInnerOp.EpetraMatrix(),ay,ax,true);

#if 0
    if (Comm().MyPID()==0)
      std::cout << std::scientific << ta.ElapsedTime() << std::flush;
#endif
  }

  {
    // Solve fluid equations for fy with the rhs fx - F(I,Gamma) sy - F(Mesh) ay

#if 0
    if (Comm().MyPID()==0)
      std::cout << "    fluid stime: " << std::flush;

    Epetra_Time tf(Comm());
#endif

    Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));

    fluidBoundOp.Multiply(false,*sy,*tmpfx);
    fx->Update(-1.0,*tmpfx,1.0);
    fluidMeshOp.Multiply(false,*ay,*tmpfx);
    fx->Update(-1.0,*tmpfx,1.0);
    fluidsolver_->Solve(fluidInnerOp.EpetraMatrix(),fy,fx,true);

    // do Richardson iteration
    LocalBlockRichardson(fluidsolver_,fluidInnerOp,fx,fy,tmpfx,fiterations_,fomega_,err_,Comm());

#if 0
    if (Comm().MyPID()==0)
      std::cout << std::scientific << tf.ElapsedTime() << std::flush;
#endif
  }

#if 0
  if (Comm().MyPID()==0)
    std::cout << "\n";
#endif

  // build solution vector

  RangeExtractor().InsertVector(*sy,0,y);
  RangeExtractor().InsertVector(*fy,1,y);
  RangeExtractor().InsertVector(*ay,2,y);

#if 0

  double sn,an,fn;

  Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));
  structInnerOp.EpetraMatrix()->Multiply(false,*sy,*tmpsx);

  Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));
  fluidInnerOp.EpetraMatrix()->Multiply(false,*fy,*tmpfx);

  Teuchos::RCP<Epetra_Vector> tmpax = Teuchos::rcp(new Epetra_Vector(DomainMap(2)));
  aleInnerOp.EpetraMatrix()->Multiply(false,*ay,*tmpax);

  tmpsx->Update(-1.0,*sx,1.0); tmpsx->Norm2(&sn);
  tmpfx->Update(-1.0,*fx,1.0); tmpfx->Norm2(&fn);
  tmpax->Update(-1.0,*ax,1.0); tmpax->Norm2(&an);

  if (Comm().MyPID()==0)
  {
    std::cout << "    structural |res|: " << std::scientific << sn
              << "    ale |res|: " << std::scientific << an
              << "    fluid |res|: " << std::scientific << fn
              << "\n";
    if (err_!=NULL)
      fprintf(err_,"    structural |res|: %e    ale |res|: %e    fluid |res|: %e\n",sn,an,fn);
  }
#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrix::SGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  const LINALG::SparseMatrix& structInnerOp = Matrix(0,0);
  const LINALG::SparseMatrix& fluidInnerOp  = Matrix(1,1);
  const LINALG::SparseMatrix& fluidMeshOp   = Matrix(1,2);
  const LINALG::SparseMatrix& fluidBoundOp  = Matrix(1,0);
  const LINALG::SparseMatrix& aleInnerOp    = Matrix(2,2);

  // Extract vector blocks
  // RHS

  const Epetra_Vector &x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  // initial guess

  Epetra_Vector &y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  Teuchos::RCP<Epetra_Vector> sy = RangeExtractor().ExtractVector(y,0);
  Teuchos::RCP<Epetra_Vector> fy = RangeExtractor().ExtractVector(y,1);
  Teuchos::RCP<Epetra_Vector> ay = RangeExtractor().ExtractVector(y,2);

  Teuchos::RCP<Epetra_Vector> sz = Teuchos::rcp(new Epetra_Vector(sy->Map()));
  Teuchos::RCP<Epetra_Vector> fz = Teuchos::rcp(new Epetra_Vector(fy->Map()));
  Teuchos::RCP<Epetra_Vector> az = Teuchos::rcp(new Epetra_Vector(ay->Map()));

  Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));
  Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));
  Teuchos::RCP<Epetra_Vector> tmpax = Teuchos::rcp(new Epetra_Vector(DomainMap(2)));

  // outer Richardson loop
  for (int run=0; run<iterations_; ++run)
  {
    Teuchos::RCP<Epetra_Vector> sx = DomainExtractor().ExtractVector(x,0);
    Teuchos::RCP<Epetra_Vector> fx = DomainExtractor().ExtractVector(x,1);
    Teuchos::RCP<Epetra_Vector> ax = DomainExtractor().ExtractVector(x,2);

    // ----------------------------------------------------------------
    // lower GS

    {
      if (run>0)
      {
        const LINALG::SparseMatrix& structBoundOp  = Matrix(0,1);

        structInnerOp.Multiply(false,*sy,*tmpsx);
        sx->Update(-1.0,*tmpsx,1.0);
        structBoundOp.Multiply(false,*fy,*tmpsx);
        sx->Update(-1.0,*tmpsx,1.0);
      }

      // Solve structure equations for sy with the rhs sx
      structuresolver_->Solve(structInnerOp.EpetraMatrix(),sz,sx,true);

      // do Richardson iteration
      LocalBlockRichardson(structuresolver_,structInnerOp,sx,sz,tmpsx,siterations_,somega_,err_,Comm());

      if (run>0)
      {
        sy->Update(omega_,*sz,1.0);
      }
      else
      {
        sy->Update(omega_,*sz,0.0);
      }
    }

    {
      // Solve ale equations for ay with the rhs ax - A(I,Gamma) sy

      if (run>0)
      {
        aleInnerOp.Multiply(false,*ay,*tmpax);
        ax->Update(-1.0,*tmpax,1.0);
      }

      if (structuresplit_)
      {
        if (run>0)
        {
          const LINALG::SparseMatrix& aleBoundOp    = Matrix(2,1);
          aleBoundOp.Multiply(false,*fy,*tmpax);
          ax->Update(-1.0,*tmpax,1.0);
        }
      }
      else
      {
        const LINALG::SparseMatrix& aleBoundOp    = Matrix(2,0);
        aleBoundOp.Multiply(false,*sy,*tmpax);
        ax->Update(-1.0,*tmpax,1.0);
      }

      alesolver_->Solve(aleInnerOp.EpetraMatrix(),az,ax,true);

      if (run>0)
      {
        ay->Update(omega_,*az,1.0);
      }
      else
      {
        ay->Update(omega_,*az,0.0);
      }
    }

    {
      // Solve fluid equations for fy with the rhs fx - F(I,Gamma) sy - F(Mesh) ay

      if (run>0)
      {
        fluidInnerOp.Multiply(false,*fy,*tmpfx);
        fx->Update(-1.0,*tmpfx,1.0);
      }

      fluidBoundOp.Multiply(false,*sy,*tmpfx);
      fx->Update(-1.0,*tmpfx,1.0);
      fluidMeshOp.Multiply(false,*ay,*tmpfx);
      fx->Update(-1.0,*tmpfx,1.0);
      fluidsolver_->Solve(fluidInnerOp.EpetraMatrix(),fz,fx,true);

      LocalBlockRichardson(fluidsolver_,fluidInnerOp,fx,fz,tmpfx,fiterations_,fomega_,err_,Comm());

      if (run>0)
      {
        fy->Update(omega_,*fz,1.0);
      }
      else
      {
        fy->Update(omega_,*fz,0.0);
      }
    }

    // ----------------------------------------------------------------
    // the symmetric part of the pc can be skipped

    if (symmetric_)
    {

      sx = DomainExtractor().ExtractVector(x,0);
      fx = DomainExtractor().ExtractVector(x,1);
      ax = DomainExtractor().ExtractVector(x,2);

      // ----------------------------------------------------------------
      // upper GS

      {
        fluidInnerOp.Multiply(false,*fy,*tmpfx);
        fx->Update(-1.0,*tmpfx,1.0);
        fluidBoundOp.Multiply(false,*sy,*tmpfx);
        fx->Update(-1.0,*tmpfx,1.0);
        fluidMeshOp.Multiply(false,*ay,*tmpfx);
        fx->Update(-1.0,*tmpfx,1.0);

        fluidsolver_->Solve(fluidInnerOp.EpetraMatrix(),fz,fx,true);

        LocalBlockRichardson(fluidsolver_,fluidInnerOp,fx,fz,tmpfx,fiterations_,fomega_,err_,Comm());
        fy->Update(omega_,*fz,1.0);
      }

      {
        aleInnerOp.Multiply(false,*ay,*tmpax);
        ax->Update(-1.0,*tmpax,1.0);
        if (structuresplit_)
        {
          const LINALG::SparseMatrix& aleBoundOp    = Matrix(2,1);
          aleBoundOp.Multiply(false,*fy,*tmpax);
          ax->Update(-1.0,*tmpax,1.0);
        }
        else
        {
          const LINALG::SparseMatrix& aleBoundOp    = Matrix(2,0);
          aleBoundOp.Multiply(false,*sy,*tmpax);
          ax->Update(-1.0,*tmpax,1.0);
        }
        alesolver_->Solve(aleInnerOp.EpetraMatrix(),az,ax,true);
        ay->Update(omega_,*az,1.0);
      }

      {
        const LINALG::SparseMatrix& structBoundOp  = Matrix(0,1);

        structInnerOp.Multiply(false,*sy,*tmpsx);
        sx->Update(-1.0,*tmpsx,1.0);
        structBoundOp.Multiply(false,*fy,*tmpsx);
        sx->Update(-1.0,*tmpsx,1.0);

        structuresolver_->Solve(structInnerOp.EpetraMatrix(),sz,sx,true);

        LocalBlockRichardson(structuresolver_,structInnerOp,sx,sz,tmpsx,siterations_,somega_,err_,Comm());
        sy->Update(omega_,*sz,1.0);
      }
    }
  }

  RangeExtractor().InsertVector(*sy,0,y);
  RangeExtractor().InsertVector(*fy,1,y);
  RangeExtractor().InsertVector(*ay,2,y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* FSI::OverlappingBlockMatrix::Label() const
{
  return "FSI::OverlappingBlockMatrix";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrix::SetupPreconditioner()
{
#ifdef BLOCKMATRIXMERGE

  // this is really evil :)
  sparse_ = Merge();
  fluidsolver_    ->Setup(sparse_->EpetraMatrix());

#if 0
  Matrix(0,0).Dump("dump-struct");
  Matrix(1,1).Dump("dump-fluid");
  Matrix(2,2).Dump("dump-ale");

  static int count;
  count++;
  std::stringstream s;
  s << "dump-" << count;
  cout << "write: " << s.str() << "\n";
  sparse_->Dump(s.str());

//   Epetra_Vector diagonal(sparse_->RowMap());
//   int err = sparse_->ExtractDiagonalCopy(diagonal);
//   diagonal.Print(cout);
#endif

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
