#ifdef CCADISCRET

#include "fsi_overlapprec.H"
#include <Epetra_Time.h>



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::OverlappingBlockMatrix::OverlappingBlockMatrix(const LINALG::MultiMapExtractor& maps,
                                                    ADAPTER::Structure& structure,
                                                    ADAPTER::Fluid& fluid,
                                                    ADAPTER::Ale& ale,
                                                    bool structuresplit,
                                                    double somega,
                                                    int siterations,
                                                    double fomega,
                                                    int fiterations,
                                                    FILE* err)
  : LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(maps,maps,81,false,true),
    structuresplit_(structuresplit),
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
#ifdef BLOCKMATRIXMERGE
  MergeSolve(X, Y);
#else
  SAFLowerGS(X, Y);
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
    if (Comm().MyPID()==0)
      std::cout << "    structural stime: " << std::flush;

    Epetra_Time ts(Comm());

    structuresolver_->Solve(structInnerOp.EpetraMatrix(),sy,sx,true);

    // do Richardson iteration
    if (siterations_ > 0)
    {
      sy->Scale(somega_);
      Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));
      Teuchos::RCP<Epetra_Vector> tmpsy = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));
      if (err_!=NULL)
        if (Comm().MyPID()==0)
          fprintf(err_,"    structure richardson (%d,%f):",siterations_,somega_);
      for (int i=0; i<siterations_; ++i)
      {
        structInnerOp.EpetraMatrix()->Multiply(false,*sy,*tmpsx);
        tmpsx->Update(1.0,*sx,-1.0);

        if (err_!=NULL)
        {
          double n;
          tmpsx->Norm2(&n);
          if (Comm().MyPID()==0)
            fprintf(err_," %e",n);
        }

        structuresolver_->Solve(structInnerOp.EpetraMatrix(),tmpsy,tmpsx,false);
        sy->Update(somega_,*tmpsy,1.0);
      }
    }

    if (Comm().MyPID()==0)
      std::cout << std::scientific << ts.ElapsedTime() << std::flush;
  }

  {
    // Solve ale equations for ay with the rhs ax - A(I,Gamma) sy

    if (Comm().MyPID()==0)
      std::cout << "    ale stime: " << std::flush;

    Epetra_Time ta(Comm());

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

    if (Comm().MyPID()==0)
      std::cout << std::scientific << ta.ElapsedTime() << std::flush;
  }

  {
    // Solve fluid equations for fy with the rhs fx - F(I,Gamma) sy - F(Mesh) ay

    if (Comm().MyPID()==0)
      std::cout << "    fluid stime: " << std::flush;

    Epetra_Time tf(Comm());

    Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));

    fluidBoundOp.Multiply(false,*sy,*tmpfx);
    fx->Update(-1.0,*tmpfx,1.0);
    fluidMeshOp.Multiply(false,*ay,*tmpfx);
    fx->Update(-1.0,*tmpfx,1.0);
    fluidsolver_->Solve(fluidInnerOp.EpetraMatrix(),fy,fx,true);

    // do Richardson iteration
    if (fiterations_ > 0)
    {
      fy->Scale(fomega_);
      Teuchos::RCP<Epetra_Vector> tmpfy = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));
      if (err_!=NULL)
        if (Comm().MyPID()==0)
          fprintf(err_,"    fluid richardson (%d,%f):",fiterations_,fomega_);
      for (int i=0; i<fiterations_; ++i)
      {
        fluidInnerOp.EpetraMatrix()->Multiply(false,*fy,*tmpfx);
        tmpfx->Update(1.0,*fx,-1.0);

        if (err_!=NULL)
        {
          double n;
          tmpfx->Norm2(&n);
          if (Comm().MyPID()==0)
            fprintf(err_," %e",n);
        }

        fluidsolver_->Solve(fluidInnerOp.EpetraMatrix(),tmpfy,tmpfx,false);
        fy->Update(fomega_,*tmpfy,1.0);
      }
      if (err_!=NULL)
        if (Comm().MyPID()==0)
          fprintf(err_,"\n");
    }

    if (Comm().MyPID()==0)
      std::cout << std::scientific << tf.ElapsedTime() << std::flush;
  }

  if (Comm().MyPID()==0)
    std::cout << "\n";

  // build solution vector

  RangeExtractor().InsertVector(*sy,0,y);
  RangeExtractor().InsertVector(*fy,1,y);
  RangeExtractor().InsertVector(*ay,2,y);

#if 1

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
void FSI::OverlappingBlockMatrix::FSALowerGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  const LINALG::SparseMatrix& structInnerOp = Matrix(0,0);
  const LINALG::SparseMatrix& fluidInnerOp  = Matrix(1,1);
  //const LINALG::SparseMatrix& fluidMeshOp   = Matrix(1,2);
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
    // Solve fluid equations for fy with the rhs fx - F(I,Gamma) sy - F(Mesh) ay

    if (Comm().MyPID()==0)
      std::cout << "    fluid stime: " << std::flush;

    Epetra_Time tf(Comm());

    fluidsolver_->Solve(fluidInnerOp.EpetraMatrix(),fy,fx,true);

    if (Comm().MyPID()==0)
      std::cout << std::scientific << tf.ElapsedTime() << std::flush;
  }

  {
    // Solve structure equations for sy with the rhs sx
    if (Comm().MyPID()==0)
      std::cout << "    structural stime: " << std::flush;

    Epetra_Time ts(Comm());

    const LINALG::SparseMatrix& structBoundOp = Matrix(0,1);

    Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));
    structBoundOp.Multiply(false,*fy,*tmpsx);
    sx->Update(-1.0,*tmpsx,1.0);
    structuresolver_->Solve(structInnerOp.EpetraMatrix(),sy,sx,true);

    if (Comm().MyPID()==0)
      std::cout << std::scientific << ts.ElapsedTime() << std::flush;
  }

  {
    // Solve ale equations for ay with the rhs ax - A(I,Gamma) sy

    if (Comm().MyPID()==0)
      std::cout << "    ale stime: " << std::flush;

    Epetra_Time ta(Comm());

    Teuchos::RCP<Epetra_Vector> tmpax = Teuchos::rcp(new Epetra_Vector(DomainMap(2)));
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
    alesolver_->Solve(aleInnerOp.EpetraMatrix(),ay,ax,true);

    if (Comm().MyPID()==0)
      std::cout << std::scientific << ta.ElapsedTime() << std::flush;
  }

  if (Comm().MyPID()==0)
    std::cout << "\n";

  // build solution vector

  RangeExtractor().InsertVector(*sy,0,y);
  RangeExtractor().InsertVector(*fy,1,y);
  RangeExtractor().InsertVector(*ay,2,y);
}


#if 0

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrix::BlockRichardson(const LINALG::SparseMatrix& Aii,
                                                  const LINALG::SparseMatrix& Aig,
                                                  const LINALG::SparseMatrix& Agi,
                                                  const LINALG::SparseMatrix& Agg,
                                                  Teuchos::RCP<LINALG::Preconditioner> Sii,
                                                  Teuchos::RCP<LINALG::Preconditioner> Sgg,
                                                  Teuchos::RCP<Epetra_Vector> xi,
                                                  Teuchos::RCP<Epetra_Vector> xg,
                                                  Teuchos::RCP<Epetra_Vector> bi,
                                                  Teuchos::RCP<Epetra_Vector> bg,
                                                  double omega,
                                                  int itenum) const
{
  Epetra_Time tr(Comm());
  if (Comm().MyPID()==0)
    std::cout << "    richardson: " << std::flush;

  Teuchos::RCP<Epetra_Vector> ti1 = Teuchos::rcp(new Epetra_Vector(bi->Map()));
  Teuchos::RCP<Epetra_Vector> ti2 = Teuchos::rcp(new Epetra_Vector(bi->Map()));
  Teuchos::RCP<Epetra_Vector> tg1 = Teuchos::rcp(new Epetra_Vector(bg->Map()));
  Teuchos::RCP<Epetra_Vector> tg2 = Teuchos::rcp(new Epetra_Vector(bg->Map()));

  for (int i=0; i<itenum; ++i)
  {
    Aii.Multiply(false,*xi,*ti1);
    Aig.Multiply(false,*xg,*ti2);
    Agi.Multiply(false,*xi,*tg1);
    Agg.Multiply(false,*xg,*tg2);

    ti1->Update(-1.0,*ti2,1.0,*bi,-1.0);
    tg1->Update(-1.0,*tg2,1.0,*bg,-1.0);

    if (i%10==0)
    {
      double n1,n2;
      ti1->Norm2(&n1);
      tg1->Norm2(&n2);
      if (Comm().MyPID()==0)
        std::cout << " " << sqrt(n1*n1+n2*n2);
    }

    Sii->Solve(Aii.EpetraMatrix(),ti2,ti1,true);
    Agi.Multiply(false,*ti2,*tg2);
    tg1->Update(-1.0,*tg2,1.0);
    Sgg->Solve(Agg.EpetraMatrix(),tg2,tg1,true);

    xi->Update(omega,*ti2,1.0);
    xg->Update(omega,*tg2,1.0);

    if (Comm().MyPID()==0)
      std::cout << "." << std::flush;
  }

  if (Comm().MyPID()==0)
    std::cout << " " << tr.ElapsedTime() << std::flush;
}

#else

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrix::BlockRichardson(const LINALG::SparseMatrix& Aii,
                                                  const LINALG::SparseMatrix& Aig,
                                                  const LINALG::SparseMatrix& Agi,
                                                  const LINALG::SparseMatrix& Agg,
                                                  Teuchos::RCP<LINALG::Preconditioner> Sii,
                                                  Teuchos::RCP<LINALG::Preconditioner> Sgg,
                                                  Teuchos::RCP<Epetra_Vector> xi,
                                                  Teuchos::RCP<Epetra_Vector> xg,
                                                  Teuchos::RCP<Epetra_Vector> bi,
                                                  Teuchos::RCP<Epetra_Vector> bg,
                                                  double omega,
                                                  int itenum) const
{
  Epetra_Time tr(Comm());
  if (Comm().MyPID()==0)
    std::cout << "    richardson: " << std::flush;

  Teuchos::RCP<Epetra_Vector> ti1 = Teuchos::rcp(new Epetra_Vector(bi->Map()));
  Teuchos::RCP<Epetra_Vector> ti2 = Teuchos::rcp(new Epetra_Vector(bi->Map()));
  Teuchos::RCP<Epetra_Vector> tg1 = Teuchos::rcp(new Epetra_Vector(bg->Map()));
  Teuchos::RCP<Epetra_Vector> tg2 = Teuchos::rcp(new Epetra_Vector(bg->Map()));

  Teuchos::RCP<Epetra_Vector> b = Teuchos::rcp(new Epetra_Vector(*bg));
  Sii->Solve(Aii.EpetraMatrix(),ti1,bi,true);
  Agi.Multiply(false,*ti1,*tg1);
  b->Update(-1.0,*tg1,1.0);

  for (int i=0; i<itenum; ++i)
  {
    Agg.Multiply(false,*xg,*tg1);

    Aig.Multiply(false,*xg,*ti2);
    Sii->Solve(Aii.EpetraMatrix(),ti1,ti2,true);
    Agi.Multiply(false,*ti1,*tg2);

    tg1->Update(1.0,*tg2,1.0,*b,-1.0);

    if (i%10==0)
    {
      double n;
      tg1->Norm2(&n);
      if (Comm().MyPID()==0)
        std::cout << " " << n;
    }

    Sgg->Solve(Agg.EpetraMatrix(),tg2,tg1,true);

    xg->Update(omega,*tg2,1.0);

    if (Comm().MyPID()==0)
      std::cout << "." << std::flush;
  }

  Aig.Multiply(false,*xg,*ti2);
  ti2->Update(1.0,*bi,-1.0);
  Sii->Solve(Aii.EpetraMatrix(),xi,ti2,true);

  if (Comm().MyPID()==0)
    std::cout << " " << tr.ElapsedTime() << std::flush;
}

#endif


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrix::UpperGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  const LINALG::SparseMatrix& structInnerOp = Matrix(0,0);
  const LINALG::SparseMatrix& fluidInnerOp  = Matrix(1,1);
  const LINALG::SparseMatrix& fluidBoundOp  = Matrix(0,1);
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

  Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));

  // Solve fluid equations for fy with the rhs fx

  fluidsolver_->Solve(fluidInnerOp.EpetraMatrix(),fy,fx,true);

  // Solve structure equations for sy with the rhs sx - F(Gamma,I) fy

  fluidBoundOp.Multiply(false,*fy,*tmpsx);
  sx->Update(-1.0,*tmpsx,1.0);
  structuresolver_->Solve(structInnerOp.EpetraMatrix(),sy,sx,true);

  // Solve ale equations for ay with the rhs ax - A(I,Gamma) sy

  alesolver_->Solve(aleInnerOp.EpetraMatrix(),ay,ax,true);

  // build solution vector

  RangeExtractor().InsertVector(*sy,0,y);
  RangeExtractor().InsertVector(*fy,1,y);
  RangeExtractor().InsertVector(*ay,2,y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrix::SGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  // this implementation is not meant to be most efficient

  UpperGS(X,Y);

  Epetra_Vector &y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  Teuchos::RCP<Epetra_Vector> sy = RangeExtractor().ExtractVector(y,0);
  Teuchos::RCP<Epetra_Vector> fy = RangeExtractor().ExtractVector(y,1);
  Teuchos::RCP<Epetra_Vector> ay = RangeExtractor().ExtractVector(y,2);

  Teuchos::RCP<Epetra_Vector> sz = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));
  Teuchos::RCP<Epetra_Vector> fz = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));
  Teuchos::RCP<Epetra_Vector> az = Teuchos::rcp(new Epetra_Vector(DomainMap(2)));

  Matrix(0,0).Apply(*sy,*sz);
  Matrix(1,1).Apply(*fy,*fz);
  Matrix(2,2).Apply(*ay,*az);

  Epetra_Vector z(Y.Map());

  RangeExtractor().InsertVector(*sz,0,z);
  RangeExtractor().InsertVector(*fz,1,z);
  RangeExtractor().InsertVector(*az,2,z);

  SAFLowerGS(z,Y);
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
  static int count;
  count++;
  std::stringstream s;
  s << "dump-" << count;
  cout << "write: " << s.str() << "\n";
  sparse_->Dump(s.str());
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
