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
static void
Invert(const LINALG::SparseMatrix& diag, Epetra_Vector& x, const Epetra_Vector& b)
{
  x.PutScalar(0.);
  const Epetra_Map& rowmap = diag.RowMap();
  const Epetra_Map& colmap = diag.ColMap();
  const Epetra_BlockMap& xmap = x.Map();

  if (not rowmap.SameAs(b.Map()))
    dserror("illegal map");

  int numelements = rowmap.NumMyElements();
  for (int i=0; i<numelements; ++i)
  {
    int NumEntries;
    int* Indices;
    double* Values;
    int err = diag.EpetraMatrix()->ExtractMyRowView(i,NumEntries,Values,Indices);
    if (err)
      dserror("ExtractMyRowView err=%d", err);
    if (NumEntries>1)
      dserror("no diagonal matrix");
    if (NumEntries==1)
    {
      int gid = colmap.GID(Indices[0]);
      if (gid<0)
        dserror("no gid for lid=%d",Indices[0]);
      int lid = xmap.LID(gid);
      if (lid<0)
        dserror("no lid for gid=%d",gid);
      x[lid] = 1./Values[0] * b[i];
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int FSI::LagrangianBlockMatrix::ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
#if 1

  const LINALG::SparseMatrix& structInnerOp = Matrix(0,0);
  const LINALG::SparseMatrix& fluidInnerOp  = Matrix(1,1);
  const LINALG::SparseMatrix& aleInnerOp    = Matrix(2,2);

  const LINALG::SparseMatrix& structCoupOpT = Matrix(0,3);
  const LINALG::SparseMatrix& fluidCoupOpT  = Matrix(1,3);
  const LINALG::SparseMatrix& aleCoupOpT    = Matrix(2,4);

  const LINALG::SparseMatrix& sfCoupOp      = Matrix(3,0);
  const LINALG::SparseMatrix& fluidCoupOp   = Matrix(3,1);
  const LINALG::SparseMatrix& saCoupOp      = Matrix(4,0);
  const LINALG::SparseMatrix& aleCoupOp     = Matrix(4,2);

  // Extract vector blocks
  // RHS

  const Epetra_Vector &x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  Teuchos::RCP<Epetra_Vector> sx = DomainExtractor().ExtractVector(x,0);
  Teuchos::RCP<Epetra_Vector> fx = DomainExtractor().ExtractVector(x,1);
  Teuchos::RCP<Epetra_Vector> ax = DomainExtractor().ExtractVector(x,2);

  Teuchos::RCP<Epetra_Vector> sfx = DomainExtractor().ExtractVector(x,3);
  Teuchos::RCP<Epetra_Vector> sax = DomainExtractor().ExtractVector(x,4);

  // initial guess

  Epetra_Vector &y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  Teuchos::RCP<Epetra_Vector> sy = RangeExtractor().ExtractVector(y,0);
  Teuchos::RCP<Epetra_Vector> fy = RangeExtractor().ExtractVector(y,1);
  Teuchos::RCP<Epetra_Vector> ay = RangeExtractor().ExtractVector(y,2);

  Teuchos::RCP<Epetra_Vector> sfy = RangeExtractor().ExtractVector(y,3);
  Teuchos::RCP<Epetra_Vector> say = RangeExtractor().ExtractVector(y,4);

#if 0

  {
    // Solve structure equations for sy with the rhs sx
    if (Comm().MyPID()==0)
      std::cout << "    structural solve: " << std::flush;

    Epetra_Time ts(Comm());

#if 0
    Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));
    structCoupOpT.Multiply(false,*sfx,*tmpsx);
    sx->Update(-1.0,*tmpsx,1.0);
#endif

    structuresolver_->Solve(structInnerOp.EpetraMatrix(),sy,sx,true);

    if (Comm().MyPID()==0)
      std::cout << ts.ElapsedTime() << std::flush;
  }

  {
    if (Comm().MyPID()==0)
      std::cout << "    ale solve: " << std::flush;

    Epetra_Time ta(Comm());

#if 0
    Teuchos::RCP<Epetra_Vector> tmpax = Teuchos::rcp(new Epetra_Vector(DomainMap(2)));
    Teuchos::RCP<Epetra_Vector> tmpay = Teuchos::rcp(new Epetra_Vector(DomainMap(2)));
    Teuchos::RCP<Epetra_Vector> tmpsax = Teuchos::rcp(new Epetra_Vector(DomainMap(4)));

    saCoupOp.Multiply(false,*sy,*tmpsax);
    Invert(aleCoupOp,*tmpax,*tmpsax);
    aleInnerOp.Multiply(false,*tmpax,*tmpay);
    tmpay->Update(1.0,*ax,1.0);
    Invert(aleCoupOpT,*say,*tmpay);

    tmpax->PutScalar(0.);
    aleCoupOpT.Multiply(false,*say,*tmpax);
    ax->Update(-1.0,*tmpax,1.0);
#endif

    alesolver_->Solve(aleInnerOp.EpetraMatrix(),ay,ax,true);

    if (Comm().MyPID()==0)
      std::cout << ta.ElapsedTime() << std::flush;
  }

  {
    if (Comm().MyPID()==0)
      std::cout << "    fluid solve: " << std::flush;

    Epetra_Time tf(Comm());

#if 0
    Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));
    Teuchos::RCP<Epetra_Vector> tmpfy = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));
    Teuchos::RCP<Epetra_Vector> tmpsfx = Teuchos::rcp(new Epetra_Vector(DomainMap(3)));

    sfCoupOp.Multiply(false,*sy,*tmpsfx);
    Invert(fluidCoupOp,*tmpfx,*tmpsfx);
    fluidInnerOp.Multiply(false,*tmpfx,*tmpfy);
    tmpfy->Update(1.0,*fx,1.0);
    Invert(fluidCoupOpT,*sfy,*tmpfy);

    tmpfx->PutScalar(0.);
    fluidCoupOpT.Multiply(false,*sfy,*tmpfx);
    fx->Update(-1.0,*tmpfx,1.0);
#endif

    fluidsolver_->Solve(fluidInnerOp.EpetraMatrix(),fy,fx,true);

    if (Comm().MyPID()==0)
      std::cout << tf.ElapsedTime() << "\n";
  }

#else

  Teuchos::RCP<Epetra_Vector> ds = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));
  Teuchos::RCP<Epetra_Vector> df = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));
  Teuchos::RCP<Epetra_Vector> da = Teuchos::rcp(new Epetra_Vector(DomainMap(2)));

  int err = structInnerOp.ExtractDiagonalCopy(*ds);
  if (err)
    dserror("ExtractDiagonalCopy: err=%d", err);

  err = fluidInnerOp.ExtractDiagonalCopy(*df);
  if (err)
    dserror("ExtractDiagonalCopy: err=%d", err);

  err = aleInnerOp.ExtractDiagonalCopy(*da);
  if (err)
    dserror("ExtractDiagonalCopy: err=%d", err);

  int numelements = sx->Map().NumMyElements();
  for (int i=0; i<numelements; ++i)
    (*sy)[i] = (*sx)[i] / (*ds)[i];

  numelements = fx->Map().NumMyElements();
  for (int i=0; i<numelements; ++i)
    (*fy)[i] = (*fx)[i] / (*df)[i];

  numelements = ax->Map().NumMyElements();
  for (int i=0; i<numelements; ++i)
    (*ay)[i] = (*ax)[i] / (*da)[i];

#endif

  // build solution vector

  y.Update(1.0,x,0.);

  RangeExtractor().InsertVector(*sy,0,y);
  RangeExtractor().InsertVector(*fy,1,y);
  RangeExtractor().InsertVector(*ay,2,y);

#if 0
  RangeExtractor().InsertVector(*sfy,3,y);
  RangeExtractor().InsertVector(*say,4,y);
#endif

  double nsx,nfx,nax,nsfx,nsax;
  double nsy,nfy,nay,nsfy,nsay;

  sx->Norm2(&nsx);
  fx->Norm2(&nfx);
  ax->Norm2(&nax);
  sfx->Norm2(&nsfx);
  sax->Norm2(&nsax);

  sy->Norm2(&nsy);
  fy->Norm2(&nfy);
  ay->Norm2(&nay);
  sfy->Norm2(&nsfy);
  say->Norm2(&nsay);

  cout << "|sx|=" << nsx
       << "\t|fx|=" << nfx
       << "\t|ax|=" << nax
       << "\t|sfx|=" << nsfx
       << "\t|sax|=" << nsax
       << "\n";
  cout << "|sy|=" << nsy
       << "\t|fy|=" << nfy
       << "\t|ay|=" << nay
       << "\t|sfy|=" << nsfy
       << "\t|say|=" << nsay
       << "\n";

#else

  // this is really evil :)
  Teuchos::RCP<LINALG::SparseMatrix> sparse = Merge();

  const Epetra_Vector &x = Teuchos::dyn_cast<const Epetra_Vector>(X);
  Epetra_Vector &y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  fluidsolver_->Solve(sparse->EpetraMatrix(),
                      Teuchos::rcp(&y,false),
                      Teuchos::rcp(new Epetra_Vector(x)),
                      true);

#endif

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* FSI::LagrangianBlockMatrix::Label() const
{
  return "FSI::OverlappingBlockMatrix";
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
