#ifdef CCADISCRET

#include "fsi_overlapprec.H"
#include <Epetra_Time.h>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::OverlappingBlockMatrix::OverlappingBlockMatrix(const LINALG::MultiMapExtractor& maps,
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
int FSI::OverlappingBlockMatrix::ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  LowerGS(X, Y);
  //UpperGS(X, Y);
  //SGS(X, Y);
  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrix::LowerGS(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
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

  {
    // Solve structure equations for sy with the rhs sx
    if (Comm().MyPID()==0)
      std::cout << "    structural solve: " << std::flush;


    Epetra_Time ts(Comm());

    structuresolver_->Solve(structInnerOp.EpetraMatrix(),sy,sx,true);

    if (Comm().MyPID()==0)
      std::cout << ts.ElapsedTime() << std::flush;
  }

  {
    // Solve fluid equations for fy with the rhs fx - F(I,Gamma) sy

    //fluidInnerOp.EpetraMatrix()->Print(cout);

    if (Comm().MyPID()==0)
      std::cout << "    fluid solve: " << std::flush;

    Epetra_Time tf(Comm());

    fluidBoundOp.Multiply(false,*sy,*tmpfx);
    fx->Update(-1.0,*tmpfx,1.0);
    fluidsolver_->Solve(fluidInnerOp.EpetraMatrix(),fy,fx,true);

    if (Comm().MyPID()==0)
      std::cout << tf.ElapsedTime() << std::flush;
  }

  {
    // Solve ale equations for ay with the rhs ax - A(I,Gamma) sy

    if (Comm().MyPID()==0)
      std::cout << "    ale solve: " << std::flush;

    Epetra_Time ta(Comm());

    aleBoundOp.Multiply(false,*sy,*tmpax);
    ax->Update(-1.0,*tmpax,1.0);
    alesolver_->Solve(aleInnerOp.EpetraMatrix(),ay,ax,true);

    if (Comm().MyPID()==0)
      std::cout << ta.ElapsedTime() << "\n";
  }

  // build solution vector

  RangeExtractor().InsertVector(*sy,0,y);
  RangeExtractor().InsertVector(*fy,1,y);
  RangeExtractor().InsertVector(*ay,2,y);
}


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

  LowerGS(z,Y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* FSI::OverlappingBlockMatrix::Label() const
{
  return "FSI::OverlappingBlockMatrix";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrix::SetupBlockPrecond()
{
  const LINALG::SparseMatrix& structInnerOp = Matrix(0,0);
  const LINALG::SparseMatrix& fluidInnerOp  = Matrix(1,1);
  const LINALG::SparseMatrix& aleInnerOp    = Matrix(2,2);

  structuresolver_->Setup(structInnerOp.EpetraMatrix());
  fluidsolver_    ->Setup(fluidInnerOp .EpetraMatrix());
  alesolver_      ->Setup(aleInnerOp   .EpetraMatrix());
}


#endif
