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
  MergeSolve(X, Y);
  return 0;

  const LINALG::SparseMatrix& S    = Matrix(0,0);

  const LINALG::SparseMatrix& Fii  = Matrix(1,1);
  const LINALG::SparseMatrix& Fig  = Matrix(1,2);
  const LINALG::SparseMatrix& Fgi  = Matrix(2,1);
  const LINALG::SparseMatrix& Fgg  = Matrix(2,2);

  const LINALG::SparseMatrix& Aii  = Matrix(3,3);
  const LINALG::SparseMatrix& Aig  = Matrix(3,4);
  const LINALG::SparseMatrix& Agi  = Matrix(4,3);
  const LINALG::SparseMatrix& Agg  = Matrix(4,4);

  //const LINALG::SparseMatrix& CSFT = Matrix(0,5);
  const LINALG::SparseMatrix& CFST = Matrix(2,5);
  const LINALG::SparseMatrix& CAST = Matrix(4,6);

  const LINALG::SparseMatrix& CSF  = Matrix(5,0);
  const LINALG::SparseMatrix& CFS  = Matrix(5,2);
  const LINALG::SparseMatrix& CSA  = Matrix(6,0);
  const LINALG::SparseMatrix& CAS  = Matrix(6,4);

//   cout << CFST.NormOne() << "  "
//        << CAST.NormOne() << "  "
//        << CSF.NormOne() << "  "
//        << CFS.NormOne() << "  "
//        << CSA.NormOne() << "  "
//        << CAS.NormOne() << "\n";

  // Extract vector blocks
  // RHS

  const Epetra_Vector &x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  Teuchos::RCP<Epetra_Vector> sx  = DomainExtractor().ExtractVector(x,0);
  Teuchos::RCP<Epetra_Vector> fix = DomainExtractor().ExtractVector(x,1);
  Teuchos::RCP<Epetra_Vector> fgx = DomainExtractor().ExtractVector(x,2);
  Teuchos::RCP<Epetra_Vector> aix = DomainExtractor().ExtractVector(x,3);
  Teuchos::RCP<Epetra_Vector> agx = DomainExtractor().ExtractVector(x,4);

  Teuchos::RCP<Epetra_Vector> sfx = DomainExtractor().ExtractVector(x,5);
  Teuchos::RCP<Epetra_Vector> sax = DomainExtractor().ExtractVector(x,6);

  // initial guess

  Epetra_Vector &y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  Teuchos::RCP<Epetra_Vector> sy  = RangeExtractor().ExtractVector(y,0);
  Teuchos::RCP<Epetra_Vector> fiy = RangeExtractor().ExtractVector(y,1);
  Teuchos::RCP<Epetra_Vector> fgy = RangeExtractor().ExtractVector(y,2);
  Teuchos::RCP<Epetra_Vector> aiy = RangeExtractor().ExtractVector(y,3);
  Teuchos::RCP<Epetra_Vector> agy = RangeExtractor().ExtractVector(y,4);

  Teuchos::RCP<Epetra_Vector> sfy = RangeExtractor().ExtractVector(y,5);
  Teuchos::RCP<Epetra_Vector> say = RangeExtractor().ExtractVector(y,6);

  {
    // Solve structure equations for sy with the rhs sx
    if (Comm().MyPID()==0)
      std::cout << "    structural solve: " << std::flush;

    Epetra_Time ts(Comm());

    structuresolver_->Solve(S.EpetraMatrix(),sy,sx,true);

    if (Comm().MyPID()==0)
      std::cout << ts.ElapsedTime() << std::flush;
  }

  {
    if (Comm().MyPID()==0)
      std::cout << "    ale solve: " << std::flush;

    Epetra_Time ta(Comm());

//     CSA.Multiply(false,*sy,*say);

//     // invert diagonal matrix CAS

//     int len = say->MyLength();
//     for (int i=0; i<len; ++i)
//     {
//       int NumEntries;
//       double *Values;
//       int *Indices;
//       int err = CAS.EpetraMatrix()->ExtractMyRowView(i,NumEntries,Values,Indices);
//       if (err)
//         dserror("ExtractMyRowView: err=%d",err);
//       if (NumEntries!=1)
//         dserror("diagonal matrix expected");
//       (*agy)[Indices[0]] = -(*say)[i] / Values[0];
//     }

//     Teuchos::RCP<Epetra_Vector> tmpaix = Teuchos::rcp(new Epetra_Vector(DomainMap(3)));
//     Aig.Multiply(false,*agy,*tmpaix);
//     aix->Update(-1.0,*tmpaix,1.0);

    alesolver_->Solve(Aii.EpetraMatrix(),aiy,aix,true);

//     Teuchos::RCP<Epetra_Vector> tmpagx = Teuchos::rcp(new Epetra_Vector(DomainMap(4)));
//     Agi.Multiply(false,*aiy,*tmpagx);
//     agx->Update(-1.0,*tmpagx,1.0);
//     Agg.Multiply(false,*agy,*tmpagx);
//     agx->Update(-1.0,*tmpagx,1.0);

//     // invert diagonal matrix CAST

//     len = agx->MyLength();
//     for (int i=0; i<len; ++i)
//     {
//       int NumEntries;
//       double *Values;
//       int *Indices;
//       int err = CAST.EpetraMatrix()->ExtractMyRowView(i,NumEntries,Values,Indices);
//       if (err)
//         dserror("ExtractMyRowView: err=%d",err);
//       if (NumEntries!=1)
//         dserror("diagonal matrix expected");
//       (*say)[Indices[0]] = -(*agx)[i] / Values[0];
//     }

    if (Comm().MyPID()==0)
      std::cout << ta.ElapsedTime() << std::flush;
  }

  {
    if (Comm().MyPID()==0)
      std::cout << "    fluid solve: " << std::flush;

    Epetra_Time tf(Comm());

//     CSF.Multiply(false,*sy,*sfy);

//     // invert diagonal matrix CFS

//     int len = sfy->MyLength();
//     for (int i=0; i<len; ++i)
//     {
//       int NumEntries;
//       double *Values;
//       int *Indices;
//       int err = CFS.EpetraMatrix()->ExtractMyRowView(i,NumEntries,Values,Indices);
//       if (err)
//         dserror("ExtractMyRowView: err=%d",err);
//       if (NumEntries!=1)
//         dserror("diagonal matrix expected");
//       (*fgy)[Indices[0]] = -(*sfy)[i] / Values[0];
//     }

//     Teuchos::RCP<Epetra_Vector> tmpfix = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));
//     Fig.Multiply(false,*fgy,*tmpfix);
//     fix->Update(-1.0,*tmpfix,1.0);

    fluidsolver_->Solve(Fii.EpetraMatrix(),fiy,fix,true);

//     Teuchos::RCP<Epetra_Vector> tmpfgx = Teuchos::rcp(new Epetra_Vector(DomainMap(2)));
//     Fgi.Multiply(false,*fiy,*tmpfgx);
//     fgx->Update(-1.0,*tmpfgx,1.0);
//     Fgg.Multiply(false,*fgy,*tmpfgx);
//     fgx->Update(-1.0,*tmpfgx,1.0);

//     // invert diagonal matrix CFST

//     len = fgx->MyLength();
//     for (int i=0; i<len; ++i)
//     {
//       int NumEntries;
//       double *Values;
//       int *Indices;
//       int err = CFST.EpetraMatrix()->ExtractMyRowView(i,NumEntries,Values,Indices);
//       if (err)
//         dserror("ExtractMyRowView: err=%d",err);
//       if (NumEntries!=1)
//         dserror("diagonal matrix expected");
//       (*sfy)[Indices[0]] = -(*fgx)[i] / Values[0];
//     }

    if (Comm().MyPID()==0)
      std::cout << tf.ElapsedTime() << "\n";
  }

  // build solution vector

  RangeExtractor().InsertVector(*sy,0,y);

  RangeExtractor().InsertVector(*fiy,1,y);
//   RangeExtractor().InsertVector(*fgy,2,y);

  RangeExtractor().InsertVector(*aiy,3,y);
//   RangeExtractor().InsertVector(*agy,4,y);

//   RangeExtractor().InsertVector(*sfy,5,y);
//   RangeExtractor().InsertVector(*say,6,y);

  return 0;
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
  const LINALG::SparseMatrix& aleInnerOp    = Matrix(3,3);

  structuresolver_->Setup(structInnerOp.EpetraMatrix());
  fluidsolver_    ->Setup(fluidInnerOp .EpetraMatrix());
  alesolver_      ->Setup(aleInnerOp   .EpetraMatrix());
}


#endif
