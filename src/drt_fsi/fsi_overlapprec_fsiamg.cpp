#ifdef CCADISCRET

#include "fsi_overlapprec.H"
#include <Epetra_Time.h>



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::OverlappingBlockMatrixFSIAMG::OverlappingBlockMatrixFSIAMG(
                                                    const LINALG::MultiMapExtractor& maps,
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
  : OverlappingBlockMatrix(maps,
                           structure,
                           fluid,
                           ale,
                           structuresplit,
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
void FSI::OverlappingBlockMatrixFSIAMG::SetupPreconditioner()
{
#ifdef BLOCKMATRIXMERGE

  // this is really evil :)
  sparse_ = Merge();
  fluidsolver_->Setup(sparse_->EpetraMatrix());

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
  
  RCP<LINALG::MapExtractor> fsidofmapex = null;
  RCP<Epetra_Map>           irownodes = null;
  
  structuresolver_->Setup(structInnerOp.EpetraMatrix());
  fluidsolver_->Setup(fluidInnerOp.EpetraMatrix(),
                      fsidofmapex,
                      fluid_.Discretization(),
                      irownodes,
                      structuresplit_);
  if (constalesolver_==Teuchos::null)
    alesolver_->Setup(aleInnerOp.EpetraMatrix());
#endif
}




/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* FSI::OverlappingBlockMatrixFSIAMG::Label() const
{
  return "FSI::OverlappingBlockMatrix_FSIAMG";
}

#endif
