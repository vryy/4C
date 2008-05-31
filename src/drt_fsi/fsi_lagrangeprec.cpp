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
