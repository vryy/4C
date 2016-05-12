/*!----------------------------------------------------------------------
\file linalg_sparsematrixbase.cpp
\brief Implementation

<pre>
\level 0
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include "linalg_sparsematrixbase.H"
#include "linalg_utils.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrixBase::SetUseTranspose(bool UseTranspose)
{
  return sysmat_->SetUseTranspose(UseTranspose);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrixBase::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  return sysmat_->Apply(X,Y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrixBase::ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  return sysmat_->ApplyInverse(X,Y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool LINALG::SparseMatrixBase::UseTranspose() const
{
  return sysmat_->UseTranspose();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool LINALG::SparseMatrixBase::HasNormInf() const
{
  return sysmat_->HasNormInf();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Comm& LINALG::SparseMatrixBase::Comm() const
{
  return sysmat_->Comm();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map& LINALG::SparseMatrixBase::OperatorDomainMap() const
{
  return sysmat_->OperatorDomainMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map& LINALG::SparseMatrixBase::OperatorRangeMap() const
{
  return sysmat_->OperatorRangeMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrixBase::MaxNumEntries() const
{
  return sysmat_->MaxNumEntries();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double LINALG::SparseMatrixBase::NormInf() const
{
  return sysmat_->NormInf();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double LINALG::SparseMatrixBase::NormOne() const
{
  return sysmat_->NormOne();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double LINALG::SparseMatrixBase::NormFrobenius() const
{
  return sysmat_->NormFrobenius();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrixBase::Multiply(bool TransA, const Epetra_Vector &x, Epetra_Vector &y) const
{
  return sysmat_->Multiply(TransA,x,y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrixBase::Multiply(bool TransA, const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  return sysmat_->Multiply(TransA,X,Y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrixBase::LeftScale(const Epetra_Vector &x)
{
  return sysmat_->LeftScale(x);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrixBase::RightScale(const Epetra_Vector &x)
{
  return sysmat_->RightScale(x);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrixBase::PutScalar(double ScalarConstant)
{
  return sysmat_->PutScalar(ScalarConstant);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrixBase::Scale(double ScalarConstant)
{
  return sysmat_->Scale(ScalarConstant);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrixBase::ReplaceDiagonalValues(const Epetra_Vector &Diagonal)
{
  return sysmat_->ReplaceDiagonalValues(Diagonal);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrixBase::ExtractDiagonalCopy(Epetra_Vector &Diagonal) const
{
  return sysmat_->ExtractDiagonalCopy(Diagonal);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::SparseMatrixBase::Add(const LINALG::SparseOperator& A, const bool transposeA, const double scalarA, const double scalarB)
{
  A.AddOther(*this, transposeA, scalarA, scalarB);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::SparseMatrixBase::AddOther(LINALG::SparseMatrixBase& B, const bool transposeA, const double scalarA, const double scalarB) const
{
  //B.Add(*this, transposeA, scalarA, scalarB);
  LINALG::Add( *sysmat_, transposeA, scalarA, B, scalarB );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::SparseMatrixBase::AddOther(LINALG::BlockSparseMatrixBase& B, const bool transposeA, const double scalarA, const double scalarB) const
{
  dserror("BlockSparseMatrix and SparseMatrix cannot be added");
}
