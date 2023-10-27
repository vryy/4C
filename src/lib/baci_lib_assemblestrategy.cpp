/*----------------------------------------------------------------------*/
/*! \file

\brief Routines for the handing a collection of element matrices and vectors to
       the actual assembly calls into one global sparse matrix and global load vector

\level 0


*/
/*----------------------------------------------------------------------*/



#include "baci_lib_assemblestrategy.H"

#include "baci_lib_discret.H"
#include "baci_linalg_sparsematrix.H"
#include "baci_linalg_sparseoperator.H"
#include "baci_linalg_utils_sparse_algebra_assemble.H"

DRT::AssembleStrategy::AssembleStrategy(int firstdofset, int seconddofset,
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
    Teuchos::RCP<Epetra_Vector> systemvector3)
    : firstdofset_(firstdofset),
      seconddofset_(seconddofset),
      systemmatrix1_(systemmatrix1),
      systemmatrix2_(systemmatrix2),
      systemvector1_(systemvector1),
      systemvector2_(systemvector2),
      systemvector3_(systemvector3)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::AssembleStrategy::Zero()
{
  if (Assemblemat1())
  {
    systemmatrix1_->Zero();
  }
  if (Assemblemat2())
  {
    systemmatrix1_->Zero();
  }
  if (Assemblevec1())
  {
    systemvector1_->PutScalar(0.0);
  }
  if (Assemblevec2())
  {
    systemvector2_->PutScalar(0.0);
  }
  if (Assemblevec3())
  {
    systemvector3_->PutScalar(0.0);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::AssembleStrategy::Complete()
{
  if (Assemblemat1())
  {
    systemmatrix1_->Complete();
  }
  if (Assemblemat2())
  {
    systemmatrix2_->Complete();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::AssembleStrategy::ClearElementStorage(int rdim, int cdim)
{
  if (Assemblemat1())
  {
    if (elematrix1_.numRows() != rdim or elematrix1_.numCols() != cdim)
      elematrix1_.shape(rdim, cdim);
    else
      elematrix1_.putScalar(0.0);
  }
  if (Assemblemat2())
  {
    if (elematrix2_.numRows() != rdim or elematrix2_.numCols() != cdim)
      elematrix2_.shape(rdim, cdim);
    else
      elematrix2_.putScalar(0.0);
  }
  if (Assemblevec1())
  {
    if (elevector1_.length() != rdim)
      elevector1_.size(rdim);
    else
      elevector1_.putScalar(0.0);
  }
  if (Assemblevec2())
  {
    if (elevector2_.length() != rdim)
      elevector2_.size(rdim);
    else
      elevector2_.putScalar(0.0);
  }
  if (Assemblevec3())
  {
    if (elevector3_.length() != rdim)
      elevector3_.size(rdim);
    else
      elevector3_.putScalar(0.0);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::AssembleStrategy::Assemble(CORE::LINALG::SparseOperator& sysmat, int eid,
    const std::vector<int>& lmstride, const CORE::LINALG::SerialDenseMatrix& Aele,
    const std::vector<int>& lm, const std::vector<int>& lmowner)
{
  sysmat.Assemble(eid, lmstride, Aele, lm, lmowner);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::AssembleStrategy::Assemble(CORE::LINALG::SparseOperator& sysmat, int eid,
    const std::vector<int>& lmstride, const CORE::LINALG::SerialDenseMatrix& Aele,
    const std::vector<int>& lmrow, const std::vector<int>& lmrowowner,
    const std::vector<int>& lmcol)
{
  sysmat.Assemble(eid, lmstride, Aele, lmrow, lmrowowner, lmcol);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::AssembleStrategy::Assemble(
    CORE::LINALG::SparseOperator& sysmat, double val, int rgid, int cgid)
{
  sysmat.Assemble(val, rgid, cgid);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::AssembleStrategy::Assemble(Epetra_Vector& V, const CORE::LINALG::SerialDenseVector& Vele,
    const std::vector<int>& lm, const std::vector<int>& lmowner)
{
  CORE::LINALG::Assemble(V, Vele, lm, lmowner);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::AssembleStrategy::Assemble(Epetra_MultiVector& V, const int n,
    const CORE::LINALG::SerialDenseVector& Vele, const std::vector<int>& lm,
    const std::vector<int>& lmowner)
{
  CORE::LINALG::Assemble(V, n, Vele, lm, lmowner);
}
