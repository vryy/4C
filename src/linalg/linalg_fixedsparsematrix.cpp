/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation

\level 2

*----------------------------------------------------------------------*/

#include "linalg_fixedsparsematrix.H"
#include "../drt_lib/drt_dserror.H"


LINALG::FixedSparseMatrix::FixedSparseMatrix(Teuchos::RCP<const Epetra_Map> dbcmap)
    : dbcmap_(dbcmap)
{
}

void LINALG::FixedSparseMatrix::SetMatrix(Teuchos::RCP<Epetra_CrsGraph> graph)
{
#if 0
  // Always create a FE matrix here. It does not hurt and allows better assembly.
  sysmat_ = Teuchos::rcp( new Epetra_FECrsMatrix( Copy, *graph ) );
#endif
  sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(::Copy, *graph));
}

void LINALG::FixedSparseMatrix::Zero() { sysmat_->PutScalar(0.0); }

void LINALG::FixedSparseMatrix::Reset() { Zero(); }

void LINALG::FixedSparseMatrix::Assemble(int eid, const Epetra_SerialDenseMatrix& Aele,
    const std::vector<int>& lmrow, const std::vector<int>& lmrowowner,
    const std::vector<int>& lmcol)
{
  if (not sysmat_->Filled()) dserror("not filled");

  const unsigned lrowdim = lmrow.size();
  const unsigned lcoldim = lmcol.size();

  const int myrank = sysmat_->Comm().MyPID();
  const Epetra_Map& rowmap = sysmat_->RowMap();
  const Epetra_Map& colmap = sysmat_->ColMap();

  std::vector<double> values(lcoldim);
  std::vector<int> localcol(lcoldim);
  for (unsigned lcol = 0; lcol < lcoldim; ++lcol)
  {
    const int cgid = lmcol[lcol];
    localcol[lcol] = colmap.LID(cgid);
  }

  // loop rows of local matrix
  for (unsigned lrow = 0; lrow < lrowdim; ++lrow)
  {
    // check ownership of row
    if (lmrowowner[lrow] != myrank) continue;

    const int rgid = lmrow[lrow];

    // if we have a Dirichlet map check if this row is a Dirichlet row
    if (dbcmap_->MyGID(rgid)) continue;

    const int rlid = rowmap.LID(rgid);

    for (unsigned lcol = 0; lcol < lcoldim; ++lcol)
    {
      values[lcol] = Aele(lrow, lcol);
    }
    const int errone = sysmat_->SumIntoMyValues(rlid, lcoldim, &values[0], &localcol[0]);
    if (errone) dserror("Epetra_CrsMatrix::SumIntoMyValues returned error code %d", errone);
  }
}

void LINALG::FixedSparseMatrix::Assemble(double val, int rgid, int cgid)
{
  if (dbcmap_->MyGID(rgid)) dserror("no assembling to Dirichlet row");

  // SumIntoGlobalValues works for filled matrices as well!
  int errone = sysmat_->SumIntoGlobalValues(rgid, 1, &val, &cgid);
  if (errone) dserror("Epetra_CrsMatrix::SumIntoGlobalValues returned error code %d", errone);
}

void LINALG::FixedSparseMatrix::Complete()
{
  if (sysmat_->Filled()) return;
#if 0
  Epetra_FECrsMatrix * mat = dynamic_cast<Epetra_FECrsMatrix*>( &*sysmat_ );
  if ( mat!=NULL )
  {
    int err = mat->GlobalAssemble( false );
    if ( err )
      dserror( "Epetra_FECrsMatrix::GlobalAssemble failed: err=%d", err );
  }
#endif
  int err = sysmat_->FillComplete(true);
  if (err) dserror("Epetra_CrsMatrix::FillComplete(domain,range) returned err=%d", err);
}

void LINALG::FixedSparseMatrix::Complete(const Epetra_Map& domainmap, const Epetra_Map& rangemap)
{
  if (sysmat_->Filled()) return;
#if 0
  Epetra_FECrsMatrix * mat = dynamic_cast<Epetra_FECrsMatrix*>( &*sysmat_ );
  if ( mat!=NULL )
  {
    int err = mat->GlobalAssemble( false );
    if ( err )
      dserror( "Epetra_FECrsMatrix::GlobalAssemble failed: err=%d", err );
  }
#endif
  int err = sysmat_->FillComplete(domainmap, rangemap, true);
  if (err) dserror("Epetra_CrsMatrix::FillComplete(domain,range) returned err=%d", err);
}

void LINALG::FixedSparseMatrix::UnComplete() { dserror("UnComplete not supported"); }

void LINALG::FixedSparseMatrix::ApplyDirichlet(
    const Teuchos::RCP<const Epetra_Vector> dbctoggle, bool diagonalblock)
{
  dserror("Dirichlet map and toggle vector cannot be combined");
}

void LINALG::FixedSparseMatrix::ApplyDirichlet(const Epetra_Map& dbcmap, bool diagonalblock)
{
  if (not sysmat_->Filled()) dserror("expect filled matrix to apply dirichlet conditions");

#ifdef DEBUG
  if (not dbcmap.SameAs(*dbcmap_))
  {
    dserror("Dirichlet maps mismatch");
  }
#endif

  if (diagonalblock)
  {
    double v = 1.0;
    int numdbc = dbcmap.NumMyElements();
    int* dbc = dbcmap.MyGlobalElements();
    for (int i = 0; i < numdbc; ++i)
    {
      int row = dbc[i];
      int err = sysmat_->ReplaceGlobalValues(row, 1, &v, &row);
      if (err < 0) dserror("Epetra_CrsMatrix::ReplaceGlobalValues returned err=%d", err);
    }
  }
}
