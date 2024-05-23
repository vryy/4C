/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of general 4C sparse matrix class

\level 0

*----------------------------------------------------------------------*/
#include "4C_linalg_sparsematrix.hpp"

#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_Transpose_RowMatrix.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <fstream>
#include <iterator>
#include <sstream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::SparseMatrix::SparseMatrix(
    Teuchos::RCP<Epetra_CrsGraph> crsgraph, Teuchos::RCP<CORE::LINALG::MultiMapExtractor> dbcmaps)
    : explicitdirichlet_(true), savegraph_(true), matrixtype_(CRS_MATRIX)
{
  sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(::Copy, *crsgraph));
  graph_ = crsgraph;
  dbcmaps_ = dbcmaps;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::SparseMatrix::SparseMatrix(const Epetra_Map& rowmap, const int npr,
    bool explicitdirichlet, bool savegraph, MatrixType matrixtype)
    : graph_(Teuchos::null),
      dbcmaps_(Teuchos::null),
      explicitdirichlet_(explicitdirichlet),
      savegraph_(savegraph),
      matrixtype_(matrixtype)
{
  if (!rowmap.UniqueGIDs()) FOUR_C_THROW("Row map is not unique");

  if (matrixtype_ == CRS_MATRIX)
    sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(::Copy, rowmap, npr, false));
  else if (matrixtype_ == FE_MATRIX)
    sysmat_ = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy, rowmap, npr, false));
  else
    FOUR_C_THROW("matrix type is not correct");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::SparseMatrix::SparseMatrix(const Epetra_Map& rowmap, std::vector<int>& numentries,
    bool explicitdirichlet, bool savegraph, MatrixType matrixtype)
    : graph_(Teuchos::null),
      dbcmaps_(Teuchos::null),
      explicitdirichlet_(explicitdirichlet),
      savegraph_(savegraph),
      matrixtype_(matrixtype)
{
  if (!rowmap.UniqueGIDs()) FOUR_C_THROW("Row map is not unique");

  if ((int)(numentries.size()) != rowmap.NumMyElements())
    FOUR_C_THROW("estimate for non zero entries per row does not match the size of row map");

  if (matrixtype_ == CRS_MATRIX)
    sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(::Copy, rowmap, numentries.data(), false));
  else if (matrixtype_ == FE_MATRIX)
    sysmat_ = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy, rowmap, numentries.data(), false));
  else
    FOUR_C_THROW("matrix type is not correct");
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::SparseMatrix::SparseMatrix(Teuchos::RCP<Epetra_CrsMatrix> matrix, DataAccess access,
    bool explicitdirichlet, bool savegraph, MatrixType matrixtype)
    : graph_(Teuchos::null),
      dbcmaps_(Teuchos::null),
      explicitdirichlet_(explicitdirichlet),
      savegraph_(savegraph),
      matrixtype_(matrixtype)
{
  if (access == Copy)
  {
    if (matrixtype_ == CRS_MATRIX)
      sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(*matrix));
    else if (matrixtype_ == FE_MATRIX)
    {
      sysmat_ = Teuchos::rcp(
          new Epetra_FECrsMatrix(*(Teuchos::rcp_dynamic_cast<Epetra_FECrsMatrix>(matrix))));
    }
    else
      FOUR_C_THROW("matrix type is not correct");
  }
  else
  {
    if (matrixtype_ == CRS_MATRIX)
      sysmat_ = matrix;
    else if (matrixtype_ == FE_MATRIX)
      sysmat_ = Teuchos::rcp_dynamic_cast<Epetra_FECrsMatrix>(matrix, true);
    else
      FOUR_C_THROW("matrix type is not correct");
  }

  if (sysmat_->Filled() and savegraph_)
  {
    graph_ = Teuchos::rcp(new Epetra_CrsGraph(sysmat_->Graph()));
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::SparseMatrix::SparseMatrix(const SparseMatrix& mat, DataAccess access)
    : CORE::LINALG::SparseMatrixBase(mat),
      explicitdirichlet_(mat.explicitdirichlet_),
      savegraph_(mat.savegraph_),
      matrixtype_(mat.matrixtype_)
{
  if (access == Copy)
  {
    // We do not care for exception proved code, so this is ok.
    *this = mat;
  }
  else
  {
    sysmat_ = mat.sysmat_;
    graph_ = mat.graph_;
    matrixtype_ = mat.matrixtype_;
    dbcmaps_ = mat.dbcmaps_;
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::SparseMatrix::SparseMatrix(
    const Epetra_Vector& diag, bool explicitdirichlet, bool savegraph, MatrixType matrixtype)
    : graph_(Teuchos::null),
      dbcmaps_(Teuchos::null),
      explicitdirichlet_(explicitdirichlet),
      savegraph_(savegraph),
      matrixtype_(matrixtype)
{
  int length = diag.Map().NumMyElements();
  Epetra_Map map(-1, length, diag.Map().MyGlobalElements(), diag.Map().IndexBase(), diag.Comm());
  if (!map.UniqueGIDs()) FOUR_C_THROW("Row map is not unique");

  if (matrixtype_ == CRS_MATRIX)
    sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(::Copy, map, map, 1, false));
  else if (matrixtype_ == FE_MATRIX)
    sysmat_ = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy, map, map, 1, false));
  else
    FOUR_C_THROW("matrix type is not correct");

  // sysmat_->FillComplete();

  for (int i = 0; i < length; ++i)
  {
    int gid = map.GID(i);
    Assemble(diag[i], gid, gid);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CORE::LINALG::SparseMatrix::Destroy(bool throw_exception)
{
  // delete first the epetra matrix object
  if (throw_exception and sysmat_.strong_count() > 1)
  {
    std::stringstream msg;
    msg << "Epetra_CrsMatrix cannot be finally deleted: The strong counter "
           "is still larger than 1. ( strong_count() = "
        << sysmat_.strong_count() << " )";
    FOUR_C_THROW(msg.str());
  }
  sysmat_ = Teuchos::null;

  // delete now also the matrix' graph
  if (throw_exception and graph_.strong_count() > 1)
  {
    std::stringstream msg;
    msg << "Epetra_CrsGraph cannot be finally deleted: The strong counter is "
           "still larger than 1. ( strong_count() = "
        << graph_.strong_count() << " )";
    FOUR_C_THROW(msg.str());
  }
  graph_ = Teuchos::null;

  // delete now also the matrix' graph
  if (throw_exception and dbcmaps_.strong_count() > 1)
  {
    std::stringstream msg;
    msg << "DBCMaps cannot be finally deleted: The strong counter is still "
           "larger than 1. ( strong_count() = "
        << dbcmaps_.strong_count() << " )";
    FOUR_C_THROW(msg.str());
  }
  dbcmaps_ = Teuchos::null;

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::SparseMatrix& CORE::LINALG::SparseMatrix::operator=(const SparseMatrix& mat)
{
  explicitdirichlet_ = mat.explicitdirichlet_;
  savegraph_ = mat.savegraph_;
  matrixtype_ = mat.matrixtype_;
  dbcmaps_ = mat.dbcmaps_;

  if (not mat.Filled())
  {
    // No communication. If just one processor fails, MPI will stop the other
    // ones as well.
    int nonzeros = mat.sysmat_->NumMyNonzeros();
    if (nonzeros > 0) FOUR_C_THROW("cannot copy non-filled matrix");
  }

  if (mat.Filled())
  {
    if (matrixtype_ == CRS_MATRIX)
      sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(*mat.sysmat_));
    else if (matrixtype_ == FE_MATRIX)
    {
      sysmat_ =
          Teuchos::rcp(new Epetra_FECrsMatrix(dynamic_cast<Epetra_FECrsMatrix&>(*mat.sysmat_)));
    }
    else
      FOUR_C_THROW("matrix type is not correct");
  }
  else
  {
    if (matrixtype_ == CRS_MATRIX)
      sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(::Copy, mat.RowMap(), 0, false));
    else if (matrixtype_ == FE_MATRIX)
      sysmat_ = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy, mat.RowMap(), 0, false));
    else
      FOUR_C_THROW("matrix type is not correct");
  }

  if (mat.graph_ != Teuchos::null)
    graph_ = Teuchos::rcp(new Epetra_CrsGraph(*mat.graph_));
  else
    graph_ = Teuchos::null;

  return *this;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::Assign(DataAccess access, const SparseMatrix& mat)
{
  if (access == Copy)
  {
    // We do not care for exception proved code, so this is ok.
    *this = mat;
  }
  else
  {
    sysmat_ = mat.sysmat_;
    graph_ = mat.graph_;
    explicitdirichlet_ = mat.explicitdirichlet_;
    savegraph_ = mat.savegraph_;
    matrixtype_ = mat.matrixtype_;
    dbcmaps_ = mat.dbcmaps_;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::Zero()
{
  if (graph_ == Teuchos::null)
  {
    if (Filled() && !explicitdirichlet_)
      sysmat_->PutScalar(0.);
    else
      Reset();
  }
  else
  {
    const Epetra_Map domainmap = sysmat_->DomainMap();
    const Epetra_Map rangemap = sysmat_->RangeMap();
    // Remove old matrix before creating a new one so we do not have old and
    // new matrix in memory at the same time!
    sysmat_ = Teuchos::null;
    if (matrixtype_ == CRS_MATRIX)
      sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(::Copy, *graph_));
    else if (matrixtype_ == FE_MATRIX)
      sysmat_ = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy, *graph_));
    else
      FOUR_C_THROW("matrix type is not correct");

    sysmat_->FillComplete(domainmap, rangemap);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::Reset()
{
  const Epetra_Map rowmap = sysmat_->RowMap();
  std::vector<int> numentries(rowmap.NumMyElements());

  const Epetra_CrsGraph& graph = sysmat_->Graph();

  if (Filled())
  {
    for (std::size_t i = 0; i < numentries.size(); ++i)
    {
      int* indices;
      int err = graph.ExtractMyRowView(i, numentries[i], indices);
      if (err != 0) FOUR_C_THROW("ExtractMyRowView failed");
    }
  }
  else
  {
    // use information about number of allocated entries not to fall back to matrix with zero size
    // otherwise assembly would be extremely expensive!
    for (std::size_t i = 0; i < numentries.size(); ++i)
    {
      numentries[i] = graph.NumAllocatedMyIndices(i);
    }
  }
  // Remove old matrix before creating a new one so we do not have old and
  // new matrix in memory at the same time!
  sysmat_ = Teuchos::null;
  if (matrixtype_ == CRS_MATRIX)
    sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(::Copy, rowmap, numentries.data(), false));
  else if (matrixtype_ == FE_MATRIX)
    sysmat_ = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy, rowmap, numentries.data(), false));
  else
    FOUR_C_THROW("matrix type is not correct");

  graph_ = Teuchos::null;
  dbcmaps_ = Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::Assemble(int eid, const std::vector<int>& lmstride,
    const CORE::LINALG::SerialDenseMatrix& Aele, const std::vector<int>& lmrow,
    const std::vector<int>& lmrowowner, const std::vector<int>& lmcol)
{
  const int lrowdim = (int)lmrow.size();
  const int lcoldim = (int)lmcol.size();
  // allow Aele to provide entries past the end of lmrow and lmcol that are
  // not used here, therefore check only for ">" rather than "!="
  if (lrowdim != (int)lmrowowner.size() || lrowdim > Aele.numRows() || lcoldim > Aele.numCols())
    FOUR_C_THROW("Mismatch in dimensions");

  const int myrank = sysmat_->Comm().MyPID();
  const Epetra_Map& rowmap = sysmat_->RowMap();
  const Epetra_Map& colmap = sysmat_->ColMap();

  //-----------------------------------------------------------------------------------
  auto& A = (CORE::LINALG::SerialDenseMatrix&)Aele;
  if (sysmat_->Filled())  // assembly in local indices
  {
#ifdef FOUR_C_ENABLE_ASSERTIONS
    // There is the case of nodes without dofs (XFEM).
    // If no row dofs are present on this proc, there is nothing to assemble.
    // However, the subsequent check for coldofs (in DEBUG mode) would incorrectly fail.
    bool doit = false;
    for (int lrow = 0; lrow < lrowdim; ++lrow)
      if (lmrowowner[lrow] == myrank)
      {
        doit = true;
        break;
      }
    if (!doit) return;
#endif

    std::vector<int> localcol(lcoldim);
    for (int lcol = 0; lcol < lcoldim; ++lcol)
    {
      const int cgid = lmcol[lcol];
      localcol[lcol] = colmap.LID(cgid);
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (localcol[lcol] < 0) FOUR_C_THROW("Sparse matrix A does not have global column %d", cgid);
#endif
    }

    // loop rows of local matrix
    for (int lrow = 0; lrow < lrowdim; ++lrow)
    {
      // check ownership of row
      if (lmrowowner[lrow] != myrank) continue;

      // check whether I have that global row
      const int rgid = lmrow[lrow];

      // if we have a Dirichlet map check if this row is a Dirichlet row
      if (dbcmaps_ != Teuchos::null and dbcmaps_->Map(1)->MyGID(rgid)) continue;

      const int rlid = rowmap.LID(rgid);

#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (rlid < 0) FOUR_C_THROW("Sparse matrix A does not have global row %d", rgid);
#endif
      int length;
      double* valview;
      int* indices;
#ifdef FOUR_C_ENABLE_ASSERTIONS
      int err =
#endif
          sysmat_->ExtractMyRowView(rlid, length, valview, indices);
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (err) FOUR_C_THROW("Epetra_CrsMatrix::ExtractMyRowView returned error code %d", err);
#endif
      const int numnode = (int)lmstride.size();
      int dofcount = 0;
      int pos = 0;
      for (int node = 0; node < numnode; ++node)
      {
        // check if 'pos' already points to the correct location before the binary search
        if (pos >= length || indices[pos] != localcol[dofcount])
        {
          int* loc = std::lower_bound(indices, indices + length, localcol[dofcount]);
#ifdef FOUR_C_ENABLE_ASSERTIONS
          if (*loc != localcol[dofcount])
            FOUR_C_THROW("Cannot find local column entry %d", localcol[dofcount]);
#endif
          pos = loc - indices;
        }
        const int stride = lmstride[node];
        // test continuity of data in sparsematrix
        bool reachedlength = false;
        bool continuous = true;
        if (stride + pos > length)
          continuous = false;
        else
        {
          for (int j = 1; j < stride; ++j)
          {
            if (indices[pos + j] == localcol[dofcount + j])
              continue;
            else
            {
              continuous = false;
              break;
            }
          }
        }

        if (continuous)
        {
          for (int j = 0; j < stride; ++j)
          {
            valview[pos++] += Aele(lrow, dofcount++);
            if (dofcount == lcoldim)
            {
              reachedlength = true;
              break;
            }
          }
        }
        else
        {
          for (int j = 0; j < stride; ++j)
          {
#ifdef FOUR_C_ENABLE_ASSERTIONS
            const int errone =
#endif
                sysmat_->SumIntoMyValues(rlid, 1, &A(lrow, dofcount), &localcol[dofcount]);
#ifdef FOUR_C_ENABLE_ASSERTIONS
            if (errone)
            {
              printf("Dimensions of A: %d x %d\n", A.numRows(), A.numCols());
              for (unsigned k = 0; k < lmstride.size(); ++k)
                printf("lmstride[%d] %d\n", k, lmstride[k]);
              FOUR_C_THROW(
                  "Epetra_CrsMatrix::SumIntoMyValues returned error code %d\nrlid %d localcol[%d] "
                  "%d dofcount %d length %d stride %d j %d node %d numnode %d",
                  errone, rlid, dofcount, localcol[dofcount], dofcount, length, stride, j, node,
                  numnode);
            }
#endif
            dofcount++;
            if (dofcount == lcoldim)
            {
              reachedlength = true;
              break;
            }
          }
        }
        if (reachedlength) break;
      }  // for (int node=0; node<numnode; ++node)
    }    // for (int lrow=0; lrow<ldim; ++lrow)
  }
  //-----------------------------------------------------------------------------------
  else  // assembly in global indices
  {
    // loop rows of local matrix
    for (int lrow = 0; lrow < lrowdim; ++lrow)
    {
      // check ownership of row
      if (lmrowowner[lrow] != myrank) continue;

      // check whether I have that global row
      const int rgid = lmrow[lrow];
      // #ifdef FOUR_C_ENABLE_ASSERTIONS
      if (!rowmap.MyGID(rgid)) FOUR_C_THROW("Proc %d does not have global row %d", myrank, rgid);
      // #endif

      // if we have a Dirichlet map check if this row is a Dirichlet row
      if (dbcmaps_ != Teuchos::null and dbcmaps_->Map(1)->MyGID(rgid)) continue;

      for (int lcol = 0; lcol < lcoldim; ++lcol)
      {
        int cgid = lmcol[lcol];
        // Now that we do not rebuild the sparse mask in each step, we
        // are bound to assemble the whole thing. Zeros included.
        const int errone = sysmat_->SumIntoGlobalValues(rgid, 1, &A(lrow, lcol), &cgid);
        if (errone > 0)
        {
          const int errtwo = sysmat_->InsertGlobalValues(rgid, 1, &A(lrow, lcol), &cgid);
          if (errtwo < 0)
            FOUR_C_THROW("Epetra_CrsMatrix::InsertGlobalValues returned error code %d", errtwo);
        }
        else if (errone)
          FOUR_C_THROW("Epetra_CrsMatrix::SumIntoGlobalValues returned error code %d", errone);
      }  // for (int lcol=0; lcol<lcoldim; ++lcol)
    }    // for (int lrow=0; lrow<lrowdim; ++lrow)
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::Assemble(int eid, const CORE::LINALG::SerialDenseMatrix& Aele,
    const std::vector<int>& lmrow, const std::vector<int>& lmrowowner,
    const std::vector<int>& lmcol)
{
  const int lrowdim = (int)lmrow.size();
  const int lcoldim = (int)lmcol.size();
  // allow Aele to provide entries past the end of lmrow and lmcol that are
  // not used here, therefore check only for ">" rather than "!="
  if (lrowdim != (int)lmrowowner.size() || lrowdim > Aele.numRows() || lcoldim > Aele.numCols())
    FOUR_C_THROW("Mismatch in dimensions");

  const int myrank = sysmat_->Comm().MyPID();
  const Epetra_Map& rowmap = sysmat_->RowMap();
  const Epetra_Map& colmap = sysmat_->ColMap();

  if (sysmat_->Filled())
  {
#ifdef FOUR_C_ENABLE_ASSERTIONS
    // There is the case of nodes without dofs (XFEM).
    // If no row dofs are present on this proc, their is nothing to assemble.
    // However, the subsequent check for coldofs (in DEBUG mode) would incorrectly fail.
    bool doit = false;
    for (int lrow = 0; lrow < lrowdim; ++lrow)
      if (lmrowowner[lrow] == myrank)
      {
        doit = true;
        break;
      }
    if (!doit) return;
#endif

    std::vector<double> values(lcoldim);
    std::vector<int> localcol(lcoldim);
    for (int lcol = 0; lcol < lcoldim; ++lcol)
    {
      const int cgid = lmcol[lcol];
      localcol[lcol] = colmap.LID(cgid);
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (localcol[lcol] < 0) FOUR_C_THROW("Sparse matrix A does not have global column %d", cgid);
#endif
    }

    // loop rows of local matrix
    for (int lrow = 0; lrow < lrowdim; ++lrow)
    {
      // check ownership of row
      if (lmrowowner[lrow] != myrank) continue;

      // check whether I have that global row
      const int rgid = lmrow[lrow];

      // if we have a Dirichlet map check if this row is a Dirichlet row
      if (dbcmaps_ != Teuchos::null and dbcmaps_->Map(1)->MyGID(rgid)) continue;

      const int rlid = rowmap.LID(rgid);
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (rlid < 0) FOUR_C_THROW("Sparse matrix A does not have global row %d", rgid);
#endif

      for (int lcol = 0; lcol < lcoldim; ++lcol)
      {
        values[lcol] = Aele(lrow, lcol);
      }
      const int errone = sysmat_->SumIntoMyValues(rlid, lcoldim, values.data(), localcol.data());
      if (errone) FOUR_C_THROW("Epetra_CrsMatrix::SumIntoMyValues returned error code %d", errone);
    }  // for (int lrow=0; lrow<ldim; ++lrow)
  }
  else
  {
    // loop rows of local matrix
    for (int lrow = 0; lrow < lrowdim; ++lrow)
    {
      // check ownership of row
      if (lmrowowner[lrow] != myrank) continue;

      // check whether I have that global row
      const int rgid = lmrow[lrow];
#ifdef FOUR_C_ENABLE_ASSERTIONS
      if (!rowmap.MyGID(rgid)) FOUR_C_THROW("Proc %d does not have global row %d", myrank, rgid);
#endif

      // if we have a Dirichlet map check if this row is a Dirichlet row
      if (dbcmaps_ != Teuchos::null and dbcmaps_->Map(1)->MyGID(rgid)) continue;

      for (int lcol = 0; lcol < lcoldim; ++lcol)
      {
        double val = Aele(lrow, lcol);
        int cgid = lmcol[lcol];

        // Now that we do not rebuild the sparse mask in each step, we
        // are bound to assemble the whole thing. Zeros included.
        const int errone = sysmat_->SumIntoGlobalValues(rgid, 1, &val, &cgid);
        if (errone > 0)
        {
          const int errtwo = sysmat_->InsertGlobalValues(rgid, 1, &val, &cgid);
          if (errtwo < 0)
            FOUR_C_THROW("Epetra_CrsMatrix::InsertGlobalValues returned error code %d", errtwo);
        }
        else if (errone)
          FOUR_C_THROW("Epetra_CrsMatrix::SumIntoGlobalValues returned error code %d", errone);
      }  // for (int lcol=0; lcol<lcoldim; ++lcol)
    }    // for (int lrow=0; lrow<lrowdim; ++lrow)
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::FEAssemble(const CORE::LINALG::SerialDenseMatrix& Aele,
    const std::vector<int>& lmrow, const std::vector<int>& lmrowowner,
    const std::vector<int>& lmcol)
{
  const int lrowdim = static_cast<int>(lmrow.size());
  const int lcoldim = static_cast<int>(lmcol.size());

  // allow Aele to provide entries past the end of lmrow and lmcol that are
  // not used here, therefore check only for ">" rather than "!="
  if (lrowdim != (int)lmrowowner.size() || lrowdim > Aele.numRows() || lcoldim > Aele.numCols())
    FOUR_C_THROW("Mismatch in dimensions");

  Teuchos::RCP<Epetra_FECrsMatrix> fe_mat =
      Teuchos::rcp_dynamic_cast<Epetra_FECrsMatrix>(sysmat_, true);
  const int myrank = fe_mat->Comm().MyPID();

  // loop rows of local matrix
  for (int lrow = 0; lrow < lrowdim; ++lrow)
  {
    // check ownership of row
    if (lmrowowner[lrow] != myrank) continue;

    const int rgid = lmrow[lrow];

    for (int lcol = 0; lcol < lcoldim; ++lcol)
    {
      double val = Aele(lrow, lcol);
      const int cgid = lmcol[lcol];
      FEAssemble(val, rgid, cgid);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::FEAssemble(const CORE::LINALG::SerialDenseMatrix& Aele,
    const std::vector<int>& lmrow, const std::vector<int>& lmcol)
{
  const int lrowdim = static_cast<int>(lmrow.size());
  const int lcoldim = static_cast<int>(lmcol.size());
  // allow Aele to provide entries past the end of lmrow and lmcol that are
  // not used here, therefore check only for ">" rather than "!="
  if (lrowdim > Aele.numRows() || lcoldim > Aele.numCols()) FOUR_C_THROW("Mismatch in dimensions");

  Teuchos::rcp_dynamic_cast<Epetra_FECrsMatrix>(sysmat_, true);

  // loop rows of local matrix
  for (int lrow = 0; lrow < lrowdim; ++lrow)
  {
    const int rgid = lmrow[lrow];

    for (int lcol = 0; lcol < lcoldim; ++lcol)
    {
      double val = Aele(lrow, lcol);
      const int cgid = lmcol[lcol];
      FEAssemble(val, rgid, cgid);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::Assemble(double val, int rgid, int cgid)
{
  if (dbcmaps_ != Teuchos::null and dbcmaps_->Map(1)->MyGID(rgid))
    FOUR_C_THROW("no assembling to Dirichlet row");

  // SumIntoGlobalValues works for filled matrices as well!
  int errone = sysmat_->SumIntoGlobalValues(rgid, 1, &val, &cgid);
  if (errone > 0)
  {
    int errtwo = sysmat_->InsertGlobalValues(rgid, 1, &val, &cgid);
    if (errtwo < 0)
      FOUR_C_THROW("Epetra_CrsMatrix::InsertGlobalValues returned error code %d", errtwo);
  }
  else if (errone)
    FOUR_C_THROW("Epetra_CrsMatrix::SumIntoGlobalValues returned error code %d", errone);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::SetValue(double val, int rgid, int cgid)
{
  if (dbcmaps_ != Teuchos::null and dbcmaps_->Map(1)->MyGID(rgid))
    FOUR_C_THROW("no assembling to Dirichlet row");

  int errone = sysmat_->ReplaceGlobalValues(rgid, 1, &val, &cgid);
  if (errone > 0)
  {
    int errtwo = sysmat_->InsertGlobalValues(rgid, 1, &val, &cgid);
    if (errtwo > 0)
      FOUR_C_THROW("Epetra_CrsMatrix::InsertGlobalValues returned error code %d", errtwo);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::FEAssemble(double val, int rgid, int cgid)
{
  // SumIntoGlobalValues works for filled matrices as well!
  int errone = (Teuchos::rcp_dynamic_cast<Epetra_FECrsMatrix>(sysmat_, true))
                   ->SumIntoGlobalValues(1, &rgid, 1, &cgid, &val);
  // if value not already present , error > 0 then use insert method
  if (errone > 0 and not Filled())
  {
    int errtwo = (Teuchos::rcp_dynamic_cast<Epetra_FECrsMatrix>(sysmat_, true))
                     ->InsertGlobalValues(1, &rgid, 1, &cgid, &val);
    if (errtwo < 0)
      FOUR_C_THROW("Epetra_FECrsMatrix::InsertGlobalValues returned error code %d", errtwo);
  }
  else if (errone < 0)
    FOUR_C_THROW("Epetra_FECrsMatrix::SumIntoGlobalValues returned error code %d", errone);
}


/*----------------------------------------------------------------------*
 |  FillComplete a matrix  (public)                          mwgee 12/06|
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::Complete(bool enforce_complete)
{
  TEUCHOS_FUNC_TIME_MONITOR("CORE::LINALG::SparseMatrix::Complete");

  // for FE_Matrix we need to gather non-local entries, independent whether matrix is filled or not
  if (matrixtype_ == FE_MATRIX)
  {
    // false indicates here that FillComplete() is not called within GlobalAssemble()
    int err = (Teuchos::rcp_dynamic_cast<Epetra_FECrsMatrix>(sysmat_))->GlobalAssemble(false);
    if (err) FOUR_C_THROW("Epetra_FECrsMatrix::GlobalAssemble() returned err=%d", err);
  }

  if (sysmat_->Filled() and not enforce_complete) return;

  int err = sysmat_->FillComplete(true);
  if (err) FOUR_C_THROW("Epetra_CrsMatrix::FillComplete(domain,range) returned err=%d", err);

  // keep mask for further use
  if (savegraph_ and graph_ == Teuchos::null)
  {
    graph_ = Teuchos::rcp(new Epetra_CrsGraph(sysmat_->Graph()));
  }
}


/*----------------------------------------------------------------------*
 |  FillComplete a matrix  (public)                          mwgee 01/08|
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::Complete(
    const Epetra_Map& domainmap, const Epetra_Map& rangemap, bool enforce_complete)
{
  TEUCHOS_FUNC_TIME_MONITOR("CORE::LINALG::SparseMatrix::Complete(domain,range)");

  // for FE_Matrix we need to gather non-local entries, independent whether matrix is filled or not
  if (matrixtype_ == FE_MATRIX)
  {
    // false indicates here that FillComplete() is not called within GlobalAssemble()
    int err = (Teuchos::rcp_dynamic_cast<Epetra_FECrsMatrix>(sysmat_))
                  ->GlobalAssemble(domainmap, rangemap, false);
    if (err) FOUR_C_THROW("Epetra_FECrsMatrix::GlobalAssemble() returned err=%d", err);
  }

  if (sysmat_->Filled() and not enforce_complete) return;

  int err = 1;
  if (enforce_complete and sysmat_->Filled())
    err = sysmat_->ExpertStaticFillComplete(domainmap, rangemap);
  else
    err = sysmat_->FillComplete(domainmap, rangemap, true);

  if (err) FOUR_C_THROW("Epetra_CrsMatrix::FillComplete(domain,range) returned err=%d", err);

  // keep mask for further use
  if (savegraph_ and graph_ == Teuchos::null)
  {
    graph_ = Teuchos::rcp(new Epetra_CrsGraph(sysmat_->Graph()));
  }
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::UnComplete()
{
  TEUCHOS_FUNC_TIME_MONITOR("CORE::LINALG::SparseMatrix::UnComplete");

  if (not Filled()) return;

  const Epetra_CrsGraph& graph = sysmat_->Graph();
  std::vector<int> nonzeros(graph.NumMyRows());
  for (std::size_t i = 0; i < nonzeros.size(); ++i)
  {
    nonzeros[i] = graph.NumMyIndices(i);
  }

  const Epetra_Map& rowmap = sysmat_->RowMap();
  const Epetra_Map& colmap = sysmat_->ColMap();
  int elements = rowmap.NumMyElements();

  Teuchos::RCP<Epetra_CrsMatrix> mat = Teuchos::null;
  if (matrixtype_ == CRS_MATRIX)
    mat = Teuchos::rcp(new Epetra_CrsMatrix(::Copy, rowmap, nonzeros.data(), false));
  else if (matrixtype_ == FE_MATRIX)
    mat = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy, rowmap, nonzeros.data(), false));
  else
    FOUR_C_THROW("matrix type is not correct");

  nonzeros.clear();
  for (int i = 0; i < elements; ++i)
  {
    int NumEntries;
    double* Values;
    int* Indices;
    // if matrix is filled, global assembly was called already and all nonlocal values are
    // distributed
    int err = sysmat_->ExtractMyRowView(i, NumEntries, Values, Indices);
    if (err) FOUR_C_THROW("ExtractMyRowView err=%d", err);
    std::vector<int> idx(NumEntries);
    for (int c = 0; c < NumEntries; ++c)
    {
      idx[c] = colmap.GID(Indices[c]);
      FOUR_C_ASSERT(idx[c] != -1, "illegal gid");
    }
    int rowgid = rowmap.GID(i);
    CORE::LINALG::InsertGlobalValues(mat, rowgid, NumEntries, Values, idx.data());
  }
  sysmat_ = mat;
  graph_ = Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  Apply dirichlet conditions  (public)                     mwgee 02/07|
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::ApplyDirichlet(const Epetra_Vector& dbctoggle, bool diagonalblock)
{
  // if matrix is filled, global assembly was called already and all nonlocal values are
  // distributed
  if (not Filled()) FOUR_C_THROW("expect filled matrix to apply dirichlet conditions");

  if (dbcmaps_ != Teuchos::null)
  {
    FOUR_C_THROW("Dirichlet map and toggle vector cannot be combined");
  }

  if (explicitdirichlet_)
  {
    // Save graph of original matrix if not done already.
    // This will never happen as the matrix is guaranteed to be filled. But to
    // make the code more explicit...
    if (savegraph_ and graph_ == Teuchos::null)
    {
      graph_ = Teuchos::rcp(new Epetra_CrsGraph(sysmat_->Graph()));
      if (not graph_->Filled()) FOUR_C_THROW("got unfilled graph from filled matrix");
    }

    // allocate a new matrix and copy all rows that are not dirichlet
    const Epetra_Map& rowmap = sysmat_->RowMap();
    const int nummyrows = sysmat_->NumMyRows();
    const int maxnumentries = sysmat_->MaxNumEntries();

    Teuchos::RCP<Epetra_CrsMatrix> Anew = Teuchos::null;
    if (matrixtype_ == CRS_MATRIX)
      Anew = Teuchos::rcp(new Epetra_CrsMatrix(::Copy, rowmap, maxnumentries, false));
    else if (matrixtype_ == FE_MATRIX)
      Anew = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy, rowmap, maxnumentries, false));
    else
      FOUR_C_THROW("matrix type is not correct");

    std::vector<int> indices(maxnumentries, 0);
    std::vector<double> values(maxnumentries, 0.0);
    for (int i = 0; i < nummyrows; ++i)
    {
      int row = sysmat_->GRID(i);
      if (dbctoggle[i] != 1.0)
      {
        int numentries;
#ifdef FOUR_C_ENABLE_ASSERTIONS
        int err = sysmat_->ExtractGlobalRowCopy(
            row, maxnumentries, numentries, values.data(), indices.data());
        if (err < 0) FOUR_C_THROW("Epetra_CrsMatrix::ExtractGlobalRowCopy returned err=%d", err);
#else
        sysmat_->ExtractGlobalRowCopy(
            row, maxnumentries, numentries, values.data(), indices.data());
#endif
          // this is also ok for FE matrices, because fill complete was called on sysmat and the
          // globalAssemble method was called already
#ifdef FOUR_C_ENABLE_ASSERTIONS
        err = Anew->InsertGlobalValues(row, numentries, values.data(), indices.data());
        if (err < 0) FOUR_C_THROW("Epetra_CrsMatrix::InsertGlobalValues returned err=%d", err);
#else
        Anew->InsertGlobalValues(row, numentries, values.data(), indices.data());
#endif
      }
      else
      {
        double v;
        if (diagonalblock)
          v = 1.0;
        else
          v = 0.0;
#ifdef FOUR_C_ENABLE_ASSERTIONS
        int err = Anew->InsertGlobalValues(row, 1, &v, &row);
        if (err < 0) FOUR_C_THROW("Epetra_CrsMatrix::InsertGlobalValues returned err=%d", err);
#else
        Anew->InsertGlobalValues(row, 1, &v, &row);
#endif
      }
    }
    sysmat_ = Anew;
    Complete();
  }
  else
  {
    const int nummyrows = sysmat_->NumMyRows();
    for (int i = 0; i < nummyrows; ++i)
    {
      if (dbctoggle[i] == 1.0)
      {
        int* indexOffset;
        int* indices;
        double* values;
#ifdef FOUR_C_ENABLE_ASSERTIONS
        int err = sysmat_->ExtractCrsDataPointers(indexOffset, indices, values);
        if (err < 0) FOUR_C_THROW("Epetra_CrsMatrix::ExtractCrsDataPointers returned err=%d", err);
#else
        sysmat_->ExtractCrsDataPointers(indexOffset, indices, values);
#endif
        // zero row
        memset(&values[indexOffset[i]], 0, (indexOffset[i + 1] - indexOffset[i]) * sizeof(double));

        if (diagonalblock)
        {
          double one = 1.0;
#ifdef FOUR_C_ENABLE_ASSERTIONS
          int err = sysmat_->SumIntoMyValues(i, 1, &one, &i);
          if (err < 0) FOUR_C_THROW("Epetra_CrsMatrix::SumIntoMyValues returned err=%d", err);
#else
          sysmat_->SumIntoMyValues(i, 1, &one, &i);
#endif
        }
      }
    }
  }
}


/*----------------------------------------------------------------------*
 |  Apply dirichlet conditions  (public)                     mwgee 02/07|
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::ApplyDirichlet(const Epetra_Map& dbctoggle, bool diagonalblock)
{
  if (not Filled()) FOUR_C_THROW("expect filled matrix to apply dirichlet conditions");

  if (dbcmaps_ != Teuchos::null)
  {
#ifdef FOUR_C_ENABLE_ASSERTIONS
    if (not dbctoggle.SameAs(*dbcmaps_->Map(1)))
    {
      FOUR_C_THROW("Dirichlet maps mismatch");
    }
#endif
    if (diagonalblock)
    {
      double v = 1.0;
      int numdbc = dbctoggle.NumMyElements();
      int* dbc = dbctoggle.MyGlobalElements();
      for (int i = 0; i < numdbc; ++i)
      {
        int row = dbc[i];
        int err = sysmat_->ReplaceGlobalValues(row, 1, &v, &row);
        if (err < 0) FOUR_C_THROW("Epetra_CrsMatrix::ReplaceGlobalValues returned err=%d", err);
      }
    }
    return;
  }

  if (explicitdirichlet_)
  {
    // Save graph of original matrix if not done already.
    // This will never happen as the matrix is guaranteed to be filled. But to
    // make the code more explicit...
    if (savegraph_ and graph_ == Teuchos::null)
    {
      graph_ = Teuchos::rcp(new Epetra_CrsGraph(sysmat_->Graph()));
      if (not graph_->Filled()) FOUR_C_THROW("got unfilled graph from filled matrix");
    }

    // allocate a new matrix and copy all rows that are not dirichlet
    const Epetra_Map& rowmap = sysmat_->RowMap();
    const int nummyrows = sysmat_->NumMyRows();
    const int maxnumentries = sysmat_->MaxNumEntries();

    // Teuchos::RCP<Epetra_CrsMatrix> Anew = Teuchos::rcp(new
    // Epetra_CrsMatrix(Copy,rowmap,maxnumentries,false));

    Teuchos::RCP<Epetra_CrsMatrix> Anew = Teuchos::null;
    if (matrixtype_ == CRS_MATRIX)
      Anew = Teuchos::rcp(new Epetra_CrsMatrix(::Copy, rowmap, maxnumentries, false));
    else if (matrixtype_ == FE_MATRIX)
      Anew = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy, rowmap, maxnumentries, false));
    else
      FOUR_C_THROW("matrix type is not correct");

    std::vector<int> indices(maxnumentries, 0);
    std::vector<double> values(maxnumentries, 0.0);
    for (int i = 0; i < nummyrows; ++i)
    {
      int row = sysmat_->GRID(i);
      if (not dbctoggle.MyGID(row))
      {
        int numentries;
#ifdef FOUR_C_ENABLE_ASSERTIONS
        int err = sysmat_->ExtractGlobalRowCopy(
            row, maxnumentries, numentries, values.data(), indices.data());
        if (err < 0) FOUR_C_THROW("Epetra_CrsMatrix::ExtractGlobalRowCopy returned err=%d", err);
#else
        sysmat_->ExtractGlobalRowCopy(
            row, maxnumentries, numentries, values.data(), indices.data());
#endif
          // this is also ok for FE matrices, because fill complete was called on sysmat and the
          // globalAssemble method was called already
#ifdef FOUR_C_ENABLE_ASSERTIONS
        err = Anew->InsertGlobalValues(row, numentries, values.data(), indices.data());
        if (err < 0) FOUR_C_THROW("Epetra_CrsMatrix::InsertGlobalValues returned err=%d", err);
#else
        Anew->InsertGlobalValues(row, numentries, values.data(), indices.data());
#endif
      }
      else
      {
        if (diagonalblock)
        {
          double v = 1.0;
#ifdef FOUR_C_ENABLE_ASSERTIONS
          int err = Anew->InsertGlobalValues(row, 1, &v, &row);
          if (err < 0) FOUR_C_THROW("Epetra_CrsMatrix::InsertGlobalValues returned err=%d", err);
#else
          Anew->InsertGlobalValues(row, 1, &v, &row);
#endif
        }
      }
    }
    Epetra_Map rangemap = sysmat_->RangeMap();
    Epetra_Map domainmap = sysmat_->DomainMap();
    sysmat_ = Anew;
    Complete(domainmap, rangemap);
  }
  else
  {
    const int nummyrows = sysmat_->NumMyRows();
    for (int i = 0; i < nummyrows; ++i)
    {
      int row = sysmat_->GRID(i);
      if (dbctoggle.MyGID(row))
      {
        int* indexOffset;
        int* indices;
        double* values;
#ifdef FOUR_C_ENABLE_ASSERTIONS
        int err = sysmat_->ExtractCrsDataPointers(indexOffset, indices, values);
        if (err < 0) FOUR_C_THROW("Epetra_CrsMatrix::ExtractCrsDataPointers returned err=%d", err);
#else
        sysmat_->ExtractCrsDataPointers(indexOffset, indices, values);
#endif
        // zero row
        memset(&values[indexOffset[i]], 0, (indexOffset[i + 1] - indexOffset[i]) * sizeof(double));

        if (diagonalblock)
        {
          double one = 1.0;
#ifdef FOUR_C_ENABLE_ASSERTIONS
          err = sysmat_->SumIntoMyValues(i, 1, &one, &i);
          if (err < 0) FOUR_C_THROW("Epetra_CrsMatrix::SumIntoMyValues returned err=%d", err);
#else
          sysmat_->SumIntoMyValues(i, 1, &one, &i);
#endif
        }
      }
    }
  }
}


/*----------------------------------------------------------------------*
 |  Apply dirichlet conditions  (public)                     mwgee 02/07|
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::apply_dirichlet_with_trafo(const CORE::LINALG::SparseMatrix& trafo,
    const Epetra_Map& dbctoggle, bool diagonalblock, bool complete)
{
  if (not Filled()) FOUR_C_THROW("expect filled matrix to apply dirichlet conditions");

  if (dbcmaps_ != Teuchos::null)
  {
    FOUR_C_THROW("Dirichlet map and transformations cannot be combined");
  }

  if (explicitdirichlet_)
  {
    // Save graph of original matrix if not done already.
    // This will never happen as the matrix is guaranteed to be filled. But to
    // make the code more explicit...
    if (savegraph_ and graph_ == Teuchos::null)
    {
      graph_ = Teuchos::rcp(new Epetra_CrsGraph(sysmat_->Graph()));
      if (not graph_->Filled()) FOUR_C_THROW("got unfilled graph from filled matrix");
    }

    // allocate a new matrix and copy all rows that are not dirichlet
    const Epetra_Map& rowmap = sysmat_->RowMap();
    const Epetra_Map& colmap = sysmat_->ColMap();
    const int nummyrows = sysmat_->NumMyRows();
    const int maxnumentries = sysmat_->MaxNumEntries();

    // prepare working arrays for extracting rows in trafo matrix
    const int trafomaxnumentries = trafo.MaxNumEntries();
    int trafonumentries = 0;
    std::vector<int> trafoindices(trafomaxnumentries, 0);
    std::vector<double> trafovalues(trafomaxnumentries, 0.0);

    // initialise matrix Anew with general size (rowmap x colmap)
    // in case of coupled problem (e.g. TSI) transform the rectangular off-diagonal block k_Td
    Teuchos::RCP<Epetra_CrsMatrix> Anew =
        Teuchos::rcp(new Epetra_CrsMatrix(::Copy, rowmap, colmap, maxnumentries, false));
    std::vector<int> indices(maxnumentries, 0);
    std::vector<double> values(maxnumentries, 0.0);
    for (int i = 0; i < nummyrows; ++i)
    {
      int row = sysmat_->GRID(i);
      if (not dbctoggle.MyGID(row))  // dof is not a Dirichlet dof
      {
        int numentries;
#ifdef FOUR_C_ENABLE_ASSERTIONS
        int err = sysmat_->ExtractGlobalRowCopy(
            row, maxnumentries, numentries, values.data(), indices.data());
        if (err < 0) FOUR_C_THROW("Epetra_CrsMatrix::ExtractGlobalRowCopy returned err=%d", err);
#else
        sysmat_->ExtractGlobalRowCopy(
            row, maxnumentries, numentries, values.data(), indices.data());
#endif

#ifdef FOUR_C_ENABLE_ASSERTIONS
        err = Anew->InsertGlobalValues(row, numentries, values.data(), indices.data());
        if (err < 0) FOUR_C_THROW("Epetra_CrsMatrix::InsertGlobalValues returned err=%d", err);
#else
        Anew->InsertGlobalValues(row, numentries, values.data(), indices.data());
#endif
      }
      else  // dof is an inclined Dirichlet dof
      {
        // diagonal block of dof with INCLINED Dirichlet boundary condition
        if (diagonalblock)
        {
          // extract values of trafo at the inclined dbc dof
#ifdef FOUR_C_ENABLE_ASSERTIONS
          int err = trafo.EpetraMatrix()->ExtractGlobalRowCopy(
              row, trafomaxnumentries, trafonumentries, trafovalues.data(), trafoindices.data());
          if (err < 0) FOUR_C_THROW("Epetra_CrsMatrix::ExtractGlobalRowCopy returned err=%d", err);
#else
          trafo.EpetraMatrix()->ExtractGlobalRowCopy(
              row, trafomaxnumentries, trafonumentries, trafovalues.data(), trafoindices.data());
#endif
        }
        // if entry of dof with inclined dbc is not a diagonal block, set zero
        // at this position
        else
        {
          trafonumentries = 1;
          trafovalues[0] = 0.0;
          trafoindices[0] = row;
        }
        // insert all these entries in transformed sysmat, i.e. in Anew
#ifdef FOUR_C_ENABLE_ASSERTIONS
        {
          int err = Anew->InsertGlobalValues(
              row, trafonumentries, trafovalues.data(), trafoindices.data());
          if (err < 0) FOUR_C_THROW("Epetra_CrsMatrix::InsertGlobalValues returned err=%d", err);
        }
#else
        Anew->InsertGlobalValues(row, trafonumentries, trafovalues.data(), trafoindices.data());
#endif
      }
    }
    // Updated sysmat_
    // normal DBC dof: '1.0' at diagonal, rest of row is blanked --> row remains
    //                 the same
    // inclined DBC: (in) rotated matrix k^{~}, i.e. '1.0' at diagonal, rest of
    //                    row is blanked for n/t/b-direction
    //               (out) matrix in global system, i.e. k: for a node with 3
    //                     dofs in x/y/z-direction, trafo block is put at the
    //                     position of the dofs of this node, rest of row is blanked
    sysmat_ = Anew;
    if (complete) Complete();
  }
  else
  {
    const int nummyrows = sysmat_->NumMyRows();

    // prepare working arrays for extracting rows in trafo matrix
    const int trafomaxnumentries = trafo.MaxNumEntries();
    int trafonumentries = 0;
    std::vector<int> trafoindices(trafomaxnumentries, 0);
    std::vector<double> trafovalues(trafomaxnumentries, 0.0);

    for (int i = 0; i < nummyrows; ++i)
    {
      int row = sysmat_->GRID(i);
      if (dbctoggle.MyGID(row))
      {
        int* indexOffset;
        int* indices;
        double* values;
#ifdef FOUR_C_ENABLE_ASSERTIONS
        int err = sysmat_->ExtractCrsDataPointers(indexOffset, indices, values);
        if (err) FOUR_C_THROW("Epetra_CrsMatrix::ExtractCrsDataPointers returned err=%d", err);
#else
        sysmat_->ExtractCrsDataPointers(indexOffset, indices, values);
#endif
        // zero row
        memset(&values[indexOffset[i]], 0, (indexOffset[i + 1] - indexOffset[i]) * sizeof(double));

        if (diagonalblock)
        {
#ifdef FOUR_C_ENABLE_ASSERTIONS
          err = trafo.EpetraMatrix()->ExtractMyRowCopy(
              i, trafomaxnumentries, trafonumentries, trafovalues.data(), trafoindices.data());
          if (err < 0) FOUR_C_THROW("Epetra_CrsMatrix::ExtractGlobalRowCopy returned err=%d", err);
#else
          trafo.EpetraMatrix()->ExtractMyRowCopy(
              i, trafomaxnumentries, trafonumentries, trafovalues.data(), trafoindices.data());
#endif

#ifdef FOUR_C_ENABLE_ASSERTIONS
          err =
              sysmat_->SumIntoMyValues(i, trafonumentries, trafovalues.data(), trafoindices.data());
          if (err < 0) FOUR_C_THROW("Epetra_CrsMatrix::SumIntoMyValues returned err=%d", err);
#else
          sysmat_->SumIntoMyValues(i, trafonumentries, trafovalues.data(), trafoindices.data());
#endif
        }
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> CORE::LINALG::SparseMatrix::extract_dirichlet_rows(
    const Teuchos::RCP<Epetra_Vector> dbctoggle)
{
  if (not Filled()) FOUR_C_THROW("expect filled matrix to extract dirichlet lines");

  Teuchos::RCP<SparseMatrix> dl =
      Teuchos::rcp(new SparseMatrix(RowMap(), MaxNumEntries(), ExplicitDirichlet(), SaveGraph()));

  const Epetra_Map& rowmap = sysmat_->RowMap();
  const Epetra_Map& colmap = sysmat_->ColMap();
  const int nummyrows = sysmat_->NumMyRows();

  const Epetra_Vector& dbct = *dbctoggle;

  std::vector<int> idx(MaxNumEntries());

  for (int i = 0; i < nummyrows; ++i)
  {
    if (dbct[i] == 1.0)
    {
      int NumEntries;
      double* Values;
      int* Indices;

      int err = sysmat_->ExtractMyRowView(i, NumEntries, Values, Indices);
      if (err) FOUR_C_THROW("ExtractMyRowView: err=%d", err);
      for (int j = 0; j < NumEntries; ++j) idx[j] = colmap.GID(Indices[j]);

      err = dl->sysmat_->InsertGlobalValues(rowmap.GID(i), NumEntries, Values, idx.data());
      if (err) FOUR_C_THROW("InsertGlobalValues: err=%d", err);
    }
  }

  dl->Complete(sysmat_->DomainMap(), RangeMap());
  return dl;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> CORE::LINALG::SparseMatrix::extract_dirichlet_rows(
    const Epetra_Map& dbctoggle)
{
  if (not Filled()) FOUR_C_THROW("expect filled matrix to extract dirichlet lines");
  if (not dbctoggle.UniqueGIDs()) FOUR_C_THROW("unique map required");

  Teuchos::RCP<SparseMatrix> dl =
      Teuchos::rcp(new SparseMatrix(RowMap(), MaxNumEntries(), ExplicitDirichlet(), SaveGraph()));

  const Epetra_Map& rowmap = sysmat_->RowMap();
  const Epetra_Map& colmap = sysmat_->ColMap();
  // const int nummyrows      = sysmat_->NumMyRows();

  std::vector<int> idx(MaxNumEntries());

  const int mylength = dbctoggle.NumMyElements();
  const int* mygids = dbctoggle.MyGlobalElements();
  for (int i = 0; i < mylength; ++i)
  {
    int gid = mygids[i];
    int lid = rowmap.LID(gid);

    if (lid < 0) FOUR_C_THROW("illegal Dirichlet map");

    int NumEntries;
    double* Values;
    int* Indices;

    int err = sysmat_->ExtractMyRowView(lid, NumEntries, Values, Indices);
    if (err) FOUR_C_THROW("ExtractMyRowView: err=%d", err);
    for (int j = 0; j < NumEntries; ++j) idx[j] = colmap.GID(Indices[j]);

    err = dl->sysmat_->InsertGlobalValues(gid, NumEntries, Values, idx.data());
    if (err) FOUR_C_THROW("InsertGlobalValues: err=%d", err);
  }

  dl->Complete(sysmat_->DomainMap(), RangeMap());
  return dl;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* CORE::LINALG::SparseMatrix::Label() const { return "CORE::LINALG::SparseMatrix"; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> CORE::LINALG::SparseMatrix::Transpose()
{
  if (not Filled()) FOUR_C_THROW("FillComplete was not called on matrix");

  EpetraExt::RowMatrix_Transpose trans;
  Teuchos::RCP<CORE::LINALG::SparseMatrix> matrix = Teuchos::null;

  if (matrixtype_ == CRS_MATRIX)
  {
    Epetra_CrsMatrix* Aprime = &(dynamic_cast<Epetra_CrsMatrix&>(trans(*sysmat_)));
    matrix = Teuchos::rcp(new SparseMatrix(
        Teuchos::rcp(Aprime, false), CORE::LINALG::Copy, explicitdirichlet_, savegraph_));
  }
  else if (matrixtype_ == FE_MATRIX)
  {
    Epetra_CrsMatrix* Aprime = &(dynamic_cast<Epetra_CrsMatrix&>(trans(*sysmat_)));
    matrix = Teuchos::rcp(new SparseMatrix(Teuchos::rcp(Aprime, false), CORE::LINALG::Copy,
        explicitdirichlet_, savegraph_, FE_MATRIX));
  }
  else
    FOUR_C_THROW("matrix type is not correct");

  return matrix;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int CORE::LINALG::SparseMatrix::ReplaceRowMap(const Epetra_BlockMap& newmap)
{
  const int err = CORE::LINALG::SparseMatrixBase::ReplaceRowMap(newmap);

  // change mask as well
  if ((not err) and savegraph_ and graph_ != Teuchos::null)
  {
    if (graph_.strong_count() > 1)
      FOUR_C_THROW("The graph_ can not be changed! (strong_count = %d)", graph_.strong_count());

    graph_ = Teuchos::rcp(new Epetra_CrsGraph(sysmat_->Graph()));
  }

  return err;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::Add(const CORE::LINALG::SparseMatrixBase& A, const bool transposeA,
    const double scalarA, const double scalarB)
{
  CORE::LINALG::Add(*A.EpetraMatrix(), transposeA, scalarA, *this, scalarB);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::Add(
    const Epetra_CrsMatrix& A, const bool transposeA, const double scalarA, const double scalarB)
{
  CORE::LINALG::Add(A, transposeA, scalarA, *this, scalarB);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::Put(const CORE::LINALG::SparseMatrix& A, const double scalarA,
    Teuchos::RCP<const Epetra_Map> rowmap)
{
  // put values onto sysmat
  if (A.GetMatrixtype() != CORE::LINALG::SparseMatrix::CRS_MATRIX)
    FOUR_C_THROW("Please check code and see wether it is save to apply it to matrix type %d",
        A.GetMatrixtype());
  Epetra_CrsMatrix* Aprime = const_cast<Epetra_CrsMatrix*>(&(*(A.EpetraMatrix())));
  if (Aprime == nullptr) FOUR_C_THROW("Cast failed");

  // Loop over Aprime's rows, extract row content and replace respective row in sysmat
  const int MaxNumEntries = EPETRA_MAX(Aprime->MaxNumEntries(), sysmat_->MaxNumEntries());

  // define row map to tackle
  // if #rowmap (a subset of #RowMap()) is provided, a selective replacing is perfomed
  const Epetra_Map* tomap = nullptr;
  if (rowmap != Teuchos::null)
    tomap = &(*rowmap);
  else
    tomap = &(RowMap());

  // working variables
  int NumEntries;
  std::vector<int> Indices(MaxNumEntries);
  std::vector<double> Values(MaxNumEntries);
  int err;

  // loop rows in #tomap and replace the rows of #this->sysmat_ with provided input matrix #A
  for (int lid = 0; lid < tomap->NumMyElements(); ++lid)
  {
    const int Row = tomap->GID(lid);
    if (Row < 0) FOUR_C_THROW("DOF not found on processor.");
    err =
        Aprime->ExtractGlobalRowCopy(Row, MaxNumEntries, NumEntries, Values.data(), Indices.data());
    if (err) FOUR_C_THROW("Epetra_CrsMatrix::ExtractGlobalRowCopy returned err=%d", err);
    if (scalarA != 1.0)
      for (int j = 0; j < NumEntries; ++j) Values[j] *= scalarA;
    err = sysmat_->ReplaceGlobalValues(Row, NumEntries, Values.data(), Indices.data());
    if (err) FOUR_C_THROW("Epetra_CrsMatrix::ReplaceGlobalValues returned err=%d", err);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::Dump(const std::string& filename)
{
  int MyRow;
  int NumEntries;
  std::stringstream rowsetname;
  rowsetname << filename << "." << Comm().MyPID() << ".row";
  std::stringstream offsetname;
  offsetname << filename << "." << Comm().MyPID() << ".off";
  std::stringstream indicesname;
  indicesname << filename << "." << Comm().MyPID() << ".idx";
  std::stringstream valuesname;
  valuesname << filename << "." << Comm().MyPID() << ".val";

  std::ofstream row(rowsetname.str().c_str());
  std::ofstream off(offsetname.str().c_str());
  std::ofstream idx(indicesname.str().c_str());
  std::ofstream val(valuesname.str().c_str());

  const Epetra_Map& rowmap = RowMap();
  const Epetra_Map& colmap = ColMap();

  if (sysmat_->Filled())
  {
    for (MyRow = 0; MyRow < sysmat_->NumMyRows(); ++MyRow)
    {
      double* Values;
      int* Indices;

      int err = sysmat_->ExtractMyRowView(MyRow, NumEntries, Values, Indices);
      if (err) FOUR_C_THROW("ExtractMyRowView failed: err=%d", err);
      row << rowmap.GID(MyRow) << "\n";
      off << NumEntries << "\n";
      // std::copy(Indices,Indices+NumEntries, std::ostream_iterator<int>(idx," "));
      for (int i = 0; i < NumEntries; ++i)
      {
        idx << colmap.GID(Indices[i]) << " ";
      }
      idx << "\n";
      val.write(reinterpret_cast<char*>(Values), NumEntries * sizeof(double));
    }
  }
  else
  {
    // Warning: does not write nonlocal values for Epetra_FECrsMatrices
    for (MyRow = 0; MyRow < sysmat_->NumMyRows(); ++MyRow)
    {
      std::vector<double> Values(sysmat_->MaxNumEntries());
      std::vector<int> Indices(sysmat_->MaxNumEntries());

      int err = sysmat_->ExtractGlobalRowCopy(
          rowmap.GID(MyRow), sysmat_->MaxNumEntries(), NumEntries, Values.data(), Indices.data());
      if (err) FOUR_C_THROW("ExtractGlobalRowCopy failed: err=%d", err);
      row << rowmap.GID(MyRow) << "\n";
      off << NumEntries << "\n";
      std::copy(Indices.data(), Indices.data() + NumEntries, std::ostream_iterator<int>(idx, " "));
      idx << "\n";
      val.write(reinterpret_cast<char*>(Values.data()), NumEntries * sizeof(double));
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::Load(const Epetra_Comm& comm, std::string& filename)
{
  std::stringstream rowsetname;
  rowsetname << filename << "." << comm.MyPID() << ".row";
  std::stringstream offsetname;
  offsetname << filename << "." << comm.MyPID() << ".off";
  std::stringstream indicesname;
  indicesname << filename << "." << comm.MyPID() << ".idx";
  std::stringstream valuesname;
  valuesname << filename << "." << comm.MyPID() << ".val";

  std::ifstream row(rowsetname.str().c_str());
  std::ifstream off(offsetname.str().c_str());
  std::ifstream idx(indicesname.str().c_str());
  std::ifstream val(valuesname.str().c_str());

  std::vector<int> rowids;
  for (;;)
  {
    int r;
    row >> r;
    if (not row) break;
    rowids.push_back(r);
  }

  Epetra_Map rowmap(-1, rowids.size(), rowids.data(), 0, comm);

  if (!rowmap.UniqueGIDs()) FOUR_C_THROW("Row map is not unique");

  std::vector<int> offnum;
  for (;;)
  {
    int o;
    off >> o;
    if (not off) break;
    offnum.push_back(o);
  }

  if (rowids.size() != offnum.size()) FOUR_C_THROW("failed to read rows and column counts");

  if (matrixtype_ == CRS_MATRIX)
    sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(::Copy, rowmap, offnum.data(), false));
  else if (matrixtype_ == FE_MATRIX)
    sysmat_ = Teuchos::rcp(new Epetra_FECrsMatrix(::Copy, rowmap, offnum.data(), false));
  else
    FOUR_C_THROW("matrix type is not correct");

  for (unsigned i = 0; i < rowids.size(); ++i)
  {
    std::vector<int> indices(offnum[i]);
    std::vector<double> values(offnum[i]);

    for (int j = 0; j < offnum[i]; ++j)
    {
      idx >> indices[j];
      if (not idx) FOUR_C_THROW("failed to read indices");
    }
    val.read(reinterpret_cast<char*>(values.data()), offnum[i] * sizeof(double));
    int err = sysmat_->InsertGlobalValues(rowids[i], offnum[i], values.data(), indices.data());
    if (err) FOUR_C_THROW("InsertGlobalValues failed: err=%d", err);
  }

  graph_ = Teuchos::null;
}


/*----------------------------------------------------------------------*
  (private)
 *----------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::Split2x2(BlockSparseMatrixBase& Abase) const
{
  // for timing of this method
  TEUCHOS_FUNC_TIME_MONITOR("CORE::LINALG::SparseMatrix::Split2x2");

  if (Abase.Rows() != 2 || Abase.Cols() != 2) FOUR_C_THROW("Can only split in 2x2 system");
  if (!Filled()) FOUR_C_THROW("SparseMatrix must be filled");
  Teuchos::RCP<Epetra_CrsMatrix> A = EpetraMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> A11 = Abase(0, 0).EpetraMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> A12 = Abase(0, 1).EpetraMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> A21 = Abase(1, 0).EpetraMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> A22 = Abase(1, 1).EpetraMatrix();
  if (A11->Filled() || A12->Filled() || A21->Filled() || A22->Filled())
    FOUR_C_THROW("Block matrix may not be filled on input");
  const Epetra_Map& A11rmap = Abase.RangeMap(0);
  const Epetra_Map& A11dmap = Abase.DomainMap(0);
  const Epetra_Map& A22rmap = Abase.RangeMap(1);
  const Epetra_Map& A22dmap = Abase.DomainMap(1);

  // find out about how the column map is linked to the individual processors.
  // this is done by filling the information about the rowmap into a vector that
  // is then exported to the column map
  Epetra_Vector dselector(A->DomainMap());
  for (int i = 0; i < dselector.MyLength(); ++i)
  {
    const int gid = A->DomainMap().GID(i);
    if (A11dmap.MyGID(gid))
      dselector[i] = 0.;
    else if (A22dmap.MyGID(gid))
      dselector[i] = 1.;
    else
      dselector[i] = -1.;
  }
  Epetra_Vector selector(A->ColMap());
  CORE::LINALG::Export(dselector, selector);

  std::vector<int> gcindices1(A->MaxNumEntries());
  std::vector<double> gvalues1(A->MaxNumEntries());
  std::vector<int> gcindices2(A->MaxNumEntries());
  std::vector<double> gvalues2(A->MaxNumEntries());
  //-------------------------------------------------- create block matrices
  const int length = A->NumMyRows();
  for (int i = 0; i < length; ++i)
  {
    int err1 = 0;
    int err2 = 0;
    int count1 = 0;
    int count2 = 0;
    const int grid = A->GRID(i);
    if (!A11rmap.MyGID(grid) && !A22rmap.MyGID(grid)) continue;
    int numentries;
    double* values;
    int* cindices;
    int err = A->ExtractMyRowView(i, numentries, values, cindices);
    if (err) FOUR_C_THROW("SparseMatrix::Split2x2: A->ExtractMyRowView returned %d", err);
    for (int j = 0; j < numentries; ++j)
    {
      const int gcid = A->ColMap().GID(cindices[j]);
      FOUR_C_ASSERT(cindices[j] < selector.MyLength(), "Internal error");
      // column is in A*1
      if (selector[cindices[j]] == 0.)
      {
        gcindices1[count1] = gcid;
        gvalues1[count1++] = values[j];
      }
      // column is in A*2
      else if (selector[cindices[j]] == 1.)
      {
        gcindices2[count2] = gcid;
        gvalues2[count2++] = values[j];
      }
      else
        FOUR_C_THROW("Could not identify column index with block, internal error.");
    }
    //======================== row belongs to A11 and A12
    if (A11rmap.MyGID(grid))
    {
      if (count1) err1 = A11->InsertGlobalValues(grid, count1, gvalues1.data(), gcindices1.data());
      if (count2) err2 = A12->InsertGlobalValues(grid, count2, gvalues2.data(), gcindices2.data());
    }
    //======================= row belongs to A21 and A22
    else
    {
      if (count1) err1 = A21->InsertGlobalValues(grid, count1, gvalues1.data(), gcindices1.data());
      if (count2) err2 = A22->InsertGlobalValues(grid, count2, gvalues2.data(), gcindices2.data());
    }
    // #ifdef FOUR_C_ENABLE_ASSERTIONS
    if (err1 < 0 || err2 < 0)
    {
      FOUR_C_THROW(
          "SparseMatrix::Split2x2: Epetra_CrsMatrix::InsertGlobalValues returned err1=%d / err2=%d",
          err1, err2);
    }
    // #endif
  }  // for (int i=0; i<A->NumMyRows(); ++i)
  // Do not complete BlockMatrix
}


/*--------------------------------------------------------------------------*
 | split SparseMatrix into an MxN BlockSparseMatrixBase          fang 02/17 |
 *--------------------------------------------------------------------------*/
void CORE::LINALG::SparseMatrix::SplitMxN(BlockSparseMatrixBase& ABlock) const
{
  // timing
  TEUCHOS_FUNC_TIME_MONITOR("CORE::LINALG::SparseMatrix::SplitMxN");

  // extract number of row/column blocks
  const int M = ABlock.Rows();
  const int N = ABlock.Cols();

  // safety checks
  if (!Filled()) FOUR_C_THROW("SparseMatrix must be filled before splitting!");
  for (int m = 0; m < M; ++m)
  {
    for (int n = 0; n < N; ++n)
      if (ABlock(m, n).EpetraMatrix()->Filled())
        FOUR_C_THROW("BlockSparseMatrixBase must not be filled before splitting!");
  }

  // extract EpetraMatrix
  const Epetra_CrsMatrix& A = *EpetraMatrix();

  // associate each global column ID of SparseMatrix with corresponding block ID of
  // BlockSparseMatrixBase this is done via an Epetra_Vector which is filled using domain map
  // information and then exported to column map
  Epetra_Vector dselector(DomainMap());
  for (int collid = 0; collid < dselector.MyLength(); ++collid)
  {
    // extract global ID of current column
    const int colgid = DomainMap().GID(collid);
    if (colgid < 0) FOUR_C_THROW("Couldn't find local column ID %d in domain map!", collid);

    // determine block ID of BlockSparseMatrixBase associated with current column
    int n(0);
    for (n = 0; n < N; ++n)
    {
      if (ABlock.DomainMap(n).MyGID(colgid))
      {
        dselector[collid] = n;
        break;
      }
    }

    // safety check
    if (n == N) FOUR_C_THROW("Matrix column was not found in BlockSparseMatrixBase!");
  }
  Epetra_Vector selector(A.ColMap());
  CORE::LINALG::Export(dselector, selector);

  // allocate vectors storing global column indexes and values of matrix entries in a given row,
  // separated by blocks allocation is done outside loop over all rows for efficiency to be on the
  // safe side, we allocate more memory than we need for most rows
  std::vector<std::vector<int>> colgids(N, std::vector<int>(A.MaxNumEntries(), -1));
  std::vector<std::vector<double>> rowvalues(N, std::vector<double>(A.MaxNumEntries(), 0.));

  // fill blocks of BlockSparseMatrixBase
  for (int rowlid = 0; rowlid < A.NumMyRows(); ++rowlid)
  {
    // extract current row of SparseMatrix
    int numentries(0);
    double* values(nullptr);
    int* indices(nullptr);
    if (A.ExtractMyRowView(rowlid, numentries, values, indices))
      FOUR_C_THROW("Row of SparseMatrix couldn't be extracted during splitting!");

    // initialize counters for number of matrix entries in current row, separated by blocks
    std::vector<unsigned> counters(N, 0);

    // assign matrix entries in current row to associated blocks of BlockSparseMatrixBase
    for (int j = 0; j < numentries; ++j)
    {
      // extract local column ID of current matrix entry
      const int collid = indices[j];
      if (collid >= selector.MyLength()) FOUR_C_THROW("Invalid local column ID %d!", collid);

      // assign current matrix entry to associated block of BlockSparseMatrixBase
      const int blockid = static_cast<int>(selector[collid]);
      colgids[blockid][counters[blockid]] = A.ColMap().GID(collid);
      rowvalues[blockid][counters[blockid]++] = values[j];
    }

    // extract global index of current matrix row
    const int rowgid = A.GRID(rowlid);

    // fill current row of BlockSparseMatrixBase by copying row entries of SparseMatrix
    int m(0);
    for (m = 0; m < M; ++m)
    {
      if (ABlock.RangeMap(m).MyGID(rowgid))
      {
        for (int n = 0; n < N; ++n)
        {
          if (counters[n])
          {
            if (ABlock(m, n).EpetraMatrix()->InsertGlobalValues(
                    rowgid, counters[n], rowvalues[n].data(), colgids[n].data()))
              FOUR_C_THROW("Couldn't insert matrix entries into BlockSparseMatrixBase!");
          }
        }
        break;
      }
    }

    // safety check
    if (m == M) FOUR_C_THROW("Matrix row was not found in BlockSparseMatrixBase!");
  }  // loop over matrix rows
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& CORE::LINALG::operator<<(std::ostream& os, const CORE::LINALG::SparseMatrix& mat)
{
  if (mat.GetMatrixtype() == SparseMatrix::CRS_MATRIX)
    os << *(const_cast<CORE::LINALG::SparseMatrix&>(mat).EpetraMatrix());
  else if (mat.GetMatrixtype() == SparseMatrix::FE_MATRIX)
  {
    os << *(Teuchos::rcp_dynamic_cast<Epetra_FECrsMatrix>(
        const_cast<CORE::LINALG::SparseMatrix&>(mat).EpetraMatrix()));
  }
  else
    FOUR_C_THROW("matrixtype does not exist");
  return os;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> CORE::LINALG::Multiply(const CORE::LINALG::SparseMatrix& A,
    bool transA, const CORE::LINALG::SparseMatrix& B, bool transB, bool completeoutput)
{
  // make sure FillComplete was called on the matrices
  if (!A.Filled()) FOUR_C_THROW("A has to be FillComplete");
  if (!B.Filled()) FOUR_C_THROW("B has to be FillComplete");

  // const int npr = A.EpetraMatrix()->MaxNumEntries()*B.EpetraMatrix()->MaxNumEntries();
  // a first guess for the bandwidth of C leading to much less memory consumption:
  const int npr = std::max(A.MaxNumEntries(), B.MaxNumEntries());

  // now create resultmatrix with correct rowmap
  Teuchos::RCP<CORE::LINALG::SparseMatrix> C;
  if (!transA)
    C = Teuchos::rcp(new SparseMatrix(A.RangeMap(), npr, A.explicitdirichlet_, A.savegraph_));
  else
    C = Teuchos::rcp(new SparseMatrix(A.DomainMap(), npr, A.explicitdirichlet_, A.savegraph_));

  int err = EpetraExt::MatrixMatrix::Multiply(
      *A.sysmat_, transA, *B.sysmat_, transB, *C->sysmat_, completeoutput);
  if (err) FOUR_C_THROW("EpetraExt::MatrixMatrix::Multiply returned err = %d", err);

  return C;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> CORE::LINALG::Multiply(const CORE::LINALG::SparseMatrix& A,
    bool transA, const CORE::LINALG::SparseMatrix& B, bool transB, bool explicitdirichlet,
    bool savegraph, bool completeoutput)
{
  // make sure FillComplete was called on the matrices
  if (!A.Filled()) FOUR_C_THROW("A has to be FillComplete");
  if (!B.Filled()) FOUR_C_THROW("B has to be FillComplete");

  // const int npr = A.EpetraMatrix()->MaxNumEntries()*B.EpetraMatrix()->MaxNumEntries();
  // a first guess for the bandwidth of C leading to much less memory consumption:
  const int npr = std::max(A.MaxNumEntries(), B.MaxNumEntries());

  // now create resultmatrix C with correct rowmap
  Teuchos::RCP<CORE::LINALG::SparseMatrix> C;
  if (!transA)
    C = Teuchos::rcp(new SparseMatrix(A.RangeMap(), npr, explicitdirichlet, savegraph));
  else
    C = Teuchos::rcp(new SparseMatrix(A.DomainMap(), npr, explicitdirichlet, savegraph));

  int err = EpetraExt::MatrixMatrix::Multiply(
      *A.sysmat_, transA, *B.sysmat_, transB, *C->sysmat_, completeoutput);
  if (err) FOUR_C_THROW("EpetraExt::MatrixMatrix::Multiply returned err = %d", err);

  return C;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> CORE::LINALG::Merge(const CORE::LINALG::SparseMatrix& Aii,
    const CORE::LINALG::SparseMatrix& Aig, const CORE::LINALG::SparseMatrix& Agi,
    const CORE::LINALG::SparseMatrix& Agg)
{
  if (not Aii.RowMap().SameAs(Aig.RowMap()) or not Agi.RowMap().SameAs(Agg.RowMap()))
    FOUR_C_THROW("row maps mismatch");

  Teuchos::RCP<Epetra_Map> rowmap = MergeMap(Aii.RowMap(), Agi.RowMap(), false);
  Teuchos::RCP<CORE::LINALG::SparseMatrix> mat =
      Teuchos::rcp(new SparseMatrix(*rowmap, std::max(Aii.MaxNumEntries() + Aig.MaxNumEntries(),
                                                 Agi.MaxNumEntries() + Agg.MaxNumEntries())));


  mat->Add(Aii, false, 1.0, 1.0);
  mat->Add(Aig, false, 1.0, 1.0);
  mat->Add(Agi, false, 1.0, 1.0);
  mat->Add(Agg, false, 1.0, 1.0);

  return mat;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> CORE::LINALG::Eye(const Epetra_Map& map)
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> eye = Teuchos::rcp(new SparseMatrix(map, 1));
  int numelements = map.NumMyElements();
  int* gids = map.MyGlobalElements();
  for (int i = 0; i < numelements; ++i)
  {
    int gid = gids[i];
    eye->Assemble(1., gid, gid);
  }
  eye->Complete();
  return eye;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> CORE::LINALG::CastToSparseMatrixAndCheckSuccess(
    Teuchos::RCP<CORE::LINALG::SparseOperator> input_matrix)
{
  auto sparse_matrix = Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(input_matrix);
  FOUR_C_ASSERT(sparse_matrix != Teuchos::null, "Matrix is not a sparse matrix!");

  return sparse_matrix;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::SparseMatrix> CORE::LINALG::CastToConstSparseMatrixAndCheckSuccess(
    Teuchos::RCP<const CORE::LINALG::SparseOperator> input_matrix)
{
  auto sparse_matrix = Teuchos::rcp_dynamic_cast<const CORE::LINALG::SparseMatrix>(input_matrix);
  FOUR_C_ASSERT(sparse_matrix != Teuchos::null, "Matrix is not a sparse matrix!");

  return sparse_matrix;
}

FOUR_C_NAMESPACE_CLOSE
