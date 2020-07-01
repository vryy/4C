/*----------------------------------------------------------------------*/
/*! \file

\brief Utilities for matrix transformations

\level 1

*/
/*---------------------------------------------------------------------*/

#include <vector>
#include <iterator>
#include "linalg_matrixtransform.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_exporter.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool LINALG::MatrixLogicalSplitAndTransform::operator()(const LINALG::SparseMatrix& src,
    const Epetra_Map& logical_range_map, const Epetra_Map& logical_domain_map, double scale,
    const ADAPTER::CouplingConverter* row_converter,
    const ADAPTER::CouplingConverter* col_converter, LINALG::SparseMatrix& dst, bool exactmatch,
    bool addmatrix)
{
  Teuchos::RCP<Epetra_CrsMatrix> esrc = src.EpetraMatrix();
  const Epetra_Map* final_range_map = &logical_range_map;
  const Epetra_Map* matching_dst_rows = &logical_range_map;

  if (row_converter)
  {
    const Epetra_Map& permsrcmap = *row_converter->PermSrcMap();

    // check if the permuted map is simply a subset of the current rowmap (no communication)
    int subset = 1;
    for (int i = 0; i < permsrcmap.NumMyElements(); ++i)
      if (!src.RowMap().MyGID(permsrcmap.GID(i)))
      {
        subset = 0;
        break;
      }

    int gsubset = 0;
    logical_range_map.Comm().MinAll(&subset, &gsubset, 1);

    // need communication -> call import on permuted map
    if (!gsubset)
    {
      if (exporter_ == Teuchos::null)
      {
        exporter_ = Teuchos::rcp(new Epetra_Export(permsrcmap, src.RowMap()));
      }

      Teuchos::RCP<Epetra_CrsMatrix> permsrc =
          Teuchos::rcp(new Epetra_CrsMatrix(::Copy, permsrcmap, 0));
      int err = permsrc->Import(*src.EpetraMatrix(), *exporter_, Insert);
      if (err) dserror("Import failed with err=%d", err);

      permsrc->FillComplete(src.DomainMap(), permsrcmap);
      esrc = permsrc;
    }

    final_range_map = &permsrcmap;
    matching_dst_rows = row_converter->DstMap().get();
  }

  SetupGidMap(col_converter ? *col_converter->SrcMap() : esrc->RowMap(), esrc->ColMap(),
      col_converter, src.Comm());

  if (!addmatrix) dst.Zero();

  InternalAdd(esrc, *final_range_map, col_converter ? *col_converter->SrcMap() : logical_domain_map,
      *matching_dst_rows, dst.EpetraMatrix(), exactmatch, scale);

  return true;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MatrixLogicalSplitAndTransform::SetupGidMap(const Epetra_Map& rowmap,
    const Epetra_Map& colmap, const ADAPTER::CouplingConverter* converter, const Epetra_Comm& comm)
{
  if (not havegidmap_)
  {
    if (converter != NULL)
    {
      DRT::Exporter ex(rowmap, colmap, comm);
      converter->FillSrcToDstMap(gidmap_);
      ex.Export(gidmap_);
    }
    else
      for (int i = 0; i < colmap.NumMyElements(); ++i) gidmap_[colmap.GID(i)] = colmap.GID(i);
    havegidmap_ = true;
  }
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MatrixLogicalSplitAndTransform::InternalAdd(Teuchos::RCP<Epetra_CrsMatrix> esrc,
    const Epetra_Map& logical_range_map, const Epetra_Map& logical_domain_map,
    const Epetra_Map& matching_dst_rows, Teuchos::RCP<Epetra_CrsMatrix> edst, bool exactmatch,
    double scale)
{
  if (not esrc->Filled()) dserror("filled source matrix expected");

  Epetra_Vector dselector(esrc->DomainMap());
  for (int i = 0; i < dselector.MyLength(); ++i)
  {
    const int gid = esrc->DomainMap().GID(i);
    if (logical_domain_map.MyGID(gid))
      dselector[i] = 1.;
    else
      dselector[i] = 0.;
  }
  Epetra_Vector selector(esrc->ColMap());
  LINALG::Export(dselector, selector);

  if (edst->Filled())
    AddIntoFilled(esrc, logical_range_map, logical_domain_map, selector, matching_dst_rows, edst,
        exactmatch, scale);
  else
    AddIntoUnfilled(esrc, logical_range_map, logical_domain_map, selector, matching_dst_rows, edst,
        exactmatch, scale);
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MatrixLogicalSplitAndTransform::AddIntoFilled(Teuchos::RCP<Epetra_CrsMatrix> esrc,
    const Epetra_Map& logical_range_map, const Epetra_Map& logical_domain_map,
    const Epetra_Vector& selector, const Epetra_Map& matching_dst_rows,
    Teuchos::RCP<Epetra_CrsMatrix> edst, bool exactmatch, double scale)
{
  const Epetra_Map& srccolmap = esrc->ColMap();
  const Epetra_Map& dstrowmap = edst->RowMap();

  // If the destination matrix is filled, we can add in local indices. This code is similar
  // to what is done in LINALG::Add(SparseMatrix, SparseMatrix) for the filled case.
  // We perform four steps:
  // 1. Identify the local column index mapping from the source to the destination matrix from
  //    on the global IDs
  // 2. Loop over the input matrix rows, extract row view in two matrices
  // 3. Match columns of row i in source matrix (called A) to the columns in the destination
  //    matrix (called B)
  // 4. Perform addition
  if (int(lidvector_.size()) != srccolmap.NumMyElements())
  {
    lidvector_.clear();
    lidvector_.resize(srccolmap.NumMyElements(), -1);
    for (std::map<int, int>::const_iterator iter = gidmap_.begin(); iter != gidmap_.end(); ++iter)
    {
      const int lid = srccolmap.LID(iter->first);
      if (lid != -1) lidvector_[lid] = edst->ColMap().LID(iter->second);
    }
  }

  int rows = logical_range_map.NumMyElements();
  for (int i = 0; i < rows; ++i)
  {
    int NumEntriesA, NumEntriesB;
    double *ValuesA, *ValuesB;
    int *IndicesA, *IndicesB;
    const int rowA = esrc->RowMap().LID(logical_range_map.GID(i));
    if (rowA == -1) dserror("Internal error");
    int err = esrc->ExtractMyRowView(rowA, NumEntriesA, ValuesA, IndicesA);
    if (err != 0) dserror("ExtractMyRowView error: %d", err);

    // identify the local row index in the destination matrix corresponding to i
    const int rowB = dstrowmap.LID(matching_dst_rows.GID(i));
    err = edst->ExtractMyRowView(rowB, NumEntriesB, ValuesB, IndicesB);
    if (err != 0) dserror("ExtractMyRowView error: %d", err);

    // loop through the columns in source matrix and find respective place in destination
    for (int jA = 0, jB = 0; jA < NumEntriesA; ++jA)
    {
      // skip entries belonging to a different block of the logical block matrix
      if (selector[IndicesA[jA]] == 0.) continue;

      const int col = lidvector_[IndicesA[jA]];
      if (col == -1)
      {
        if (exactmatch)
          dserror("gid %d not found in map for lid %d at %d", srccolmap.GID(IndicesA[jA]),
              IndicesA[jA], jA);
        else
          continue;
      }

      // try linear search in B
      while (jB < NumEntriesB && IndicesB[jB] < col) ++jB;

      // did not find index in linear search (re-indexing from A.ColMap() to B.ColMap()
      // might pass through the indices differently), try binary search
      if (jB == NumEntriesB || IndicesB[jB] != col)
        jB = std::lower_bound(&IndicesB[0], &IndicesB[0] + NumEntriesB, col) - &IndicesB[0];

      // not found, sparsity pattern of B does not contain the index from A -> terminate
      if (jB == NumEntriesB || IndicesB[jB] != col)
      {
        dserror(
            "Source matrix entry with global row ID %d and global column ID %d couldn't be added to"
            " destination matrix entry with global row ID %d and unknown global column ID %d!",
            esrc->RowMap().GID(i), srccolmap.GID(IndicesA[jA]), matching_dst_rows.GID(i),
            edst->ColMap().GID(col));
      }

      ValuesB[jB] += ValuesA[jA] * scale;
    }
  }
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void LINALG::MatrixLogicalSplitAndTransform::AddIntoUnfilled(Teuchos::RCP<Epetra_CrsMatrix> esrc,
    const Epetra_Map& logical_range_map, const Epetra_Map& logical_domain_map,
    const Epetra_Vector& selector, const Epetra_Map& matching_dst_rows,
    Teuchos::RCP<Epetra_CrsMatrix> edst, bool exactmatch, double scale)
{
  const Epetra_Map& srccolmap = esrc->ColMap();

  // standard code for the unfilled case
  std::vector<int> idx;
  std::vector<double> vals;
  int rows = logical_range_map.NumMyElements();
  for (int i = 0; i < rows; ++i)
  {
    int NumEntries;
    double* Values;
    int* Indices;
    int err = esrc->ExtractMyRowView(
        esrc->RowMap().LID(logical_range_map.GID(i)), NumEntries, Values, Indices);
    if (err != 0) dserror("ExtractMyRowView error: %d", err);

    idx.clear();
    vals.clear();

    for (int j = 0; j < NumEntries; ++j)
    {
      // skip entries belonging to a different block of the logical block matrix
      if (selector[Indices[j]] == 0.) continue;

      int gid = srccolmap.GID(Indices[j]);
      std::map<int, int>::const_iterator iter = gidmap_.find(gid);
      if (iter != gidmap_.end())
      {
        idx.push_back(iter->second);
        vals.push_back(Values[j] * scale);
      }
      else
      {
        // only complain if an exact match is demanded
        if (exactmatch) dserror("gid %d not found in map for lid %d at %d", gid, Indices[j], j);
      }
    }

    NumEntries = vals.size();
    const int globalRow = matching_dst_rows.GID(i);

    // put row into matrix
    //
    // We might want to preserve a Dirichlet row in our destination matrix
    // here as well. Skip for now.

    if (edst->NumAllocatedGlobalEntries(globalRow) == 0)
    {
      int err = edst->InsertGlobalValues(globalRow, NumEntries, &vals[0], &idx[0]);
      if (err < 0) dserror("InsertGlobalValues error: %d", err);
    }
    else
      for (int j = 0; j < NumEntries; ++j)
      {
        // add all values, including zeros, as we need a proper matrix graph
        int err = edst->SumIntoGlobalValues(globalRow, 1, &vals[j], &idx[j]);
        if (err > 0)
        {
          err = edst->InsertGlobalValues(globalRow, 1, &vals[j], &idx[j]);
          if (err < 0) dserror("InsertGlobalValues error: %d", err);
        }
        else if (err < 0)
          dserror("SumIntoGlobalValues error: %d", err);
      }
  }
}



bool LINALG::MatrixRowTransform::operator()(const LINALG::SparseMatrix& src, double scale,
    const ::ADAPTER::CouplingConverter& converter, LINALG::SparseMatrix& dst, bool addmatrix)
{
  return transformer(
      src, src.RangeMap(), src.DomainMap(), scale, &converter, NULL, dst, false, addmatrix);
}



bool LINALG::MatrixColTransform::operator()(const Epetra_Map&, const Epetra_Map&,
    const LINALG::SparseMatrix& src, double scale, const ::ADAPTER::CouplingConverter& converter,
    LINALG::SparseMatrix& dst, bool exactmatch, bool addmatrix)
{
  return transformer(
      src, src.RangeMap(), src.DomainMap(), scale, NULL, &converter, dst, exactmatch, addmatrix);
}



bool LINALG::MatrixRowColTransform::operator()(const LINALG::SparseMatrix& src, double scale,
    const ::ADAPTER::CouplingConverter& rowconverter,
    const ::ADAPTER::CouplingConverter& colconverter, LINALG::SparseMatrix& dst, bool exactmatch,
    bool addmatrix)
{
  return transformer(src, src.RangeMap(), src.DomainMap(), scale, &rowconverter, &colconverter, dst,
      exactmatch, addmatrix);
}
