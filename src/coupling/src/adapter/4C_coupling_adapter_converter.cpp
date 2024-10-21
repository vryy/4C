#include "4C_coupling_adapter_converter.hpp"

#include "4C_comm_exporter.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linalg_vector.hpp"

#include <Epetra_Export.h>
#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/

Teuchos::RCP<Core::LinAlg::Vector<double>> Coupling::Adapter::CouplingMasterConverter::src_to_dst(
    Teuchos::RCP<const Core::LinAlg::Vector<double>> source_vector) const
{
  return coup_.master_to_slave(*source_vector);
}

Teuchos::RCP<Core::LinAlg::Vector<double>> Coupling::Adapter::CouplingMasterConverter::dst_to_src(
    Teuchos::RCP<const Core::LinAlg::Vector<double>> destination_vector) const
{
  return coup_.slave_to_master(*destination_vector);
}

Teuchos::RCP<const Epetra_Map> Coupling::Adapter::CouplingMasterConverter::src_map() const
{
  return coup_.master_dof_map();
}

Teuchos::RCP<const Epetra_Map> Coupling::Adapter::CouplingMasterConverter::dst_map() const
{
  return coup_.slave_dof_map();
}

Teuchos::RCP<const Epetra_Map> Coupling::Adapter::CouplingMasterConverter::perm_src_map() const
{
  return coup_.perm_master_dof_map();
}

Teuchos::RCP<const Epetra_Map> Coupling::Adapter::CouplingMasterConverter::perm_dst_map() const
{
  return coup_.perm_slave_dof_map();
}

void Coupling::Adapter::CouplingMasterConverter::fill_src_to_dst_map(
    std::map<int, int>& rowmap) const
{
  coup_.fill_master_to_slave_map(rowmap);
}


Teuchos::RCP<Core::LinAlg::Vector<double>> Coupling::Adapter::CouplingSlaveConverter::src_to_dst(
    Teuchos::RCP<const Core::LinAlg::Vector<double>> source_vector) const
{
  return coup_.slave_to_master(*source_vector);
}

Teuchos::RCP<Core::LinAlg::Vector<double>> Coupling::Adapter::CouplingSlaveConverter::dst_to_src(
    Teuchos::RCP<const Core::LinAlg::Vector<double>> destination_vector) const
{
  return coup_.master_to_slave(*destination_vector);
}

Teuchos::RCP<const Epetra_Map> Coupling::Adapter::CouplingSlaveConverter::src_map() const
{
  return coup_.slave_dof_map();
}

Teuchos::RCP<const Epetra_Map> Coupling::Adapter::CouplingSlaveConverter::dst_map() const
{
  return coup_.master_dof_map();
}

Teuchos::RCP<const Epetra_Map> Coupling::Adapter::CouplingSlaveConverter::perm_src_map() const
{
  return coup_.perm_slave_dof_map();
}

Teuchos::RCP<const Epetra_Map> Coupling::Adapter::CouplingSlaveConverter::perm_dst_map() const
{
  return coup_.perm_master_dof_map();
}

void Coupling::Adapter::CouplingSlaveConverter::fill_src_to_dst_map(
    std::map<int, int>& rowmap) const
{
  coup_.fill_slave_to_master_map(rowmap);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Coupling::Adapter::MatrixLogicalSplitAndTransform::operator()(
    const Core::LinAlg::SparseMatrix& src, const Epetra_Map& logical_range_map,
    const Epetra_Map& logical_domain_map, double scale, const CouplingConverter* row_converter,
    const CouplingConverter* col_converter, Core::LinAlg::SparseMatrix& dst, bool exactmatch,
    bool addmatrix)
{
  Teuchos::RCP<Epetra_CrsMatrix> esrc = src.epetra_matrix();
  const Epetra_Map* final_range_map = &logical_range_map;
  const Epetra_Map* matching_dst_rows = &logical_range_map;

  if (row_converter)
  {
    const Epetra_Map& permsrcmap = *row_converter->perm_src_map();

    // check if the permuted map is simply a subset of the current rowmap (no communication)
    int subset = 1;
    for (int i = 0; i < permsrcmap.NumMyElements(); ++i)
      if (!src.row_map().MyGID(permsrcmap.GID(i)))
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
        exporter_ = Teuchos::make_rcp<Epetra_Export>(permsrcmap, src.row_map());
      }

      Teuchos::RCP<Epetra_CrsMatrix> permsrc =
          Teuchos::make_rcp<Epetra_CrsMatrix>(::Copy, permsrcmap, 0);
      int err = permsrc->Import(*src.epetra_matrix(), *exporter_, Insert);
      if (err) FOUR_C_THROW("Import failed with err=%d", err);

      permsrc->FillComplete(src.domain_map(), permsrcmap);
      esrc = permsrc;
    }

    final_range_map = &permsrcmap;
    matching_dst_rows = row_converter->dst_map().get();
  }

  setup_gid_map(col_converter ? *col_converter->src_map() : esrc->RowMap(), esrc->ColMap(),
      col_converter, src.Comm());

  if (!addmatrix) dst.zero();

  internal_add(*esrc, *final_range_map,
      col_converter ? *col_converter->src_map() : logical_domain_map, *matching_dst_rows,
      *dst.epetra_matrix(), exactmatch, scale);

  return true;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::MatrixLogicalSplitAndTransform::setup_gid_map(const Epetra_Map& rowmap,
    const Epetra_Map& colmap, const CouplingConverter* converter, const Epetra_Comm& comm)
{
  if (not havegidmap_)
  {
    if (converter != nullptr)
    {
      Core::Communication::Exporter ex(rowmap, colmap, comm);
      converter->fill_src_to_dst_map(gidmap_);
      ex.do_export(gidmap_);
    }
    else
      for (int i = 0; i < colmap.NumMyElements(); ++i) gidmap_[colmap.GID(i)] = colmap.GID(i);
    havegidmap_ = true;
  }
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::MatrixLogicalSplitAndTransform::internal_add(Epetra_CrsMatrix& esrc,
    const Epetra_Map& logical_range_map, const Epetra_Map& logical_domain_map,
    const Epetra_Map& matching_dst_rows, Epetra_CrsMatrix& edst, bool exactmatch, double scale)
{
  if (not esrc.Filled()) FOUR_C_THROW("filled source matrix expected");

  Core::LinAlg::Vector<double> dselector(esrc.DomainMap());
  for (int i = 0; i < dselector.MyLength(); ++i)
  {
    const int gid = esrc.DomainMap().GID(i);
    if (logical_domain_map.MyGID(gid))
      dselector[i] = 1.;
    else
      dselector[i] = 0.;
  }
  Core::LinAlg::Vector<double> selector(esrc.ColMap());
  Core::LinAlg::export_to(dselector, selector);

  if (edst.Filled())
    add_into_filled(esrc, logical_range_map, logical_domain_map, selector, matching_dst_rows, edst,
        exactmatch, scale);
  else
    add_into_unfilled(esrc, logical_range_map, logical_domain_map, selector, matching_dst_rows,
        edst, exactmatch, scale);
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::MatrixLogicalSplitAndTransform::add_into_filled(Epetra_CrsMatrix& esrc,
    const Epetra_Map& logical_range_map, const Epetra_Map& logical_domain_map,
    const Core::LinAlg::Vector<double>& selector, const Epetra_Map& matching_dst_rows,
    Epetra_CrsMatrix& edst, bool exactmatch, double scale)
{
  const Epetra_Map& srccolmap = esrc.ColMap();
  const Epetra_Map& dstrowmap = edst.RowMap();

  // If the destination matrix is filled, we can add in local indices. This code is similar
  // to what is done in Core::LinAlg::Add(SparseMatrix, SparseMatrix) for the filled case.
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
      if (lid != -1) lidvector_[lid] = edst.ColMap().LID(iter->second);
    }
  }

  int rows = logical_range_map.NumMyElements();
  for (int i = 0; i < rows; ++i)
  {
    int NumEntriesA, NumEntriesB;
    double *ValuesA, *ValuesB;
    int *IndicesA, *IndicesB;
    const int rowA = esrc.RowMap().LID(logical_range_map.GID(i));
    if (rowA == -1) FOUR_C_THROW("Internal error");
    int err = esrc.ExtractMyRowView(rowA, NumEntriesA, ValuesA, IndicesA);
    if (err != 0) FOUR_C_THROW("ExtractMyRowView error: %d", err);

    // identify the local row index in the destination matrix corresponding to i
    const int rowB = dstrowmap.LID(matching_dst_rows.GID(i));
    err = edst.ExtractMyRowView(rowB, NumEntriesB, ValuesB, IndicesB);
    if (err != 0) FOUR_C_THROW("ExtractMyRowView error: %d", err);

    // loop through the columns in source matrix and find respective place in destination
    for (int jA = 0, jB = 0; jA < NumEntriesA; ++jA)
    {
      // skip entries belonging to a different block of the logical block matrix
      if (selector[IndicesA[jA]] == 0.) continue;

      const int col = lidvector_[IndicesA[jA]];
      if (col == -1)
      {
        if (exactmatch)
          FOUR_C_THROW("gid %d not found in map for lid %d at %d", srccolmap.GID(IndicesA[jA]),
              IndicesA[jA], jA);
        else
          continue;
      }

      // try linear search in B
      while (jB < NumEntriesB && IndicesB[jB] < col) ++jB;

      // did not find index in linear search (re-indexing from A.ColMap() to B.ColMap()
      // might pass through the indices differently), try binary search
      if (jB == NumEntriesB || IndicesB[jB] != col)
        jB = std::lower_bound(IndicesB, IndicesB + NumEntriesB, col) - IndicesB;

      // not found, sparsity pattern of B does not contain the index from A -> terminate
      if (jB == NumEntriesB || IndicesB[jB] != col)
      {
        FOUR_C_THROW(
            "Source matrix entry with global row ID %d and global column ID %d couldn't be added to"
            " destination matrix entry with global row ID %d and unknown global column ID %d!",
            esrc.RowMap().GID(i), srccolmap.GID(IndicesA[jA]), matching_dst_rows.GID(i),
            edst.ColMap().GID(col));
      }

      ValuesB[jB] += ValuesA[jA] * scale;
    }
  }
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Coupling::Adapter::MatrixLogicalSplitAndTransform::add_into_unfilled(Epetra_CrsMatrix& esrc,
    const Epetra_Map& logical_range_map, const Epetra_Map& logical_domain_map,
    const Core::LinAlg::Vector<double>& selector, const Epetra_Map& matching_dst_rows,
    Epetra_CrsMatrix& edst, bool exactmatch, double scale)
{
  const Epetra_Map& srccolmap = esrc.ColMap();

  // standard code for the unfilled case
  std::vector<int> idx;
  std::vector<double> vals;
  int rows = logical_range_map.NumMyElements();
  for (int i = 0; i < rows; ++i)
  {
    int NumEntries;
    double* Values;
    int* Indices;
    int err = esrc.ExtractMyRowView(
        esrc.RowMap().LID(logical_range_map.GID(i)), NumEntries, Values, Indices);
    if (err != 0) FOUR_C_THROW("ExtractMyRowView error: %d", err);

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
        if (exactmatch)
          FOUR_C_THROW("gid %d not found in map for lid %d at %d", gid, Indices[j], j);
      }
    }

    NumEntries = vals.size();
    const int globalRow = matching_dst_rows.GID(i);

    // put row into matrix
    //
    // We might want to preserve a Dirichlet row in our destination matrix
    // here as well. Skip for now.

    if (edst.NumAllocatedGlobalEntries(globalRow) == 0)
    {
      int err = edst.InsertGlobalValues(globalRow, NumEntries, vals.data(), idx.data());
      if (err < 0) FOUR_C_THROW("InsertGlobalValues error: %d", err);
    }
    else
      for (int j = 0; j < NumEntries; ++j)
      {
        // add all values, including zeros, as we need a proper matrix graph
        int err = edst.SumIntoGlobalValues(globalRow, 1, &vals[j], &idx[j]);
        if (err > 0)
        {
          err = edst.InsertGlobalValues(globalRow, 1, &vals[j], &idx[j]);
          if (err < 0) FOUR_C_THROW("InsertGlobalValues error: %d", err);
        }
        else if (err < 0)
          FOUR_C_THROW("SumIntoGlobalValues error: %d", err);
      }
  }
}



bool Coupling::Adapter::MatrixRowTransform::operator()(const Core::LinAlg::SparseMatrix& src,
    double scale, const CouplingConverter& converter, Core::LinAlg::SparseMatrix& dst,
    bool addmatrix)
{
  return transformer_(
      src, src.range_map(), src.domain_map(), scale, &converter, nullptr, dst, false, addmatrix);
}



bool Coupling::Adapter::MatrixColTransform::operator()(const Epetra_Map&, const Epetra_Map&,
    const Core::LinAlg::SparseMatrix& src, double scale, const CouplingConverter& converter,
    Core::LinAlg::SparseMatrix& dst, bool exactmatch, bool addmatrix)
{
  return transformer_(src, src.range_map(), src.domain_map(), scale, nullptr, &converter, dst,
      exactmatch, addmatrix);
}



bool Coupling::Adapter::MatrixRowColTransform::operator()(const Core::LinAlg::SparseMatrix& src,
    double scale, const CouplingConverter& rowconverter, const CouplingConverter& colconverter,
    Core::LinAlg::SparseMatrix& dst, bool exactmatch, bool addmatrix)
{
  return transformer_(src, src.range_map(), src.domain_map(), scale, &rowconverter, &colconverter,
      dst, exactmatch, addmatrix);
}

FOUR_C_NAMESPACE_CLOSE
