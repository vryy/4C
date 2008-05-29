
#ifdef CCADISCRET

#include <vector>

#include "fsi_matrixtransform.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_exporter.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MatrixTransform::InsertValues(Teuchos::RCP<Epetra_CrsMatrix> edst,
                                        const Epetra_Map& dstrowmap,
                                        const Epetra_Map& dstmap,
                                        int row,
                                        int NumEntries,
                                        double *Values,
                                        int *Indices)
{
  if (not edst->Filled())
  {
    // put row into matrix
    int err = edst->InsertGlobalValues(dstmap.GID(row), NumEntries, Values, Indices);
    if (err)
      dserror("InsertGlobalValues error: %d", err);
  }
  else
  {
    const Epetra_Map& dstcolmap = edst->ColMap();
    for (int j=0; j<NumEntries; ++j)
    {
      Indices[j] = dstcolmap.LID(Indices[j]);
    }

    int lid = dstrowmap.LID(dstmap.GID(row));
    int err = edst->ReplaceMyValues(lid, NumEntries, Values, Indices);
    if (err)
      dserror("ReplaceMyValues error: %d", err);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix>
FSI::MatrixRowTransform::Redistribute(const LINALG::SparseMatrix& src,
                                      const Epetra_Map& permsrcmap)
{
  if (exporter_==Teuchos::null)
  {
    exporter_ = Teuchos::rcp(new Epetra_Export(permsrcmap, src.RowMap()));
  }

  Teuchos::RCP<Epetra_CrsMatrix> permsrc = Teuchos::rcp(new Epetra_CrsMatrix(Copy,permsrcmap,src.MaxNumEntries()));
  int err = permsrc->Import(*src.EpetraMatrix(),*exporter_,Insert);
  if (err)
    dserror("Import failed with err=%d",err);

  permsrc->FillComplete(src.DomainMap(),permsrcmap);
  return permsrc;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void
FSI::MatrixRowTransform::MatrixInsert(Teuchos::RCP<Epetra_CrsMatrix> permsrc,
                                      const Epetra_Map& dstmap,
                                      Teuchos::RCP<Epetra_CrsMatrix> edst)
{
  Epetra_Map dstrowmap = edst->RowMap();
  Epetra_Map srccolmap = permsrc->ColMap();

  int rows = permsrc->NumMyRows();
  for (int i=0; i<rows; ++i)
  {
    int NumEntries;
    double *Values;
    int *Indices;
    int err = permsrc->ExtractMyRowView(i, NumEntries, Values, Indices);
    if (err!=0)
      dserror("ExtractMyRowView error: %d", err);

    // pull indices back to global
    std::vector<int> idx(NumEntries);
    for (int j=0; j<NumEntries; ++j)
    {
      idx[j] = srccolmap.GID(Indices[j]);
    }

    InsertValues(edst,dstrowmap,dstmap,i,NumEntries,Values,&idx[0]);

#if 0
    if (not edst->Filled())
    {
      // put row into matrix
      err = edst->InsertGlobalValues(dstmap.GID(i), NumEntries, Values, &idx[0]);
      if (err)
        dserror("InsertGlobalValues error: %d", err);
    }
    else
    {
      const Epetra_Map& dstcolmap = edst->ColMap();
      for (int j=0; j<NumEntries; ++j)
      {
        idx[j] = dstcolmap.LID(idx[j]);
      }

      int lid = dstrowmap.LID(dstmap.GID(i));
      err = edst->ReplaceMyValues(lid, NumEntries, Values, &idx[0]);
      if (err)
        dserror("ReplaceMyValues error: %d", err);
    }
#endif
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool
FSI::MatrixRowTransform::operator()(const LINALG::SparseMatrix& src,
                                    double scale,
                                    const ADAPTER::Coupling::Converter& converter,
                                    LINALG::SparseMatrix& dst)
{
  const Epetra_Map& permsrcmap = *converter.PermSrcMap();
  const Epetra_Map& dstmap     = *converter.DstMap();

  Teuchos::RCP<Epetra_CrsMatrix> permsrc = Redistribute(src,permsrcmap);
  permsrc->Scale(scale);

  dst.Zero();
  Teuchos::RCP<Epetra_CrsMatrix> edst = dst.EpetraMatrix();

  MatrixInsert(permsrc,dstmap,edst);

  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool
FSI::MatrixColTransform::operator()(const LINALG::BlockSparseMatrixBase& fullsrc,
                                    const LINALG::SparseMatrix& src,
                                    double scale,
                                    const ADAPTER::Coupling::Converter& converter,
                                    LINALG::SparseMatrix& dst,
                                    bool exactmatch)
{
  if (not havegidmap_)
  {
    DRT::Exporter ex(fullsrc.FullRowMap(),fullsrc.FullColMap(),src.Comm());
    converter.FillSrcToDstMap(gidmap_);
    ex.Export(gidmap_);
    havegidmap_ = true;
  }

  const Epetra_Map& rowmap = src.RowMap();
  const Epetra_Map& colmap = src.ColMap();

  dst.Zero();

  Teuchos::RCP<Epetra_CrsMatrix> esrc = src.EpetraMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> edst = dst.EpetraMatrix();

  int rows = esrc->NumMyRows();
  for (int i=0; i<rows; ++i)
  {
    int NumEntries;
    double *Values;
    int *Indices;
    int err = esrc->ExtractMyRowView(i, NumEntries, Values, Indices);
    if (err!=0)
      dserror("ExtractMyRowView error: %d", err);

    std::vector<int> idx;
    std::vector<double> vals;
    idx.reserve(NumEntries);
    if (not exactmatch)
      vals.reserve(NumEntries);

    for (int j=0; j<NumEntries; ++j)
    {
      int gid = colmap.GID(Indices[j]);
      std::map<int,int>::const_iterator iter = gidmap_.find(gid);
      if (iter!=gidmap_.end())
      {
        idx.push_back(iter->second);
        if (not exactmatch)
          vals.push_back(Values[j]);
      }
      else
      {
        // only complain if an exact match is demanded
        if (exactmatch)
          dserror("gid %d not found in map for lid %d at %d", gid, Indices[j], j);
      }
    }

    // If we do not demand an exact match, some values of this row might have
    // been skipped. Thus we redirect the value pointer and reset the number
    // of entries on this row.
    if (not exactmatch)
    {
      Values = &vals[0];
      NumEntries = vals.size();
    }

    InsertValues(edst,rowmap,rowmap,i,NumEntries,Values,&idx[0]);
  }

  if (scale!=1.)
    edst->Scale(scale);

  return true;
}


#endif
