
#ifdef CCADISCRET

#include <vector>

#include "fsi_matrixtransform.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_exporter.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
static void
InsertValues(Teuchos::RCP<Epetra_CrsMatrix> edst,
             const Epetra_Map& dstrowmap,
             const Epetra_Map& dstmap,
             int row,
             int NumEntries,
             const double *Values,
             int *Indices)
{
  if (not edst->Filled())
  {
    // put row into matrix
    int err = edst->InsertGlobalValues(dstmap.GID(row), NumEntries, const_cast<double*>(Values), Indices);
    if (err<0)
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
    int err = edst->ReplaceMyValues(lid, NumEntries, const_cast<double*>(Values), Indices);
    if (err)
      dserror("ReplaceMyValues error: %d", err);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
static void
AddValues(Teuchos::RCP<Epetra_CrsMatrix> edst,
          const Epetra_Map& dstrowmap,
          const Epetra_Map& dstmap,
          int row,
          int NumEntries,
          const double *Values,
          int *Indices)
{
  if (not edst->Filled())
  {
    row = dstmap.GID(row);
    // put row into matrix
    for (int j=0; j<NumEntries; ++j)
    {
      if (Values[j]!=0.0)
      {
        int err = edst->SumIntoGlobalValues(row, 1, const_cast<double*>(&Values[j]), &Indices[j]);
        if (err>0)
        {
          err = edst->InsertGlobalValues(row, 1, const_cast<double*>(&Values[j]), &Indices[j]);
          if (err<0)
            dserror("InsertGlobalValues error: %d", err);
        }
        else if (err<0)
          dserror("SumIntoGlobalValues error: %d", err);
      }
    }
  }
  else
  {
    const Epetra_Map& dstcolmap = edst->ColMap();
    for (int j=0; j<NumEntries; ++j)
    {
      Indices[j] = dstcolmap.LID(Indices[j]);
    }

    int lid = dstrowmap.LID(dstmap.GID(row));
    int err = edst->SumIntoMyValues(lid, NumEntries, const_cast<double*>(Values), Indices);
    if (err)
      dserror("SumIntoMyValues error: %d", err);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix>
FSI::UTILS::MatrixRowTransform::Redistribute(const LINALG::SparseMatrix& src,
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
FSI::UTILS::MatrixRowTransform::MatrixInsert(Teuchos::RCP<Epetra_CrsMatrix> esrc,
                                             const Epetra_Map& dstmap,
                                             Teuchos::RCP<Epetra_CrsMatrix> edst,
                                             bool addmatrix)
{
  if (not esrc->Filled())
    dserror("filled source matrix expected");

  Epetra_Map dstrowmap = edst->RowMap();
  Epetra_Map srccolmap = esrc->ColMap();

  int rows = esrc->NumMyRows();
  for (int i=0; i<rows; ++i)
  {
    int NumEntries;
    double *Values;
    int *Indices;
    int err = esrc->ExtractMyRowView(i, NumEntries, Values, Indices);
    if (err!=0)
      dserror("ExtractMyRowView error: %d", err);

    // pull indices back to global
    std::vector<int> idx(NumEntries);
    for (int j=0; j<NumEntries; ++j)
    {
      idx[j] = srccolmap.GID(Indices[j]);
    }

    if (addmatrix)
      AddValues(edst,dstrowmap,dstmap,i,NumEntries,Values,&idx[0]);
    else
      InsertValues(edst,dstrowmap,dstmap,i,NumEntries,Values,&idx[0]);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool
FSI::UTILS::MatrixRowTransform::operator()(const LINALG::SparseMatrix& src,
                                           double scale,
                                           const ADAPTER::Coupling::Converter& converter,
                                           LINALG::SparseMatrix& dst,
                                           bool addmatrix)
{
  const Epetra_Map& permsrcmap = *converter.PermSrcMap();
  const Epetra_Map& dstmap     = *converter.DstMap();

  Teuchos::RCP<Epetra_CrsMatrix> permsrc = Redistribute(src,permsrcmap);
  permsrc->Scale(scale);

  if (not addmatrix)
    dst.Zero();

  Teuchos::RCP<Epetra_CrsMatrix> edst = dst.EpetraMatrix();
  MatrixInsert(permsrc,dstmap,edst,addmatrix);

  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void
FSI::UTILS::MatrixColTransform::SetupGidMap(const Epetra_Map& rowmap,
                                            const Epetra_Map& colmap,
                                            const ADAPTER::Coupling::Converter& converter,
                                            const Epetra_Comm& comm)
{
  if (not havegidmap_)
  {
    DRT::Exporter ex(rowmap,colmap,comm);
    converter.FillSrcToDstMap(gidmap_);
    ex.Export(gidmap_);
    havegidmap_ = true;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void
FSI::UTILS::MatrixColTransform::MatrixInsert(Teuchos::RCP<Epetra_CrsMatrix> esrc,
                                             const Epetra_Map& dstmap,
                                             Teuchos::RCP<Epetra_CrsMatrix> edst,
                                             bool exactmatch,
                                             bool addmatrix)
{
  if (not esrc->Filled())
    dserror("filled source matrix expected");

  const Epetra_Map& dstrowmap = edst->RowMap();
  const Epetra_Map& srccolmap = esrc->ColMap();

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
      int gid = srccolmap.GID(Indices[j]);
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

    if (addmatrix)
      AddValues(edst,dstrowmap,dstmap,i,NumEntries,Values,&idx[0]);
    else
      InsertValues(edst,dstrowmap,dstmap,i,NumEntries,Values,&idx[0]);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool
FSI::UTILS::MatrixColTransform::operator()(const LINALG::BlockSparseMatrixBase& fullsrc,
                                           const LINALG::SparseMatrix& src,
                                           double scale,
                                           const ADAPTER::Coupling::Converter& converter,
                                           LINALG::SparseMatrix& dst,
                                           bool exactmatch,
                                           bool addmatrix)
{
  SetupGidMap(fullsrc.FullRowMap(),fullsrc.FullColMap(),converter,src.Comm());

  if (not addmatrix)
    dst.Zero();

  Teuchos::RCP<Epetra_CrsMatrix> esrc = src.EpetraMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> edst = dst.EpetraMatrix();

  MatrixInsert(esrc,esrc->RowMap(),edst,exactmatch,addmatrix);

  if (scale!=1.)
    edst->Scale(scale);

  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool
FSI::UTILS::MatrixRowColTransform::operator()(const LINALG::SparseMatrix& src,
                                              double scale,
                                              const ADAPTER::Coupling::Converter& rowconverter,
                                              const ADAPTER::Coupling::Converter& colconverter,
                                              LINALG::SparseMatrix& dst,
                                              bool exactmatch,
                                              bool addmatrix)
{
  const Epetra_Map& permsrcmap = *rowconverter.PermSrcMap();
  const Epetra_Map& dstmap     = *rowconverter.DstMap();

  Teuchos::RCP<Epetra_CrsMatrix> permsrc = rowtrans_.Redistribute(src,permsrcmap);
  if (scale!=1.)
    permsrc->Scale(scale);

  if (not addmatrix)
    dst.Zero();

  Teuchos::RCP<Epetra_CrsMatrix> edst = dst.EpetraMatrix();

  coltrans_.SetupGidMap(*rowconverter.SrcMap(),permsrc->ColMap(),colconverter,src.Comm());
  coltrans_.MatrixInsert(permsrc,dstmap,edst,exactmatch,addmatrix);

  return true;
}


#endif
