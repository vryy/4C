


#include <vector>
#include <iterator>
#include "fsi_matrixtransform.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_exporter.H"
#include "../linalg/linalg_blocksparsematrix.H"


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
    //
    // We might want to preserve a Dirichlet row in our destination matrix
    // here as well. Skip for now.

    for (int j=0; j<NumEntries; ++j)
    {
      // add all values, including zeros, as we need a proper matrix graph
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
  else
  {
    const Epetra_Map& dstcolmap = edst->ColMap();
    for (int j=0; j<NumEntries; ++j)
    {
      int gid = Indices[j];
      int lid = dstcolmap.LID(gid);
      if ( lid < 0 )
      {
        dserror("illegal local id: lid=%d, gid=%d",lid,gid);
      }
      Indices[j] = lid;
    }

    // We have to care for Dirichlet conditions in the filled destination
    // matrix. If there is just the diagonal entry we assume the row to be
    // Dirichlet and force the source matrix to have a Dirichlet row as
    // well. If this is not desired another flag would be needed.

    int lid = dstrowmap.LID(dstmap.GID(row));

    int myNumEntries;
    double *myValues;
    int *myIndices;
    int err = edst->ExtractMyRowView(lid, myNumEntries, myValues, myIndices);
    if (err)
      dserror("ExtractMyRowView error: %d on row lid=%d.\nI'm totally lost here.",err,lid);

    if (myNumEntries>=NumEntries)
    {
      // The normal case. This has to match.

      err = edst->SumIntoMyValues(lid, NumEntries, const_cast<double*>(Values), Indices);
      if (err)
      {
        const Epetra_Comm& comm = dstrowmap.Comm();
        for (int i=0; i<comm.NumProc(); ++i)
        {
          if (i==comm.MyPID())
          {
            std::cout << "PROC " << i << ":\n";
            std::cout << "actual line: ";
            std::copy(myIndices, myIndices+myNumEntries, std::ostream_iterator<int>(std::cout, " "));
            std::cout << "\ngiven  line: ";
            std::copy(Indices, Indices+NumEntries, std::ostream_iterator<int>(std::cout, " "));
            std::cout << "\n";
          }
          comm.Barrier();
        }

        dserror("SumIntoMyValues error: %d on row lid=%d.\nMaybe unfill of matrix block is needed.", err, lid);
      }
    }
    else if (myNumEntries==1)
    {
      // we have a dirichlet line in our destination matrix
      if (myIndices[0]!=lid)
        dserror("Single entry row without diagonal value. Confused.");
      for (int j=0; j<NumEntries; ++j)
      {
        if (Indices[j]==lid)
        {
          err = edst->SumIntoMyValues(lid, 1, const_cast<double*>(&Values[j]), &Indices[j]);
          if (err)
            dserror("SumIntoMyValues error: %d on row lid=%d.\nThis is not supposed to happen.", err, lid);
        }
        else
        {
          if (Values[j]!=0.0)
          {
            dserror("Attempt to add a non-Dirichlet row to a filled Dirichlet row.\n"
                    "Did you specify all required boundary conditions?");
          }
        }
      }
    }
    else
    {
      dserror("Filled destination matrix row shorter than src row: %d, %d. Panic.", myNumEntries, NumEntries);
    }
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
                                           const ADAPTER::CouplingConverter& converter,
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
                                            const ADAPTER::CouplingConverter& converter,
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
                                             bool addmatrix,
                                             double scale)
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
    vals.reserve(NumEntries);

    for (int j=0; j<NumEntries; ++j)
    {
      int gid = srccolmap.GID(Indices[j]);
      std::map<int,int>::const_iterator iter = gidmap_.find(gid);
      if (iter!=gidmap_.end())
      {
        idx.push_back(iter->second);
        vals.push_back(Values[j]*scale);
      }
      else
      {
        // only complain if an exact match is demanded
        if (exactmatch)
          dserror("gid %d not found in map for lid %d at %d", gid, Indices[j], j);
      }
    }

    Values = &vals[0];
    NumEntries = vals.size();

    if (addmatrix)
      AddValues(edst,dstrowmap,dstmap,i,NumEntries,Values,&idx[0]);
    else
      InsertValues(edst,dstrowmap,dstmap,i,NumEntries,Values,&idx[0]);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool
FSI::UTILS::MatrixColTransform::operator()(const Epetra_Map& rowmap,
                                           const Epetra_Map& colmap,
                                           const LINALG::SparseMatrix& src,
                                           double scale,
                                           const ADAPTER::CouplingConverter& converter,
                                           LINALG::SparseMatrix& dst,
                                           bool exactmatch,
                                           bool addmatrix)
{
  SetupGidMap(rowmap,colmap,converter,src.Comm());

  if (not addmatrix)
    dst.Zero();

  Teuchos::RCP<Epetra_CrsMatrix> esrc = src.EpetraMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> edst = dst.EpetraMatrix();

  MatrixInsert(esrc,esrc->RowMap(),edst,exactmatch,addmatrix,scale);

  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool
FSI::UTILS::MatrixRowColTransform::operator()(const LINALG::SparseMatrix& src,
                                              double scale,
                                              const ADAPTER::CouplingConverter& rowconverter,
                                              const ADAPTER::CouplingConverter& colconverter,
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

  coltrans_.SetupGidMap(*colconverter.SrcMap(),permsrc->ColMap(),colconverter,src.Comm());
  coltrans_.MatrixInsert(permsrc,dstmap,edst,exactmatch,addmatrix,1.0);

  return true;
}


