#ifdef CCADISCRET

#include "linalg_systemmatrix.H"
#include "linalg_utils.H"
#include "drt_dserror.H"

#include <EpetraExt_Transpose_RowMatrix.h>
#include <EpetraExt_MatrixMatrix.h>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::SparseMatrix::SparseMatrix(bool explicitdirichlet, bool savegraph)
  : explicitdirichlet_(explicitdirichlet),
    savegraph_(savegraph)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::SparseMatrix::SparseMatrix(const Epetra_Map& rowmap, const int npr, bool explicitdirichlet, bool savegraph)
  : explicitdirichlet_(explicitdirichlet),
    savegraph_(savegraph)
{
  Setup(rowmap,npr);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::SparseMatrix::SparseMatrix(const Epetra_CrsMatrix& matrix, bool explicitdirichlet, bool savegraph)
  : explicitdirichlet_(explicitdirichlet),
    savegraph_(savegraph)
{
  sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(matrix));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::SparseMatrix::SparseMatrix(Teuchos::RCP<Epetra_CrsMatrix> matrix, bool explicitdirichlet, bool savegraph)
  : sysmat_(matrix),
    explicitdirichlet_(explicitdirichlet),
    savegraph_(savegraph)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::SparseMatrix::SparseMatrix(const SparseMatrix& mat)
  : explicitdirichlet_(mat.explicitdirichlet_),
    savegraph_(mat.savegraph_)
{
  if (mat.sysmat_!=Teuchos::null)
    sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(*mat.sysmat_));

  if (mat.graph_!=Teuchos::null)
    graph_ = Teuchos::rcp(new Epetra_CrsGraph(*mat.graph_));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::SparseMatrix::SparseMatrix(const Epetra_Vector& diag, bool explicitdirichlet, bool savegraph)
  : explicitdirichlet_(explicitdirichlet),
    savegraph_(savegraph)
{
  int length = diag.Map().NumMyElements();
  Epetra_Map map(-1,length,diag.Map().MyGlobalElements(),
                 diag.Map().IndexBase(),diag.Comm());
  Setup(map,1);
  for (int i=0; i<length; ++i)
  {
    int gid = diag.Map().GID(i);
    Assemble(diag[i],gid,gid);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::SparseMatrix::~SparseMatrix()
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::SparseMatrix::Setup(const Epetra_Map& rowmap, const int npr)
{
  if (!rowmap.UniqueGIDs())
    dserror("Row map is not unique");
  sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,rowmap,npr,false));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::SparseMatrix::Zero()
{
  // graph_!=Teuchos::null if savegraph_==false only
  if (graph_==Teuchos::null)
  {
    const Epetra_Map& rowmap = sysmat_->RowMap();
    int mne = sysmat_->MaxNumEntries();
    sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,rowmap,mne,false));
  }
  else
  {
    const Epetra_Map domainmap = sysmat_->DomainMap();
    const Epetra_Map rangemap = sysmat_->RangeMap();
    sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph_));
    sysmat_->FillComplete(domainmap,rangemap);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::SparseMatrix::Reset()
{
  Epetra_Map rowmap = sysmat_->RowMap();
  int maxnumentries = sysmat_->MaxNumEntries();
  sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,rowmap,maxnumentries,false));
  graph_ = Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  assemble a matrix  (public)                               popp 01/08|
 *----------------------------------------------------------------------*/
void LINALG::SparseMatrix::Assemble(const Epetra_SerialDenseMatrix& Aele,
                                          const std::vector<int>& lmrow,
                                          const std::vector<int>& lmrowowner,
                                          const std::vector<int>& lmcol)
{
  const int lrowdim = (int)lmrow.size();
  const int lcoldim = (int)lmcol.size();
  if (lrowdim!=(int)lmrowowner.size() || lrowdim!=Aele.M() || lcoldim!=Aele.N())
    dserror("Mismatch in dimensions");

  const int myrank = sysmat_->Comm().MyPID();
  const Epetra_Map& rowmap = sysmat_->RowMap();
  const Epetra_Map& colmap = sysmat_->ColMap();

  if (sysmat_->Filled())
  {
    // loop rows of local matrix
    for (int lrow=0; lrow<lrowdim; ++lrow)
    {
      // check ownership of row
      if (lmrowowner[lrow] != myrank) continue;

      // check whether I have that global row
      int rgid = lmrow[lrow];
      int rlid = rowmap.LID(rgid);
#ifdef DEBUG
      if (rlid<0) dserror("Sparse matrix A does not have global row %d",rgid);
#endif

      for (int lcol=0; lcol<lcoldim; ++lcol)
      {
        double val = Aele(lrow,lcol);
        int cgid = lmcol[lcol];
        int clid = colmap.LID(cgid);
#ifdef DEBUG
        if (clid<0) dserror("Sparse matrix A does not have global column %d",cgid);
#endif
        int errone = sysmat_->SumIntoMyValues(rlid,1,&val,&clid);
        if (errone)
          dserror("Epetra_CrsMatrix::SumIntoMyValues returned error code %d",errone);
      } // for (int lcol=0; lcol<ldim; ++lcol)
    } // for (int lrow=0; lrow<ldim; ++lrow)
  }
  else
  {
    // loop rows of local matrix
    for (int lrow=0; lrow<lrowdim; ++lrow)
    {
      // check ownership of row
      if (lmrowowner[lrow] != myrank) continue;

      // check whether I have that global row
      int rgid = lmrow[lrow];
#ifdef DEBUG
      if (!(rowmap.MyGID(rgid))) dserror("Sparse matrix A does not have global row %d",rgid);
#endif

      for (int lcol=0; lcol<lcoldim; ++lcol)
      {
        double val = Aele(lrow,lcol);
        int cgid = lmcol[lcol];

        // Now that we do not rebuild the sparse mask in each step, we
        // are bound to assemble the whole thing. Zeros included.
        int errone = sysmat_->SumIntoGlobalValues(rgid,1,&val,&cgid);
        if (errone>0)
        {
          int errtwo = sysmat_->InsertGlobalValues(rgid,1,&val,&cgid);
          if (errtwo<0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned error code %d",errtwo);
        }
        else if (errone)
          dserror("Epetra_CrsMatrix::SumIntoGlobalValues returned error code %d",errone);
      } // for (int lcol=0; lcol<lcoldim; ++lcol)
    } // for (int lrow=0; lrow<lrowdim; ++lrow)
  }
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::SparseMatrix::Assemble(double val, int rgid, int cgid)
{
  // SumIntoGlobalValues works for filled matrices as well!
  int errone = sysmat_->SumIntoGlobalValues(rgid,1,&val,&cgid);
  if (errone>0)
  {
    int errtwo = sysmat_->InsertGlobalValues(rgid,1,&val,&cgid);
    if (errtwo<0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned error code %d",errtwo);
  }
  else if (errone)
    dserror("Epetra_CrsMatrix::SumIntoGlobalValues returned error code %d",errone);
}


/*----------------------------------------------------------------------*
 |  FillComplete a matrix  (public)                          mwgee 12/06|
 *----------------------------------------------------------------------*/
void LINALG::SparseMatrix::Complete()
{
  if (sysmat_->Filled()) return;

  int err = sysmat_->FillComplete(true);
  if (err) dserror("Epetra_CrsMatrix::FillComplete() returned err=%d",err);

  // keep mask for further use
  if (savegraph_ and graph_==Teuchos::null)
  {
    graph_ = Teuchos::rcp(new Epetra_CrsGraph(sysmat_->Graph()));
  }
}


/*----------------------------------------------------------------------*
 |  FillComplete a matrix  (public)                          mwgee 01/08|
 *----------------------------------------------------------------------*/
void  LINALG::SparseMatrix::Complete(const Epetra_Map& domainmap, const Epetra_Map& rangemap)
{
  if (sysmat_->Filled()) return;

  int err = sysmat_->FillComplete(domainmap,rangemap,true);
  if (err) dserror("Epetra_CrsMatrix::FillComplete(domain,range) returned err=%d",err);

  // keep mask for further use
  if (savegraph_ and graph_==Teuchos::null)
  {
    graph_ = Teuchos::rcp(new Epetra_CrsGraph(sysmat_->Graph()));
  }
}


/*----------------------------------------------------------------------*
 |  Apply dirichlet conditions  (public)                     mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::SparseMatrix::ApplyDirichlet(const Teuchos::RCP<Epetra_Vector> dbctoggle)
{
  if (not Filled())
    dserror("expect filled matrix to apply dirichlet conditions");

  const Epetra_Vector& dbct = *dbctoggle;

  if (explicitdirichlet_)
  {
    // Save graph of original matrix if not done already.
    // This will never happen as the matrix is guaranteed to be filled. But to
    // make the code more explicit...
    if (savegraph_ and graph_==Teuchos::null)
    {
      graph_ = Teuchos::rcp(new Epetra_CrsGraph(sysmat_->Graph()));
      if (not graph_->Filled())
        dserror("got unfilled graph from filled matrix");
    }

    // allocate a new matrix and copy all rows that are not dirichlet
    const Epetra_Map& rowmap = sysmat_->RowMap();
    const int nummyrows      = sysmat_->NumMyRows();
    const int maxnumentries  = sysmat_->MaxNumEntries();

    Teuchos::RCP<Epetra_CrsMatrix> Anew = Teuchos::rcp(new Epetra_CrsMatrix(Copy,rowmap,maxnumentries,false));
    vector<int> indices(maxnumentries,0);
    vector<double> values(maxnumentries,0.0);
    for (int i=0; i<nummyrows; ++i)
    {
      int row = sysmat_->GRID(i);
      if (dbct[i]!=1.0)
      {
        int numentries;
        int err = sysmat_->ExtractGlobalRowCopy(row,maxnumentries,numentries,&values[0],&indices[0]);
#ifdef DEBUG
        if (err) dserror("Epetra_CrsMatrix::ExtractGlobalRowCopy returned err=%d",err);
#endif
        err = Anew->InsertGlobalValues(row,numentries,&values[0],&indices[0]);
#ifdef DEBUG
        if (err<0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned err=%d",err);
#endif
      }
      else
      {
        double one = 1.0;
        int err = Anew->InsertGlobalValues(row,1,&one,&row);
#ifdef DEBUG
        if (err<0) dserror("Epetra_CrsMatrix::InsertGlobalValues returned err=%d",err);
#endif
      }
    }
    sysmat_ = Anew;
    Complete();
  }
  else
  {
    const int nummyrows = sysmat_->NumMyRows();
    for (int i=0; i<nummyrows; ++i)
    {
      if (dbct[i]==1.0)
      {
        int *indexOffset;
        int *indices;
        double *values;
        int err = sysmat_->ExtractCrsDataPointers(indexOffset, indices, values);
#ifdef DEBUG
        if (err) dserror("Epetra_CrsMatrix::ExtractCrsDataPointers returned err=%d",err);
#endif
        // zero row
        memset(&values[indexOffset[i]], 0,
               (indexOffset[i+1]-indexOffset[i])*sizeof(double));

        double one = 1.0;
        err = sysmat_->SumIntoMyValues(i,1,&one,&i);
#ifdef DEBUG
        if (err<0) dserror("Epetra_CrsMatrix::SumIntoMyValues returned err=%d",err);
#endif
      }
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrix::SetUseTranspose(bool UseTranspose)
{
  return sysmat_->SetUseTranspose(UseTranspose);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrix::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  return sysmat_->Apply(X,Y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrix::ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  return sysmat_->ApplyInverse(X,Y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* LINALG::SparseMatrix::Label() const
{
  return "LINALG::SparseMatrix";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool LINALG::SparseMatrix::UseTranspose() const
{
  return sysmat_->UseTranspose();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool LINALG::SparseMatrix::HasNormInf() const
{
  return sysmat_->HasNormInf();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Comm& LINALG::SparseMatrix::Comm() const
{
  return sysmat_->Comm();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map& LINALG::SparseMatrix::OperatorDomainMap() const
{
  return sysmat_->OperatorDomainMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map& LINALG::SparseMatrix::OperatorRangeMap() const
{
  return sysmat_->OperatorRangeMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrix::MaxNumEntries() const
{
  return sysmat_->MaxNumEntries();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double LINALG::SparseMatrix::NormInf() const
{
  return sysmat_->NormInf();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double LINALG::SparseMatrix::NormOne() const
{
  return sysmat_->NormOne();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double LINALG::SparseMatrix::NormFrobenius() const
{
  return sysmat_->NormFrobenius();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrix::Multiply(bool TransA, const Epetra_Vector &x, Epetra_Vector &y) const
{
  return sysmat_->Multiply(TransA,x,y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrix::Multiply(bool TransA, const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  return sysmat_->Multiply(TransA,X,Y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrix::LeftScale(const Epetra_Vector &x)
{
  return sysmat_->LeftScale(x);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrix::RightScale(const Epetra_Vector &x)
{
  return sysmat_->RightScale(x);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrix::PutScalar(double ScalarConstant)
{
  return sysmat_->PutScalar(ScalarConstant);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrix::Scale(double ScalarConstant)
{
  return sysmat_->Scale(ScalarConstant);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrix::ReplaceDiagonalValues(const Epetra_Vector &Diagonal)
{
  return sysmat_->ReplaceDiagonalValues(Diagonal);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::SparseMatrix::ExtractDiagonalCopy(Epetra_Vector &Diagonal) const
{
  return sysmat_->ExtractDiagonalCopy(Diagonal);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> LINALG::SparseMatrix::Transpose()
{
  if (not Filled()) dserror("FillComplete was not called on matrix");

  EpetraExt::RowMatrix_Transpose trans;
  Epetra_CrsMatrix* Aprime = &(dynamic_cast<Epetra_CrsMatrix&>(trans(*sysmat_)));
  return Teuchos::rcp(new SparseMatrix(*Aprime,explicitdirichlet_,savegraph_));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::SparseMatrix::Add(const LINALG::SparseMatrix& A,
                               const bool transposeA,
                               const double scalarA,
                               const double scalarB)
{
  if (!A.Filled()) dserror("FillComplete was not called on A");
  if (Filled()) dserror("FillComplete was called on me before");

  Epetra_CrsMatrix* Aprime = NULL;
  RCP<EpetraExt::RowMatrix_Transpose> Atrans;
  if (transposeA)
  {
    Atrans = rcp(new EpetraExt::RowMatrix_Transpose(false,NULL,false));
    Aprime = &(dynamic_cast<Epetra_CrsMatrix&>(((*Atrans)(const_cast<Epetra_CrsMatrix&>(*A.sysmat_)))));
  }
  else
  {
    Aprime = const_cast<Epetra_CrsMatrix*>(&*A.sysmat_);
  }

  if (scalarB == 0.0)
    sysmat_->PutScalar(0.0);
  else if (scalarB != 1.0)
    sysmat_->Scale(scalarB);

  //Loop over Aprime's rows and sum into
  int MaxNumEntries = EPETRA_MAX( Aprime->MaxNumEntries(), sysmat_->MaxNumEntries() );
  int NumEntries;
  vector<int>    Indices(MaxNumEntries);
  vector<double> Values(MaxNumEntries);

  const int NumMyRows = Aprime->NumMyRows();
  int Row, err;
  if (scalarA)
  {
    for( int i = 0; i < NumMyRows; ++i )
    {
      Row = Aprime->GRID(i);
      int ierr = Aprime->ExtractGlobalRowCopy(Row,MaxNumEntries,NumEntries,&Values[0],&Indices[0]);
      if (ierr) dserror("Epetra_CrsMatrix::ExtractGlobalRowCopy returned err=%d",ierr);
      if (scalarA != 1.0)
        for (int j = 0; j < NumEntries; ++j) Values[j] *= scalarA;
      for (int j=0; j<NumEntries; ++j)
      {
        err = sysmat_->SumIntoGlobalValues(Row,1,&Values[j],&Indices[j]);
        if (err<0 || err==2)
          err = sysmat_->InsertGlobalValues(Row,1,&Values[j],&Indices[j]);
        if (err < 0)
          dserror("Epetra_CrsMatrix::InsertGlobalValues returned err=%d",err);
      }
    }
  }
}


/*----------------------------------------------------------------------*
  (private)
 *----------------------------------------------------------------------*/
void LINALG::SparseMatrix::Split2x2(BlockSparseMatrixBase& Abase)
{
  if (Abase.Rows() != 2 || Abase.Cols() != 2) dserror("Can only split in 2x2 system");
  if (!Filled()) dserror("SparsMatrix must be filled");
  Teuchos::RCP<Epetra_CrsMatrix> A   = EpetraMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> A11 = Abase(0,0).EpetraMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> A12 = Abase(0,1).EpetraMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> A21 = Abase(1,0).EpetraMatrix();
  Teuchos::RCP<Epetra_CrsMatrix> A22 = Abase(1,1).EpetraMatrix();
  if (A11->Filled() || A12->Filled() || A21->Filled() || A22->Filled())
    dserror("Block matrix may not be filled on input");
  const Epetra_Comm& Comm    = Abase.Comm();
  const Epetra_Map&  A11rmap = A11->RangeMap();
  const Epetra_Map&  A11dmap = A11->DomainMap();
  const Epetra_Map&  A22rmap = A22->RangeMap();
  const Epetra_Map&  A22dmap = A22->DomainMap();

  //----------------------------- create a parallel redundant map of A11domainmap
  map<int,int> a11gmap;
  {
    vector<int> a11global(A11dmap.NumGlobalElements());
    int count=0;
    for (int proc=0; proc<Comm.NumProc(); ++proc)
    {
      int length = 0;
      if (proc==Comm.MyPID())
      {
        for (int i=0; i<A11dmap.NumMyElements(); ++i)
        {
          a11global[count+length] = A11dmap.GID(i);
          ++length;
        }
      }
      Comm.Broadcast(&length,1,proc);
      Comm.Broadcast(&a11global[count],length,proc);
      count += length;
    }
    if (count != A11dmap.NumGlobalElements())
      dserror("SparseMatrix::Split2x2: mismatch in dimensions");

    // create the map
    for (int i=0; i<count; ++i)
      a11gmap[a11global[i]] = 1;
    a11global.clear();
  }

  vector<int>    gcindices(A->MaxNumEntries());
  vector<double> gvalues(A->MaxNumEntries());

  //----------------------------------------------------- create matrix A11
  if (A11rmap.NumGlobalElements()>0 && A11dmap.NumGlobalElements()>0)
  {
    {
      for (int i=0; i<A->NumMyRows(); ++i)
      {
        const int grid = A->GRID(i);
        if (A11rmap.MyGID(grid)==false) continue;
        int     numentries;
        double* values;
        int*    cindices;
        int err = A->ExtractMyRowView(i,numentries,values,cindices);
        if (err) dserror("SparseMatrix::Split2x2: A->ExtractMyRowView returned %i",err);
        int count=0;
        for (int j=0; j<numentries; ++j)
        {
          const int gcid = A->ColMap().GID(cindices[j]);
          // see whether we have gcid as part of a11gmap
          map<int,int>::iterator curr = a11gmap.find(gcid);
          if (curr==a11gmap.end()) continue;
          gcindices[count] = gcid;
          gvalues[count] = values[j];
          ++count;
        }
        err = A11->InsertGlobalValues(grid,count,&gvalues[0],&gcindices[0]);
        if (err<0) dserror("SparseMatrix::Split2x2: A->InsertGlobalValues returned %i",err);
      }
    }
  }

  //--------------------------------------------------- create matrix A22
  if (A22rmap.NumGlobalElements()>0 && A22dmap.NumGlobalElements()>0)
  {
    for (int i=0; i<A->NumMyRows(); ++i)
    {
      const int grid = A->GRID(i);
      if (A22rmap.MyGID(grid)==false) continue;
      int     numentries;
      double* values;
      int*    cindices;
      int err = A->ExtractMyRowView(i,numentries,values,cindices);
      if (err) dserror("SparseMatrix::Split2x2: A->ExtractMyRowView returned %i",err);
      int count=0;
      for (int j=0; j<numentries; ++j)
      {
        const int gcid = A->ColMap().GID(cindices[j]);
        // see whether we have gcid as part of a11gmap
        map<int,int>::iterator curr = a11gmap.find(gcid);
        if (curr!=a11gmap.end()) continue;
        gcindices[count] = gcid;
        gvalues[count]    = values[j];
        ++count;
      }
      err = A22->InsertGlobalValues(grid,count,&gvalues[0],&gcindices[0]);
      if (err<0) dserror("SparseMatrix::Split2x2: A->InsertGlobalValues returned %i",err);
    }
  }


  //---------------------------------------------------- create matrix A12
  if (A11rmap.NumGlobalElements()>0 && A22dmap.NumGlobalElements()>0)
  {
    for (int i=0; i<A->NumMyRows(); ++i)
    {
      const int grid = A->GRID(i);
      if (A11rmap.MyGID(grid)==false) continue;
      int     numentries;
      double* values;
      int*    cindices;
      int err = A->ExtractMyRowView(i,numentries,values,cindices);
      if (err) dserror("SparseMatrix::Split2x2: A->ExtractMyRowView returned %i",err);
      int count=0;
      for (int j=0; j<numentries; ++j)
      {
        const int gcid = A->ColMap().GID(cindices[j]);
        // see whether we have gcid as part of a11gmap
        map<int,int>::iterator curr = a11gmap.find(gcid);
        if (curr!=a11gmap.end()) continue;
        gcindices[count] = gcid;
        gvalues[count] = values[j];
        ++count;
      }
      err = A12->InsertGlobalValues(grid,count,&gvalues[0],&gcindices[0]);
      if (err<0) dserror("SparseMatrix::Split2x2: A->InsertGlobalValues returned %i",err);
    }
  }

  //----------------------------------------------------------- create A21
  if (A22rmap.NumGlobalElements()>0 && A11dmap.NumGlobalElements()>0)
  {
    for (int i=0; i<A->NumMyRows(); ++i)
    {
      const int grid = A->GRID(i);
      if (A22rmap.MyGID(grid)==false) continue;
      int     numentries;
      double* values;
      int*    cindices;
      int err = A->ExtractMyRowView(i,numentries,values,cindices);
      if (err) dserror("SparseMatrix::Split2x2: A->ExtractMyRowView returned %i",err);
      int count=0;
      for (int j=0; j<numentries; ++j)
      {
        const int gcid = A->ColMap().GID(cindices[j]);
        // see whether we have gcid as part of a11gmap
        map<int,int>::iterator curr = a11gmap.find(gcid);
        if (curr==a11gmap.end()) continue;
        gcindices[count] = gcid;
        gvalues[count] = values[j];
        ++count;
      }
      err = A21->InsertGlobalValues(grid,count,&gvalues[0],&gcindices[0]);
      if (err<0) dserror("SparseMatrix::Split2x2: A->InsertGlobalValues returned %i",err);
    }
  }

  // Do not complete BlockMatrix
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const LINALG::SparseMatrix& mat)
{
  os << *(const_cast<LINALG::SparseMatrix&>(mat).EpetraMatrix());
  return os;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> LINALG::Multiply(const LINALG::SparseMatrix& A,
                                                    bool transA,
                                                    const LINALG::SparseMatrix& B,
                                                    bool transB,
                                                    bool completeoutput)
{
  // make sure FillComplete was called on the matrices
  if (!A.Filled()) dserror("A has to be FillComplete");
  if (!B.Filled()) dserror("B has to be FillComplete");

  // create resultmatrix with correct rowmap
  Teuchos::RCP<LINALG::SparseMatrix> C;
  if (!transA)
    C = Teuchos::rcp(new SparseMatrix(A.RangeMap(),20,A.explicitdirichlet_,A.savegraph_));
  else
    C = Teuchos::rcp(new SparseMatrix(A.DomainMap(),20,A.explicitdirichlet_,A.savegraph_));

  int err = EpetraExt::MatrixMatrix::Multiply(*A.sysmat_,transA,*B.sysmat_,transB,*C->sysmat_,completeoutput);
  if (err) dserror("EpetraExt::MatrixMatrix::Multiply returned err = &d",err);

  return C;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::BlockSparseMatrixBase::BlockSparseMatrixBase(const MultiMapExtractor& domainmaps,
                                                     const MultiMapExtractor& rangemaps)

  : domainmaps_(domainmaps),
    rangemaps_(rangemaps)
{
  blocks_.resize(Rows()*Cols());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::BlockSparseMatrixBase::Zero()
{
  for (unsigned i=0; i<blocks_.size(); ++i)
    blocks_[i].Zero();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::BlockSparseMatrixBase::Complete()
{
  for (int r=0; r<Rows(); ++r)
  {
    for (int c=0; c<Cols(); ++c)
    {
      Matrix(r,c).Complete(DomainMap(c),RangeMap(r));
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::BlockSparseMatrixBase::Complete(const Epetra_Map& domainmap, const Epetra_Map& rangemap)
{
  dserror("Complete with arguments not supported for block matrices");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool LINALG::BlockSparseMatrixBase::Filled() const
{
  for (unsigned i=0; i<blocks_.size(); ++i)
    if (not blocks_[i].Filled())
      return false;
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::BlockSparseMatrixBase::ApplyDirichlet(const Teuchos::RCP<Epetra_Vector> dbctoggle)
{
  int rows = Rows();
  int cols = Cols();
  for (int rblock=0; rblock<rows; ++rblock)
  {
    Teuchos::RCP<Epetra_Vector> rowtoggle = rangemaps_.ExtractVector(dbctoggle,rblock);
    for (int cblock=0; cblock<cols; ++cblock)
    {
      LINALG::SparseMatrix& bmat = Matrix(rblock,cblock);
      bmat.ApplyDirichlet(rowtoggle);
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::BlockSparseMatrixBase::SetUseTranspose(bool UseTranspose)
{
  if (UseTranspose)
    dserror("transposed block matrix not implemented");
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::BlockSparseMatrixBase::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  int rows = Rows();
  int cols = Cols();
  Y.PutScalar(0.0);

  if (not UseTranspose())
  {
    for (int rblock=0; rblock<rows; ++rblock)
    {
      Teuchos::RCP<Epetra_MultiVector> rowresult = rangemaps_.Vector(rblock,Y.NumVectors());
      Teuchos::RCP<Epetra_MultiVector> rowy      = rangemaps_.Vector(rblock,Y.NumVectors());
      for (int cblock=0; cblock<cols; ++cblock)
      {
        Teuchos::RCP<Epetra_MultiVector> colx = domainmaps_.ExtractVector(X,cblock);
        const LINALG::SparseMatrix& bmat = Matrix(rblock,cblock);
        int err = bmat.Apply(*colx,*rowy);
        if (err!=0)
          dserror("failed to apply vector to matrix: err=%d",err);
        rowresult->Update(1.0,*rowy,1.0);
      }
      rangemaps_.InsertVector(*rowy,rblock,Y);
    }
  }
  else
  {
    dserror("transposed block matrices not supported");
  }

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::BlockSparseMatrixBase::ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  dserror("LINALG::BlockSparseMatrixBase::ApplyInverse not implemented");
  return -1;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double LINALG::BlockSparseMatrixBase::NormInf() const
{
  return -1;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* LINALG::BlockSparseMatrixBase::Label() const
{
  return "LINALG::BlockSparseMatrixBase";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool LINALG::BlockSparseMatrixBase::UseTranspose() const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool LINALG::BlockSparseMatrixBase::HasNormInf() const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Comm& LINALG::BlockSparseMatrixBase::Comm() const
{
  return FullDomainMap().Comm();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map& LINALG::BlockSparseMatrixBase::OperatorDomainMap() const
{
  return FullDomainMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map& LINALG::BlockSparseMatrixBase::OperatorRangeMap() const
{
  return FullRangeMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const LINALG::BlockSparseMatrixBase& mat)
{
  for (int i=0; i<mat.Rows(); ++i)
    for (int j=0; j<mat.Cols(); ++j)
    {
      if (mat.Comm().MyPID()==0)
        os << "====================================Matrix block (" << i << "," << j << "):" << endl;
      fflush(stdout);
      os << mat(i,j);
    }
  return os;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::DefaultBlockMatrixStrategy::DefaultBlockMatrixStrategy(BlockSparseMatrixBase& mat)
  : mat_(mat)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::DefaultBlockMatrixStrategy::RowBlock(int lrow, int rgid)
{
  int rows = mat_.Rows();
  for (int rblock=0; rblock<rows; ++rblock)
  {
    if (mat_.RangeMap(rblock).MyGID(rgid))
    {
      return rblock;
    }
  }
  return -1;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::DefaultBlockMatrixStrategy::ColBlock(int rblock, int lcol, int cgid)
{
  int cols = mat_.Cols();
  for (int cblock = 0; cblock<cols; ++cblock)
  {
    SparseMatrix& matrix = mat_.Matrix(rblock,cblock);

    // If we have a filled matrix we know the column map already.
    if (matrix.Filled())
    {
      if (matrix.ColMap().MyGID(cgid))
      {
        return cblock;
      }
    }

    // otherwise we can get just the non-ghost entries right now
    else if (matrix.DomainMap().MyGID(cgid))
    {
      return cblock;
    }
  }

  // ghost entries in a non-filled matrix will have to be done later

  return -1;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::DefaultBlockMatrixStrategy::Assemble(double val,
                                                   int lrow, int rgid, int rblock,
                                                   int lcol, int cgid, int cblock)
{
#ifdef DEBUG
  if (rblock==-1)
    dserror("no block entry found for row gid=%d",rgid);
#endif

  if (cblock>-1)
  {
    SparseMatrix& matrix = mat_.Matrix(rblock,cblock);
    matrix.Assemble(val,rgid,cgid);
  }
  else
  {
    // ghost entry in non-filled matrix. Save for later insertion.
    ghost_[rgid][cgid] += val;
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::DefaultBlockMatrixStrategy::Complete()
{
  if (mat_.Filled())
  {
    if (ghost_.size()!=0)
    {
      dserror("no unresolved ghost entries in a filled block matrix allowed");
    }
    return;
  }

  // finish ghost entries

  int rows = mat_.Rows();
  int cols = mat_.Cols();

  std::set<int> cgids;

  // get the list of all ghost entries gids
  for (int rblock=0; rblock<rows; ++rblock)
  {
    const Epetra_Map& rowmap = mat_.RangeMap(rblock);

    for (int rlid=0; rlid<rowmap.NumMyElements(); ++rlid)
    {
      int rgid = rowmap.GID(rlid);
      std::transform(ghost_[rgid].begin(),
                     ghost_[rgid].end(),
                     std::inserter(cgids,cgids.begin()),
                     select1st<std::map<int,double>::value_type>());
    }
  }

  std::vector<int> cgidlist;
  cgidlist.reserve(cgids.size());
  cgidlist.assign(cgids.begin(),cgids.end());
  cgids.clear();

  // get to know the native processors of each ghost entry
  // this is expensive!

  std::vector<int> cpidlist(cgidlist.size());
  std::vector<int> clidlist(cgidlist.size());

  int err = mat_.FullDomainMap().RemoteIDList(cgidlist.size(),&cgidlist[0],&cpidlist[0],&clidlist[0]);
  if (err!=0)
    dserror("RemoteIDList failed");

  // never mind the lids
  clidlist.clear();

  const Epetra_Comm& comm = mat_.FullRangeMap().Comm();
  const int numproc = comm.NumProc();

  // Send the ghost gids to their respective processor to ask for the domain
  // map the gids belong to.

  std::vector<std::vector<int> > ghostgids(comm.NumProc());
  for (unsigned i=0; i<cgidlist.size(); ++i)
  {
    ghostgids[cpidlist[i]].push_back(cgidlist[i]);
  }

  cpidlist.clear();
  cgidlist.clear();

  std::vector<std::vector<int> > requests;
  AllToAllCommunication(comm, ghostgids, requests);

  // Now all gids are at the processors that own them. Lets find the owning
  // block for each of them.

  std::vector<std::vector<int> > block(comm.NumProc());

  for (int proc=0; proc<numproc; ++proc)
  {
    for (unsigned i=0; i<requests[proc].size(); ++i)
    {
      int gid = requests[proc][i];
      for (int cblock=0; cblock<cols; ++cblock)
      {
        // assume row and range equal domain
        const Epetra_Map& domainmap = mat_.DomainMap(cblock);
        if (domainmap.MyGID(gid))
        {
          block[proc].push_back(cblock);
          break;
        }
      }

      if (block[proc].size()!=i+1)
      {
        dserror("gid %d not owned by any domain map",gid);
      }
    }
  }

  // communicate our findings back
  requests.clear();
  AllToAllCommunication(comm, block, requests);
  block.clear();

  // store domain block number for each ghost gid

  std::map<int,int> ghostmap;
  for (int proc=0; proc<numproc; ++proc)
  {
    if (requests[proc].size()!=ghostgids[proc].size())
    {
      dserror("size mismatch panic");
    }

    for (unsigned i=0; i<requests[proc].size(); ++i)
    {
      int cblock = requests[proc][i];
      int cgid = ghostgids[proc][i];

      if (ghostmap.find(cgid)!=ghostmap.end())
        dserror("column gid %d defined more often that once",cgid);

      ghostmap[cgid] = cblock;
    }
  }

  requests.clear();
  ghostgids.clear();

  // and finally do the assembly of ghost entries

  for (std::map<int,std::map<int,double> >::iterator irow=ghost_.begin();
       irow!=ghost_.end();
       ++irow)
  {
    // most stupid way to find the right row
    int rgid = irow->first;
    int rblock = RowBlock(0, rgid);
    if (rblock==-1)
      dserror("row finding panic");

    for (std::map<int,double>::iterator icol=irow->second.begin();
         icol!=irow->second.end();
         ++icol)
    {
      int cgid = icol->first;
      if (ghostmap.find(cgid)==ghostmap.end())
        dserror("unknown ghost gid %d",cgid);

      int cblock = ghostmap[cgid];
      double val = icol->second;

      SparseMatrix& matrix = mat_.Matrix(rblock,cblock);
      matrix.Assemble(val,rgid,cgid);
    }
  }

  ghost_.clear();
}


#endif
