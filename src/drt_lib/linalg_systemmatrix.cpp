#ifdef CCADISCRET

#include "linalg_systemmatrix.H"
#include "linalg_utils.H"
#include "drt_dserror.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::SystemMatrix::~SystemMatrix()
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::SingleSystemMatrix::SingleSystemMatrix(bool realdirichlet)
  : realdirichlet_(realdirichlet)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::SingleSystemMatrix::~SingleSystemMatrix()
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::SingleSystemMatrix::Setup(const Epetra_Map& rowmap, const int npr)
{
  sysmat_ = CreateMatrix(rowmap, npr);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::SingleSystemMatrix::Zero()
{
  if (mask_==Teuchos::null)
  {
    const Epetra_Map& rowmap = sysmat_->RowMap();
    int mne = sysmat_->MaxNumEntries();
    sysmat_ = CreateMatrix(rowmap, mne);
  }
  else
  {
    const Epetra_Map domainmap = sysmat_->DomainMap();
    const Epetra_Map rangemap = sysmat_->RangeMap();
    sysmat_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *mask_));
    sysmat_->FillComplete(domainmap,rangemap);
  }
}


/*----------------------------------------------------------------------*
 |  assemble a matrix  (public)                               popp 01/08|
 *----------------------------------------------------------------------*/
void LINALG::SingleSystemMatrix::Assemble(const Epetra_SerialDenseMatrix& Aele,
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

  // this 'Assemble' is not implemented for a Filled() matrix A
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
      if (!(rowmap.MyGID(rgid))) dserror("Sparse matrix A does not have global row %d",rgid);

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
void LINALG::SingleSystemMatrix::Assemble(double val, int rgid, int cgid)
{
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
void LINALG::SingleSystemMatrix::Complete()
{
  Complete(sysmat_->OperatorDomainMap(),sysmat_->OperatorRangeMap());
}


/*----------------------------------------------------------------------*
 |  FillComplete a matrix  (public)                          mwgee 01/08|
 *----------------------------------------------------------------------*/
void  LINALG::SingleSystemMatrix::Complete(const Epetra_Map& domainmap, const Epetra_Map& rangemap)
{
  if (sysmat_->Filled()) return;

  // keep mask for further use
  if (mask_==Teuchos::null)
  {
    mask_ = Teuchos::rcp(new Epetra_CrsGraph(sysmat_->Graph()));
  }

  int err = sysmat_->FillComplete(domainmap,rangemap,true);
  if (err) dserror("Epetra_CrsMatrix::FillComplete(domain,range) returned err=%d",err);
  return;
}


/*----------------------------------------------------------------------*
 |  Apply dirichlet conditions  (public)                     mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::SingleSystemMatrix::ApplyDirichlet(const Teuchos::RCP<Epetra_Vector> dbctoggle)
{
  const Epetra_Vector& dbct = *dbctoggle;

  if (realdirichlet_)
  {
    // allocate a new matrix and copy all rows that are not dirichlet
    const Epetra_Map& rowmap = sysmat_->RowMap();
    const int nummyrows      = sysmat_->NumMyRows();
    const int maxnumentries  = sysmat_->MaxNumEntries();

    Teuchos::RCP<Epetra_CrsMatrix> Anew = LINALG::CreateMatrix(rowmap,maxnumentries);
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
    const int nummyrows      = sysmat_->NumMyRows();
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
LINALG::BlockSystemMatrixBase::BlockSystemMatrixBase(const Epetra_Map& fullrangemap,
                                                     const Epetra_Map& fulldomainmap,
                                                     int rows,
                                                     int cols,
                                                     std::vector<Epetra_Map> rangemaps,
                                                     std::vector<Epetra_Map> domainmaps)
  : fullrangemap_(fullrangemap),
    fulldomainmap_(fulldomainmap),
    rows_(rows),
    cols_(cols),
    rangemaps_(rangemaps),
    domainmaps_(domainmaps)
{
  blocks_.resize(rows*cols);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::BlockSystemMatrixBase::Zero()
{
  for (unsigned i=0; i<blocks_.size(); ++i)
    blocks_[i].Zero();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::BlockSystemMatrixBase::Complete()
{
  for (unsigned i=0; i<blocks_.size(); ++i)
    blocks_[i].Complete();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::BlockSystemMatrixBase::Complete(const Epetra_Map& domainmap, const Epetra_Map& rangemap)
{
  dserror("Complete with arguments not supported for block matrices");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool LINALG::BlockSystemMatrixBase::Filled() const
{
  for (unsigned i=0; i<blocks_.size(); ++i)
    if (not blocks_[i].Filled())
      return false;
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::DefaultBlockMatrixCondition::DefaultBlockMatrixCondition(BlockSystemMatrixBase* mat)
  : mat_(mat)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::DefaultBlockMatrixCondition::RowBlock(int lrow, int rgid)
{
  int rows = mat_->Rows();
  for (int rblock=0; rblock<rows; ++rblock)
  {
    if (mat_->RangeMap(rblock).MyGID(rgid))
    {
      return rblock;
    }
  }
  return -1;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::DefaultBlockMatrixCondition::ColBlock(int rblock, int lcol, int cgid)
{
  int cols = mat_->Cols();
  for (int cblock = 0; cblock<cols; ++cblock)
  {
    SingleSystemMatrix& matrix = mat_->Matrix(rblock,cblock);

    // If we have a filled matrix we know the column map already.
    if (matrix.Filled())
    {
      if (matrix.ColMap().MyGID(cgid))
      {
        return cblock;
      }
    }

    // otherwise we can get just the non-ghost entries right now
    else if (matrix.RowMap().MyGID(cgid))
    {
      return cblock;
    }
  }

  // ghost entries in a non-filled matrix will have to be done later

  return -1;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::DefaultBlockMatrixCondition::Assemble(double val,
                                                   int lrow, int rgid, int rblock,
                                                   int lcol, int cgid, int cblock)
{
#ifdef DEBUG
  if (rblock==-1)
    dserror("no block entry found for row gid=%d",rgid);
#endif

  if (cblock>-1)
  {
    SingleSystemMatrix& matrix = mat_->Matrix(rblock,cblock);
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
void LINALG::DefaultBlockMatrixCondition::Complete()
{
  if (mat_->Filled())
  {
    if (ghost_.size()!=0)
    {
      dserror("no unresolved ghost entries in a filled block matrix allowed");
    }
    return;
  }

  // finish ghost entries

  int rows = mat_->Rows();
  int cols = mat_->Cols();

  std::set<int> cgids;

  // get the list of all ghost entries gids
  for (int rblock=0; rblock<rows; ++rblock)
  {
    const Epetra_Map& rowmap = mat_->RangeMap(rblock);

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

  int err = mat_->FullDomainMap().RemoteIDList(cgidlist.size(),&cgidlist[0],&cpidlist[0],&clidlist[0]);
  if (err!=0)
    dserror("RemoteIDList failed");

  // never mind the lids
  clidlist.clear();

  const Epetra_Comm& comm = mat_->FullRangeMap().Comm();
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
        const Epetra_Map& domainmap = mat_->DomainMap(cblock);
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

      SingleSystemMatrix& matrix = mat_->Matrix(rblock,cblock);
      matrix.Assemble(val,rgid,cgid);
    }
  }

  ghost_.clear();
}


#endif
