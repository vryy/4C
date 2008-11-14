/*!----------------------------------------------------------------------
\file linalg_blocksparsematrix.cpp

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich
              
Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed, 
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de) 
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de                   

-------------------------------------------------------------------------
</pre>
<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "linalg_blocksparsematrix.H"
#include "linalg_utils.H"
#include "drt_dserror.H"

#include <EpetraExt_Transpose_RowMatrix.h>
#include <EpetraExt_MatrixMatrix.h>
#include <Teuchos_TimeMonitor.hpp>



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::BlockSparseMatrixBase::BlockSparseMatrixBase(const MultiMapExtractor& domainmaps,
                                                     const MultiMapExtractor& rangemaps,
                                                     int npr,
                                                     bool explicitdirichlet,
                                                     bool savegraph)
  : domainmaps_(domainmaps),
    rangemaps_(rangemaps)
{
  blocks_.reserve(Rows()*Cols());

  // add sparse matrices in row,column order
  for (int r=0; r<Rows(); ++r)
  {
    for (int c=0; c<Cols(); ++c)
    {
      blocks_.push_back(SparseMatrix(RangeMap(r),npr,explicitdirichlet,savegraph));
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> LINALG::BlockSparseMatrixBase::Merge() const
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::BlockSparseMatrixBase::Merge");

  const SparseMatrix& m00 = Matrix(0,0);
  Teuchos::RCP<SparseMatrix> sparse = Teuchos::rcp(new SparseMatrix(*fullrowmap_,
                                                                    m00.MaxNumEntries(),
                                                                    m00.ExplicitDirichlet(),
                                                                    m00.SaveGraph()));
  for (unsigned i=0; i<blocks_.size(); ++i)
  {
    sparse->Add(blocks_[i],false,1.0,1.0);
  }
  if (Filled())
  {
    sparse->Complete(FullDomainMap(),FullRangeMap());
  }
  return sparse;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::BlockSparseMatrixBase::Assign(int r, int c, Epetra_DataAccess access, SparseMatrix& mat)
{
#ifdef DEBUG
  if (not Matrix(r,c).RowMap().SameAs(mat.RowMap()))
    dserror("cannot assign nonmatching matrices");
#endif
  Matrix(r,c).Assign(access,mat);
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

  // FIXME: We already know the full range map from the MapExtractor. Here
  // we build it again and call it row map. Maybe we need the distinction
  // between row and range maps for block matrices one day. Right now it is
  // not supported but still there different maps for a start... (The worst of
  // both worlds.)

  if (fullrowmap_==Teuchos::null)
  {
    // build full row map
    int rowmaplength = 0;
    for (int r=0; r<Rows(); ++r)
    {
      rowmaplength += Matrix(r,0).RowMap().NumMyElements();
    }
    std::vector<int> rowmapentries;
    rowmapentries.reserve(rowmaplength);
    for (int r=0; r<Rows(); ++r)
    {
      const Epetra_Map& rowmap = Matrix(r,0).RowMap();
      copy(rowmap.MyGlobalElements(),
           rowmap.MyGlobalElements()+rowmap.NumMyElements(),
           back_inserter(rowmapentries));
    }
    fullrowmap_ = Teuchos::rcp(new Epetra_Map(-1,rowmapentries.size(),&rowmapentries[0],0,Comm()));
  }

  if (fullcolmap_==Teuchos::null)
  {
    // build full col map
    int colmaplength = 0;
    for (int c=0; c<Cols(); ++c)
    {
      colmaplength += Matrix(0,c).ColMap().NumMyElements();
    }
    std::vector<int> colmapentries;
    colmapentries.reserve(colmaplength);
    for (int c=0; c<Cols(); ++c)
    {
      const Epetra_Map& colmap = Matrix(0,c).ColMap();
      copy(colmap.MyGlobalElements(),
           colmap.MyGlobalElements()+colmap.NumMyElements(),
           back_inserter(colmapentries));
    }
    fullcolmap_ = Teuchos::rcp(new Epetra_Map(-1,colmapentries.size(),&colmapentries[0],0,Comm()));
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
void LINALG::BlockSparseMatrixBase::UnComplete()
{
  for (unsigned i=0; i<blocks_.size(); ++i)
    blocks_[i].UnComplete();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::BlockSparseMatrixBase::ApplyDirichlet(const Teuchos::RCP<Epetra_Vector> dbctoggle, bool diagonalblock)
{
  int rows = Rows();
  int cols = Cols();
  for (int rblock=0; rblock<rows; ++rblock)
  {
    Teuchos::RCP<Epetra_Vector> rowtoggle = rangemaps_.ExtractVector(dbctoggle,rblock);
    for (int cblock=0; cblock<cols; ++cblock)
    {
      LINALG::SparseMatrix& bmat = Matrix(rblock,cblock);
      bmat.ApplyDirichlet(rowtoggle,diagonalblock and rblock==cblock);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::BlockSparseMatrixBase::ApplyDirichlet(const Epetra_Map& dbcmap, bool diagonalblock)
{
  const int rows = Rows();
  const int cols = Cols();
  for (int rblock=0; rblock<rows; ++rblock)
  {
    for (int cblock=0; cblock<cols; ++cblock)
    {
      LINALG::SparseMatrix& bmat = Matrix(rblock,cblock);
      bmat.ApplyDirichlet(dbcmap,diagonalblock and rblock==cblock);
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
      rangemaps_.InsertVector(*rowresult,rblock,Y);
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
ostream& LINALG::operator << (ostream& os, const LINALG::BlockSparseMatrixBase& mat)
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
void LINALG::DefaultBlockMatrixStrategy::Complete()
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::DefaultBlockMatrixStrategy::Complete");

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


