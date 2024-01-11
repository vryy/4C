/*----------------------------------------------------------------------*/
/*! \file

\brief block sparse matrix implementation

\level 0


*----------------------------------------------------------------------*/

#include "baci_linalg_blocksparsematrix.H"

#include "baci_linalg_utils_densematrix_communication.H"
#include "baci_linalg_utils_sparse_algebra_manipulation.H"

#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_Transpose_RowMatrix.h>
#include <Teuchos_TimeMonitor.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::BlockSparseMatrixBase::BlockSparseMatrixBase(const MultiMapExtractor& domainmaps,
    const MultiMapExtractor& rangemaps, int npr, bool explicitdirichlet, bool savegraph)
    : domainmaps_(domainmaps), rangemaps_(rangemaps), usetranspose_(false)
{
  blocks_.reserve(Rows() * Cols());

  // add sparse matrices in row,column order
  for (int r = 0; r < Rows(); ++r)
  {
    for (int c = 0; c < Cols(); ++c)
    {
      blocks_.emplace_back(RangeMap(r), npr, explicitdirichlet, savegraph);
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CORE::LINALG::BlockSparseMatrixBase::Destroy(bool throw_exception_for_blocks)
{
  /// destroy matrix blocks
  for (auto& block : blocks_)
  {
    block.Destroy(throw_exception_for_blocks);
  }
  /// destroy full matrix row map
  if (fullrowmap_.strong_count() > 1)
  {
    dserror("fullrowmap_ cannot be finally deleted - any RCP (%i>1) still points to it",
        fullrowmap_.strong_count());
  }
  fullrowmap_ = Teuchos::null;

  /// destroy full matrix column map
  if (fullcolmap_.strong_count() > 1)
  {
    dserror("fullrowmap_ cannot be finally deleted - any RCP (%i>1) still points to it",
        fullrowmap_.strong_count());
  }
  fullcolmap_ = Teuchos::null;

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::SparseMatrix> CORE::LINALG::BlockSparseMatrixBase::Merge(
    bool explicitdirichlet) const
{
  TEUCHOS_FUNC_TIME_MONITOR("CORE::LINALG::BlockSparseMatrixBase::Merge");

  const SparseMatrix& m00 = Matrix(0, 0);

  Teuchos::RCP<SparseMatrix> sparse =
      Teuchos::rcp(new SparseMatrix(*fullrowmap_, m00.MaxNumEntries(), explicitdirichlet));
  for (const auto& block : blocks_)
  {
    sparse->Add(block, false, 1.0, 1.0);
  }
  if (Filled())
  {
    sparse->Complete(FullDomainMap(), FullRangeMap());
  }
  return sparse;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::BlockSparseMatrixBase::Assign(
    int r, int c, DataAccess access, const SparseMatrix& mat)
{
#ifdef BACI_DEBUG
  if (not Matrix(r, c).RowMap().SameAs(mat.RowMap())) dserror("cannot assign nonmatching matrices");
#endif
  Matrix(r, c).Assign(access, mat);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::BlockSparseMatrixBase::Zero()
{
  for (auto& block : blocks_) block.Zero();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::BlockSparseMatrixBase::Reset()
{
  for (int i = 0; i < Rows(); ++i)
  {
    for (int j = 0; j < Cols(); ++j)
    {
      Matrix(i, j).Reset();
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::BlockSparseMatrixBase::Complete(bool enforce_complete)
{
  for (int r = 0; r < Rows(); ++r)
  {
    for (int c = 0; c < Cols(); ++c)
    {
      Matrix(r, c).Complete(DomainMap(c), RangeMap(r), enforce_complete);
    }
  }

  fullrowmap_ = Teuchos::rcp(new Epetra_Map(*(rangemaps_.FullMap())));

  if (fullcolmap_ == Teuchos::null)
  {
    // build full col map
    std::vector<int> colmapentries;
    for (int c = 0; c < Cols(); ++c)
    {
      for (int r = 0; r < Rows(); ++r)
      {
        const Epetra_Map& colmap = Matrix(r, c).ColMap();
        colmapentries.insert(colmapentries.end(), colmap.MyGlobalElements(),
            colmap.MyGlobalElements() + colmap.NumMyElements());
      }
    }
    std::sort(colmapentries.begin(), colmapentries.end());
    colmapentries.erase(
        std::unique(colmapentries.begin(), colmapentries.end()), colmapentries.end());
    fullcolmap_ =
        Teuchos::rcp(new Epetra_Map(-1, colmapentries.size(), colmapentries.data(), 0, Comm()));
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::BlockSparseMatrixBase::Complete(
    const Epetra_Map& domainmap, const Epetra_Map& rangemap, bool enforce_complete)
{
  dserror("Complete with arguments not supported for block matrices");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CORE::LINALG::BlockSparseMatrixBase::Filled() const
{
  for (const auto& block : blocks_)
    if (not block.Filled()) return false;
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::BlockSparseMatrixBase::UnComplete()
{
  for (auto& block : blocks_) block.UnComplete();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::BlockSparseMatrixBase::ApplyDirichlet(
    const Epetra_Vector& dbctoggle, bool diagonalblock)
{
  int rows = Rows();
  int cols = Cols();
  for (int rblock = 0; rblock < rows; ++rblock)
  {
    Teuchos::RCP<Epetra_Vector> rowtoggle = rangemaps_.ExtractVector(dbctoggle, rblock);
    for (int cblock = 0; cblock < cols; ++cblock)
    {
      CORE::LINALG::SparseMatrix& bmat = Matrix(rblock, cblock);
      bmat.ApplyDirichlet(*rowtoggle, diagonalblock and rblock == cblock);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::BlockSparseMatrixBase::ApplyDirichlet(
    const Epetra_Map& dbcmap, bool diagonalblock)
{
  const int rows = Rows();
  const int cols = Cols();
  for (int rblock = 0; rblock < rows; ++rblock)
  {
    for (int cblock = 0; cblock < cols; ++cblock)
    {
      CORE::LINALG::SparseMatrix& bmat = Matrix(rblock, cblock);
      bmat.ApplyDirichlet(dbcmap, diagonalblock and rblock == cblock);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CORE::LINALG::BlockSparseMatrixBase::IsDbcApplied(
    const Epetra_Map& dbcmap, bool diagonalblock, const CORE::LINALG::SparseMatrix* trafo) const
{
  const int rows = Rows();
  const int cols = Cols();

  for (int rblock = 0; rblock < rows; ++rblock)
  {
    for (int cblock = 0; cblock < cols; ++cblock)
    {
      if (not Matrix(rblock, cblock)
                  .IsDbcApplied(dbcmap, diagonalblock and (rblock == cblock), trafo))
        return false;
    }
  }
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int CORE::LINALG::BlockSparseMatrixBase::SetUseTranspose(bool UseTranspose)
{
  for (auto& block : blocks_) block.SetUseTranspose(UseTranspose);
  usetranspose_ = UseTranspose;
  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int CORE::LINALG::BlockSparseMatrixBase::Apply(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  int rows = Rows();
  int cols = Cols();
  Y.PutScalar(0.0);

  if (not UseTranspose())
  {
    for (int rblock = 0; rblock < rows; ++rblock)
    {
      Teuchos::RCP<Epetra_MultiVector> rowresult = rangemaps_.Vector(rblock, Y.NumVectors());
      Teuchos::RCP<Epetra_MultiVector> rowy = rangemaps_.Vector(rblock, Y.NumVectors());
      for (int cblock = 0; cblock < cols; ++cblock)
      {
        Teuchos::RCP<Epetra_MultiVector> colx = domainmaps_.ExtractVector(X, cblock);
        const CORE::LINALG::SparseMatrix& bmat = Matrix(rblock, cblock);
        int err = bmat.Apply(*colx, *rowy);
        if (err != 0)
          dserror("failed to apply vector to matrix block (%d,%d): err=%d", rblock, cblock, err);
        rowresult->Update(1.0, *rowy, 1.0);
      }
      rangemaps_.InsertVector(*rowresult, rblock, Y);
    }
  }
  else
  {
    for (int rblock = 0; rblock < cols; ++rblock)
    {
      Teuchos::RCP<Epetra_MultiVector> rowresult = rangemaps_.Vector(rblock, Y.NumVectors());
      Teuchos::RCP<Epetra_MultiVector> rowy = rangemaps_.Vector(rblock, Y.NumVectors());
      for (int cblock = 0; cblock < rows; ++cblock)
      {
        Teuchos::RCP<Epetra_MultiVector> colx = domainmaps_.ExtractVector(X, cblock);
        const CORE::LINALG::SparseMatrix& bmat = Matrix(cblock, rblock);
        int err = bmat.Apply(*colx, *rowy);
        if (err != 0) dserror("failed to apply vector to matrix: err=%d", err);
        rowresult->Update(1.0, *rowy, 1.0);
      }
      rangemaps_.InsertVector(*rowresult, rblock, Y);
    }
  }

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int CORE::LINALG::BlockSparseMatrixBase::ApplyInverse(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  dserror("CORE::LINALG::BlockSparseMatrixBase::ApplyInverse not implemented");
  return -1;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::BlockSparseMatrixBase::Add(const CORE::LINALG::SparseOperator& A,
    const bool transposeA, const double scalarA, const double scalarB)
{
  A.AddOther(*this, transposeA, scalarA, scalarB);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::BlockSparseMatrixBase::AddOther(CORE::LINALG::BlockSparseMatrixBase& A,
    const bool transposeA, const double scalarA, const double scalarB) const
{
  A.Add(*this, transposeA, scalarA, scalarB);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::BlockSparseMatrixBase::AddOther(CORE::LINALG::SparseMatrixBase& A,
    const bool transposeA, const double scalarA, const double scalarB) const
{
  dserror("BlockSparseMatrix and SparseMatrix cannot be added");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::BlockSparseMatrixBase::Add(const CORE::LINALG::BlockSparseMatrixBase& A,
    const bool transposeA, const double scalarA, const double scalarB)
{
  for (int i = 0; i < Rows(); i++)
  {
    for (int j = 0; j < Cols(); j++)
    {
      if (transposeA)
        Matrix(i, j).Add(A.Matrix(j, i), transposeA, scalarA, scalarB);
      else
        Matrix(i, j).Add(A.Matrix(i, j), transposeA, scalarA, scalarB);
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int CORE::LINALG::BlockSparseMatrixBase::Scale(double ScalarConstant)
{
  for (int i = 0; i < Rows(); i++)
  {
    for (int j = 0; j < Cols(); j++)
    {
      int err = Matrix(i, j).Scale(ScalarConstant);
      if (err != 0) dserror("Scaling of matrix block (%d,%d) failed", i, j);
    }
  }
  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int CORE::LINALG::BlockSparseMatrixBase::Multiply(
    bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (TransA) dserror("transpose multiply not implemented for BlockSparseMatrix");
  return Apply(X, Y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double CORE::LINALG::BlockSparseMatrixBase::NormInf() const { return -1; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* CORE::LINALG::BlockSparseMatrixBase::Label() const
{
  return "CORE::LINALG::BlockSparseMatrixBase";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CORE::LINALG::BlockSparseMatrixBase::UseTranspose() const { return usetranspose_; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool CORE::LINALG::BlockSparseMatrixBase::HasNormInf() const { return false; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Comm& CORE::LINALG::BlockSparseMatrixBase::Comm() const
{
  return FullDomainMap().Comm();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map& CORE::LINALG::BlockSparseMatrixBase::OperatorDomainMap() const
{
  return FullDomainMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map& CORE::LINALG::BlockSparseMatrixBase::OperatorRangeMap() const
{
  return FullRangeMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::BlockSparseMatrixBase::GetPartialExtractor(
    const MultiMapExtractor& full_extractor, const std::vector<unsigned>& block_ids,
    MultiMapExtractor& partial_extractor) const
{
  const unsigned num_blocks = block_ids.size();

  Teuchos::RCP<Epetra_Map> full_map = Teuchos::null;

  std::vector<Teuchos::RCP<const Epetra_Map>> p_block_maps;
  p_block_maps.reserve(num_blocks);

  for (const int id : block_ids)
  {
    p_block_maps.push_back(full_extractor.Map(id));

    full_map = MergeMap(full_map, full_extractor.Map(id), false);
  }

  partial_extractor.Setup(*full_map, p_block_maps);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> CORE::LINALG::Multiply(
    const BlockSparseMatrixBase& A, bool transA, const BlockSparseMatrixBase& B, bool transB,
    bool explicitdirichlet, bool savegraph, bool completeoutput)
{
  if (!A.Filled() || !B.Filled())
    dserror("CORE::LINALG::BlockSparseMatrixBase::MatrixMultiyply: we expect A and B to be filled");

  if (A.Cols() != B.Rows() /*|| !A.FullDomainMap().SameAs(B.FullRowMap())*/)
  {
    dserror("CORE::LINALG::BlockSparseMatrixBase::MatrixMultiply: A and B not compatible");
  }


  int npr = 81;  // estimated number of entries per row in each block

  // generate result matrix
  Teuchos::RCP<CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>> C =
      Teuchos::rcp(new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
          B.DomainExtractor(), A.RangeExtractor(), npr, explicitdirichlet, savegraph));

  // nested loop over all blocks
  for (int i = 0; i < C->Rows(); i++)
  {
    for (int j = 0; j < C->Cols(); j++)
    {
      // build block C(i,j)
      Teuchos::RCP<SparseMatrix> Cij = Teuchos::rcp(
          new CORE::LINALG::SparseMatrix(C->RangeMap(i), npr, explicitdirichlet, savegraph));

      for (int l = 0; l < C->Cols(); l++)
      {
        // build submatrices for row i and j
        Teuchos::RCP<SparseMatrix> tmpij =
            CORE::LINALG::Multiply(A.Matrix(i, l), false, B.Matrix(l, j), false, true);
        Cij->Add(*tmpij, false, 1.0, 1.0);
      }

      // complete Cij with correct range and domain map
      if (completeoutput) Cij->Complete(C->DomainMap(j), C->RangeMap(i));

      // assign Cij block
      C->Assign(i, j, CORE::LINALG::View, *Cij);
    }
  }

  if (completeoutput) C->Complete();

  return Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(C);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>>
CORE::LINALG::BlockMatrix2x2(CORE::LINALG::SparseMatrix& A00, CORE::LINALG::SparseMatrix& A01,
    CORE::LINALG::SparseMatrix& A10, CORE::LINALG::SparseMatrix& A11)
{
  if (!A00.RangeMap().SameAs(A01.RangeMap()) || !A00.DomainMap().SameAs(A10.DomainMap()) ||
      !A01.DomainMap().SameAs(A11.DomainMap()) || !A10.RangeMap().SameAs(A11.RangeMap()))
    dserror("CORE::LINALG::BlockMatrix2x2: block maps are not compatible.");


  // generate range map
  std::vector<Teuchos::RCP<const Epetra_Map>> range_maps;
  range_maps.reserve(2);

  range_maps.emplace_back(Teuchos::rcp(new Epetra_Map(A00.RangeMap())));
  range_maps.emplace_back(Teuchos::rcp(new Epetra_Map(A10.RangeMap())));
  Teuchos::RCP<const Epetra_Map> range_map = MultiMapExtractor::MergeMaps(range_maps);
  Teuchos::RCP<MultiMapExtractor> rangeMMex =
      Teuchos::rcp(new MultiMapExtractor(*range_map, range_maps));

  // generate domain map
  std::vector<Teuchos::RCP<const Epetra_Map>> domain_maps;
  domain_maps.reserve(2);

  domain_maps.emplace_back(Teuchos::rcp(new Epetra_Map(A00.DomainMap())));
  domain_maps.emplace_back(Teuchos::rcp(new Epetra_Map(A01.DomainMap())));
  Teuchos::RCP<const Epetra_Map> domain_map = MultiMapExtractor::MergeMaps(domain_maps);
  Teuchos::RCP<MultiMapExtractor> domainMMex =
      Teuchos::rcp(new MultiMapExtractor(*domain_map, domain_maps));

  // generate result matrix
  Teuchos::RCP<CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>> C =
      Teuchos::rcp(new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
          *domainMMex, *rangeMMex));
  // Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> Cb =
  // Teuchos::rcp_dynamic_cast<CORE::LINALG::BlockSparseMatrixBase>(C);
  // assign matrices
  C->Assign(0, 0, CORE::LINALG::View, A00);
  C->Assign(0, 1, CORE::LINALG::View, A01);
  C->Assign(1, 0, CORE::LINALG::View, A10);
  C->Assign(1, 1, CORE::LINALG::View, A11);

  C->Complete();

  return C;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& CORE::LINALG::operator<<(
    std::ostream& os, const CORE::LINALG::BlockSparseMatrixBase& mat)
{
  for (int i = 0; i < mat.Rows(); ++i)
  {
    for (int j = 0; j < mat.Cols(); ++j)
    {
      if (mat.Comm().MyPID() == 0)
        os << "====================================Matrix block (" << i << "," << j
           << "):" << std::endl;
      fflush(stdout);
      os << mat(i, j);
    }
  }
  return os;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::DefaultBlockMatrixStrategy::DefaultBlockMatrixStrategy(BlockSparseMatrixBase& mat)
    : mat_(mat), scratch_lcols_(mat_.Rows())
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::DefaultBlockMatrixStrategy::Complete(bool enforce_complete)
{
  TEUCHOS_FUNC_TIME_MONITOR("CORE::LINALG::DefaultBlockMatrixStrategy::Complete");

  if (mat_.Filled() and not enforce_complete)
  {
    if (ghost_.size() != 0)
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
  for (int rblock = 0; rblock < rows; ++rblock)
  {
    const Epetra_Map& rowmap = mat_.RangeMap(rblock);

    for (int rlid = 0; rlid < rowmap.NumMyElements(); ++rlid)
    {
      int rgid = rowmap.GID(rlid);
      std::transform(ghost_[rgid].begin(), ghost_[rgid].end(), std::inserter(cgids, cgids.begin()),
          select1st<std::map<int, double>::value_type>());
    }
  }

  std::vector<int> cgidlist;
  cgidlist.reserve(cgids.size());
  cgidlist.assign(cgids.begin(), cgids.end());
  cgids.clear();

  // get to know the native processors of each ghost entry
  // this is expensive!

  std::vector<int> cpidlist(cgidlist.size());

  int err =
      mat_.FullDomainMap().RemoteIDList(cgidlist.size(), cgidlist.data(), cpidlist.data(), nullptr);
  if (err != 0) dserror("RemoteIDList failed");

  const Epetra_Comm& comm = mat_.FullRangeMap().Comm();
  const int numproc = comm.NumProc();

  // Send the ghost gids to their respective processor to ask for the domain
  // map the gids belong to.

  std::vector<std::vector<int>> ghostgids(comm.NumProc());
  for (unsigned i = 0; i < cgidlist.size(); ++i)
  {
    ghostgids[cpidlist[i]].push_back(cgidlist[i]);
  }

  cpidlist.clear();
  cgidlist.clear();

  std::vector<std::vector<int>> requests;
  AllToAllCommunication(comm, ghostgids, requests);

  // Now all gids are at the processors that own them. Lets find the owning
  // block for each of them.

  std::vector<std::vector<int>> block(comm.NumProc());

  for (int proc = 0; proc < numproc; ++proc)
  {
    for (unsigned i = 0; i < requests[proc].size(); ++i)
    {
      int gid = requests[proc][i];
      for (int cblock = 0; cblock < cols; ++cblock)
      {
        // assume row and range equal domain
        const Epetra_Map& domainmap = mat_.DomainMap(cblock);
        if (domainmap.MyGID(gid))
        {
          block[proc].push_back(cblock);
          break;
        }
      }

      if (block[proc].size() != i + 1)
      {
        dserror("gid %d not owned by any domain map", gid);
      }
    }
  }

  // communicate our findings back
  requests.clear();
  AllToAllCommunication(comm, block, requests);
  block.clear();

  // store domain block number for each ghost gid

  std::map<int, int> ghostmap;
  for (int proc = 0; proc < numproc; ++proc)
  {
    if (requests[proc].size() != ghostgids[proc].size())
    {
      dserror("size mismatch panic");
    }

    for (unsigned i = 0; i < requests[proc].size(); ++i)
    {
      int cblock = requests[proc][i];
      int cgid = ghostgids[proc][i];

      if (ghostmap.find(cgid) != ghostmap.end())
        dserror("column gid %d defined more often that once", cgid);

      ghostmap[cgid] = cblock;
    }
  }

  requests.clear();
  ghostgids.clear();

  // and finally do the assembly of ghost entries

  for (auto& irow : ghost_)
  {
    // most stupid way to find the right row
    int rgid = irow.first;
    int rblock = RowBlock(rgid);
    if (rblock == -1) dserror("row finding panic");

    for (auto& icol : irow.second)
    {
      int cgid = icol.first;
      if (ghostmap.find(cgid) == ghostmap.end()) dserror("unknown ghost gid %d", cgid);

      int cblock = ghostmap[cgid];
      double val = icol.second;

      SparseMatrix& matrix = mat_.Matrix(rblock, cblock);
      matrix.Assemble(val, rgid, cgid);
    }
  }

  ghost_.clear();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase>
CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(
    Teuchos::RCP<CORE::LINALG::SparseOperator> input_matrix)
{
  auto block_matrix = Teuchos::rcp_dynamic_cast<CORE::LINALG::BlockSparseMatrixBase>(input_matrix);
  dsassert(block_matrix != Teuchos::null, "Matrix is not a block matrix!");

  return block_matrix;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const CORE::LINALG::BlockSparseMatrixBase>
CORE::LINALG::CastToConstBlockSparseMatrixBaseAndCheckSuccess(
    Teuchos::RCP<const CORE::LINALG::SparseOperator> input_matrix)
{
  auto block_matrix =
      Teuchos::rcp_dynamic_cast<const CORE::LINALG::BlockSparseMatrixBase>(input_matrix);
  dsassert(block_matrix != Teuchos::null, "Matrix is not a block matrix!");

  return block_matrix;
}

BACI_NAMESPACE_CLOSE
