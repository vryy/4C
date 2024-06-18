/*----------------------------------------------------------------------*/
/*! \file

\brief block sparse matrix implementation

\level 0


*----------------------------------------------------------------------*/

#include "4C_linalg_blocksparsematrix.hpp"

#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"

#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_Transpose_RowMatrix.h>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::LinAlg::BlockSparseMatrixBase::BlockSparseMatrixBase(const MultiMapExtractor& domainmaps,
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
bool Core::LinAlg::BlockSparseMatrixBase::Destroy(bool throw_exception_for_blocks)
{
  /// destroy matrix blocks
  for (auto& block : blocks_)
  {
    block.Destroy(throw_exception_for_blocks);
  }
  /// destroy full matrix row map
  if (fullrowmap_.strong_count() > 1)
  {
    FOUR_C_THROW("fullrowmap_ cannot be finally deleted - any RCP (%i>1) still points to it",
        fullrowmap_.strong_count());
  }
  fullrowmap_ = Teuchos::null;

  /// destroy full matrix column map
  if (fullcolmap_.strong_count() > 1)
  {
    FOUR_C_THROW("fullrowmap_ cannot be finally deleted - any RCP (%i>1) still points to it",
        fullrowmap_.strong_count());
  }
  fullcolmap_ = Teuchos::null;

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseMatrix> Core::LinAlg::BlockSparseMatrixBase::Merge(
    bool explicitdirichlet) const
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::LinAlg::BlockSparseMatrixBase::Merge");

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
void Core::LinAlg::BlockSparseMatrixBase::Assign(
    int r, int c, DataAccess access, const SparseMatrix& mat)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (not Matrix(r, c).RowMap().SameAs(mat.RowMap()))
    FOUR_C_THROW("cannot assign nonmatching matrices");
#endif
  Matrix(r, c).Assign(access, mat);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::Zero()
{
  for (auto& block : blocks_) block.Zero();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::reset()
{
  for (int i = 0; i < Rows(); ++i)
  {
    for (int j = 0; j < Cols(); ++j)
    {
      Matrix(i, j).reset();
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::Complete(bool enforce_complete)
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
void Core::LinAlg::BlockSparseMatrixBase::Complete(
    const Epetra_Map& domainmap, const Epetra_Map& rangemap, bool enforce_complete)
{
  FOUR_C_THROW("Complete with arguments not supported for block matrices");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Core::LinAlg::BlockSparseMatrixBase::Filled() const
{
  for (const auto& block : blocks_)
    if (not block.Filled()) return false;
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::UnComplete()
{
  for (auto& block : blocks_) block.UnComplete();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::ApplyDirichlet(
    const Epetra_Vector& dbctoggle, bool diagonalblock)
{
  int rows = Rows();
  int cols = Cols();
  for (int rblock = 0; rblock < rows; ++rblock)
  {
    Teuchos::RCP<Epetra_Vector> rowtoggle = rangemaps_.ExtractVector(dbctoggle, rblock);
    for (int cblock = 0; cblock < cols; ++cblock)
    {
      Core::LinAlg::SparseMatrix& bmat = Matrix(rblock, cblock);
      bmat.ApplyDirichlet(*rowtoggle, diagonalblock and rblock == cblock);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::ApplyDirichlet(
    const Epetra_Map& dbcmap, bool diagonalblock)
{
  const int rows = Rows();
  const int cols = Cols();
  for (int rblock = 0; rblock < rows; ++rblock)
  {
    for (int cblock = 0; cblock < cols; ++cblock)
    {
      Core::LinAlg::SparseMatrix& bmat = Matrix(rblock, cblock);
      bmat.ApplyDirichlet(dbcmap, diagonalblock and rblock == cblock);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Core::LinAlg::BlockSparseMatrixBase::IsDbcApplied(
    const Epetra_Map& dbcmap, bool diagonalblock, const Core::LinAlg::SparseMatrix* trafo) const
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
int Core::LinAlg::BlockSparseMatrixBase::SetUseTranspose(bool UseTranspose)
{
  for (auto& block : blocks_) block.SetUseTranspose(UseTranspose);
  usetranspose_ = UseTranspose;
  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::LinAlg::BlockSparseMatrixBase::Apply(
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
        const Core::LinAlg::SparseMatrix& bmat = Matrix(rblock, cblock);
        int err = bmat.Apply(*colx, *rowy);
        if (err != 0)
          FOUR_C_THROW(
              "failed to apply vector to matrix block (%d,%d): err=%d", rblock, cblock, err);
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
        const Core::LinAlg::SparseMatrix& bmat = Matrix(cblock, rblock);
        int err = bmat.Apply(*colx, *rowy);
        if (err != 0) FOUR_C_THROW("failed to apply vector to matrix: err=%d", err);
        rowresult->Update(1.0, *rowy, 1.0);
      }
      rangemaps_.InsertVector(*rowresult, rblock, Y);
    }
  }

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::LinAlg::BlockSparseMatrixBase::ApplyInverse(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  FOUR_C_THROW("Core::LinAlg::BlockSparseMatrixBase::ApplyInverse not implemented");
  return -1;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::Add(const Core::LinAlg::SparseOperator& A,
    const bool transposeA, const double scalarA, const double scalarB)
{
  A.AddOther(*this, transposeA, scalarA, scalarB);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::AddOther(Core::LinAlg::BlockSparseMatrixBase& A,
    const bool transposeA, const double scalarA, const double scalarB) const
{
  A.Add(*this, transposeA, scalarA, scalarB);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::AddOther(Core::LinAlg::SparseMatrixBase& A,
    const bool transposeA, const double scalarA, const double scalarB) const
{
  FOUR_C_THROW("BlockSparseMatrix and SparseMatrix cannot be added");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::Add(const Core::LinAlg::BlockSparseMatrixBase& A,
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
int Core::LinAlg::BlockSparseMatrixBase::Scale(double ScalarConstant)
{
  for (int i = 0; i < Rows(); i++)
  {
    for (int j = 0; j < Cols(); j++)
    {
      int err = Matrix(i, j).Scale(ScalarConstant);
      if (err != 0) FOUR_C_THROW("Scaling of matrix block (%d,%d) failed", i, j);
    }
  }
  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::LinAlg::BlockSparseMatrixBase::Multiply(
    bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (TransA) FOUR_C_THROW("transpose multiply not implemented for BlockSparseMatrix");
  return Apply(X, Y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double Core::LinAlg::BlockSparseMatrixBase::NormInf() const { return -1; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* Core::LinAlg::BlockSparseMatrixBase::Label() const
{
  return "Core::LinAlg::BlockSparseMatrixBase";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Core::LinAlg::BlockSparseMatrixBase::UseTranspose() const { return usetranspose_; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Core::LinAlg::BlockSparseMatrixBase::HasNormInf() const { return false; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Comm& Core::LinAlg::BlockSparseMatrixBase::Comm() const
{
  return FullDomainMap().Comm();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map& Core::LinAlg::BlockSparseMatrixBase::OperatorDomainMap() const
{
  return FullDomainMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map& Core::LinAlg::BlockSparseMatrixBase::OperatorRangeMap() const
{
  return FullRangeMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::get_partial_extractor(
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

  partial_extractor.setup(*full_map, p_block_maps);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> Core::LinAlg::Multiply(
    const BlockSparseMatrixBase& A, bool transA, const BlockSparseMatrixBase& B, bool transB,
    bool explicitdirichlet, bool savegraph, bool completeoutput)
{
  if (!A.Filled() || !B.Filled())
    FOUR_C_THROW(
        "Core::LinAlg::BlockSparseMatrixBase::MatrixMultiyply: we expect A and B to be filled");

  if (A.Cols() != B.Rows() /*|| !A.FullDomainMap().SameAs(B.FullRowMap())*/)
  {
    FOUR_C_THROW("Core::LinAlg::BlockSparseMatrixBase::MatrixMultiply: A and B not compatible");
  }


  int npr = 81;  // estimated number of entries per row in each block

  // generate result matrix
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>> C =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          B.DomainExtractor(), A.RangeExtractor(), npr, explicitdirichlet, savegraph));

  // nested loop over all blocks
  for (int i = 0; i < C->Rows(); i++)
  {
    for (int j = 0; j < C->Cols(); j++)
    {
      // build block C(i,j)
      Teuchos::RCP<SparseMatrix> Cij = Teuchos::rcp(
          new Core::LinAlg::SparseMatrix(C->RangeMap(i), npr, explicitdirichlet, savegraph));

      for (int l = 0; l < C->Cols(); l++)
      {
        // build submatrices for row i and j
        Teuchos::RCP<SparseMatrix> tmpij =
            Core::LinAlg::Multiply(A.Matrix(i, l), false, B.Matrix(l, j), false, true);
        Cij->Add(*tmpij, false, 1.0, 1.0);
      }

      // complete Cij with correct range and domain map
      if (completeoutput) Cij->Complete(C->DomainMap(j), C->RangeMap(i));

      // assign Cij block
      C->Assign(i, j, Core::LinAlg::View, *Cij);
    }
  }

  if (completeoutput) C->Complete();

  return Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(C);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>>
Core::LinAlg::BlockMatrix2x2(Core::LinAlg::SparseMatrix& A00, Core::LinAlg::SparseMatrix& A01,
    Core::LinAlg::SparseMatrix& A10, Core::LinAlg::SparseMatrix& A11)
{
  if (!A00.RangeMap().SameAs(A01.RangeMap()) || !A00.DomainMap().SameAs(A10.DomainMap()) ||
      !A01.DomainMap().SameAs(A11.DomainMap()) || !A10.RangeMap().SameAs(A11.RangeMap()))
    FOUR_C_THROW("Core::LinAlg::BlockMatrix2x2: block maps are not compatible.");


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
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>> C =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *domainMMex, *rangeMMex));
  // Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> Cb =
  // Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(C);
  // assign matrices
  C->Assign(0, 0, Core::LinAlg::View, A00);
  C->Assign(0, 1, Core::LinAlg::View, A01);
  C->Assign(1, 0, Core::LinAlg::View, A10);
  C->Assign(1, 1, Core::LinAlg::View, A11);

  C->Complete();

  return C;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& Core::LinAlg::operator<<(
    std::ostream& os, const Core::LinAlg::BlockSparseMatrixBase& mat)
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
Core::LinAlg::DefaultBlockMatrixStrategy::DefaultBlockMatrixStrategy(BlockSparseMatrixBase& mat)
    : mat_(mat), scratch_lcols_(mat_.Rows())
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::DefaultBlockMatrixStrategy::Complete(bool enforce_complete)
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::LinAlg::DefaultBlockMatrixStrategy::Complete");

  if (mat_.Filled() and not enforce_complete)
  {
    if (ghost_.size() != 0)
    {
      FOUR_C_THROW("no unresolved ghost entries in a filled block matrix allowed");
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
          Select1st<std::map<int, double>::value_type>());
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
  if (err != 0) FOUR_C_THROW("RemoteIDList failed");

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
        FOUR_C_THROW("gid %d not owned by any domain map", gid);
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
      FOUR_C_THROW("size mismatch panic");
    }

    for (unsigned i = 0; i < requests[proc].size(); ++i)
    {
      int cblock = requests[proc][i];
      int cgid = ghostgids[proc][i];

      if (ghostmap.find(cgid) != ghostmap.end())
        FOUR_C_THROW("column gid %d defined more often that once", cgid);

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
    int rblock = row_block(rgid);
    if (rblock == -1) FOUR_C_THROW("row finding panic");

    for (auto& icol : irow.second)
    {
      int cgid = icol.first;
      if (ghostmap.find(cgid) == ghostmap.end()) FOUR_C_THROW("unknown ghost gid %d", cgid);

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
Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>
Core::LinAlg::CastToBlockSparseMatrixBaseAndCheckSuccess(
    Teuchos::RCP<Core::LinAlg::SparseOperator> input_matrix)
{
  auto block_matrix = Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(input_matrix);
  FOUR_C_ASSERT(block_matrix != Teuchos::null, "Matrix is not a block matrix!");

  return block_matrix;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Core::LinAlg::BlockSparseMatrixBase>
Core::LinAlg::CastToConstBlockSparseMatrixBaseAndCheckSuccess(
    Teuchos::RCP<const Core::LinAlg::SparseOperator> input_matrix)
{
  auto block_matrix =
      Teuchos::rcp_dynamic_cast<const Core::LinAlg::BlockSparseMatrixBase>(input_matrix);
  FOUR_C_ASSERT(block_matrix != Teuchos::null, "Matrix is not a block matrix!");

  return block_matrix;
}

FOUR_C_NAMESPACE_CLOSE
