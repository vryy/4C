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
  blocks_.reserve(rows() * cols());

  // add sparse matrices in row,column order
  for (int r = 0; r < rows(); ++r)
  {
    for (int c = 0; c < cols(); ++c)
    {
      blocks_.emplace_back(range_map(r), npr, explicitdirichlet, savegraph);
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Core::LinAlg::BlockSparseMatrixBase::destroy(bool throw_exception_for_blocks)
{
  /// destroy matrix blocks
  for (auto& block : blocks_)
  {
    block.destroy(throw_exception_for_blocks);
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
Teuchos::RCP<Core::LinAlg::SparseMatrix> Core::LinAlg::BlockSparseMatrixBase::merge(
    bool explicitdirichlet) const
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::LinAlg::BlockSparseMatrixBase::Merge");

  const SparseMatrix& m00 = matrix(0, 0);

  Teuchos::RCP<SparseMatrix> sparse =
      Teuchos::rcp(new SparseMatrix(*fullrowmap_, m00.max_num_entries(), explicitdirichlet));
  for (const auto& block : blocks_)
  {
    sparse->add(block, false, 1.0, 1.0);
  }
  if (filled())
  {
    sparse->complete(full_domain_map(), full_range_map());
  }
  return sparse;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::assign(
    int r, int c, DataAccess access, const SparseMatrix& mat)
{
#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (not matrix(r, c).row_map().SameAs(mat.row_map()))
    FOUR_C_THROW("cannot assign nonmatching matrices");
#endif
  matrix(r, c).assign(access, mat);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::zero()
{
  for (auto& block : blocks_) block.zero();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::reset()
{
  for (int i = 0; i < rows(); ++i)
  {
    for (int j = 0; j < cols(); ++j)
    {
      matrix(i, j).reset();
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::complete(bool enforce_complete)
{
  for (int r = 0; r < rows(); ++r)
  {
    for (int c = 0; c < cols(); ++c)
    {
      matrix(r, c).complete(domain_map(c), range_map(r), enforce_complete);
    }
  }

  fullrowmap_ = Teuchos::rcp(new Epetra_Map(*(rangemaps_.full_map())));

  if (fullcolmap_ == Teuchos::null)
  {
    // build full col map
    std::vector<int> colmapentries;
    for (int c = 0; c < cols(); ++c)
    {
      for (int r = 0; r < rows(); ++r)
      {
        const Epetra_Map& colmap = matrix(r, c).col_map();
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
void Core::LinAlg::BlockSparseMatrixBase::complete(
    const Epetra_Map& domainmap, const Epetra_Map& rangemap, bool enforce_complete)
{
  FOUR_C_THROW("Complete with arguments not supported for block matrices");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Core::LinAlg::BlockSparseMatrixBase::filled() const
{
  for (const auto& block : blocks_)
    if (not block.filled()) return false;
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::un_complete()
{
  for (auto& block : blocks_) block.un_complete();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::apply_dirichlet(
    const Epetra_Vector& dbctoggle, bool diagonalblock)
{
  for (int rblock = 0; rblock < rows(); ++rblock)
  {
    Teuchos::RCP<Epetra_Vector> rowtoggle = rangemaps_.extract_vector(dbctoggle, rblock);
    for (int cblock = 0; cblock < cols(); ++cblock)
    {
      Core::LinAlg::SparseMatrix& bmat = matrix(rblock, cblock);
      bmat.apply_dirichlet(*rowtoggle, diagonalblock and rblock == cblock);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::apply_dirichlet(
    const Epetra_Map& dbcmap, bool diagonalblock)
{
  for (int rblock = 0; rblock < rows(); ++rblock)
  {
    for (int cblock = 0; cblock < cols(); ++cblock)
    {
      Core::LinAlg::SparseMatrix& bmat = matrix(rblock, cblock);
      bmat.apply_dirichlet(dbcmap, diagonalblock and rblock == cblock);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Core::LinAlg::BlockSparseMatrixBase::is_dbc_applied(
    const Epetra_Map& dbcmap, bool diagonalblock, const Core::LinAlg::SparseMatrix* trafo) const
{
  for (int rblock = 0; rblock < rows(); ++rblock)
  {
    for (int cblock = 0; cblock < cols(); ++cblock)
    {
      if (not matrix(rblock, cblock)
                  .is_dbc_applied(dbcmap, diagonalblock and (rblock == cblock), trafo))
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
  Y.PutScalar(0.0);

  if (not UseTranspose())
  {
    for (int rblock = 0; rblock < rows(); ++rblock)
    {
      Teuchos::RCP<Epetra_MultiVector> rowresult = rangemaps_.vector(rblock, Y.NumVectors());
      Teuchos::RCP<Epetra_MultiVector> rowy = rangemaps_.vector(rblock, Y.NumVectors());
      for (int cblock = 0; cblock < cols(); ++cblock)
      {
        Teuchos::RCP<Epetra_MultiVector> colx = domainmaps_.extract_vector(X, cblock);
        const Core::LinAlg::SparseMatrix& bmat = matrix(rblock, cblock);
        int err = bmat.Apply(*colx, *rowy);
        if (err != 0)
          FOUR_C_THROW(
              "failed to apply vector to matrix block (%d,%d): err=%d", rblock, cblock, err);
        rowresult->Update(1.0, *rowy, 1.0);
      }
      rangemaps_.insert_vector(*rowresult, rblock, Y);
    }
  }
  else
  {
    for (int rblock = 0; rblock < cols(); ++rblock)
    {
      Teuchos::RCP<Epetra_MultiVector> rowresult = rangemaps_.vector(rblock, Y.NumVectors());
      Teuchos::RCP<Epetra_MultiVector> rowy = rangemaps_.vector(rblock, Y.NumVectors());
      for (int cblock = 0; cblock < rows(); ++cblock)
      {
        Teuchos::RCP<Epetra_MultiVector> colx = domainmaps_.extract_vector(X, cblock);
        const Core::LinAlg::SparseMatrix& bmat = matrix(cblock, rblock);
        int err = bmat.Apply(*colx, *rowy);
        if (err != 0) FOUR_C_THROW("failed to apply vector to matrix: err=%d", err);
        rowresult->Update(1.0, *rowy, 1.0);
      }
      rangemaps_.insert_vector(*rowresult, rblock, Y);
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
void Core::LinAlg::BlockSparseMatrixBase::add(const Core::LinAlg::SparseOperator& A,
    const bool transposeA, const double scalarA, const double scalarB)
{
  A.add_other(*this, transposeA, scalarA, scalarB);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::add_other(Core::LinAlg::BlockSparseMatrixBase& A,
    const bool transposeA, const double scalarA, const double scalarB) const
{
  A.add(*this, transposeA, scalarA, scalarB);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::add_other(Core::LinAlg::SparseMatrixBase& A,
    const bool transposeA, const double scalarA, const double scalarB) const
{
  FOUR_C_THROW("BlockSparseMatrix and SparseMatrix cannot be added");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::BlockSparseMatrixBase::add(const Core::LinAlg::BlockSparseMatrixBase& A,
    const bool transposeA, const double scalarA, const double scalarB)
{
  for (int i = 0; i < rows(); i++)
  {
    for (int j = 0; j < cols(); j++)
    {
      if (transposeA)
        matrix(i, j).add(A.matrix(j, i), transposeA, scalarA, scalarB);
      else
        matrix(i, j).add(A.matrix(i, j), transposeA, scalarA, scalarB);
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::LinAlg::BlockSparseMatrixBase::scale(double ScalarConstant)
{
  for (int i = 0; i < rows(); i++)
  {
    for (int j = 0; j < cols(); j++)
    {
      int err = matrix(i, j).scale(ScalarConstant);
      if (err != 0) FOUR_C_THROW("Scaling of matrix block (%d,%d) failed", i, j);
    }
  }
  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::LinAlg::BlockSparseMatrixBase::multiply(
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
  return full_domain_map().Comm();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map& Core::LinAlg::BlockSparseMatrixBase::OperatorDomainMap() const
{
  return full_domain_map();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map& Core::LinAlg::BlockSparseMatrixBase::OperatorRangeMap() const
{
  return full_range_map();
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
  if (!A.filled() || !B.filled())
    FOUR_C_THROW(
        "Core::LinAlg::BlockSparseMatrixBase::MatrixMultiyply: we expect A and B to be filled");

  if (A.cols() != B.rows() /*|| !A.FullDomainMap().SameAs(B.FullRowMap())*/)
  {
    FOUR_C_THROW("Core::LinAlg::BlockSparseMatrixBase::MatrixMultiply: A and B not compatible");
  }


  int npr = 81;  // estimated number of entries per row in each block

  // generate result matrix
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>> C =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          B.domain_extractor(), A.range_extractor(), npr, explicitdirichlet, savegraph));

  // nested loop over all blocks
  for (int i = 0; i < C->rows(); i++)
  {
    for (int j = 0; j < C->cols(); j++)
    {
      // build block C(i,j)
      Teuchos::RCP<SparseMatrix> Cij = Teuchos::rcp(
          new Core::LinAlg::SparseMatrix(C->range_map(i), npr, explicitdirichlet, savegraph));

      for (int l = 0; l < C->cols(); l++)
      {
        // build submatrices for row i and j
        Teuchos::RCP<SparseMatrix> tmpij =
            Core::LinAlg::Multiply(A.matrix(i, l), false, B.matrix(l, j), false, true);
        Cij->add(*tmpij, false, 1.0, 1.0);
      }

      // complete Cij with correct range and domain map
      if (completeoutput) Cij->complete(C->domain_map(j), C->range_map(i));

      // assign Cij block
      C->assign(i, j, Core::LinAlg::View, *Cij);
    }
  }

  if (completeoutput) C->complete();

  return Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(C);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>>
Core::LinAlg::BlockMatrix2x2(Core::LinAlg::SparseMatrix& A00, Core::LinAlg::SparseMatrix& A01,
    Core::LinAlg::SparseMatrix& A10, Core::LinAlg::SparseMatrix& A11)
{
  if (!A00.range_map().SameAs(A01.range_map()) || !A00.domain_map().SameAs(A10.domain_map()) ||
      !A01.domain_map().SameAs(A11.domain_map()) || !A10.range_map().SameAs(A11.range_map()))
    FOUR_C_THROW("Core::LinAlg::BlockMatrix2x2: block maps are not compatible.");


  // generate range map
  std::vector<Teuchos::RCP<const Epetra_Map>> range_maps;
  range_maps.reserve(2);

  range_maps.emplace_back(Teuchos::rcp(new Epetra_Map(A00.range_map())));
  range_maps.emplace_back(Teuchos::rcp(new Epetra_Map(A10.range_map())));
  Teuchos::RCP<const Epetra_Map> range_map = MultiMapExtractor::merge_maps(range_maps);
  Teuchos::RCP<MultiMapExtractor> rangeMMex =
      Teuchos::rcp(new MultiMapExtractor(*range_map, range_maps));

  // generate domain map
  std::vector<Teuchos::RCP<const Epetra_Map>> domain_maps;
  domain_maps.reserve(2);

  domain_maps.emplace_back(Teuchos::rcp(new Epetra_Map(A00.domain_map())));
  domain_maps.emplace_back(Teuchos::rcp(new Epetra_Map(A01.domain_map())));
  Teuchos::RCP<const Epetra_Map> domain_map = MultiMapExtractor::merge_maps(domain_maps);
  Teuchos::RCP<MultiMapExtractor> domainMMex =
      Teuchos::rcp(new MultiMapExtractor(*domain_map, domain_maps));

  // generate result matrix
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>> C =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *domainMMex, *rangeMMex));
  // Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> Cb =
  // Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(C);
  // assign matrices
  C->assign(0, 0, Core::LinAlg::View, A00);
  C->assign(0, 1, Core::LinAlg::View, A01);
  C->assign(1, 0, Core::LinAlg::View, A10);
  C->assign(1, 1, Core::LinAlg::View, A11);

  C->complete();

  return C;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& Core::LinAlg::operator<<(
    std::ostream& os, const Core::LinAlg::BlockSparseMatrixBase& mat)
{
  for (int i = 0; i < mat.rows(); ++i)
  {
    for (int j = 0; j < mat.cols(); ++j)
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
    : mat_(mat), scratch_lcols_(mat_.rows())
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::DefaultBlockMatrixStrategy::complete(bool enforce_complete)
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::LinAlg::DefaultBlockMatrixStrategy::Complete");

  if (mat_.filled() and not enforce_complete)
  {
    if (ghost_.size() != 0)
    {
      FOUR_C_THROW("no unresolved ghost entries in a filled block matrix allowed");
    }
    return;
  }

  // finish ghost entries

  int rows = mat_.rows();
  int cols = mat_.cols();

  std::set<int> cgids;

  // get the list of all ghost entries gids
  for (int rblock = 0; rblock < rows; ++rblock)
  {
    const Epetra_Map& rowmap = mat_.range_map(rblock);

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

  int err = mat_.full_domain_map().RemoteIDList(
      cgidlist.size(), cgidlist.data(), cpidlist.data(), nullptr);
  if (err != 0) FOUR_C_THROW("RemoteIDList failed");

  const Epetra_Comm& comm = mat_.full_range_map().Comm();
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
        const Epetra_Map& domainmap = mat_.domain_map(cblock);
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

      SparseMatrix& matrix = mat_.matrix(rblock, cblock);
      matrix.assemble(val, rgid, cgid);
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
