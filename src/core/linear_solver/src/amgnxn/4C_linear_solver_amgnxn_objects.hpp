/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_AMGNXN_OBJECTS_HPP
#define FOUR_C_LINEAR_SOLVER_AMGNXN_OBJECTS_HPP

#include "4C_config.hpp"

#include "4C_linalg_blocksparsematrix.hpp"

#include <Epetra_MultiVector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


namespace Core::LinearSolver::AMGNxN
{
  class BlockedVector
  {
   public:
    BlockedVector(int size) : vectors_(size, Teuchos::null) {}

    bool HasOnlyOneBlock() const { return vectors_.size() == 1; }

    int GetNumBlocks() const { return vectors_.size(); }

    Teuchos::RCP<Epetra_MultiVector> GetVector(int i) const { return vectors_[i]; }

    void SetVector(Teuchos::RCP<Epetra_MultiVector> V, int i)
    {
      vectors_[i] = V;
      return;
    }

    BlockedVector GetBlockedVector(const std::vector<int>& blocks) const;

    Teuchos::RCP<BlockedVector> GetBlockedVectorRCP(const std::vector<int>& blocks) const;

    void Update(double a_V, const BlockedVector& V, double a_this);

    void PutScalar(double a);

    BlockedVector DeepCopy() const;

    Teuchos::RCP<BlockedVector> DeepCopyRCP() const;

    Teuchos::RCP<BlockedVector> NewRCP(bool ZeroIt = false) const;

   private:
    std::vector<Teuchos::RCP<Epetra_MultiVector>> vectors_;
  };

  class BlockedMatrix
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~BlockedMatrix() = default;

    BlockedMatrix() = default;
    BlockedMatrix(int rows, int cols)
        : matrices_(rows * cols, Teuchos::null), rows_(rows), cols_(cols)
    {
    }

    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> GetMatrix(int i, int j) const
    {
      return matrices_[i * GetNumCols() + j];
    }  // Row major order

    virtual void SetMatrix(Teuchos::RCP<Core::LinAlg::SparseMatrix> A, int i, int j)
    {
      matrices_[i * GetNumCols() + j] = A;
      return;
    }

    Teuchos::RCP<BlockedMatrix> GetBlockedMatrixRCP(
        const std::vector<int>& row_blocks, const std::vector<int>& col_blocks) const
    {
      return Teuchos::rcp(new BlockedMatrix(GetBlockedMatrix(row_blocks, col_blocks)));
    }

    virtual BlockedMatrix GetBlockedMatrix(
        const std::vector<int>& row_blocks, const std::vector<int>& col_blocks) const;


    virtual void Apply(const BlockedVector& in, BlockedVector& out) const;

    bool HasOnlyOneBlock() const { return GetNumRows() * GetNumCols() == 1; }

    int GetNumRows() const { return rows_; }

    int GetNumCols() const { return cols_; }

    virtual Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> get_block_sparse_matrix(
        Core::LinAlg::DataAccess access);

    void ParseBlocks(const std::string& block_string, const std::vector<int>& blocks,
        std::vector<std::vector<int>>& superblocks_to_blocks,
        std::vector<std::vector<int>>& superblocks_to_blocks_local);

    virtual Teuchos::RCP<BlockedVector> new_domain_blocked_vector(
        int NV, bool ZeroIt = false) const;

    virtual Teuchos::RCP<BlockedVector> new_range_blocked_vector(int NV, bool ZeroIt = false) const;

   protected:
    std::vector<Teuchos::RCP<Core::LinAlg::SparseMatrix>> matrices_;  // Row major order
    int rows_;
    int cols_;
  };

  class DiagonalBlockedMatrix : public BlockedMatrix
  {
   public:
    DiagonalBlockedMatrix(int rows)
    {
      matrices_.assign(rows, Teuchos::null);
      rows_ = rows;
      cols_ = rows;
    }

    void SetMatrix(Teuchos::RCP<Core::LinAlg::SparseMatrix> A, int i, int j) override
    {
      if (i != j) FOUR_C_THROW("You can only set diagonal blocks");
      matrices_[i] = A;
    }

    Teuchos::RCP<Core::LinAlg::SparseMatrix> GetMatrix(int i, int j) const override
    {
      if (i != j) FOUR_C_THROW("You can only get diagonal blocks");
      return matrices_[i];
    }

    BlockedMatrix GetBlockedMatrix(
        const std::vector<int>& row_blocks, const std::vector<int>& col_blocks) const override
    {
      FOUR_C_THROW("Function not implemented yet.");
      return *this;
    }

    virtual Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> get_block_sparse_matrix()
    {
      FOUR_C_THROW("Function not implemented yet.");
      return Teuchos::null;
    }

    void Apply(const BlockedVector& in, BlockedVector& out) const override;

    Teuchos::RCP<BlockedVector> new_domain_blocked_vector(
        int NV, bool ZeroIt = false) const override;

    Teuchos::RCP<BlockedVector> new_range_blocked_vector(
        int NV, bool ZeroIt = false) const override;
  };

  // class BlockAggrupator
  //{
  //
  //  public:
  //    std::vector<int> SuperBlock2Blocks
  //
  //};

  // class GenericSmoother
  //{
  //  public:
  //    virtual void Solve(const BlockedVector& X, BlockedVector& Y, bool InitialGuessIsZero) =
  //    0;
  //};

}  // namespace Core::LinearSolver::AMGNxN

FOUR_C_NAMESPACE_CLOSE

#endif
