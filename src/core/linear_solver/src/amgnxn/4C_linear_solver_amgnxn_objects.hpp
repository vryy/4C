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

    bool has_only_one_block() const { return vectors_.size() == 1; }

    int get_num_blocks() const { return vectors_.size(); }

    Teuchos::RCP<Epetra_MultiVector> get_vector(int i) const { return vectors_[i]; }

    void set_vector(Teuchos::RCP<Epetra_MultiVector> V, int i)
    {
      vectors_[i] = V;
      return;
    }

    BlockedVector get_blocked_vector(const std::vector<int>& blocks) const;

    Teuchos::RCP<BlockedVector> get_blocked_vector_rcp(const std::vector<int>& blocks) const;

    void update(double a_V, const BlockedVector& V, double a_this);

    void put_scalar(double a);

    BlockedVector deep_copy() const;

    Teuchos::RCP<BlockedVector> deep_copy_rcp() const;

    Teuchos::RCP<BlockedVector> new_rcp(bool ZeroIt = false) const;

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

    virtual Teuchos::RCP<Core::LinAlg::SparseMatrix> get_matrix(int i, int j) const
    {
      return matrices_[i * get_num_cols() + j];
    }  // Row major order

    virtual void set_matrix(Teuchos::RCP<Core::LinAlg::SparseMatrix> A, int i, int j)
    {
      matrices_[i * get_num_cols() + j] = A;
      return;
    }

    Teuchos::RCP<BlockedMatrix> get_blocked_matrix_rcp(
        const std::vector<int>& row_blocks, const std::vector<int>& col_blocks) const
    {
      return Teuchos::rcp(new BlockedMatrix(get_blocked_matrix(row_blocks, col_blocks)));
    }

    virtual BlockedMatrix get_blocked_matrix(
        const std::vector<int>& row_blocks, const std::vector<int>& col_blocks) const;


    virtual void apply(const BlockedVector& in, BlockedVector& out) const;

    bool has_only_one_block() const { return get_num_rows() * get_num_cols() == 1; }

    int get_num_rows() const { return rows_; }

    int get_num_cols() const { return cols_; }

    virtual Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> get_block_sparse_matrix(
        Core::LinAlg::DataAccess access);

    void parse_blocks(const std::string& block_string, const std::vector<int>& blocks,
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

    void set_matrix(Teuchos::RCP<Core::LinAlg::SparseMatrix> A, int i, int j) override
    {
      if (i != j) FOUR_C_THROW("You can only set diagonal blocks");
      matrices_[i] = A;
    }

    Teuchos::RCP<Core::LinAlg::SparseMatrix> get_matrix(int i, int j) const override
    {
      if (i != j) FOUR_C_THROW("You can only get diagonal blocks");
      return matrices_[i];
    }

    BlockedMatrix get_blocked_matrix(
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

    void apply(const BlockedVector& in, BlockedVector& out) const override;

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
