/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_AMGNXN_OBJECTS_HPP
#define FOUR_C_LINEAR_SOLVER_AMGNXN_OBJECTS_HPP

// Trilinos includes
#include "baci_config.hpp"

#include <Epetra_MultiVector.h>
#include <Teuchos_RCP.hpp>

// Baci includes
#include "baci_linalg_blocksparsematrix.hpp"

BACI_NAMESPACE_OPEN


namespace CORE::LINEAR_SOLVER::AMGNXN
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

    virtual Teuchos::RCP<CORE::LINALG::SparseMatrix> GetMatrix(int i, int j) const
    {
      return matrices_[i * GetNumCols() + j];
    }  // Row major order

    virtual void SetMatrix(Teuchos::RCP<CORE::LINALG::SparseMatrix> A, int i, int j)
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

    virtual Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> GetBlockSparseMatrix(
        CORE::LINALG::DataAccess access);

    void ParseBlocks(const std::string& block_string, const std::vector<int>& blocks,
        std::vector<std::vector<int>>& superblocks_to_blocks,
        std::vector<std::vector<int>>& superblocks_to_blocks_local);

    virtual Teuchos::RCP<BlockedVector> NewDomainBlockedVector(int NV, bool ZeroIt = false) const;

    virtual Teuchos::RCP<BlockedVector> NewRangeBlockedVector(int NV, bool ZeroIt = false) const;

   protected:
    std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> matrices_;  // Row major order
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

    void SetMatrix(Teuchos::RCP<CORE::LINALG::SparseMatrix> A, int i, int j) override
    {
      if (i != j) dserror("You can only set diagonal blocks");
      matrices_[i] = A;
    }

    Teuchos::RCP<CORE::LINALG::SparseMatrix> GetMatrix(int i, int j) const override
    {
      if (i != j) dserror("You can only get diagonal blocks");
      return matrices_[i];
    }

    BlockedMatrix GetBlockedMatrix(
        const std::vector<int>& row_blocks, const std::vector<int>& col_blocks) const override
    {
      dserror("Function not implemented yet.");
      return *this;
    }

    virtual Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> GetBlockSparseMatrix()
    {
      dserror("Function not implemented yet.");
      return Teuchos::null;
    }

    void Apply(const BlockedVector& in, BlockedVector& out) const override;

    Teuchos::RCP<BlockedVector> NewDomainBlockedVector(int NV, bool ZeroIt = false) const override;

    Teuchos::RCP<BlockedVector> NewRangeBlockedVector(int NV, bool ZeroIt = false) const override;
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

}  // namespace CORE::LINEAR_SOLVER::AMGNXN

BACI_NAMESPACE_CLOSE

#endif  // SOLVER_AMGNXN_OBJECTS_H
