/*!----------------------------------------------------------------------
\file solver_amgnxn_objects.H

<pre>
\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15262
Created on: Feb 27, 2014
</pre>
*----------------------------------------------------------------------*/
#ifndef AMGNXN_OBJECTS_H
#define AMGNXN_OBJECTS_H
#ifdef HAVE_MueLu

// Trilinos includes
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

// Baci includes
#include "../linalg/linalg_blocksparsematrix.H"


namespace LINALG
{
  namespace SOLVER
  {
    namespace AMGNXN
    {
      class BlockedVector
      {
       public:
        BlockedVector(int size) : vectors_(size, Teuchos::null), size_(size) {}

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
        int size_;
      };


      class BlockedMatrix
      {
       public:
        BlockedMatrix() {}
        BlockedMatrix(int rows, int cols)
            : matrices_(rows * cols, Teuchos::null), rows_(rows), cols_(cols)
        {
        }

        virtual Teuchos::RCP<SparseMatrix> GetMatrix(int i, int j) const
        {
          return matrices_[i * GetNumCols() + j];
        }  // Row major order

        virtual void SetMatrix(Teuchos::RCP<SparseMatrix> A, int i, int j)
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

        virtual Teuchos::RCP<BlockSparseMatrixBase> GetBlockSparseMatrix(DataAccess access);

        void ParseBlocks(const std::string& block_string, const std::vector<int>& blocks,
            std::vector<std::vector<int>>& superblocks_to_blocks,
            std::vector<std::vector<int>>& superblocks_to_blocks_local);

        virtual Teuchos::RCP<BlockedVector> NewDomainBlockedVector(
            int NV, bool ZeroIt = false) const;

        virtual Teuchos::RCP<BlockedVector> NewRangeBlockedVector(
            int NV, bool ZeroIt = false) const;

       protected:
        std::vector<Teuchos::RCP<SparseMatrix>> matrices_;  // Row major order
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

        virtual void SetMatrix(Teuchos::RCP<SparseMatrix> A, int i, int j)
        {
          if (i != j) dserror("You can only set diagonal blocks");
          matrices_[i] = A;
        }

        virtual Teuchos::RCP<SparseMatrix> GetMatrix(int i, int j) const
        {
          if (i != j) dserror("You can only get diagonal blocks");
          return matrices_[i];
        }

        virtual BlockedMatrix GetBlockedMatrix(
            const std::vector<int>& row_blocks, const std::vector<int>& col_blocks) const
        {
          dserror("Function not implemented yet.");
          return *this;
        }

        virtual Teuchos::RCP<BlockSparseMatrixBase> GetBlockSparseMatrix()
        {
          dserror("Function not implemented yet.");
          return Teuchos::null;
        }

        virtual void Apply(const BlockedVector& in, BlockedVector& out) const;

        virtual Teuchos::RCP<BlockedVector> NewDomainBlockedVector(
            int NV, bool ZeroIt = false) const;

        virtual Teuchos::RCP<BlockedVector> NewRangeBlockedVector(
            int NV, bool ZeroIt = false) const;
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


    }  // namespace AMGNXN
  }    // namespace SOLVER
}  // namespace LINALG

#endif  // HAVE_MueLu
#endif  // AMGNXN_OBJECTS_H
