/*----------------------------------------------------------------------*/
/*! \file

\brief  special assemble/split strategy for blockmatrices arising
        when simulating electrochemical problems with ion transport

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_UTILS_SPLITSTRATEGY_H
#define SCATRA_UTILS_SPLITSTRATEGY_H

#include "../linalg/linalg_blocksparsematrix.H"

namespace SCATRA
{
  /*!
  \brief concentrations-el.potential split strategy

   */
  class SplitStrategy
  {
   public:
    /// construct with a block matrix base
    explicit SplitStrategy(LINALG::BlockSparseMatrixBase& mat)
        : mat_(mat),
          matrix00_(mat_.Matrix(0, 0)),
          matrix01_(mat_.Matrix(0, 1)),
          matrix10_(mat_.Matrix(1, 0)),
          matrix11_(mat_.Matrix(1, 1)),
          numscal_(-1),
          numdofpernode_(-1)
    {
    }

    /// find row block to a given row gid
    int RowBlock(int lrow, int rgid)
    {
      if ((lrow % numdofpernode_) < numscal_) return 0;
      return 1;
    }

    /// find column block to a given column gid
    int ColBlock(int rblock, int lcol, int cgid)
    {
      if ((lcol % numdofpernode_) < numscal_) return 0;
      return 1;
    }

    /// assemble into the given block
    /* For Electrochemistry-applications the upper-left block matrix
     * A00 obeys a sparse block-diagonal substructure. We implement here a special
     * assembly strategy in order to avoid assembling zeros that make up at
     * least 22% of the entries in the element matrix.
     * The well-known VelPresSplitStrategy formerly used for this purpose
     * shows this drawback.
     * Effects: Faster assembly and a more sparse global matrix (and graph)
     */
    void Assemble(int eid, int myrank, const std::vector<int>& lmstride,
        const Epetra_SerialDenseMatrix& Aele, const std::vector<int>& lmrow,
        const std::vector<int>& lmrowowner, const std::vector<int>& lmcol)
    {
      const int lrowdim = (int)lmrow.size();
      const int lcoldim = (int)lmcol.size();

      if (mat_.Filled())
      {
        // We use the maps of the matrix to gain fast access to the LID's.
        // Assembling with SumIntoMyValues based on LID's is two times faster
        // than inserting single values based on the global row and column ids.

        // There is the case of nodes without dofs (XFEM).
        // If no row dofs are present on this proc, their is nothing to assemble.
        // However, the subsequent check for coldofs (in DEBUG mode) would incorrectly fail.
        bool doit = false;
        for (int lrow = 0; lrow < lrowdim; ++lrow)
          if (lmrowowner[lrow] == myrank)
          {
            doit = true;  // This proc owns at least one row of interest
            break;
          }
        if (!doit) return;

        // get the maps
        const Epetra_Map& colmap00 = mat_.Matrix(0, 0).ColMap();
        const Epetra_Map& colmap01 = mat_.Matrix(0, 1).ColMap();
        const Epetra_Map& colmap10 = mat_.Matrix(1, 0).ColMap();
        const Epetra_Map& colmap11 = mat_.Matrix(1, 1).ColMap();
        const Epetra_Map& rowmap00 = mat_.Matrix(0, 0).RowMap();
        const Epetra_Map& rowmap01 = mat_.Matrix(0, 1).RowMap();
        const Epetra_Map& rowmap10 = mat_.Matrix(1, 0).RowMap();
        const Epetra_Map& rowmap11 = mat_.Matrix(1, 1).RowMap();

        // prepare vectors for holding column local ids and the values to be assembled
        const int nnode = lcoldim / numdofpernode_;
        std::vector<double> values0(numscal_ * nnode);
        std::vector<double> values1(nnode);
        std::vector<double> values00(nnode);
        std::map<int, std::vector<int>> localcol00map;
        for (int i = 0; i < numscal_; i++)
        {
          std::vector<int> vec(nnode);
          localcol00map.insert(std::pair<int, std::vector<int>>(i, vec));
        }
        std::vector<int> localcol01(nnode);
        std::vector<int> localcol10(numscal_ * nnode);
        std::vector<int> localcol11(nnode);

        // fill vectors with the column map LID's
        for (int inode = 0; inode < nnode; ++inode)
        {
          const int p = inode * numdofpernode_;
          // the concentrations
          for (int sclarid = 0; sclarid < numscal_; ++sclarid)
          {
            const int lcol = p + sclarid;
            const int cgid = lmcol[lcol];
            (localcol00map[sclarid])[inode] = (colmap00.LID(cgid));
            const int pos = inode * numscal_ + sclarid;
            localcol10[pos] = (colmap10.LID(cgid));
          }
          // the el. potential
          const int lcol = p + numscal_;
          const int cgid = lmcol[lcol];
          localcol01[inode] = (colmap01.LID(cgid));
          localcol11[inode] = (colmap11.LID(cgid));
        }

        int rlid0;
        int rlid1;
        int errone;

        // loop rows of local matrix and assemble each of them
        for (int lrow = 0; lrow < lrowdim; ++lrow)
        {
          // check ownership of row
          if (lmrowowner[lrow] != myrank) continue;

          // get global id and block number for current row
          const int rgid = lmrow[lrow];
          const int scalarid = lrow % numdofpernode_;

          if (scalarid < numscal_)
          {
            // the current row belongs to the transport equation of species 'scalarid'
            rlid0 = rowmap00.LID(rgid);
            rlid1 = rowmap01.LID(rgid);
#ifdef DEBUG
            if (rlid0 < 0) dserror("Sparse matrix A00 does not have global row %d", rgid);
            if (rlid1 < 0) dserror("Sparse matrix A01 does not have global row %d", rgid);
#endif

            // separate the (non-zero!) values of the current row
            for (int j = 0; j < nnode; ++j)
            {
              int lcol = j * numdofpernode_ + scalarid;
              double val = Aele(lrow, lcol);
              values00[j] = val;

              lcol = j * numdofpernode_ + numscal_;
              val = Aele(lrow, lcol);
              values1[j] = val;
            }

            // assemble
            errone = matrix00_.EpetraMatrix()->SumIntoMyValues(
                rlid0, nnode, &values00[0], &(localcol00map[scalarid])[0]);
            if (errone)
              dserror("Epetra_CrsMatrix::SumIntoMyValues returned error code %d for A00", errone);
            errone = matrix01_.EpetraMatrix()->SumIntoMyValues(
                rlid1, nnode, &values1[0], &localcol01[0]);
            if (errone)
              dserror("Epetra_CrsMatrix::SumIntoMyValues returned error code %d for A01", errone);
          }
          else
          {
            // the current row belongs to the equation for the electric potential
            rlid0 = rowmap10.LID(rgid);
            rlid1 = rowmap11.LID(rgid);
#ifdef DEBUG
            if (rlid0 < 0) dserror("Sparse matrix A10 does not have global row %d", rgid);
            if (rlid1 < 0) dserror("Sparse matrix A11 does not have global row %d", rgid);
#endif
            // separate the values of the current row
            int nodespassed = 0;
            for (int lcol = 0; lcol < lcoldim; ++lcol)
            {
              double val = Aele(lrow, lcol);
              const int rest = lcol % numdofpernode_;
              if (rest < numscal_)
              {
                int pos = nodespassed * numscal_ + rest;
                values0[pos] = val;
              }
              else
              {
                values1[nodespassed] = val;
                nodespassed++;
              }
            }

            // assemble
            errone = matrix10_.EpetraMatrix()->SumIntoMyValues(
                rlid0, nnode * numscal_, &values0[0], &localcol10[0]);
            if (errone)
              dserror("Epetra_CrsMatrix::SumIntoMyValues returned error code %d for A10", errone);
            errone = matrix11_.EpetraMatrix()->SumIntoMyValues(
                rlid1, nnode, &values1[0], &localcol11[0]);
            if (errone)
              dserror("Epetra_CrsMatrix::SumIntoMyValues returned error code %d for A11", errone);
          }
        }  // for (int lrow=0; lrow<ldim; ++lrow)
      }
      else
      {
        // the maps of the blockmatrices are not yet available; use global ids instead

        // loop rows of local matrix
        for (int lrow = 0; lrow < lrowdim; ++lrow)
        {
          // check ownership of row
          if (lmrowowner[lrow] != myrank) continue;

          int rgid = lmrow[lrow];
          int rblock = RowBlock(lrow, rgid);
          const int scalarid = lrow % numdofpernode_;

          if (scalarid < numscal_)
          {
            // special treatment for block matrix A00:
            for (int lcol = scalarid; lcol < lcoldim; lcol += numdofpernode_)
            {
              double val = Aele(lrow, lcol);
              int cgid = lmcol[lcol];
              int cblock = ColBlock(rblock, lcol, cgid);

              LINALG::SparseMatrix& matrix = mat_.Matrix(rblock, cblock);
              matrix.Assemble(val, rgid, cgid);
            }
            // values for block matrix A01:
            for (int lcol = numscal_; lcol < lcoldim; lcol += numdofpernode_)
            {
              double val = Aele(lrow, lcol);
              int cgid = lmcol[lcol];
              int cblock = ColBlock(rblock, lcol, cgid);

              LINALG::SparseMatrix& matrix = mat_.Matrix(rblock, cblock);
              matrix.Assemble(val, rgid, cgid);
            }
          }
          else  // rblock == 1, assemble everything
          {
            for (int lcol = 0; lcol < lcoldim; ++lcol)
            {
              double val = Aele(lrow, lcol);
              int cgid = lmcol[lcol];
              int cblock = ColBlock(rblock, lcol, cgid);

              LINALG::SparseMatrix& matrix = mat_.Matrix(rblock, cblock);
              matrix.Assemble(val, rgid, cgid);
            }
          }
        }
      }
    }

    /// assemble into the given block
    void Assemble(double val, int rgid, int cgid)
    {
      int rblock = RowBlock(0, rgid);
      int cblock = ColBlock(rblock, 0, cgid);
      LINALG::SparseMatrix& matrix = mat_.Matrix(rblock, cblock);
      matrix.Assemble(val, rgid, cgid);
    }

    /// assemble the remaining ghost entries
    void Complete() {}

    /// set number of concentration dofs
    void SetNumScal(const int numscal)
    {
      numscal_ = numscal;
      numdofpernode_ = numscal + 1;
    }

   private:
    /// my block matrix base
    LINALG::BlockSparseMatrixBase& mat_;

    // the four sub-matrices of the whole matrix
    LINALG::SparseMatrix& matrix00_;
    LINALG::SparseMatrix& matrix01_;
    LINALG::SparseMatrix& matrix10_;
    LINALG::SparseMatrix& matrix11_;

    /// number of concentration dofs
    int numscal_;

    /// number of dofs per node (= numdim_ +1)
    int numdofpernode_;

  };  // class SplitStrategy

}  // namespace SCATRA

#endif
