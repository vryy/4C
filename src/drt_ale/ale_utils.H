/*--------------------------------------------------------------------------*/
/*! \file

\brief Utility routines for ale mesh tying

\maintainer Matthias Mayr

\level 2
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_ALE_ALE_UTILS_H_
#define SRC_DRT_ALE_ALE_UTILS_H_


#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"


namespace DRT
{
  class Discretization;
}

namespace ALE
{
  namespace UTILS
  {
    /// (FSI) interface block matrix split strategy
    class InterfaceSplitStrategy : public LINALG::DefaultBlockMatrixStrategy
    {
     public:
      explicit InterfaceSplitStrategy(LINALG::BlockSparseMatrixBase& mat)
          : LINALG::DefaultBlockMatrixStrategy(mat)
      {
      }

      /// assemble into the given block
      void Assemble(int eid, int myrank, const std::vector<int>& lmstride,
          const Epetra_SerialDenseMatrix& Aele, const std::vector<int>& lmrow,
          const std::vector<int>& lmrowowner, const std::vector<int>& lmcol)
      {
        if (condelements_->find(eid) != condelements_->end())
        {
          // if we have an element with conditioned nodes, we have to do the
          // default assembling
          LINALG::DefaultBlockMatrixStrategy::Assemble(
              eid, myrank, lmstride, Aele, lmrow, lmrowowner, lmcol);
        }
        else
        {
          // if there are no conditioned nodes we can simply assemble to the
          // internal matrix
          LINALG::SparseMatrix& matrix = Mat().Matrix(0, 0);
          matrix.Assemble(eid, lmstride, Aele, lmrow, lmrowowner, lmcol);
        }
      }

      void Assemble(double val, int rgid, int cgid)
      {
        // forward single value assembling
        LINALG::DefaultBlockMatrixStrategy::Assemble(val, rgid, cgid);
      }

      void SetCondElements(Teuchos::RCP<std::set<int>> condelements)
      {
        condelements_ = condelements;
      }

     private:
      Teuchos::RCP<std::set<int>> condelements_;
    };
  }  // namespace UTILS
}  // namespace ALE


#endif /* SRC_DRT_ALE_ALE_UTILS_H_ */
