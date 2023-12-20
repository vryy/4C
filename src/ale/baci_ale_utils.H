/*--------------------------------------------------------------------------*/
/*! \file

\brief Utility routines for ale mesh tying


\level 2
*/
/*--------------------------------------------------------------------------*/
#ifndef BACI_ALE_UTILS_H
#define BACI_ALE_UTILS_H


#include "baci_config.H"

#include "baci_lib_discret.H"
#include "baci_linalg_utils_sparse_algebra_math.H"

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}

namespace ALE
{
  namespace UTILS
  {
    /// (FSI) interface block matrix split strategy
    class InterfaceSplitStrategy : public CORE::LINALG::DefaultBlockMatrixStrategy
    {
     public:
      explicit InterfaceSplitStrategy(CORE::LINALG::BlockSparseMatrixBase& mat)
          : CORE::LINALG::DefaultBlockMatrixStrategy(mat)
      {
      }

      /// assemble into the given block
      void Assemble(int eid, int myrank, const std::vector<int>& lmstride,
          const CORE::LINALG::SerialDenseMatrix& Aele, const std::vector<int>& lmrow,
          const std::vector<int>& lmrowowner, const std::vector<int>& lmcol)
      {
        if (condelements_->find(eid) != condelements_->end())
        {
          // if we have an element with conditioned nodes, we have to do the
          // default assembling
          CORE::LINALG::DefaultBlockMatrixStrategy::Assemble(
              eid, myrank, lmstride, Aele, lmrow, lmrowowner, lmcol);
        }
        else
        {
          // if there are no conditioned nodes we can simply assemble to the
          // internal matrix
          CORE::LINALG::SparseMatrix& matrix = Mat().Matrix(0, 0);
          matrix.Assemble(eid, lmstride, Aele, lmrow, lmrowowner, lmcol);
        }
      }

      void Assemble(double val, int rgid, int cgid)
      {
        // forward single value assembling
        CORE::LINALG::DefaultBlockMatrixStrategy::Assemble(val, rgid, cgid);
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


BACI_NAMESPACE_CLOSE

#endif  // ALE_UTILS_H
