/*--------------------------------------------------------------------------*/
/*! \file

\brief Utility routines for ale mesh tying


\level 2
*/
/*--------------------------------------------------------------------------*/
#ifndef FOUR_C_ALE_UTILS_HPP
#define FOUR_C_ALE_UTILS_HPP


#include "4C_config.hpp"

#include "4C_lib_discret.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Discret
{
  class Discretization;
}

namespace ALE
{
  namespace UTILS
  {
    /// (FSI) interface block matrix split strategy
    class InterfaceSplitStrategy : public Core::LinAlg::DefaultBlockMatrixStrategy
    {
     public:
      explicit InterfaceSplitStrategy(Core::LinAlg::BlockSparseMatrixBase& mat)
          : Core::LinAlg::DefaultBlockMatrixStrategy(mat)
      {
      }

      /// assemble into the given block
      void Assemble(int eid, int myrank, const std::vector<int>& lmstride,
          const Core::LinAlg::SerialDenseMatrix& Aele, const std::vector<int>& lmrow,
          const std::vector<int>& lmrowowner, const std::vector<int>& lmcol)
      {
        if (condelements_->find(eid) != condelements_->end())
        {
          // if we have an element with conditioned nodes, we have to do the
          // default assembling
          Core::LinAlg::DefaultBlockMatrixStrategy::Assemble(
              eid, myrank, lmstride, Aele, lmrow, lmrowowner, lmcol);
        }
        else
        {
          // if there are no conditioned nodes we can simply assemble to the
          // internal matrix
          Core::LinAlg::SparseMatrix& matrix = mat().Matrix(0, 0);
          matrix.Assemble(eid, lmstride, Aele, lmrow, lmrowowner, lmcol);
        }
      }

      void Assemble(double val, int rgid, int cgid)
      {
        // forward single value assembling
        Core::LinAlg::DefaultBlockMatrixStrategy::Assemble(val, rgid, cgid);
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


FOUR_C_NAMESPACE_CLOSE

#endif
