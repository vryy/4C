/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble the fbi coupling contributions


\level 1

*/
/*-----------------------------------------------------------*/


#ifndef FLUID_ASSEMBLYBLOCK_STRATEGY_H_
#define FLUID_ASSEMBLYBLOCK_STRATEGY_H_

#include <vector>

#include "dserror.H"
#include "Teuchos_RCP.hpp"
#include <Epetra_FEVector.h>
#include "fluid_assembly_strategy.H"


namespace FBI
{
  namespace UTILS
  {
    /**
     * \brief This class assembles the contributions of fluid beam mesh tying pairs into the global
     * matrices in the case that a SparseBlockMatrix is used in the fluid problem.
     *
     * The form of the fluid matrix and in an extension the required assembly method
     * depend on the fluid problem, particularly if mesh tying is used.
     */
    class FBIBlockAssemblyStrategy : public FBIAssemblyStrategy
    {
     public:
      /**
       * \brief Destructor.
       */
      ~FBIBlockAssemblyStrategy() override = default;

      /**
       * \brief Calls the correct assembly method for the used global fluid matrix depending on the
       * fluid problem
       *
       * \param[in, out] cff fluid coupling matrix
       * \param[in] eid element gid
       * \param[in] Aele dense matrix to be assembled
       * \param[in] lmrow vector with row gids
       * \param[in] lmrowowner vector with owner procs of row gids
       * \param[in] lmcol vector with column gids
       */
      void AssembleFluidMatrix(Teuchos::RCP<LINALG::SparseOperator> cff, int eid,
          const std::vector<int>& lmstride, const Epetra_SerialDenseMatrix& Aele,
          const std::vector<int>& lmrow, const std::vector<int>& lmrowowner,
          const std::vector<int>& lmcol) override;
    };
  }  // namespace UTILS
}  // namespace FBI

#endif
