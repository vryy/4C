/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble the fbi coupling contributions


\level 1

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_FBI_FLUIDBLOCKMATRIX_ASSEMBLY_STRATEGY_HPP
#define FOUR_C_FBI_FLUIDBLOCKMATRIX_ASSEMBLY_STRATEGY_HPP

#include "4C_config.hpp"

#include "4C_fbi_fluid_assembly_strategy.hpp"
#include "4C_utils_exceptions.hpp"

#include <Teuchos_RCP.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN


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
       * \brief Calls the correct assembly method for the used global fluid matrix depending on the
       * fluid problem
       *
       * \param[in, out] cff fluid coupling matrix
       * \param[in] elegid element gid
       * \param[in] elemat dense matrix to be assembled
       * \param[in] lmrow vector with row gids
       * \param[in] lmrowowner vector with owner procs of row gids
       * \param[in] lmcol vector with column gids
       */
      void assemble_fluid_matrix(Teuchos::RCP<Core::LinAlg::SparseOperator> cff, int elegid,
          const std::vector<int>& lmstride, const Core::LinAlg::SerialDenseMatrix& elemat,
          const std::vector<int>& lmrow, const std::vector<int>& lmrowowner,
          const std::vector<int>& lmcol) override;
    };
  }  // namespace UTILS
}  // namespace FBI

FOUR_C_NAMESPACE_CLOSE

#endif
