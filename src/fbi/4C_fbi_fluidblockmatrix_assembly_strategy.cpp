/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble the fbi coupling contributions


\level 1

*/
/*-----------------------------------------------------------*/

#include "4C_fbi_fluidblockmatrix_assembly_strategy.hpp"

#include "4C_beam3_base.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_serialdensematrix.hpp"
#include "4C_linalg_serialdensevector.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/

void FBI::UTILS::FBIBlockAssemblyStrategy::assemble_fluid_matrix(
    Teuchos::RCP<Core::LinAlg::SparseOperator> cff, int elegid, const std::vector<int>& lmstride,
    const Core::LinAlg::SerialDenseMatrix& elemat, const std::vector<int>& lmrow,
    const std::vector<int>& lmrowowner, const std::vector<int>& lmcol)
{
  cff->assemble(elegid, lmstride, elemat, lmrow, lmrowowner, lmcol);
}

FOUR_C_NAMESPACE_CLOSE
