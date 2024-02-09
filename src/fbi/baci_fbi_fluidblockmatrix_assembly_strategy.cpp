/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble the fbi coupling contributions


\level 1

*/
/*-----------------------------------------------------------*/

#include "baci_fbi_fluidblockmatrix_assembly_strategy.hpp"

#include "baci_beam3_base.hpp"
#include "baci_lib_discret.hpp"
#include "baci_lib_element.hpp"
#include "baci_linalg_blocksparsematrix.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"
#include "baci_linalg_sparsematrix.hpp"
#include "baci_utils_exceptions.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/

void FBI::UTILS::FBIBlockAssemblyStrategy::AssembleFluidMatrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> cff, int elegid, const std::vector<int>& lmstride,
    const CORE::LINALG::SerialDenseMatrix& elemat, const std::vector<int>& lmrow,
    const std::vector<int>& lmrowowner, const std::vector<int>& lmcol)
{
  cff->Assemble(elegid, lmstride, elemat, lmrow, lmrowowner, lmcol);
}

BACI_NAMESPACE_CLOSE
