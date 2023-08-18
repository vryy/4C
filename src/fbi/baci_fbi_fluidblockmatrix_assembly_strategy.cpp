/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble the fbi coupling contributions


\level 1

*/
/*-----------------------------------------------------------*/

#include "baci_fbi_fluidblockmatrix_assembly_strategy.H"

#include "baci_beam3_base.H"
#include "baci_lib_discret.H"
#include "baci_lib_element.H"
#include "baci_linalg_blocksparsematrix.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_linalg_serialdensevector.H"
#include "baci_linalg_sparsematrix.H"
#include "baci_utils_exceptions.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/

void FBI::UTILS::FBIBlockAssemblyStrategy::AssembleFluidMatrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> cff, int elegid, const std::vector<int>& lmstride,
    const CORE::LINALG::SerialDenseMatrix& elemat, const std::vector<int>& lmrow,
    const std::vector<int>& lmrowowner, const std::vector<int>& lmcol)
{
  cff->Assemble(elegid, lmstride, elemat, lmrow, lmrowowner, lmcol);
}
