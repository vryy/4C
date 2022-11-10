/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble the fbi coupling contributions


\level 1

*/
/*-----------------------------------------------------------*/

#include "./fluidblockmatrix_assembly_strategy.H"
#include "drt_discret.H"
#include "drt_dserror.H"
#include "drt_element.H"
#include "beam3_base.H"
#include "linalg_serialdensevector.H"
#include "linalg_serialdensematrix.H"
#include "linalg_sparsematrix.H"
#include "linalg_blocksparsematrix.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/

void FBI::UTILS::FBIBlockAssemblyStrategy::AssembleFluidMatrix(
    Teuchos::RCP<LINALG::SparseOperator> cff, int elegid, const std::vector<int>& lmstride,
    const Epetra_SerialDenseMatrix& elemat, const std::vector<int>& lmrow,
    const std::vector<int>& lmrowowner, const std::vector<int>& lmcol)
{
  cff->Assemble(elegid, lmstride, elemat, lmrow, lmrowowner, lmcol);
}
