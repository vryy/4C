/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble the fbi coupling contributions


\level 1

*/
/*-----------------------------------------------------------*/

#include "fbi_fluid_assembly_strategy.H"
#include "lib_discret.H"
#include "lib_dserror.H"
#include "lib_element.H"
#include "beam3_base.H"
#include "linalg_serialdensevector.H"
#include "linalg_serialdensematrix.H"
#include "linalg_sparsematrix.H"
#include "linalg_blocksparsematrix.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/

void FBI::UTILS::FBIAssemblyStrategy::Assemble(const DRT::Discretization& discretization1,
    const DRT::Discretization& discretization2, std::vector<int> const& elegid,
    std::vector<LINALG::SerialDenseVector> const& elevec,
    std::vector<std::vector<LINALG::SerialDenseMatrix>> const& elemat,
    Teuchos::RCP<Epetra_FEVector>& f1, Teuchos::RCP<Epetra_FEVector>& f2,
    Teuchos::RCP<LINALG::SparseMatrix>& c11, Teuchos::RCP<LINALG::SparseOperator> c22,
    Teuchos::RCP<LINALG::SparseMatrix>& c12, Teuchos::RCP<LINALG::SparseMatrix>& c21)
{
  // the entries of elevecX  belong to the Dofs of the element with GID elegidX
  // the rows    of elematXY belong to the Dofs of the element with GID elegidX
  // the columns of elematXY belong to the Dofs of the element with GID elegidY
  const DRT::Element* ele1 = discretization1.gElement(elegid[0]);
  const DRT::Element* ele2 = discretization2.gElement(elegid[1]);

  // get element location vector and ownerships
  std::vector<int> lmrow1;
  std::vector<int> lmrow2;
  std::vector<int> lmrowowner1;
  std::vector<int> lmrowowner2;
  std::vector<int> lmstride;

  ele1->LocationVector(discretization1, lmrow1, lmrowowner1, lmstride);
  ele2->LocationVector(discretization2, lmrow2, lmrowowner2, lmstride);

  // assemble both element vectors into global system vector
  if (f1 != Teuchos::null)
  {
    f1->SumIntoGlobalValues(elevec[0].Length(), &lmrow1[0], elevec[0].Values());
  }
  if (f2 != Teuchos::null)
  {
    f2->SumIntoGlobalValues(elevec[1].Length(), &lmrow2[0], elevec[1].Values());
  }

  // and finally also assemble stiffness contributions
  if (c11 != Teuchos::null)
  {
    c11->FEAssemble(elemat[0][0], lmrow1, lmrow1);
  }
  if (c12 != Teuchos::null)
  {
    c12->FEAssemble(elemat[0][1], lmrow1, lmrow2);
  }
  if (c21 != Teuchos::null)
  {
    c21->FEAssemble(elemat[1][0], lmrow2, lmrow1);
  }
  if (c22 != Teuchos::null)
  {
    AssembleFluidMatrix(c22, elegid[1], lmstride, elemat[1][1], lmrow2, lmrowowner2, lmrow2);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/

void FBI::UTILS::FBIAssemblyStrategy::AssembleFluidMatrix(Teuchos::RCP<LINALG::SparseOperator> cff,
    int elegid, const std::vector<int>& lmstride, const Epetra_SerialDenseMatrix& elemat,
    const std::vector<int>& lmrow, const std::vector<int>& lmrowowner,
    const std::vector<int>& lmcol)
{
  Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(cff, true)->FEAssemble(elemat, lmrow, lmcol);
}
