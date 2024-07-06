/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble the fbi coupling contributions


\level 1

*/
/*-----------------------------------------------------------*/

#include "4C_fbi_fluid_assembly_strategy.hpp"

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

void FBI::UTILS::FBIAssemblyStrategy::assemble(const Core::FE::Discretization& discretization1,
    const Core::FE::Discretization& discretization2, std::vector<int> const& elegid,
    std::vector<Core::LinAlg::SerialDenseVector> const& elevec,
    std::vector<std::vector<Core::LinAlg::SerialDenseMatrix>> const& elemat,
    Teuchos::RCP<Epetra_FEVector>& f1, Teuchos::RCP<Epetra_FEVector>& f2,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& c11, Teuchos::RCP<Core::LinAlg::SparseOperator> c22,
    Teuchos::RCP<Core::LinAlg::SparseMatrix>& c12, Teuchos::RCP<Core::LinAlg::SparseMatrix>& c21)
{
  // the entries of elevecX  belong to the Dofs of the element with GID elegidX
  // the rows    of elematXY belong to the Dofs of the element with GID elegidX
  // the columns of elematXY belong to the Dofs of the element with GID elegidY
  const Core::Elements::Element* ele1 = discretization1.g_element(elegid[0]);
  const Core::Elements::Element* ele2 = discretization2.g_element(elegid[1]);

  // get element location vector and ownerships
  std::vector<int> lmrow1;
  std::vector<int> lmrow2;
  std::vector<int> lmrowowner1;
  std::vector<int> lmrowowner2;
  std::vector<int> lmstride;

  ele1->location_vector(discretization1, lmrow1, lmrowowner1, lmstride);
  ele2->location_vector(discretization2, lmrow2, lmrowowner2, lmstride);

  // assemble both element vectors into global system vector
  if (f1 != Teuchos::null)
  {
    f1->SumIntoGlobalValues(elevec[0].length(), lmrow1.data(), elevec[0].values());
  }
  if (f2 != Teuchos::null)
  {
    f2->SumIntoGlobalValues(elevec[1].length(), lmrow2.data(), elevec[1].values());
  }

  // and finally also assemble stiffness contributions
  if (c11 != Teuchos::null)
  {
    c11->fe_assemble(elemat[0][0], lmrow1, lmrow1);
  }
  if (c12 != Teuchos::null)
  {
    c12->fe_assemble(elemat[0][1], lmrow1, lmrow2);
  }
  if (c21 != Teuchos::null)
  {
    c21->fe_assemble(elemat[1][0], lmrow2, lmrow1);
  }
  if (c22 != Teuchos::null)
  {
    assemble_fluid_matrix(c22, elegid[1], lmstride, elemat[1][1], lmrow2, lmrowowner2, lmrow2);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/

void FBI::UTILS::FBIAssemblyStrategy::assemble_fluid_matrix(
    Teuchos::RCP<Core::LinAlg::SparseOperator> cff, int elegid, const std::vector<int>& lmstride,
    const Core::LinAlg::SerialDenseMatrix& elemat, const std::vector<int>& lmrow,
    const std::vector<int>& lmrowowner, const std::vector<int>& lmcol)
{
  Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(cff, true)->fe_assemble(
      elemat, lmrow, lmcol);
}

FOUR_C_NAMESPACE_CLOSE
