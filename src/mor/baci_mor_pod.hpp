/*----------------------------------------------------------------------*/
/*! \file

\brief Model Order Reduction (MOR) using Proper Orthogonal Decomposition (POD)

\level 2


 *----------------------------------------------------------------------*/

#ifndef BACI_MOR_POD_HPP
#define BACI_MOR_POD_HPP

#include "baci_config.hpp"

#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace CORE::LINALG
{
  class SparseMatrix;
  //  class DefaultBlockMatrixStrategy;
  //  class BlockSparseMatrix;
}  // namespace CORE::LINALG
namespace MOR
{
  class ProperOrthogonalDecomposition
  {
   public:
    /*!
        \brief Constructor
     */
    ProperOrthogonalDecomposition(Teuchos::RCP<DRT::Discretization> discr);


    /*!
        \brief Read in binary matrix from file
     */
    /* Read matrix from specified binary file
     * The binary file has to be formatted like this:
     * Number of Rows: int
     * Number of Columns: int
     * Values (row-wise): float
     */
    void ReadMatrix(std::string filename, Teuchos::RCP<Epetra_MultiVector>& projmatrix);

    Teuchos::RCP<CORE::LINALG::SparseMatrix> ReduceDiagnoal(
        Teuchos::RCP<CORE::LINALG::SparseMatrix> M);

    Teuchos::RCP<CORE::LINALG::SparseMatrix> ReduceOffDiagonal(
        Teuchos::RCP<CORE::LINALG::SparseMatrix> M);

    Teuchos::RCP<Epetra_MultiVector> ReduceRHS(Teuchos::RCP<Epetra_MultiVector> v);

    Teuchos::RCP<Epetra_Vector> ReduceResidual(Teuchos::RCP<Epetra_Vector> v);

    Teuchos::RCP<Epetra_Vector> ExtendSolution(Teuchos::RCP<Epetra_Vector> v);

    bool HaveMOR() { return havemor_; };

    Teuchos::RCP<Epetra_MultiVector> GetPODMatrix() { return projmatrix_; };

    int GetRedDim() { return projmatrix_->NumVectors(); };

   private:
    void MultiplyEpetraMultiVectors(Teuchos::RCP<Epetra_MultiVector>, char,
        Teuchos::RCP<Epetra_MultiVector>, char, Teuchos::RCP<Epetra_Map>,
        Teuchos::RCP<Epetra_Import>, Teuchos::RCP<Epetra_MultiVector>);

    void EpetraMultiVectorToLINALGSparseMatrix(Teuchos::RCP<Epetra_MultiVector> multivect,
        Teuchos::RCP<Epetra_Map> rangemap, Teuchos::RCP<Epetra_Map> domainmap,
        Teuchos::RCP<CORE::LINALG::SparseMatrix> sparsemat);

    bool IsOrthogonal(Teuchos::RCP<Epetra_MultiVector> M);

    Teuchos::RCP<DRT::Discretization> actdisc_;
    int myrank_;  //!< ID of actual processor in parallel
    Teuchos::ParameterList morparams_;
    bool havemor_;
    Teuchos::RCP<Epetra_MultiVector> projmatrix_;  //!< projection matrix for POD


    Teuchos::RCP<Epetra_Map> structmapr_;  ///< unique map of structure dofs after POD-MOR
    Teuchos::RCP<Epetra_Map>
        redstructmapr_;  ///< full redundant map of structure dofs after POD-MOR
    Teuchos::RCP<Epetra_Import> structrimpo_;     ///< importer for fully redundant map of structure
                                                  ///< dofs after POD-MOR into distributed one
    Teuchos::RCP<Epetra_Import> structrinvimpo_;  ///< importer for distributed map of structure
                                                  ///< dofs after POD-MOR into fully redundant one

  };  // class
}  // namespace MOR
BACI_NAMESPACE_CLOSE

#endif
