/*----------------------------------------------------------------------*/
/*! \file

\brief Model Order Reduction (MOR) using Proper Orthogonal Decomposition (POD)

\level 2


 *----------------------------------------------------------------------*/

#ifndef FOUR_C_MOR_POD_HPP
#define FOUR_C_MOR_POD_HPP

#include "4C_config.hpp"

#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::LinAlg
{
  class SparseMatrix;
  //  class DefaultBlockMatrixStrategy;
  //  class BlockSparseMatrix;
}  // namespace Core::LinAlg
namespace ModelOrderRed
{
  class ProperOrthogonalDecomposition
  {
   public:
    /*!
        \brief Constructor
     */
    ProperOrthogonalDecomposition(Teuchos::RCP<Core::FE::Discretization> discr);


    /*!
        \brief Read in binary matrix from file
     */
    /* Read matrix from specified binary file
     * The binary file has to be formatted like this:
     * Number of Rows: int
     * Number of Columns: int
     * Values (row-wise): float
     */
    void read_matrix(std::string filename, Teuchos::RCP<Epetra_MultiVector>& projmatrix);

    Teuchos::RCP<Core::LinAlg::SparseMatrix> reduce_diagnoal(
        Teuchos::RCP<Core::LinAlg::SparseMatrix> M);

    Teuchos::RCP<Core::LinAlg::SparseMatrix> reduce_off_diagonal(
        Teuchos::RCP<Core::LinAlg::SparseMatrix> M);

    Teuchos::RCP<Epetra_MultiVector> reduce_rhs(Teuchos::RCP<Epetra_MultiVector> v);

    Teuchos::RCP<Epetra_Vector> reduce_residual(Teuchos::RCP<Epetra_Vector> v);

    Teuchos::RCP<Epetra_Vector> extend_solution(Teuchos::RCP<Epetra_Vector> v);

    bool have_mor() { return havemor_; };

    Teuchos::RCP<Epetra_MultiVector> get_pod_matrix() { return projmatrix_; };

    int get_red_dim() { return projmatrix_->NumVectors(); };

   private:
    void multiply_epetra_multi_vectors(Teuchos::RCP<Epetra_MultiVector>, char,
        Teuchos::RCP<Epetra_MultiVector>, char, Teuchos::RCP<Epetra_Map>,
        Teuchos::RCP<Epetra_Import>, Teuchos::RCP<Epetra_MultiVector>);

    void epetra_multi_vector_to_linalg_sparse_matrix(Teuchos::RCP<Epetra_MultiVector> multivect,
        Teuchos::RCP<Epetra_Map> rangemap, Teuchos::RCP<Epetra_Map> domainmap,
        Teuchos::RCP<Core::LinAlg::SparseMatrix> sparsemat);

    bool is_orthogonal(Teuchos::RCP<Epetra_MultiVector> M);

    Teuchos::RCP<Core::FE::Discretization> actdisc_;
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
}  // namespace ModelOrderRed
FOUR_C_NAMESPACE_CLOSE

#endif
