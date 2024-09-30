/*----------------------------------------------------------------------*/
/*! \file

\brief Model Order Reduction (MOR) using Proper Orthogonal Decomposition (POD)

\level 2


 *----------------------------------------------------------------------*/

#ifndef FOUR_C_CARDIOVASCULAR0D_MOR_POD_HPP
#define FOUR_C_CARDIOVASCULAR0D_MOR_POD_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_MultiVector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::LinAlg
{
  class SparseMatrix;
}  // namespace Core::LinAlg
namespace Cardiovascular0D
{
  class ProperOrthogonalDecomposition
  {
   public:
    /*!
        \brief Constructor
     */
    ProperOrthogonalDecomposition(Teuchos::RCP<const Epetra_Map> full_model_dof_row_map_,
        const std::string& pod_matrix_file_name, const std::string& absolute_path_to_input_file);

    //! M_red = V^T * M * V
    Teuchos::RCP<Core::LinAlg::SparseMatrix> reduce_diagnoal(
        Teuchos::RCP<Core::LinAlg::SparseMatrix> M);

    //! M_red = V^T * M
    Teuchos::RCP<Core::LinAlg::SparseMatrix> reduce_off_diagonal(
        Teuchos::RCP<Core::LinAlg::SparseMatrix> M);

    //! v_red = V^T * v
    Teuchos::RCP<Epetra_MultiVector> reduce_rhs(Teuchos::RCP<Epetra_MultiVector> v);

    //! v_red = V^T * v
    Teuchos::RCP<Core::LinAlg::Vector<double>> reduce_residual(
        Teuchos::RCP<Core::LinAlg::Vector<double>> v);

    //! v = V * v_red
    Teuchos::RCP<Core::LinAlg::Vector<double>> extend_solution(
        Teuchos::RCP<Core::LinAlg::Vector<double>> v);

    bool have_mor() { return havemor_; };

    int get_red_dim() { return projmatrix_->NumVectors(); };

   private:
    /*! \brief Read POD basis vectors from file
     *
     * Read matrix from specified binary file
     * The binary file has to be formatted like this:
     * Number of Rows: int
     * Number of Columns: int
     * Values (row-wise): float
     */
    void read_pod_basis_vectors_from_file(
        const std::string& absolute_path_to_pod_file, Teuchos::RCP<Epetra_MultiVector>& projmatrix);

    //! Multiply two Epetra MultiVectors
    void multiply_epetra_multi_vectors(Teuchos::RCP<Epetra_MultiVector>, char,
        Teuchos::RCP<Epetra_MultiVector>, char, Teuchos::RCP<Epetra_Map>,
        Teuchos::RCP<Epetra_Import>, Teuchos::RCP<Epetra_MultiVector>);

    //! Epetra_MultiVector to Core::LinAlg::SparseMatrix
    void epetra_multi_vector_to_linalg_sparse_matrix(Teuchos::RCP<Epetra_MultiVector> multivect,
        Teuchos::RCP<Epetra_Map> rangemap, Teuchos::RCP<Epetra_Map> domainmap,
        Teuchos::RCP<Core::LinAlg::SparseMatrix> sparsemat);

    //! Check orthogonality of POD basis vectors with M^T * M - I == 0
    bool is_pod_basis_orthogonal(const Epetra_MultiVector& M);

    /// DOF row map of the full model, i.e. map of POD basis vectors
    Teuchos::RCP<const Epetra_Map> full_model_dof_row_map_;

    //! Flag to indicate usage of model order reduction
    bool havemor_ = false;

    //! Projection matrix for POD
    Teuchos::RCP<Epetra_MultiVector> projmatrix_;

    //! Unique map of structure dofs after POD-MOR
    Teuchos::RCP<Epetra_Map> structmapr_;

    //! Full redundant map of structure dofs after POD-MOR
    Teuchos::RCP<Epetra_Map> redstructmapr_;

    //! Importer for fully redundant map of structure dofs after POD-MOR into distributed one
    Teuchos::RCP<Epetra_Import> structrimpo_;

    //! Importer for distributed map of structure dofs after POD-MOR into fully redundant one
    Teuchos::RCP<Epetra_Import> structrinvimpo_;

  };  // class
}  // namespace Cardiovascular0D
FOUR_C_NAMESPACE_CLOSE

#endif
