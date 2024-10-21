#ifndef FOUR_C_CARDIOVASCULAR0D_MOR_POD_HPP
#define FOUR_C_CARDIOVASCULAR0D_MOR_POD_HPP

#include "4C_config.hpp"

#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

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
    Teuchos::RCP<Core::LinAlg::SparseMatrix> reduce_diagnoal(Core::LinAlg::SparseMatrix& M);

    //! M_red = V^T * M
    Teuchos::RCP<Core::LinAlg::SparseMatrix> reduce_off_diagonal(Core::LinAlg::SparseMatrix& M);

    //! v_red = V^T * v
    Teuchos::RCP<Core::LinAlg::MultiVector<double>> reduce_rhs(
        Core::LinAlg::MultiVector<double>& v);

    //! v_red = V^T * v
    Teuchos::RCP<Core::LinAlg::Vector<double>> reduce_residual(Core::LinAlg::Vector<double>& v);

    //! v = V * v_red
    Teuchos::RCP<Core::LinAlg::Vector<double>> extend_solution(Core::LinAlg::Vector<double>& v);

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
    void read_pod_basis_vectors_from_file(const std::string& absolute_path_to_pod_file,
        Teuchos::RCP<Core::LinAlg::MultiVector<double>>& projmatrix);

    //! Multiply two Epetra MultiVectors
    void multiply_epetra_multi_vectors(Core::LinAlg::MultiVector<double>&, char,
        Core::LinAlg::MultiVector<double>&, char, Epetra_Map&, Epetra_Import&,
        Core::LinAlg::MultiVector<double>&);

    //! Core::LinAlg::MultiVector<double> to Core::LinAlg::SparseMatrix
    void epetra_multi_vector_to_linalg_sparse_matrix(Core::LinAlg::MultiVector<double>& multivect,
        Epetra_Map& rangemap, Teuchos::RCP<Epetra_Map> domainmap,
        Core::LinAlg::SparseMatrix& sparsemat);

    //! Check orthogonality of POD basis vectors with M^T * M - I == 0
    bool is_pod_basis_orthogonal(const Core::LinAlg::MultiVector<double>& M);

    /// DOF row map of the full model, i.e. map of POD basis vectors
    Teuchos::RCP<const Epetra_Map> full_model_dof_row_map_;

    //! Flag to indicate usage of model order reduction
    bool havemor_ = false;

    //! Projection matrix for POD
    Teuchos::RCP<Core::LinAlg::MultiVector<double>> projmatrix_;

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
