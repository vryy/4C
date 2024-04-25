/*----------------------------------------------------------------------*/
/*! \file

\brief Utilities for matrix equilibration

\level 1
*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_LINALG_EQUILIBRATE_HPP
#define FOUR_C_LINALG_EQUILIBRATE_HPP

#include "4C_config.hpp"

#include <Epetra_CrsMatrix.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINALG
{
  enum class MatrixType;
  class MultiMapExtractor;
  class SparseMatrix;
  class SparseOperator;

  //! method, how a system of linear equations is equilibrated.
  enum class EquilibrationMethod
  {
    none,
    rows_full,
    rows_maindiag,
    columns_full,
    columns_maindiag,
    rowsandcolumns_full,
    rowsandcolumns_maindiag,
    symmetry,
    local
  };

  /*!
   * \brief Equilibrates linear system of equations (matrix and right hand side).
   *
   * Equilibration can be done by row or by column.
   * - By row means scaling each row (matrix and right hand side) with a scalar (inverse of sum of
   * matrix entries of the row).
   * - By column means scaling each column of the matrix with a scalar (inverse of sum of matrix
   * entries of the column). This could be seen as defining an individual system of units for each
   * dof. Thus, after solving the solution vector has to be rescaled (in UnequilibrateIncrement()).
   * - Combining equilibration of rows and columns is possible.
   *
   * Equilibration can be done for symmetric matrix in a way that it keeps the symmetry and has
   * maximum value of 1.0 in each column and row
   *
   * In case of block matrices the inverse of the scale factor can be calculated from the rows /
   * columns of main-diagonal block or of the entire matrix (EquilibrationMethod
   * xxx_full/maindiag).
   *
   * In case of block matrices each block can be equilibrated in with a different method
   *
   * Note: Call EquilibrateSystem() before Solve() and UnequilibrateIncrement() after Solve() within
   * specific problem. In case of linear problems, where the matrices do not change, it is possible
   * to use EquilibrateMatrix() after assembly and EquilibrateRHS() before every Solve() call.
   */
  class Equilibration
  {
   public:
    Equilibration(Teuchos::RCP<const Epetra_Map> dofrowmap);

    virtual ~Equilibration() = default;

    /*!
     * @brief equilibrate matrix if necessary
     *
     * @param[in,out] systemmatrix  system matrix
     * @param[in]     blockmaps     (block) map(s) of system matrix
     */
    virtual void EquilibrateMatrix(Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> blockmaps) const = 0;

    /*!
     * @brief equilibrate right hand side
     *
     * @param[in,out] residual  residual vector
     */
    virtual void EquilibrateRHS(Teuchos::RCP<Epetra_Vector> residual) const = 0;

    /*!
     * @brief equilibrate global system of equations if necessary
     *
     * @param[in,out] systemmatrix  system matrix
     * @param[in,out] residual      residual vector
     * @param[in]     blockmaps     (block) map(s) of system matrix
     */
    void EquilibrateSystem(Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<Epetra_Vector> residual,
        Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> blockmaps) const;

    /*!
     * @brief unequilibrate global increment vector if necessary
     *
     * @param[in,out] increment  increment vector
     */
    virtual void UnequilibrateIncrement(Teuchos::RCP<Epetra_Vector> increment) const = 0;

   protected:
    /*!
     * @brief compute inverse sums of absolute values of matrix column entries
     *
     * @param[in]  matrix      matrix
     * @param[out] invcolsums  inverse sums of absolute values of column entries in matrix
     */
    void ComputeInvColSums(const CORE::LINALG::SparseMatrix& matrix,
        Teuchos::RCP<Epetra_Vector> invcolsums, const EquilibrationMethod method) const;

    /*!
     * @brief compute inverse sums of absolute values of matrix row entries
     *
     * @param[in]  matrix      matrix
     * @param[out] invrowsums  inverse sums of absolute values of row entries in matrix
     */
    void ComputeInvRowSums(const CORE::LINALG::SparseMatrix& matrix,
        Teuchos::RCP<Epetra_Vector> invrowsums, const EquilibrationMethod method) const;

    /*!
     * @brief compute scaling of matrix to keep symmetry
     *
     * Algorithm: d_i = (max(sqrt(A_ii), p_i))^(-1), with p_i = max_{1<=j<=i-1}(d_j*A_ij)
     *
     * @param[in]  matrix      symmetric matrix A
     * @param[out] invsymmetry scale vector to keep symmetry in matrix d
     */
    void ComputeInvSymmetry(
        const CORE::LINALG::SparseMatrix& matrix, Teuchos::RCP<Epetra_Vector> invsymmetry) const;

    /*!
     * @brief equilibrate matrix columns
     *
     * @param[in,out] matrix      matrix
     * @param[in]     invcolsums  sums of absolute values of column entries in matrix
     */
    void EquilibrateMatrixColumns(
        CORE::LINALG::SparseMatrix& matrix, Teuchos::RCP<const Epetra_Vector> invcolsums) const;

    /*!
     * @brief equilibrate matrix rows
     *
     * @param[in,out] matrix      matrix
     * @param[in]     invrowsums  sums of absolute values of row entries in matrix
     */
    void EquilibrateMatrixRows(
        CORE::LINALG::SparseMatrix& matrix, Teuchos::RCP<const Epetra_Vector> invrowsums) const;

    //! inverse sums of absolute values of column entries in global system matrix
    Teuchos::RCP<Epetra_Vector> invcolsums_;

    //! inverse sums of absolute values of row entries in global system matrix
    Teuchos::RCP<Epetra_Vector> invrowsums_;
  };

  /*----------------------------------------------------------------------*
   *----------------------------------------------------------------------*/
  //! Global system matrix is equilibrated using the same method
  class EquilibrationUniversal : public Equilibration
  {
   public:
    EquilibrationUniversal(EquilibrationMethod method, Teuchos::RCP<const Epetra_Map> dofrowmap);
    //! return equilibration method
    EquilibrationMethod Method() const { return method_; }

    void EquilibrateRHS(Teuchos::RCP<Epetra_Vector> residual) const override;

    void UnequilibrateIncrement(Teuchos::RCP<Epetra_Vector> increment) const override;

   private:
    EquilibrationMethod method_;
  };

  /*----------------------------------------------------------------------*
   *----------------------------------------------------------------------*/
  //! System matrix is a sparse matrix
  class EquilibrationSparse : public EquilibrationUniversal
  {
   public:
    EquilibrationSparse(EquilibrationMethod method, Teuchos::RCP<const Epetra_Map> dofrowmap);
    void EquilibrateMatrix(Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> blockmaps) const override;

    /*!
     * @brief equilibrate matrix if necessary
     *
     * @param[in,out] systemmatrix  system matrix
     */
    void EquilibrateMatrix(Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix) const;
  };

  /*----------------------------------------------------------------------*
   *----------------------------------------------------------------------*/
  //! System matrix is a block matrix
  class EquilibrationBlock : public EquilibrationUniversal
  {
   public:
    EquilibrationBlock(EquilibrationMethod method, Teuchos::RCP<const Epetra_Map> dofrowmap);
    void EquilibrateMatrix(Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> blockmaps) const override;
  };

  /*----------------------------------------------------------------------*
   *----------------------------------------------------------------------*/
  //! System matrix is a block matrix and every block should be scaled individually
  class EquilibrationBlockSpecific : public Equilibration
  {
   public:
    EquilibrationBlockSpecific(
        const std::vector<EquilibrationMethod>& method, Teuchos::RCP<const Epetra_Map> dofrowmap);
    void EquilibrateMatrix(Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> blockmaps) const override;

    void EquilibrateRHS(const Teuchos::RCP<Epetra_Vector> residual) const override;

    void UnequilibrateIncrement(Teuchos::RCP<Epetra_Vector> increment) const override;

   private:
    std::vector<EquilibrationMethod> method_blocks_;
  };

  /*----------------------------------------------------------------------*
   *----------------------------------------------------------------------*/
  //! System matrix can be a sparse or block matrix, but no equilibration is performed
  class EquilibrationNone : public Equilibration
  {
   public:
    EquilibrationNone(Teuchos::RCP<const Epetra_Map> dofrowmap) : Equilibration(dofrowmap) {}
    void EquilibrateMatrix(Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix,
        Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> blockmaps) const override
    {
    }

    void EquilibrateRHS(Teuchos::RCP<Epetra_Vector> residual) const override {}

    void UnequilibrateIncrement(Teuchos::RCP<Epetra_Vector> increment) const override {}
  };

  /*----------------------------------------------------------------------*
   *----------------------------------------------------------------------*/
  /*! @brief build specific equilibration method
   *
   * @param[in] type   type of matrix
   * @param[in] method equilibration method
   * @param dofrowmap  degree of freedom row map
   *
   * @return equilibration method
   */
  Teuchos::RCP<CORE::LINALG::Equilibration> BuildEquilibration(MatrixType type,
      const std::vector<EquilibrationMethod>& method, Teuchos::RCP<const Epetra_Map> dofrowmap);
}  // namespace CORE::LINALG
FOUR_C_NAMESPACE_CLOSE

#endif
