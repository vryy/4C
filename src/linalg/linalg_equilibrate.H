/*----------------------------------------------------------------------*/
/*! \file

\brief Utilities for matrix equilibration

\level 1
*/
/*---------------------------------------------------------------------*/

#ifndef LINALG_EQUILIBRATE_H
#define LINALG_EQUILIBRATE_H

#include <Teuchos_RCP.hpp>
#include <Epetra_CrsMatrix.h>

namespace LINALG
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
    rowsandcolumns_maindiag
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
   * In case of block matrices the inverse of the scale factor can be calculated from the rows /
   * columns of main-diagonal block or of the entire matrix (EquilibrationMethod
   * xxx_full/main_diag).
   *
   * Note: Call EquilibrateSystem() before Solve() and UnequilibrateIncrement() after Solve() within
   * specific problem
   */
  class Equilibration
  {
   public:
    //! constructor
    Equilibration(EquilibrationMethod method, const Epetra_Map& dofrowmap);

    //! destructor
    virtual ~Equilibration() = default;

    /*!
     * @brief equilibrate global system of equations if necessary
     *
     * @param[in,out] systemmatrix  system matrix
     * @param[in,out] residual      residual vector
     * @param[in]     blockmaps     (block) map(s) of system matrix
     */
    virtual void EquilibrateSystem(const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        const Teuchos::RCP<Epetra_Vector>& residual,
        const LINALG::MultiMapExtractor& blockmaps) const = 0;

    /*!
     * @brief unequilibrate global increment vector if necessary
     *
     * @param[in,out] increment  increment vector
     */
    virtual void UnequilibrateIncrement(const Teuchos::RCP<Epetra_Vector>& increment) const;

   protected:
    /*!
     * @brief compute inverse sums of absolute values of matrix column entries
     *
     * @param[in]  matrix      matrix
     * @param[out] invcolsums  inverse sums of absolute values of column entries in matrix
     */
    void ComputeInvColSums(
        const LINALG::SparseMatrix& matrix, const Teuchos::RCP<Epetra_Vector>& invcolsums) const;

    /*!
     * @brief compute inverse sums of absolute values of matrix row entries
     *
     * @param[in]  matrix      matrix
     * @param[out] invrowsums  inverse sums of absolute values of row entries in matrix
     */
    void ComputeInvRowSums(
        const LINALG::SparseMatrix& matrix, const Teuchos::RCP<Epetra_Vector>& invrowsums) const;

    /*!
     * @brief equilibrate matrix columns
     *
     * @param[in,out] matrix      matrix
     * @param[in]     invcolsums  sums of absolute values of column entries in matrix
     */
    void EquilibrateMatrixColumns(
        LINALG::SparseMatrix& matrix, const Teuchos::RCP<Epetra_Vector>& invcolsums) const;

    /*!
     * @brief equilibrate matrix rows
     *
     * @param[in,out] matrix      matrix
     * @param[in]     invrowsums  sums of absolute values of row entries in matrix
     */
    void EquilibrateMatrixRows(
        LINALG::SparseMatrix& matrix, const Teuchos::RCP<Epetra_Vector>& invrowsums) const;

    /*!
     * @brief equilibrate right hand side
     *
     * @param[in,out] residual  residual vector
     */
    void EquilibrateRHS(const Teuchos::RCP<Epetra_Vector>& residual) const;

    //! return equilibration method
    EquilibrationMethod Method() const { return method_; }

    //! inverse sums of absolute values of column entries in global system matrix
    Teuchos::RCP<Epetra_Vector> invcolsums_;

    //! inverse sums of absolute values of row entries in global system matrix
    Teuchos::RCP<Epetra_Vector> invrowsums_;

   private:
    //! equilibration method
    EquilibrationMethod method_;
  };

  //! System matrix is a sparse matrix
  class EquilibrationSparse : public Equilibration
  {
   public:
    //! constructor
    EquilibrationSparse(EquilibrationMethod method, const Epetra_Map& dofrowmap)
        : Equilibration(method, dofrowmap)
    {
    }

    ~EquilibrationSparse() override = default;

    void EquilibrateSystem(const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        const Teuchos::RCP<Epetra_Vector>& residual,
        const LINALG::MultiMapExtractor& blockmaps) const override;
  };

  //! System matrix is a block matrix
  class EquilibrationBlock : public Equilibration
  {
   public:
    //! constructor
    EquilibrationBlock(EquilibrationMethod method, const Epetra_Map& dofrowmap)
        : Equilibration(method, dofrowmap)
    {
    }

    ~EquilibrationBlock() override = default;

    void EquilibrateSystem(const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        const Teuchos::RCP<Epetra_Vector>& residual,
        const LINALG::MultiMapExtractor& blockmaps) const override;
  };

  //! System matrix can be a sparse or block matrix
  class EquilibrationNone : public Equilibration
  {
   public:
    //! constructor
    EquilibrationNone(EquilibrationMethod method, const Epetra_Map& dofrowmap)
        : Equilibration(method, dofrowmap)
    {
    }

    ~EquilibrationNone() override = default;

    void EquilibrateSystem(const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        const Teuchos::RCP<Epetra_Vector>& residual,
        const LINALG::MultiMapExtractor& blockmaps) const override
    {
    }

    void UnequilibrateIncrement(const Teuchos::RCP<Epetra_Vector>& increment) const override {}
  };

  //! Factory for equilibration
  class EquilibrationFactory
  {
   public:
    //! constructor
    EquilibrationFactory() = default;

    /*!
     * @brief build specific equilibration method
     *
     * @param[in] type   type of matrix
     * @param[in] method equilibration method
     * @param dofrowmap  degree of freedom row map
     *
     * @return equilibration method
     */
    Teuchos::RCP<LINALG::Equilibration> BuildEquilibration(
        MatrixType type, EquilibrationMethod method, const Epetra_Map& dofrowmap);
  };

}  // namespace LINALG
#endif
