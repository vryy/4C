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
   * - Combining equilibration of rows and colums is possible.
   *
   * In case of block matrices the inverse of the scale factor can be calculated from the rows /
   * colums of main-diagonal block or of the entire matrix (EquilibrationMethod xxx_full/main_diag).
   *
   * Note: Call EquilibrateSystem() before Solve() and UnequilibrateIncrement() after Solve() within
   * specific problem
   */
  class Equilibration
  {
   public:
    //! constructor
    Equilibration(EquilibrationMethod method, Epetra_Map dofrowmap);

    //! destructor
    virtual ~Equilibration() = default;

    //! equilibrate global system of equations if necessary
    virtual void EquilibrateSystem(const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        const Teuchos::RCP<Epetra_Vector>& residual,
        const LINALG::MultiMapExtractor& blockmaps) const = 0;

    //! unequilibrate global increment vector if necessary
    virtual void UnequilibrateIncrement(const Teuchos::RCP<Epetra_Vector>& increment) const;

   protected:
    //! compute inverse sums of absolute values of matrix column entries
    void ComputeInvColSums(
        const LINALG::SparseMatrix& matrix, const Teuchos::RCP<Epetra_Vector>& invcolsums) const;

    //! compute inverse sums of absolute values of matrix row entries
    void ComputeInvRowSums(
        const LINALG::SparseMatrix& matrix, const Teuchos::RCP<Epetra_Vector>& invrowsums) const;

    //! equilibrate matrix columns
    void EquilibrateMatrixColumns(
        LINALG::SparseMatrix& matrix, const Teuchos::RCP<Epetra_Vector>& invcolsums) const;

    //! equilibrate matrix rows
    void EquilibrateMatrixRows(
        LINALG::SparseMatrix& matrix, const Teuchos::RCP<Epetra_Vector>& invrowsums) const;

    //! equilibrate right hand side
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
    EquilibrationSparse(EquilibrationMethod method, Epetra_Map dofrowmap)
        : Equilibration(method, dofrowmap)
    {
      return;
    }

    //! destructor
    ~EquilibrationSparse() override = default;

    //! equilibrate global system of equations if necessary
    void EquilibrateSystem(const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        const Teuchos::RCP<Epetra_Vector>& residual,
        const LINALG::MultiMapExtractor& blockmaps) const override;
  };

  //! System matrix is a block matrix
  class EquilibrationBlock : public Equilibration
  {
   public:
    //! constructor
    EquilibrationBlock(EquilibrationMethod method, Epetra_Map dofrowmap)
        : Equilibration(method, dofrowmap)
    {
      return;
    }

    //! destructor
    ~EquilibrationBlock() override = default;

    //! equilibrate global system of equations if necessary
    void EquilibrateSystem(const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        const Teuchos::RCP<Epetra_Vector>& residual,
        const LINALG::MultiMapExtractor& blockmaps) const override;
  };

  //! System matrix is a block matrix
  class EquilibrationNone : public Equilibration
  {
   public:
    //! constructor
    EquilibrationNone(EquilibrationMethod method, Epetra_Map dofrowmap)
        : Equilibration(method, dofrowmap)
    {
      return;
    }

    //! destructor
    ~EquilibrationNone() override = default;

    //! equilibrate global system of equations if necessary
    void EquilibrateSystem(const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        const Teuchos::RCP<Epetra_Vector>& residual,
        const LINALG::MultiMapExtractor& blockmaps) const override
    {
      return;
    }

    //! unequilibrate global increment vector if necessary
    void UnequilibrateIncrement(const Teuchos::RCP<Epetra_Vector>& increment) const override
    {
      return;
    }
  };

  //! Factory for equilibration
  class EquilibrationFactory
  {
   public:
    //! constructor
    EquilibrationFactory() { return; }

    //! build specifc equilibration type
    Teuchos::RCP<LINALG::Equilibration> BuildEquilibration(
        MatrixType type, EquilibrationMethod method, Epetra_Map dofrowmap);
  };

}  // namespace LINALG
#endif
