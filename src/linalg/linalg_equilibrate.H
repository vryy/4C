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
  class MultiMapExtractor;
  class SparseMatrix;
  class SparseOperator;

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
   * Call Solve() within Problem between EquilibrateSystem() and UnequilibrateIncrement()
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
    void UnequilibrateIncrement(const Teuchos::RCP<Epetra_Vector>& increment) const;

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
    EquilibrationMethod Method() const { return method_; };

    //! inverse sums of absolute values of column entries in global system matrix
    Teuchos::RCP<Epetra_Vector> invcolsums_;

    //! inverse sums of absolute values of row entries in global system matrix
    Teuchos::RCP<Epetra_Vector> invrowsums_;

   private:
    //! equilibration method
    EquilibrationMethod method_;
  };

  /*!
   * System matrix is a sparse matrix
   */
  class EquilibrationSparse : public Equilibration
  {
   public:
    //! constructor
    EquilibrationSparse(EquilibrationMethod method, Epetra_Map dofrowmap);

    //! destructor
    virtual ~EquilibrationSparse() = default;

    //! equilibrate global system of equations if necessary
    void EquilibrateSystem(const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        const Teuchos::RCP<Epetra_Vector>& residual,
        const LINALG::MultiMapExtractor& blockmaps) const override;
  };

  /*!
   * System matrix is a block matrix
   */
  class EquilibrationBlock : public Equilibration
  {
   public:
    //! constructor
    EquilibrationBlock(EquilibrationMethod method, Epetra_Map dofrowmap);

    //! destructor
    virtual ~EquilibrationBlock() = default;

    //! equilibrate global system of equations if necessary
    void EquilibrateSystem(const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
        const Teuchos::RCP<Epetra_Vector>& residual,
        const LINALG::MultiMapExtractor& blockmaps) const override;
  };

}  // namespace LINALG
#endif
