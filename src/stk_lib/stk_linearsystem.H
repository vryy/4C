#ifdef STKADAPTIVE

#ifndef STK_LINEARSYSTEM_H
#define STK_LINEARSYSTEM_H

#include <map>
#include <set>

#include <Epetra_FECrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_Import.h>

#include <Teuchos_RCP.hpp>

namespace STK
{
  class Discretization;

/// a linear system consists of a system matrix, a vector of unknowns and a rhs.
/*!
  Support off-row assembly such that each element is just once evaluated.

  Support hanging nodes. Hanging node constraits are implicitly removed.

  \author u.kue
 */
class LinearSystem
{
public:

  /// construct from a discretization
  LinearSystem( STK::Discretization & dis );

  /// @name Assemble

  /// assemble element matrix
  void Assemble(int eid,
                const Epetra_SerialDenseMatrix& Aele,
                const std::vector<int>& lmrow,
                const std::vector<int>& lmrowowner,
                const std::vector<int>& lmcol);

  /// assemble rhs vector
  void Assemble(int eid,
                const Epetra_SerialDenseVector& Vele,
                const std::vector<int>& lm);

  /// complete any global objects
  void Complete();

  //@}

  /// return sparse matrix
  Teuchos::RCP<Epetra_FECrsMatrix> Matrix() { return sysmat_; }

  /// return rhs
  Teuchos::RCP<Epetra_Vector> RHS() { return rowrhs_; }

  /// zero sparse matrix and rhs vector
  void Zero();

  /// put Dirichlet values according to known Dirichlet map into sparse matrix
  /// and rhs vector
  void ApplyDirichlet( Epetra_Vector & dirichlet );

private:

  /// create sparse matrix graph including Dirichlet rows and hanging node rows
  Teuchos::RCP<Epetra_CrsGraph> MatrixGraph( STK::Discretization & dis );

  /// sparse matrix
  Teuchos::RCP<Epetra_FECrsMatrix> sysmat_;

  /// non-overlapping (row) rhs vector
  Teuchos::RCP<Epetra_Vector> rowrhs_;

  /// overlapping (column) rhs vector
  Teuchos::RCP<Epetra_Vector> rhs_;

  /// row Dirichlet map
  Epetra_Map rowdbcmap_;

  /// column Dirichlet map
  Epetra_Map coldbcmap_;

  /// transfer operator from column rhs to row rhs
  Epetra_Import importer_;

  /// support data for hanging nodes
  struct DofData
  {
    /// real node ids this node hangs on
    std::set<int> realdofmap;

    /// distribution factor
    double fact;
  };

  /// hanging node dof to real dof and factor
  std::map<int, DofData> hndofmap_;
};

}

#endif
#endif
