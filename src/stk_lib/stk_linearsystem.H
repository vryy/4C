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

 */
class LinearSystem
{
public:

  LinearSystem( STK::Discretization & dis );

  /// @name Assemble

  /// assemble element matrix
  void Assemble(int eid,
                const Epetra_SerialDenseMatrix& Aele,
                const std::vector<int>& lmrow,
                const std::vector<int>& lmrowowner,
                const std::vector<int>& lmcol);

  /// assemble rhs
  void Assemble(int eid,
                const Epetra_SerialDenseVector& Vele,
                const std::vector<int>& lm);

  /// complete any global objects
  void Complete();

  //@}

  Teuchos::RCP<Epetra_FECrsMatrix> Matrix() { return sysmat_; }

  Teuchos::RCP<Epetra_Vector> RHS() { return rowrhs_; }

  void Zero();

  void ApplyDirichlet( Epetra_Vector & dirichlet );

private:

  Teuchos::RCP<Epetra_CrsGraph> MatrixGraph( STK::Discretization & dis );

  Teuchos::RCP<Epetra_FECrsMatrix> sysmat_;
  Teuchos::RCP<Epetra_Vector> rowrhs_;
  Teuchos::RCP<Epetra_Vector> rhs_;

  /// Dirichlet map
  Epetra_Map rowdbcmap_;
  Epetra_Map coldbcmap_;

  Epetra_Import importer_;

  struct DofData
  {
    std::set<int> realdofmap;
    double fact;
  };

  /// hn dof to real dof and factor
  std::map<int, DofData> hndofmap_;
};

}

#endif
#endif
