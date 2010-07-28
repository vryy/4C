#ifndef LINALG_FIXEDSPARSEMATRIX_H
#define LINALG_FIXEDSPARSEMATRIX_H

#include "linalg_sparsematrixbase.H"

#include <Epetra_CrsGraph.h>

namespace LINALG
{

  class FixedSparseMatrix : public SparseMatrixBase
  {
  public:

    FixedSparseMatrix( Teuchos::RCP<const Epetra_Map> dbcmap );

    virtual const char* Label() const { return "LINALG::FixedSparseMatrix"; }

    virtual void SetMatrix( Teuchos::RCP<Epetra_CrsGraph> graph );

    virtual void Zero();

    virtual void Reset();

    virtual void Assemble(int eid,
                          const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lmrow,
                          const std::vector<int>& lmrowowner,
                          const std::vector<int>& lmcol);

    virtual void Assemble(double val, int rgid, int cgid);

    virtual void Complete();

    virtual void Complete(const Epetra_Map& domainmap, const Epetra_Map& rangemap);

    virtual void UnComplete();

    virtual void ApplyDirichlet(const Teuchos::RCP<const Epetra_Vector> dbctoggle, bool diagonalblock=true);

    virtual void ApplyDirichlet(const Epetra_Map& dbcmap, bool diagonalblock=true);

  protected:

    Teuchos::RCP<const Epetra_Map> dbcmap_;
  };

}

#endif
