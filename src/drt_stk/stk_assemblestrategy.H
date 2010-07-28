#ifdef STKADAPTIVE

#ifndef STK_ASSEMBLESTRATEGY_H
#define STK_ASSEMBLESTRATEGY_H

#include <map>
#include <set>
#include <vector>

#include <Epetra_Map.h>

#include "../drt_lib/drt_assemblestrategy.H"

namespace DRT
{
  class Discretization;
}

namespace STK
{
  class Mesh;

  class AssembleStrategy : public DRT::AssembleStrategy
  {
  public:

    AssembleStrategy( DRT::Discretization & dis,
                      STK::Mesh & mesh,
                      Teuchos::RCP<const Epetra_Map> dbcmap,
                      Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
                      Teuchos::RCP<LINALG::SparseOperator> systemmatrix2,
                      Teuchos::RCP<Epetra_Vector>          systemvector1,
                      Teuchos::RCP<Epetra_Vector>          systemvector2,
                      Teuchos::RCP<Epetra_Vector>          systemvector3 );

    virtual Teuchos::RCP<Epetra_CrsGraph> MatrixGraph( DRT::Discretization & dis, Teuchos::RCP<const Epetra_Map> dbcmap );

    virtual void Assemble(LINALG::SparseOperator& sysmat,
                          int eid,
                          const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lm,
                          const std::vector<int>& lmowner);

    virtual void Assemble(LINALG::SparseOperator& sysmat,
                          int eid,
                          const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lmrow,
                          const std::vector<int>& lmrowowner,
                          const std::vector<int>& lmcol);

    virtual void Assemble(LINALG::SparseOperator& sysmat, double val, int rgid, int cgid);

    virtual void Assemble(Epetra_Vector& V, const Epetra_SerialDenseVector& Vele,
                          const std::vector<int>& lm, const std::vector<int>& lmowner);

    virtual void Assemble(Epetra_MultiVector& V, const int n, const Epetra_SerialDenseVector& Vele,
                          const std::vector<int>& lm, const std::vector<int>& lmowner);

  private:

    Teuchos::RCP<const Epetra_Map> dbcmap_;

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
