#ifdef STKADAPTIVE

#ifndef STK_FIXEDSPARSEMATRIX_H
#define STK_FIXEDSPARSEMATRIX_H

#include <map>
#include <set>
#include <vector>

#include "../linalg/linalg_fixedsparsematrix.H"

namespace DRT
{
  class Discretization;
}

namespace STK
{
  class Mesh;

  class FixedSparseMatrix : public LINALG::FixedSparseMatrix
  {
  public:
    FixedSparseMatrix( Teuchos::RCP<const Epetra_Map> dbcmap );

    virtual const char* Label() const { return "STK::FixedSparseMatrix"; }

#if 0
    virtual void Assemble(int eid,
                          const Epetra_SerialDenseMatrix& Aele,
                          const std::vector<int>& lmrow,
                          const std::vector<int>& lmrowowner,
                          const std::vector<int>& lmcol);
#endif

  private:

#if 0
    struct DofData
    {
      std::set<int> realdofmap;
      double fact;
    };

    /// hn dof to real dof and factor
    std::map<int, DofData> hndofmap_;
#endif
  };
}

#endif
#endif
