/*!----------------------------------------------------------------------
\file adapter_coupling_volmortar.H

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifndef ADAPTER_COUPLING_VOLMORTAR_H_
#define ADAPTER_COUPLING_VOLMORTAR_H_

#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>

namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class SparseMatrix;
}

namespace ADAPTER
{
  class MortarVolCoupl
  {
  public:

    /// empty constructor
    MortarVolCoupl();

    //setup...
    void Setup(Teuchos::RCP<DRT::Discretization> slavedis,
               Teuchos::RCP<DRT::Discretization> masterdis,
               const Epetra_Comm&                comm);

    Teuchos::RCP<const LINALG::SparseMatrix> GetPMatrixAB() const {return pmatrixA_;};
    Teuchos::RCP<const LINALG::SparseMatrix> GetPMatrixBA() const {return pmatrixB_;};

    Teuchos::RCP<const Epetra_Vector> ApplyVectorMappingAB(Teuchos::RCP<const Epetra_Vector> vec) const ;
    Teuchos::RCP<const Epetra_Vector> ApplyVectorMappingBA(Teuchos::RCP<const Epetra_Vector> vec) const ;

    Teuchos::RCP<LINALG::SparseMatrix> ApplyMatrixMappingAB(Teuchos::RCP<const LINALG::SparseMatrix> mat) const ;
    Teuchos::RCP<LINALG::SparseMatrix> ApplyMatrixMappingBA(Teuchos::RCP<const LINALG::SparseMatrix> mat) const ;

  private:
    Teuchos::RCP<LINALG::SparseMatrix> pmatrixA_;
    Teuchos::RCP<LINALG::SparseMatrix> pmatrixB_;

  };
}

#endif /* ADAPTER_COUPLING_VOLMORTAR_H_ */
