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

/*---------------------------------------------------------------------*
 | headers                                                 farah 10/13 |
 *---------------------------------------------------------------------*/
#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include "../linalg/linalg_fixedsizematrix.H"

/*---------------------------------------------------------------------*
 | forward declarations                                    farah 10/13 |
 *---------------------------------------------------------------------*/
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
  /// Class for calling volmortar coupling and proper parallel redistr.
  class MortarVolCoupl
  {
  public:

    /*!
    \brief Empty constructor

    */
    MortarVolCoupl();

    /*!
    \brief Call parallel redistr. and evaluate volmortar coupl.

    */
    void Setup(Teuchos::RCP<DRT::Discretization> slavedis,
               Teuchos::RCP<DRT::Discretization> masterdis);

    /*!
    \brief Get coupling matrices for field A and B

    */
    Teuchos::RCP<const LINALG::SparseMatrix> GetPMatrixAB() const {return pmatrixA_;};
    Teuchos::RCP<const LINALG::SparseMatrix> GetPMatrixBA() const {return pmatrixB_;};

    /*!
    \brief Mortar mapping for A to B and B to A - for vectors

    */
    Teuchos::RCP<const Epetra_Vector> ApplyVectorMappingAB(Teuchos::RCP<const Epetra_Vector> vec) const ;
    Teuchos::RCP<const Epetra_Vector> ApplyVectorMappingBA(Teuchos::RCP<const Epetra_Vector> vec) const ;

    /*!
    \brief Mortar mapping for A to B and B to A - for matrices

    */
    Teuchos::RCP<LINALG::SparseMatrix> ApplyMatrixMappingAB(Teuchos::RCP<const LINALG::SparseMatrix> mat) const ;
    Teuchos::RCP<LINALG::SparseMatrix> ApplyMatrixMappingBA(Teuchos::RCP<const LINALG::SparseMatrix> mat) const ;

    //@}
  private:
    Teuchos::RCP<LINALG::SparseMatrix> pmatrixA_;             /// P matrix for field A
    Teuchos::RCP<LINALG::SparseMatrix> pmatrixB_;             /// P matrix for field B

  };
}

#endif /* ADAPTER_COUPLING_VOLMORTAR_H_ */
