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
  class Element;
}

namespace LINALG
{
  class SparseMatrix;
}

namespace VOLMORTAR
{
  namespace UTILS
  {
    class DefaultMaterialStrategy;
  }
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
    void Setup(Teuchos::RCP<DRT::Discretization> dis1,
               Teuchos::RCP<DRT::Discretization> dis2,
               std::vector<int>* coupleddof12,
               std::vector<int>* coupleddof21,
               Teuchos::RCP<VOLMORTAR::UTILS::DefaultMaterialStrategy> materialstrategy = Teuchos::null
               );

    /*!
    \brief Get coupling matrices for field 1 and 2

    */
    Teuchos::RCP<const LINALG::SparseMatrix> GetPMatrix12() const {return P12_;};
    Teuchos::RCP<const LINALG::SparseMatrix> GetPMatrix21() const {return P21_;};

    /*!
    \brief Mortar mapping for 1 to 2 and 2 to 1 - for vectors

    */
    Teuchos::RCP<const Epetra_Vector> ApplyVectorMapping12(Teuchos::RCP<const Epetra_Vector> vec) const ;
    Teuchos::RCP<const Epetra_Vector> ApplyVectorMapping21(Teuchos::RCP<const Epetra_Vector> vec) const ;

    /*!
    \brief Mortar mapping for 1 to 2 and 2 to 1 - for matrices

    */
    Teuchos::RCP<LINALG::SparseMatrix> ApplyMatrixMapping12(Teuchos::RCP<const LINALG::SparseMatrix> mat) const ;
    Teuchos::RCP<LINALG::SparseMatrix> ApplyMatrixMapping21(Teuchos::RCP<const LINALG::SparseMatrix> mat) const ;

    //@}

    //! assign materials
    void AssignMaterials(
        Teuchos::RCP<DRT::Discretization> dis1,
        Teuchos::RCP<DRT::Discretization> dis2,
        Teuchos::RCP<VOLMORTAR::UTILS::DefaultMaterialStrategy> materialstrategy = Teuchos::null);

  private:
    // mortar matrices and projector
    // s1 = P12 * s2
    // s2 = P21 * s1
    Teuchos::RCP<LINALG::SparseMatrix> P12_;     // global Mortar projection matrix P Omega_2 -> Omega_1
    Teuchos::RCP<LINALG::SparseMatrix> P21_;     // global Mortar projection matrix P Omega_1 -> Omega_2

  };
}

#endif /* ADAPTER_COUPLING_VOLMORTAR_H_ */
