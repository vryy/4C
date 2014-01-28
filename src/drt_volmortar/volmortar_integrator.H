/*!----------------------------------------------------------------------
\file volmortar_integrator.H

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

/*---------------------------------------------------------------------*
 | definitions                                             farah 01/14 |
 *---------------------------------------------------------------------*/
#ifndef VOLMORTAR_INTEGRATOR_H_
#define VOLMORTAR_INTEGRATOR_H_

/*---------------------------------------------------------------------*
 | headers                                                 farah 01/14 |
 *---------------------------------------------------------------------*/
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

/*---------------------------------------------------------------------*
 | forward declarations                                    farah 01/14 |
 *---------------------------------------------------------------------*/
namespace LINALG
{
class SerialDenseVector;
class SparseMatrix;
}

namespace MORTAR
{
class IntCell;
}

namespace VOLMORTAR
{
class Cell;

template<DRT::Element::DiscretizationType distypeS,DRT::Element::DiscretizationType distypeM>
class VolMortarIntegrator
{
public:

  // constructor
  VolMortarIntegrator();

  // destructor
  virtual ~VolMortarIntegrator() {};

  //! ns_: number of slave element nodes
  static const int ns_ = DRT::UTILS::DisTypeToNumNodePerEle<distypeS>::numNodePerElement;

  //! nm_: number of master element nodes
  static const int nm_ = DRT::UTILS::DisTypeToNumNodePerEle<distypeM>::numNodePerElement;

  //! number of space dimensions ("+1" due to considering only interface elements)
  static const int ndim_ = DRT::UTILS::DisTypeToDim<distypeS>::dim;

  /*!
  \brief Initialize Gauss rule (points, weights)

  */
  void InitializeGP(bool integrateele=false);

  /*!
  \brief Integrate cell for 2D problems

  */
  void IntegrateCells2D(
      DRT::Element& sele,
      DRT::Element& mele,
      Teuchos::RCP<MORTAR::IntCell> cell,
      LINALG::SparseMatrix& dmatrix,
      LINALG::SparseMatrix& mmatrix,
      Teuchos::RCP<const DRT::Discretization> slavedis,
      Teuchos::RCP<const DRT::Discretization> masterdis);

  /*!
  \brief Integrate cell for 3D problems

  */
  void IntegrateCells3D(
       DRT::Element& sele,
       DRT::Element& mele,
       Teuchos::RCP<VOLMORTAR::Cell> cell,
       LINALG::SparseMatrix& dmatrix,
       LINALG::SparseMatrix& mmatrix,
       Teuchos::RCP<const DRT::Discretization> slavedis,
       Teuchos::RCP<const DRT::Discretization> masterdis);

  /*!
  \brief Integrate ele for 3D problems

  */
  void IntegrateSele3D(
       DRT::Element& sele,
       DRT::Element& mele,
       LINALG::SparseMatrix& dmatrix,
       LINALG::SparseMatrix& mmatrix,
       Teuchos::RCP<const DRT::Discretization> slavedis,
       Teuchos::RCP<const DRT::Discretization> masterdis);

  /*!
  \brief Integrate ele for 3D problems

  */
  void IntegrateMele3D(
       DRT::Element& sele,
       DRT::Element& mele,
       LINALG::SparseMatrix& dmatrix,
       LINALG::SparseMatrix& mmatrix,
       Teuchos::RCP<const DRT::Discretization> slavedis,
       Teuchos::RCP<const DRT::Discretization> masterdis);

  /*!
  \brief Initialize Gauss rule (points, weights)

  */
  bool CheckMapping2D(
      DRT::Element& sele,
      DRT::Element& mele,
      double* sxi,
      double* mxi);

  bool CheckMapping3D(
      DRT::Element& sele,
      DRT::Element& mele,
      double* sxi,
      double* mxi);

  //@}
private:
  int                           ngp_;         // number of Gauss points
  Epetra_SerialDenseMatrix      coords_;      // Gauss point coordinates
  std::vector<double>           weights_;     // Gauss point weights
};

}

#endif /* VOLMORTAR_INTEGRATOR_H_ */
