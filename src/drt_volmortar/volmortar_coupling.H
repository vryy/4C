/*!----------------------------------------------------------------------
\file volmortar_coupling.H

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef VOLMORTAR_COUPLING_H_
#define VOLMORTAR_COUPLING_H_

#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include "../drt_mortar/mortar_coupling3d_classes.H"

// forward declarations
namespace DRT
{
  class Discretization;
  class Element;
}

namespace LINALG
{
  class SparseMatrix;
}

namespace MORTAR
{
  class IntCell;
  class Vertex;
}
/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms

*/
namespace VOLMORTAR
{
class Cell;

class VolMortarCoupl
{
public:

  /*!
  \brief Constructor

  */
  VolMortarCoupl(const int dim,const Epetra_Comm& comm,
                 Teuchos::RCP<const DRT::Discretization> slavedis,
                 Teuchos::RCP<const DRT::Discretization> master);

  /*!
  \brief Destructor

  */
  virtual ~VolMortarCoupl() {}

  /*!
  \brief Get discretization of slave

  */
  virtual Teuchos::RCP<const DRT::Discretization> ADiscret() const { return Adiscret_; }

  /*!
  \brief Get discretization of master

  */
  virtual Teuchos::RCP<const DRT::Discretization> BDiscret() const { return Bdiscret_; }

  /*!
  \brief Initialize / reset volmortar coupling

  */
  virtual void Initialize();

  /*!
  \brief Evaluate volmortar coupling

  */
  virtual void Evaluate();

  virtual void DefineVerticesSlave(DRT::Element& ele,
                                   std::vector<MORTAR::Vertex>& SlaveVertices);

  virtual void DefineVerticesMaster(DRT::Element& ele,
                                    std::vector<MORTAR::Vertex>& SlaveVertices);

  virtual void Complete();

  virtual void CreateProjectionOpterator();


  virtual bool PolygonClippingConvexHull(std::vector<MORTAR::Vertex>& poly1,
                                                      std::vector<MORTAR::Vertex>& poly2,
                                                      std::vector<MORTAR::Vertex>& respoly,
                                                      DRT::Element& sele,
                                                      DRT::Element& mele,
                                                      double& tol);

  virtual bool DelaunayTriangulation(std::vector<Teuchos::RCP<MORTAR::IntCell> >&  cells,
                                     std::vector<MORTAR::Vertex> & clip,
                                     double tol);

  virtual void Integrate2D(DRT::Element& sele,
                           DRT::Element& mele,
                           std::vector<Teuchos::RCP<MORTAR::IntCell> >& cells);

  virtual void Integrate3DCell(DRT::Element& sele,
                               DRT::Element& mele,
                               std::vector<Teuchos::RCP<Cell> >& cells);

  virtual void Integrate3DSele(DRT::Element& sele,
                              DRT::Element& mele);

  virtual void Integrate3DMele(DRT::Element& sele,
                              DRT::Element& mele);

  virtual bool CheckCut(DRT::Element& sele,
                        DRT::Element& mele);

  virtual bool CheckEleIntegration(DRT::Element& sele,
                                   DRT::Element& mele);

  virtual void PerformCut(DRT::Element* sele,
                          DRT::Element* mele,
                          std::vector<Teuchos::RCP<Cell> >& cells);

  virtual double* Auxn() { return auxn_; }

  Teuchos::RCP<LINALG::SparseMatrix> GetPMatrixAB() {return pmatrixA_;};
  Teuchos::RCP<LINALG::SparseMatrix> GetPMatrixBA() {return pmatrixB_;};

protected:
  // don't want = operator and cctor
  VolMortarCoupl operator = (const VolMortarCoupl& old);
  VolMortarCoupl(const VolMortarCoupl& old);

  const int                                dim_;      // dimension of problem (2D or 3D)

  Teuchos::RCP<const DRT::Discretization>  Adiscret_; // the discretization of field A
  Teuchos::RCP<const DRT::Discretization>  Bdiscret_; // the discretization of field B

  Teuchos::RCP<LINALG::SparseMatrix> dmatrixA_;       // global Mortar matrix D for field A
  Teuchos::RCP<LINALG::SparseMatrix> dmatrixB_;       // global Mortar matrix D for field B

  Teuchos::RCP<LINALG::SparseMatrix> mmatrixA_;       // global Mortar matrix M for field A
  Teuchos::RCP<LINALG::SparseMatrix> mmatrixB_;       // global Mortar matrix M for field B

  Teuchos::RCP<LINALG::SparseMatrix> pmatrixA_;       // global Mortar projection matrix P for field A
  Teuchos::RCP<LINALG::SparseMatrix> pmatrixB_;       // global Mortar projection matrix P for field B


  double                             auxn_[3];        // normal of auxiliary plane
  int                                polygoncounter_;
  int                                cellcounter_;

  Teuchos::RCP<Epetra_Comm>          comm_;
};

}//namespace...

// << operator
std::ostream& operator << (std::ostream& os, const VOLMORTAR::VolMortarCoupl& coupl);

#endif /* VOLMORTAR_COUPLING_H_ */
