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

/*---------------------------------------------------------------------*
 | headers                                                 farah 01/14 |
 *---------------------------------------------------------------------*/
#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include "../drt_mortar/mortar_coupling3d_classes.H"

/*---------------------------------------------------------------------*
 | forward declarations                                    farah 01/14 |
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

/// Class for generating projection matrices for 3D mortar
/*!
  The idea is to clue to non-matching meshes together using the mortar method.
  In general, this works for displacement fields, as well as any other field
  (e.g. temperature field in tsi). 
  The constructor expects the two discretizations, which are filled properly. 
  Both discretization are expected to have at least two dofsets of which the first
  of one discretization is meant to be coupled with the second of the other 
  discretization. I.e. in TSI the structure must have temperature dofs as second 
  dof set and the thermo discretization must have displacement dofs as second dof set.
  When calling Evaluate() this class will identify volume cells (using polygon 
  clipping in 2D and the cut algorthm in 3D), and build a volmortar integrator class, 
  which evaluates the two projection matrixes. 

  \author vuong 01/14
*/
class VolMortarCoupl
{
public:

  /*!
  \brief Constructor

  */
  VolMortarCoupl(int dim,
                 const Epetra_Comm& comm,
                 Teuchos::RCP<DRT::Discretization> slavedis,
                 Teuchos::RCP<DRT::Discretization> masterdis);

  /*!
  \brief Destructor

  */
  virtual ~VolMortarCoupl() {}

  /*!
  \brief Get discretization of A-field

  */
  virtual Teuchos::RCP<const DRT::Discretization> ADiscret() const { return Adiscret_; }

  /*!
  \brief Get discretization of B-field

  */
  virtual Teuchos::RCP<const DRT::Discretization> BDiscret() const { return Bdiscret_; }

  /*!
  \brief get auxiliary plane normal (2D)

  */
  virtual double* Auxn() { return auxn_; }

  /*!
  \brief check if we need cut (3D)

  */
  virtual bool CheckCut(DRT::Element& sele,
                        DRT::Element& mele);

  /*!
  \brief check if we can integrate element-wise (3D)

  */
  virtual bool CheckEleIntegration(DRT::Element& sele,
                                   DRT::Element& mele);

  /*!
  \brief complete created matrices

  */
  virtual void Complete();

  /*!
  \brief compute projection matrices D^-1 * M

  */
  virtual void CreateProjectionOpterator();

  /*!
  \brief define vertices for 2D polygon clipping (master)

  */
  virtual void DefineVerticesMaster(DRT::Element& ele,
                                    std::vector<MORTAR::Vertex>& SlaveVertices);

  /*!
  \brief define vertices for 2D polygon clipping (slave)

  */
  virtual void DefineVerticesSlave(DRT::Element& ele,
                                   std::vector<MORTAR::Vertex>& SlaveVertices);

  /*!
  \brief create integration cells for 2D volmortar

  */
  virtual bool DelaunayTriangulation(std::vector<Teuchos::RCP<MORTAR::IntCell> >&  cells,
                                     std::vector<MORTAR::Vertex> & clip,
                                     double tol);

  /*!
  \brief Evaluate volmortar coupling (basic routine)

  */
  virtual void Evaluate();

  /*!
  \brief get projection matrix A

  */
  Teuchos::RCP<LINALG::SparseMatrix> GetPMatrixAB() {return pmatrixA_;};

  /*!
  \brief get projection matrix B

  */
  Teuchos::RCP<LINALG::SparseMatrix> GetPMatrixBA() {return pmatrixB_;};

  /*!
  \brief Initialize / reset volmortar coupling

  */
  virtual void Initialize();

  /*!
  \brief perform 2D integration

  */
  virtual void Integrate2D(DRT::Element& sele,
                           DRT::Element& mele,
                           std::vector<Teuchos::RCP<MORTAR::IntCell> >& cells);

  /*!
  \brief perform 3D element-wise integration

  */
  virtual void Integrate3D(DRT::Element& sele,
                           DRT::Element& mele,
                           int domain);

  /*!
  \brief perform 3D integration of created cells

  */
  virtual void Integrate3DCell(DRT::Element& sele,
                               DRT::Element& mele,
                               std::vector<Teuchos::RCP<Cell> >& cells);


  /*!
  \brief perform cut and create integration cells (3D)

  */
  virtual void PerformCut(DRT::Element* sele,
                          DRT::Element* mele,
                          std::vector<Teuchos::RCP<Cell> >& cells);

  /*!
  \brief perform 2D polygon clipping

  */
  virtual bool PolygonClippingConvexHull(std::vector<MORTAR::Vertex>& poly1,
                                                      std::vector<MORTAR::Vertex>& poly2,
                                                      std::vector<MORTAR::Vertex>& respoly,
                                                      DRT::Element& sele,
                                                      DRT::Element& mele,
                                                      double& tol);

  /*!
  \brief Output for evaluation status -- progress

  */
  virtual void PrintStatus(int& i);

  //@}
protected:
  // don't want = operator and cctor
  VolMortarCoupl operator = (const VolMortarCoupl& old);
  VolMortarCoupl(const VolMortarCoupl& old);

  const int                                dim_;       /// dimension of problem (2D or 3D)

  Teuchos::RCP<const DRT::Discretization>  Adiscret_;  /// the discretization of field A
  Teuchos::RCP<const DRT::Discretization>  Bdiscret_;  /// the discretization of field B

  Teuchos::RCP<LINALG::SparseMatrix> dmatrixA_;        /// global Mortar matrix D for field A
  Teuchos::RCP<LINALG::SparseMatrix> dmatrixB_;        /// global Mortar matrix D for field B

  Teuchos::RCP<LINALG::SparseMatrix> mmatrixA_;        /// global Mortar matrix M for field A
  Teuchos::RCP<LINALG::SparseMatrix> mmatrixB_;        /// global Mortar matrix M for field B

  Teuchos::RCP<LINALG::SparseMatrix> pmatrixA_;        /// global Mortar projection matrix P for field A
  Teuchos::RCP<LINALG::SparseMatrix> pmatrixB_;        /// global Mortar projection matrix P for field B


  double                             auxn_[3];         /// normal of auxiliary plane (2D problems)
  int                                polygoncounter_;  /// counter for created polygons/polyhedra
  int                                cellcounter_;     /// counter for created integration cells
  int                                inteles_;         /// counter for element-wise integration

  Teuchos::RCP<Epetra_Comm>          commA_;           /// communicator for field A
  Teuchos::RCP<Epetra_Comm>          commB_;           /// communicator for field B
};

}//namespace...

// << operator
std::ostream& operator << (std::ostream& os, const VOLMORTAR::VolMortarCoupl& coupl);

#endif /* VOLMORTAR_COUPLING_H_ */
