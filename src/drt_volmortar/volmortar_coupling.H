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
#include "../drt_cut/cut_utils.H"
#include "../drt_inpar/inpar_volmortar.H"

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

namespace GEO
{
class SearchTree;
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

namespace UTILS
{
class DefaultMaterialStrategy;
}

/// Class for generating projection matrices for 3D mortar
/*!
 The idea is to glue two non-matching meshes together using the mortar method.
 In general, this works for displacement fields, as well as any other field
 (e.g. temperature field in tsi).
 The constructor expects the two discretizations, which are filled properly.
 Both discretization are expected to have at least two dofsets of which the first
 of one discretization is meant to be coupled with the second of the other
 discretization. I.e. in TSI the structure must have temperature dofs as second
 dof set and the thermo discretization must have displacement dofs as second dof set.
 When calling Evaluate() this class will identify volume cells (using polygon
 clipping in 2D and the cut algorithm in 3D) OR skip this and ignore weak discontinuities,
 and build a volmortar integrator class, which evaluates the two projection matrices.

 \author vuong 01/14
 */
class VolMortarCoupl
{

public:

  /*!
   \brief Constructor

   */
  VolMortarCoupl(int dim,
                 Teuchos::RCP<DRT::Discretization> slavedis,
                 Teuchos::RCP<DRT::Discretization> masterdis,
                 Teuchos::RCP<VOLMORTAR::UTILS::DefaultMaterialStrategy> materialstrategy);

  /*!
   \brief Destructor

   */
  virtual ~VolMortarCoupl()
  {
  }

  /*!
   \brief Evaluate volmortar coupling (basic routine)

   */
  virtual void EvaluateVolmortar();

  /*!
   \brief Evaluate consistent interpolation (NO VOLMORTAR)

   */
  virtual void EvaluateConsistentInterpolation();

  /*!
   \brief get projection matrix A

   */
  Teuchos::RCP<LINALG::SparseMatrix> GetPMatrixAB()
  {
    return pmatrixA_;
  }
  ;

  /*!
   \brief get projection matrix B

   */
  Teuchos::RCP<LINALG::SparseMatrix> GetPMatrixBA()
  {
    return pmatrixB_;
  }
  ;

protected:
  /*!
   \brief Get discretization of A-field

   */
  virtual Teuchos::RCP<const DRT::Discretization> ADiscret() const
  {
    return Adiscret_;
  }

  /*!
   \brief Assemble p matrix for cons. interpolation approach

   */
  virtual void AssembleConsistentInterpolation_ADis(DRT::Node* node);

  /*!
   \brief Assemble p matrix for cons. interpolation approach

   */
  virtual void AssembleConsistentInterpolation_BDis(DRT::Node* node);

  /*!
   \brief check initial coupling constraint

   */
  virtual void AssignMaterials();

  /*!
   \brief Get discretization of B-field

   */
  virtual Teuchos::RCP<DRT::Discretization> BDiscret() const
  {
    return Bdiscret_;
  }

  /*!
   \brief get auxiliary plane normal (2D)

   */
  virtual double* Auxn()
  {
    return auxn_;
  }

  /*!
   \brief calc dops for background mesh

   */
  virtual std::map<int, LINALG::Matrix<9, 2> > CalcBackgroundDops(
      Teuchos::RCP<DRT::Discretization> searchdis);

  /*!
   \brief calc dops for one element

   */
  virtual LINALG::Matrix<9, 2> CalcDop(DRT::Element& ele);

  /*!
   \brief check if we need cut (3D)

   */
  virtual bool CheckCut(DRT::Element& sele, DRT::Element& mele);

  /*!
   \brief check if we can integrate element-wise (3D)

   */
  virtual bool CheckEleIntegration(DRT::Element& sele, DRT::Element& mele);

  /*!
   \brief check initial coupling constraint

   */
  virtual void CheckInitialResiduum();

  /*!
   \brief complete created matrices

   */
  virtual void Complete();

  /*!
   \brief compute projection matrices D^-1 * M

   */
  virtual void CreateProjectionOpterator();

  /*!
   \brief compute trafo operator

   */
  virtual void CreateTrafoOperator(DRT::Element& ele,
      Teuchos::RCP<DRT::Discretization> searchdis, bool dis,
      std::set<int>& donebefore);

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
  virtual bool DelaunayTriangulation(
      std::vector<Teuchos::RCP<MORTAR::IntCell> >& cells,
      std::vector<MORTAR::Vertex> & clip, double tol);

  /*!
   \brief Evaluate element-based

   */
  virtual void EvaluateElements();

  /*!
   \brief Evaluate segment-based

   */
  virtual void EvaluateSegments();

  /*!
   \brief Evaluate segment-based for 2D problems

   */
  virtual void EvaluateSegments2D(DRT::Element& Aele, DRT::Element& Bele);

  /*!
   \brief Evaluate segment-based for 3D problems

   */
  virtual void EvaluateSegments3D(DRT::Element* Aele, DRT::Element* Bele);

  /*!
   \brief get adjacent node ids for quadr. dual shape functions (trafo calculation)

   */
  std::vector<int> GetAdjacentNodes(DRT::Element::DiscretizationType shape,
      int& lid);

  /*!
   \brief Initialize / reset volmortar coupling

   */
  virtual void Initialize();

  /*!
   \brief Initialize DOP normals for DOP calculation (Search algorithm)

   */
  virtual void InitDopNormals();

  /*!
   \brief Initialize search tree

   */
  virtual Teuchos::RCP<GEO::SearchTree> InitSearch(
      Teuchos::RCP<DRT::Discretization> searchdis);

  /*!
   \brief perform 2D integration

   */
  virtual void Integrate2D(DRT::Element& sele, DRT::Element& mele,
      std::vector<Teuchos::RCP<MORTAR::IntCell> >& cells);

  /*!
   \brief perform 3D element-wise integration

   */
  virtual void Integrate3D(DRT::Element& sele, DRT::Element& mele, int domain);

  /*!
   \brief perform 3D element-wise integration for ADis

   */
  virtual void Integrate3DEleBased_ADis(DRT::Element& sele,
      std::vector<int>& foundeles);

  /*!
   \brief perform 3D element-wise integration for BDis

   */
  virtual void Integrate3DEleBased_BDis(DRT::Element& sele,
      std::vector<int>& foundeles);

  /*!
   \brief perform 3D integration of created cells

   */
  virtual void Integrate3DCell(DRT::Element& sele, DRT::Element& mele,
      std::vector<Teuchos::RCP<Cell> >& cells);

  /*!
   \brief perform 3D integration of created cells

   */
  virtual void Integrate3DCell_DirectDivergence(DRT::Element& sele,
      DRT::Element& mele, bool switched_conf = false);
  /*!
   \brief perform mesh init procedure

   */
  virtual void MeshInit();

  /*!
   \brief get parameter list

   */
  Teuchos::ParameterList& Params()
  {
    return params_;
  }
  ;

  /*!
   \brief perform cut and create integration cells (3D)

   */
  virtual void PerformCut(DRT::Element* sele, DRT::Element* mele,
      bool switched_conf = false);

  /*!
   \brief perform 2D polygon clipping

   */
  virtual bool PolygonClippingConvexHull(std::vector<MORTAR::Vertex>& poly1,
      std::vector<MORTAR::Vertex>& poly2, std::vector<MORTAR::Vertex>& respoly,
      DRT::Element& sele, DRT::Element& mele, double& tol);

  /*!
   \brief Output for evaluation status -- progress

   */
  virtual void PrintStatus(int& i, bool dis_switch = false);

  /*!
   \brief Get required parameters and check for validity

   */
  virtual void ReadAndCheckInput();

  /*!
   \brief search algorithm

   */
  virtual std::vector<int> Search(DRT::Element& ele,
      Teuchos::RCP<GEO::SearchTree> SearchTree,
      std::map<int, LINALG::Matrix<9, 2> >& currentKDOPs);

  //@}

  // don't want = operator and cctor
  VolMortarCoupl operator =(const VolMortarCoupl& old);
  VolMortarCoupl(const VolMortarCoupl& old);

  // global problem information
  const int                          dim_;             /// dimension of problem (2D or 3D)
  Teuchos::ParameterList             params_;          /// global parameter list for volmortar coupling
  // discretizations
  Teuchos::RCP<DRT::Discretization>  Adiscret_;        /// the discretization of field A
  Teuchos::RCP<DRT::Discretization>  Bdiscret_;        /// the discretization of field B

  // mortar matrices and projector
  Teuchos::RCP<LINALG::SparseMatrix> dmatrixA_;        /// global Mortar matrix D for field A
  Teuchos::RCP<LINALG::SparseMatrix> dmatrixB_;        /// global Mortar matrix D for field B
  Teuchos::RCP<LINALG::SparseMatrix> mmatrixA_;        /// global Mortar matrix M for field A
  Teuchos::RCP<LINALG::SparseMatrix> mmatrixB_;        /// global Mortar matrix M for field B
  Teuchos::RCP<LINALG::SparseMatrix> pmatrixA_;        /// global Mortar projection matrix P for field A
  Teuchos::RCP<LINALG::SparseMatrix> pmatrixB_;        /// global Mortar projection matrix P for field B

  // trafo matrix for quadr. elements
  Teuchos::RCP<LINALG::SparseMatrix> TA_;              /// global trafo matrix for A dis
  Teuchos::RCP<LINALG::SparseMatrix> TB_;              /// global trafo matrix for B dis

  // quantity for 2D segmentation
  double                             auxn_[3];         /// normal of auxiliary plane (2D problems)

  // counter and stat. information
  double                             volume_;          /// overall volume
  int                                polygoncounter_;  /// counter for created polygons/polyhedra
  int                                cellcounter_;     /// counter for created integration cells
  int                                inteles_;         /// counter for element-wise integration

  // communicator
  Teuchos::RCP<Epetra_Comm>          comm_;            /// communicator
  int                                myrank_;          /// my proc id

  //cut specific quantities
  GEO::CUT::plain_volumecell_set     volcell_;         /// set of volume cells for direct divergence integration

  // search algorithm
  LINALG::Matrix<9,3>                dopnormals_;      /// dop normals for seach algorithm

  //input
  INPAR::VOLMORTAR::DualQuad         dualquad_;        /// type of quadratic weighting interpolation

  Teuchos::RCP<VOLMORTAR::UTILS::DefaultMaterialStrategy> materialstrategy_;
};

}//namespace...


#endif /* VOLMORTAR_COUPLING_H_ */
