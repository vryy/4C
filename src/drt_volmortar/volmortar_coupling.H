/*----------------------------------------------------------------------*/
/*! \file

\level 1


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
}  // namespace DRT

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
}  // namespace MORTAR

namespace VOLMORTAR
{
  class Cell;

  namespace UTILS
  {
    class DefaultMaterialStrategy;
  }

  /// Class for generating projection matrices for volumetric coupling
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
    VolMortarCoupl(int dim, Teuchos::RCP<DRT::Discretization> dis1,
        Teuchos::RCP<DRT::Discretization> dis2, std::vector<int>* coupleddof12 = NULL,
        std::vector<int>* coupleddof21 = NULL, std::pair<int, int>* dofsets12 = NULL,
        std::pair<int, int>* dofsets21 = NULL,
        Teuchos::RCP<VOLMORTAR::UTILS::DefaultMaterialStrategy> materialstrategy = Teuchos::null);

    /*!
     \brief Destructor

     */
    virtual ~VolMortarCoupl() {}

    /*!
     \brief Evaluate volmortar coupling (basic routine)

     */
    virtual void EvaluateVolmortar();

    /*!
     \brief Evaluate consistent interpolation (NO VOLMORTAR)

     */
    virtual void EvaluateConsistentInterpolation();

    /*!
     \brief get projection matrix 2 --> 1

     */
    Teuchos::RCP<LINALG::SparseMatrix> GetPMatrix12() { return P12_; };

    /*!
     \brief get projection matrix 1 --> 2

     */
    Teuchos::RCP<LINALG::SparseMatrix> GetPMatrix21() { return P21_; };

    /*!
     \brief assign materials

     */
    virtual void AssignMaterials();

   private:
    /*!
     \brief Assemble p matrix for cons. interpolation approach

     */
    virtual void AssembleConsistentInterpolation_P12(DRT::Node* node, std::vector<int>& foundeles);

    /*!
     \brief Assemble p matrix for cons. interpolation approach

     */
    virtual void AssembleConsistentInterpolation_P21(DRT::Node* node, std::vector<int>& foundeles);

    /*!
     \brief get auxiliary plane normal (2D)

     */
    virtual double* Auxn() { return auxn_; }

    /*!
     \brief Build maps based n coupling dofs

     */
    virtual void BuildMaps(Teuchos::RCP<DRT::Discretization>& dis,
        Teuchos::RCP<const Epetra_Map>& dofmap, const std::vector<int>* coupleddof,
        const int* nodes, int numnode, int dofset);

    /*!
     \brief calc dops for background mesh

     */
    virtual std::map<int, LINALG::Matrix<9, 2>> CalcBackgroundDops(
        Teuchos::RCP<DRT::Discretization> searchdis);

    /*!
     \brief calc dops for one element

     */
    virtual LINALG::Matrix<9, 2> CalcDop(DRT::Element& ele);

    /*!
     \brief center triangulation (if delaunay fails)

     */
    virtual bool CenterTriangulation(std::vector<Teuchos::RCP<MORTAR::IntCell>>& cells,
        std::vector<MORTAR::Vertex>& clip, double tol);

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
    virtual void CreateProjectionOperator();

    /*!
     \brief compute trafo operator

     */
    virtual void CreateTrafoOperator(DRT::Element& ele, Teuchos::RCP<DRT::Discretization> searchdis,
        bool dis, std::set<int>& donebefore);

    /*!
     \brief define vertices for 2D polygon clipping (master)

     */
    virtual void DefineVerticesMaster(
        DRT::Element& ele, std::vector<MORTAR::Vertex>& SlaveVertices);

    /*!
     \brief define vertices for 2D polygon clipping (slave)

     */
    virtual void DefineVerticesSlave(DRT::Element& ele, std::vector<MORTAR::Vertex>& SlaveVertices);

    /*!
     \brief create integration cells for 2D volmortar

     */
    virtual bool DelaunayTriangulation(std::vector<Teuchos::RCP<MORTAR::IntCell>>& cells,
        std::vector<MORTAR::Vertex>& clip, double tol);

    /*!
     \brief Get discretization of Omega_1

     */
    virtual Teuchos::RCP<const DRT::Discretization> Discret1() const { return dis1_; }

    /*!
     \brief Get discretization of Omega_2

     */
    virtual Teuchos::RCP<DRT::Discretization> Discret2() const { return dis2_; }

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
    std::vector<int> GetAdjacentNodes(DRT::Element::DiscretizationType shape, int& lid);

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
    virtual Teuchos::RCP<GEO::SearchTree> InitSearch(Teuchos::RCP<DRT::Discretization> searchdis);

    /*!
     \brief perform 2D integration

     */
    virtual void Integrate2D(
        DRT::Element& sele, DRT::Element& mele, std::vector<Teuchos::RCP<MORTAR::IntCell>>& cells);

    /*!
     \brief perform 3D element-wise integration

     */
    virtual void Integrate3D(DRT::Element& sele, DRT::Element& mele, int domain);

    /*!
     \brief perform 3D element-wise integration for P12

     */
    virtual void Integrate3DEleBased_P12(DRT::Element& sele, std::vector<int>& foundeles);

    /*!
     \brief perform 3D element-wise integration for BDis

     */
    virtual void Integrate3DEleBased_P21(DRT::Element& sele, std::vector<int>& foundeles);

    /*!
     \brief perform 3D element-wise integration for ADis for meshinit

     */
    virtual void Integrate3DEleBased_ADis_MeshInit(
        DRT::Element& Aele, std::vector<int>& foundeles, int dofseta, int dofsetb);

    /*!
     \brief perform 3D element-wise integration for BDis for meshinit

     */
    virtual void Integrate3DEleBased_BDis_MeshInit(
        DRT::Element& Bele, std::vector<int>& foundeles, int dofseta, int dofsetb);
    /*!
     \brief perform 3D integration of created cells

     */
    virtual void Integrate3DCell(
        DRT::Element& sele, DRT::Element& mele, std::vector<Teuchos::RCP<Cell>>& cells);

    /*!
     \brief perform 3D integration of created cells

     */
    virtual void Integrate3DCell_DirectDivergence(
        DRT::Element& sele, DRT::Element& mele, bool switched_conf = false);
    /*!
     \brief perform mesh init procedure

     */
    virtual void MeshInit();

    /*!
     \brief get parameter list

     */
    Teuchos::ParameterList& Params() { return params_; };

    /*!
     \brief perform cut and create integration cells (3D)

     */
    virtual void PerformCut(DRT::Element* sele, DRT::Element* mele, bool switched_conf = false);

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
    virtual std::vector<int> Search(DRT::Element& ele, Teuchos::RCP<GEO::SearchTree> SearchTree,
        std::map<int, LINALG::Matrix<9, 2>>& currentKDOPs);

    // don't want = operator and cctor
    VolMortarCoupl operator=(const VolMortarCoupl& old);
    VolMortarCoupl(const VolMortarCoupl& old);

    //! @name global problem information
    const int dim_;                  /// dimension of problem (2D or 3D)
    Teuchos::ParameterList params_;  /// global parameter list for volmortar coupling
    std::pair<int, int>
        dofset12_;  /// dofset number dofs of Omega_2 and Omega_1 in P Omega_2 -> Omega_1
    std::pair<int, int>
        dofset21_;  /// dofset number dofs of Omega_1 and Omega_2 in P Omega_1 -> Omega_2

    Teuchos::RCP<Epetra_Comm> comm_;  /// communicator
    int myrank_;                      /// my proc id

    //@}

    //! @name discretizations
    Teuchos::RCP<DRT::Discretization> dis1_;  /// the discretization Omega_1
    Teuchos::RCP<DRT::Discretization> dis2_;  /// the discretization Omega_2
    //@}

    //! @name mortar matrices and projector
    // s1 = D1^-1 * M12 * s2  = P12 * s2
    // s2 = D2^-1 * M21 * s1  = P21 * s1
    Teuchos::RCP<LINALG::SparseMatrix> D1_;   /// global Mortar matrix D1  for Omega_1
    Teuchos::RCP<LINALG::SparseMatrix> D2_;   /// global Mortar matrix D2  for Omega_2
    Teuchos::RCP<LINALG::SparseMatrix> M12_;  /// global Mortar matrix M12 for Omega_1
    Teuchos::RCP<LINALG::SparseMatrix> M21_;  /// global Mortar matrix M21 for Omega_2
    Teuchos::RCP<LINALG::SparseMatrix>
        P12_;  /// global Mortar projection matrix P Omega_2 -> Omega_1
    Teuchos::RCP<LINALG::SparseMatrix>
        P21_;  /// global Mortar projection matrix P Omega_1 -> Omega_2
    //@}

    //! @name trafo matrices for quadr. elements
    Teuchos::RCP<LINALG::SparseMatrix> T1_;  /// global trafo matrix for Omega_1
    Teuchos::RCP<LINALG::SparseMatrix> T2_;  /// global trafo matrix for Omega_2
    //@}

    //! @name maps
    Teuchos::RCP<const Epetra_Map>
        P12_dofrowmap_;  /// row map of projection matrix P Omega_2 -> Omega_1
    Teuchos::RCP<const Epetra_Map>
        P12_dofdomainmap_;  /// domain map of projection matrix P Omega_2 -> Omega_1
    Teuchos::RCP<const Epetra_Map>
        P21_dofrowmap_;  /// row map of projection matrix P Omega_1 -> Omega_2
    Teuchos::RCP<const Epetra_Map>
        P21_dofdomainmap_;  /// domain map of projection matrix P Omega_1 -> Omega_2
    Teuchos::RCP<const Epetra_Map>
        P12_dofcolmap_;  /// column map of projection matrix P Omega_2 -> Omega_1
    Teuchos::RCP<const Epetra_Map>
        P21_dofcolmap_;  /// column map of projection matrix P Omega_1 -> Omega_2
    //@}

    // quantity for 2D segmentation
    double auxn_[3];  /// normal of auxiliary plane (2D problems)

    //! @name counter and stat. information
    double volume_;       /// overall volume
    int polygoncounter_;  /// counter for created polygons/polyhedra
    int cellcounter_;     /// counter for created integration cells
    int inteles_;         /// counter for element-wise integration
    //@}

    // cut specific quantities
    GEO::CUT::plain_volumecell_set
        volcell_;  /// set of volume cells for direct divergence integration

    // search algorithm
    LINALG::Matrix<9, 3> dopnormals_;  /// dop normals for seach algorithm

    // input
    INPAR::VOLMORTAR::DualQuad dualquad_;  /// type of quadratic weighting interpolation

    /// strategy for element information transfer (mainly material, but can be more)
    Teuchos::RCP<VOLMORTAR::UTILS::DefaultMaterialStrategy> materialstrategy_;

    //! @name mesh initialization

    // maps for mesh init
    Teuchos::RCP<Epetra_Map> Xa_;
    Teuchos::RCP<Epetra_Map> Xb_;
    Teuchos::RCP<Epetra_Map> mergedmap_;

    // mortar matrices for mesh init
    Teuchos::RCP<LINALG::SparseMatrix> dmatrixXA_;  /// global Mortar matrix D for field A
    Teuchos::RCP<LINALG::SparseMatrix> dmatrixXB_;  /// global Mortar matrix D for field B
    Teuchos::RCP<LINALG::SparseMatrix> mmatrixXA_;  /// global Mortar matrix M for field A
    Teuchos::RCP<LINALG::SparseMatrix> mmatrixXB_;  /// global Mortar matrix M for field B

    //@}
  };

}  // namespace VOLMORTAR


#endif /* VOLMORTAR_COUPLING_H_ */
