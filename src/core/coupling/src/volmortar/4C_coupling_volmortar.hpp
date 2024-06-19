/*----------------------------------------------------------------------*/
/*! \file

\level 1


 *----------------------------------------------------------------------*/

#ifndef FOUR_C_COUPLING_VOLMORTAR_HPP
#define FOUR_C_COUPLING_VOLMORTAR_HPP

/*---------------------------------------------------------------------*
 | headers                                                 farah 01/14 |
 *---------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_cut_utils.hpp"
#include "4C_mortar_coupling3d_classes.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------*
 | forward declarations                                    farah 01/14 |
 *---------------------------------------------------------------------*/
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE
namespace Core::Elements
{
  class Element;
}

namespace Core::LinAlg
{
  class SparseMatrix;
}

namespace Core::Geo
{
  class SearchTree;
}

namespace Mortar
{
  class IntCell;
  class Vertex;
}  // namespace Mortar

namespace Core::VolMortar
{
  class Cell;

  // Type of integration procedure
  enum IntType
  {
    inttype_segments,  ///< cut procedure of volume meshes
    inttype_elements   ///< fast, elementwise integration
  };

  // Type of cut procedure
  enum CutType
  {
    cuttype_directdivergence,  ///< direct divergence for integration
    cuttype_tessellation       ///< tessellation of volume meshes
  };

  // Type weighting function for quadr. problems
  enum DualQuad
  {
    dualquad_no_mod,   ///< no modification
    dualquad_lin_mod,  ///< linear modification
    dualquad_quad_mod  ///< quadr. modification
  };

  // Type of coupling
  enum CouplingType
  {
    couplingtype_volmortar,  ///< volmortar
    couplingtype_coninter    ///< consist. interpolation
  };

  // Type of coupling
  enum Shapefcn
  {
    shape_dual,  ///< dual shape functions
    shape_std    ///< std. shape functions --> lumped
  };

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
   When calling evaluate() this class will identify volume cells (using polygon
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
    VolMortarCoupl(int dim, Teuchos::RCP<Core::FE::Discretization> dis1,
        Teuchos::RCP<Core::FE::Discretization> dis2,
        const Teuchos::ParameterList& volmortar_parameters,
        std::vector<int>* coupleddof12 = nullptr, std::vector<int>* coupleddof21 = nullptr,
        std::pair<int, int>* dofset12 = nullptr, std::pair<int, int>* dofset21 = nullptr,
        Teuchos::RCP<Core::VolMortar::UTILS::DefaultMaterialStrategy> materialstrategy =
            Teuchos::null);

    /*!
     \brief Destructor

     */
    virtual ~VolMortarCoupl() = default;
    /*!
     \brief Evaluate volmortar coupling (basic routine)

     */
    virtual void EvaluateVolmortar();

    /*!
     \brief Evaluate consistent interpolation (NO Core::VOLMORTAR)

     */
    virtual void evaluate_consistent_interpolation();

    /*!
     \brief get projection matrix 2 --> 1

     */
    Teuchos::RCP<Core::LinAlg::SparseMatrix> GetPMatrix12() { return p12_; };

    /*!
     \brief get projection matrix 1 --> 2

     */
    Teuchos::RCP<Core::LinAlg::SparseMatrix> GetPMatrix21() { return p21_; };

    /*!
     \brief assign materials

     */
    virtual void AssignMaterials();

   private:
    /*!
     \brief Assemble p matrix for cons. interpolation approach

     */
    virtual void assemble_consistent_interpolation_p12(
        Core::Nodes::Node* node, std::vector<int>& foundeles);

    /*!
     \brief Assemble p matrix for cons. interpolation approach

     */
    virtual void assemble_consistent_interpolation_p21(
        Core::Nodes::Node* node, std::vector<int>& foundeles);

    /*!
     \brief get auxiliary plane normal (2D)

     */
    virtual double* auxn() { return auxn_; }

    /*!
     \brief Build maps based n coupling dofs

     */
    virtual void build_maps(Teuchos::RCP<Core::FE::Discretization>& dis,
        Teuchos::RCP<const Epetra_Map>& dofmap, const std::vector<int>* coupleddof,
        const int* nodes, int numnode, int dofset);

    /*!
     \brief calc dops for background mesh

     */
    virtual std::map<int, Core::LinAlg::Matrix<9, 2>> calc_background_dops(
        Teuchos::RCP<Core::FE::Discretization> searchdis);

    /*!
     \brief calc dops for one element

     */
    virtual Core::LinAlg::Matrix<9, 2> calc_dop(Core::Elements::Element& ele);

    /*!
     \brief center triangulation (if delaunay fails)

     */
    virtual bool center_triangulation(std::vector<Teuchos::RCP<Mortar::IntCell>>& cells,
        std::vector<Mortar::Vertex>& clip, double tol);

    /*!
     \brief check if we need cut (3D)

     */
    virtual bool check_cut(Core::Elements::Element& sele, Core::Elements::Element& mele);

    /*!
     \brief check if we can integrate element-wise (3D)

     */
    virtual bool check_ele_integration(
        Core::Elements::Element& sele, Core::Elements::Element& mele);

    /*!
     \brief check initial coupling constraint

     */
    virtual void check_initial_residuum();

    /*!
     \brief complete created matrices

     */
    virtual void complete();

    /*!
     \brief compute projection matrices D^-1 * M

     */
    virtual void create_projection_operator();

    /*!
     \brief compute trafo operator

     */
    virtual void create_trafo_operator(Core::Elements::Element& ele,
        Teuchos::RCP<Core::FE::Discretization> searchdis, bool dis, std::set<int>& donebefore);

    /*!
     \brief define vertices for 2D polygon clipping (master)

     */
    virtual void define_vertices_master(
        Core::Elements::Element& ele, std::vector<Mortar::Vertex>& slave_vertices);

    /*!
     \brief define vertices for 2D polygon clipping (slave)

     */
    virtual void define_vertices_slave(
        Core::Elements::Element& ele, std::vector<Mortar::Vertex>& slave_vertices);

    /*!
     \brief create integration cells for 2D volmortar

     */
    virtual bool delaunay_triangulation(std::vector<Teuchos::RCP<Mortar::IntCell>>& cells,
        std::vector<Mortar::Vertex>& clip, double tol);

    /*!
     \brief Get discretization of Omega_1

     */
    virtual Teuchos::RCP<const Core::FE::Discretization> discret1() const { return dis1_; }

    /*!
     \brief Get discretization of Omega_2

     */
    virtual Teuchos::RCP<Core::FE::Discretization> discret2() const { return dis2_; }

    /*!
     \brief Evaluate element-based

     */
    virtual void evaluate_elements();

    /*!
     \brief Evaluate segment-based

     */
    virtual void evaluate_segments();

    /*!
     \brief Evaluate segment-based for 2D problems

     */
    virtual void evaluate_segments2_d(Core::Elements::Element& Aele, Core::Elements::Element& Bele);

    /*!
     \brief Evaluate segment-based for 3D problems

     */
    virtual void evaluate_segments3_d(Core::Elements::Element* Aele, Core::Elements::Element* Bele);

    /*!
     \brief get adjacent node ids for quadr. dual shape functions (trafo calculation)

     */
    std::vector<int> get_adjacent_nodes(Core::FE::CellType shape, int& lid);

    /*!
     \brief Initialize / reset volmortar coupling

     */
    virtual void initialize();

    /*!
     \brief Initialize DOP normals for DOP calculation (Search algorithm)

     */
    virtual void init_dop_normals();

    /*!
     \brief Initialize search tree

     */
    virtual Teuchos::RCP<Core::Geo::SearchTree> init_search(
        Teuchos::RCP<Core::FE::Discretization> searchdis);

    /*!
     \brief perform 2D integration

     */
    virtual void integrate2_d(Core::Elements::Element& sele, Core::Elements::Element& mele,
        std::vector<Teuchos::RCP<Mortar::IntCell>>& cells);

    /*!
     \brief perform 3D element-wise integration

     */
    virtual void integrate3_d(
        Core::Elements::Element& sele, Core::Elements::Element& mele, int domain);

    /*!
     \brief perform 3D element-wise integration for P12

     */
    virtual void integrate3_d_ele_based_p12(
        Core::Elements::Element& Aele, std::vector<int>& foundeles);

    /*!
     \brief perform 3D element-wise integration for BDis

     */
    virtual void integrate3_d_ele_based_p21(
        Core::Elements::Element& Bele, std::vector<int>& foundeles);

    /*!
     \brief perform 3D element-wise integration for ADis for meshinit

     */
    virtual void integrate3_d_ele_based_a_dis_mesh_init(
        Core::Elements::Element& Aele, std::vector<int>& foundeles, int dofseta, int dofsetb);

    /*!
     \brief perform 3D element-wise integration for BDis for meshinit

     */
    virtual void integrate3_d_ele_based_b_dis_mesh_init(
        Core::Elements::Element& Bele, std::vector<int>& foundeles, int dofsetb, int dofseta);
    /*!
     \brief perform 3D integration of created cells

     */
    virtual void integrate3_d_cell(Core::Elements::Element& sele, Core::Elements::Element& mele,
        std::vector<Teuchos::RCP<Cell>>& cells);

    /*!
     \brief perform 3D integration of created cells

     */
    virtual void integrate3_d_cell_direct_divergence(
        Core::Elements::Element& sele, Core::Elements::Element& mele, bool switched_conf = false);
    /*!
     \brief perform mesh init procedure

     */
    virtual void mesh_init();

    /*!
     \brief get parameter list

     */
    Teuchos::ParameterList& params() { return params_; };

    /*!
     \brief perform cut and create integration cells (3D)

     */
    virtual void perform_cut(
        Core::Elements::Element* sele, Core::Elements::Element* mele, bool switched_conf = false);

    /*!
     \brief perform 2D polygon clipping

     */
    virtual bool polygon_clipping_convex_hull(std::vector<Mortar::Vertex>& poly1,
        std::vector<Mortar::Vertex>& poly2, std::vector<Mortar::Vertex>& respoly,
        Core::Elements::Element& sele, Core::Elements::Element& mele, double& tol);

    /*!
     \brief Output for evaluation status -- progress

     */
    virtual void print_status(int& i, bool dis_switch = false);

    /*!
     \brief Get required parameters and check for validity

     */
    virtual void read_and_check_input(const Teuchos::ParameterList& volmortar_parameters);

    /*!
     \brief search algorithm

     */
    virtual std::vector<int> search(Core::Elements::Element& ele,
        Teuchos::RCP<Core::Geo::SearchTree> SearchTree,
        std::map<int, Core::LinAlg::Matrix<9, 2>>& currentKDOPs);

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
    Teuchos::RCP<Core::FE::Discretization> dis1_;  /// the discretization Omega_1
    Teuchos::RCP<Core::FE::Discretization> dis2_;  /// the discretization Omega_2
    //@}

    //! @name mortar matrices and projector
    // s1 = D1^-1 * M12 * s2  = P12 * s2
    // s2 = D2^-1 * M21 * s1  = P21 * s1
    Teuchos::RCP<Core::LinAlg::SparseMatrix> d1_;   /// global Mortar matrix D1  for Omega_1
    Teuchos::RCP<Core::LinAlg::SparseMatrix> d2_;   /// global Mortar matrix D2  for Omega_2
    Teuchos::RCP<Core::LinAlg::SparseMatrix> m12_;  /// global Mortar matrix M12 for Omega_1
    Teuchos::RCP<Core::LinAlg::SparseMatrix> m21_;  /// global Mortar matrix M21 for Omega_2
    Teuchos::RCP<Core::LinAlg::SparseMatrix>
        p12_;  /// global Mortar projection matrix P Omega_2 -> Omega_1
    Teuchos::RCP<Core::LinAlg::SparseMatrix>
        p21_;  /// global Mortar projection matrix P Omega_1 -> Omega_2
    //@}

    //! @name trafo matrices for quadr. elements
    Teuchos::RCP<Core::LinAlg::SparseMatrix> t1_;  /// global trafo matrix for Omega_1
    Teuchos::RCP<Core::LinAlg::SparseMatrix> t2_;  /// global trafo matrix for Omega_2
    //@}

    //! @name maps
    Teuchos::RCP<const Epetra_Map>
        p12_dofrowmap_;  /// row map of projection matrix P Omega_2 -> Omega_1
    Teuchos::RCP<const Epetra_Map>
        p12_dofdomainmap_;  /// domain map of projection matrix P Omega_2 -> Omega_1
    Teuchos::RCP<const Epetra_Map>
        p21_dofrowmap_;  /// row map of projection matrix P Omega_1 -> Omega_2
    Teuchos::RCP<const Epetra_Map>
        p21_dofdomainmap_;  /// domain map of projection matrix P Omega_1 -> Omega_2
    Teuchos::RCP<const Epetra_Map>
        p12_dofcolmap_;  /// column map of projection matrix P Omega_2 -> Omega_1
    Teuchos::RCP<const Epetra_Map>
        p21_dofcolmap_;  /// column map of projection matrix P Omega_1 -> Omega_2
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
    Core::Geo::Cut::plain_volumecell_set
        volcell_;  /// set of volume cells for direct divergence integration

    // search algorithm
    Core::LinAlg::Matrix<9, 3> dopnormals_;  /// dop normals for seach algorithm

    // input
    DualQuad dualquad_;  /// type of quadratic weighting interpolation

    /// strategy for element information transfer (mainly material, but can be more)
    Teuchos::RCP<Core::VolMortar::UTILS::DefaultMaterialStrategy> materialstrategy_;

    //! @name mesh initialization

    // maps for mesh init
    Teuchos::RCP<Epetra_Map> xa_;
    Teuchos::RCP<Epetra_Map> xb_;
    Teuchos::RCP<Epetra_Map> mergedmap_;

    // mortar matrices for mesh init
    Teuchos::RCP<Core::LinAlg::SparseMatrix> dmatrix_xa_;  /// global Mortar matrix D for field A
    Teuchos::RCP<Core::LinAlg::SparseMatrix> dmatrix_xb_;  /// global Mortar matrix D for field B
    Teuchos::RCP<Core::LinAlg::SparseMatrix> mmatrix_xa_;  /// global Mortar matrix M for field A
    Teuchos::RCP<Core::LinAlg::SparseMatrix> mmatrix_xb_;  /// global Mortar matrix M for field B

    //@}
  };

}  // namespace Core::VolMortar


FOUR_C_NAMESPACE_CLOSE

#endif
