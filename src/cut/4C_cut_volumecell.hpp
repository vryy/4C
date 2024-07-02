/*---------------------------------------------------------------------*/
/*! \file

\brief Cut Volumecell (3d-Object described by 2d surfaces(facets))

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_VOLUMECELL_HPP
#define FOUR_C_CUT_VOLUMECELL_HPP

#include "4C_config.hpp"

#include "4C_cut_enum.hpp"
#include "4C_cut_facet.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
  {
    class Point;
    class Facet;
    class Element;
    class BoundaryCell;
    class IntegrationCell;
    class IntegrationCellCreator;
    class Mesh;


    /// a partial volume of a cutted element of any shape
    /*!
      \brief Volume cell is a partial (or whole) element that is defined by its
      facets.

      It knows

      - any number of integration cells (that have simple shapes)
      - any number of boundary cells (at the cut facets with outward normals)
      - its position to the cut surface
      - the dof set number of all element nodes
     */


    class VolumeCell
    {
     public:
      VolumeCell(const plain_facet_set& facets,
          const std::map<std::pair<Point*, Point*>, plain_facet_set>& volume_lines,
          Element* element);

      /*!
      \brief Check whether the volumecell has this point as one of its defining entities
      */
      bool Contains(Point* p);

      bool Contains(Core::LinAlg::Matrix<3, 1>& x);

      bool IsEqual(const plain_facet_set& vcell) const;

      void Neighbors(Point* p, const plain_volumecell_set& cells, const plain_volumecell_set& done,
          plain_volumecell_set& connected, plain_element_set& elements);

      void Neighbors(Point* p, const plain_volumecell_set& cells, const plain_volumecell_set& done,
          plain_volumecell_set& connected);

      /*!
      \brief Return the parent element from which this volumecell is generated
       */
      Element* parent_element() { return element_; }

      /*!
      \brief Return the ID of the parent element, from which this volumecell is generated
       */
      int GetParentElementId() const;

      void GetAllPoints(Mesh& mesh, PointSet& cut_points);

      void create_tet4_integration_cells(Mesh& mesh, const std::vector<std::vector<Point*>>& tets,
          const std::map<Facet*, std::vector<Point*>>& sides_xyz);

      void GetIntegrationCells(plain_integrationcell_set& cells);

      /*!
      \brief Construct the Gaussian quadrature rule by performing moment fitting over the
      volumecells
       */
      void moment_fit_gauss_weights(Element* elem, Mesh& mesh, bool include_inner,
          Core::Geo::Cut::BCellGaussPts BCellgausstype);

      /*!
      \brief Construct the Gaussian quadrature rule by triangulating the facets of volumecell and
      applying divergence theorem
       */
      void direct_divergence_gauss_rule(Element* elem, Mesh& mesh, bool include_inner,
          Core::Geo::Cut::BCellGaussPts BCellgausstype =
              Core::Geo::Cut::BCellGaussPts_Tessellation);

      /*!
      \brief Project the integration rule generated w.r to the global coordinates of the element to
      its local coordinate system
       */
      void project_gauss_points_to_local_coodinates();

      /*!
      \brief Generate boundarycells for the cut facets. This is used when volumecells are treated by
      momentfitting and boundarycells are by tessellation
       */
      void generate_boundary_cells(Mesh& mesh, const Core::Geo::Cut::Point::PointPosition posi,
          Element* elem, int BaseNos, Core::Geo::Cut::BCellGaussPts BCellgausstype);

      void generate_boundary_cells_level_set_side(Mesh& mesh,
          const Core::Geo::Cut::Point::PointPosition posi, Element* elem, Facet* facet, int BaseNos,
          Core::Geo::Cut::BCellGaussPts BCellgausstype);

      /*!
      \brief Get the boundarycells generated from this volumecell
         get a map of boundary cells for all cutting sides, key= side-Id, value= vector of boundary
      cells note that the boundary cells of subsides with the same side id are stored now in one key
       */
      void GetBoundaryCells(std::map<int, std::vector<Core::Geo::Cut::BoundaryCell*>>& bcells);
      void get_boundary_cells_to_be_integrated(
          std::map<int, std::vector<Core::Geo::Cut::BoundaryCell*>>& bcells);


      /// collect cut sides
      void CollectCutSides(plain_int_set& cutside_ids);

      ///
      void ConnectNodalDOFSets(bool include_inner);

      /*!
      \brief Return the position of volumecell whether they are in the fluid region in structure
       */
      Core::Geo::Cut::Point::PointPosition Position() const { return position_; }

      /*!
      \brief Identify the position of volumecell whether they are in the fluid region in structure
       */
      void Position(Core::Geo::Cut::Point::PointPosition position);

      inline void print() const { print(std::cout); };
      void print(std::ostream& stream) const;

      const plain_facet_set& Facets() const { return facets_; }

      bool IsCut(Side* s)
      {
        for (plain_facet_set::iterator fac = facets_.begin(); fac != facets_.end(); fac++)
        {
          if ((*fac)->IsCutSide(s)) return true;
        }
        return false;
      }

      /*!
      \brief Return the integrationcells generated from this volumecell
       */
      const plain_integrationcell_set& IntegrationCells() const { return integrationcells_; }

      /*!
      \brief Return the boundarycells generated from this volumecell
       */
      const plain_boundarycell_set& BoundaryCells() const { return bcells_; }

      bool Empty() const { return integrationcells_.size() == 0 and bcells_.size() == 0; }

      void Disconnect();

      const std::vector<int>& NodalDofSet() const { return nodaldofset_; }

      void SetNodalDofSet(const std::vector<int>& nds)
      {
        nodaldofset_.clear();
        std::copy(nds.begin(), nds.end(), std::inserter(nodaldofset_, nodaldofset_.end()));
      }


      /*!
      \brief Create integration cell of specified shape
       */
      void NewIntegrationCell(Mesh& mesh, Core::FE::CellType shape, const std::vector<Point*>& x);

      /// create line2 volume cell
      void NewLine2Cell(Mesh& mesh, const std::vector<Point*>& points);

      /// create tri3 volume cell
      void NewTri3Cell(Mesh& mesh, const std::vector<Point*>& points);

      /// create quad4 volume cell
      void NewQuad4Cell(Mesh& mesh, const std::vector<Point*>& points);

      /// create hex8 volume cell
      void NewHex8Cell(Mesh& mesh, const std::vector<Point*>& points);

      /// create tet4 volume cell
      IntegrationCell* NewTet4Cell(Mesh& mesh, const std::vector<Point*>& points);

      /// create wedge6 volume cell
      void NewWedge6Cell(Mesh& mesh, const std::vector<Point*>& points);

      /// create pryramid5 volume cell
      void NewPyramid5Cell(Mesh& mesh, const std::vector<Point*>& points);

      /*!
      \brief Create boundary cells of specified shape
      */
      void NewBoundaryCell(
          Mesh& mesh, Core::FE::CellType shape, Facet* f, const std::vector<Point*>& x);

      /*!
      \brief Create Point1 boundary cell
      */
      void NewPoint1Cell(Mesh& mesh, Facet* f, const std::vector<Point*>& x);

      /// Create Line2 boundary cell
      void NewLine2Cell(Mesh& mesh, Facet* f, const std::vector<Point*>& x);

      /*!
      \brief Create Tri3 boundary cell
      */
      void NewTri3Cell(Mesh& mesh, Facet* f, const std::vector<Point*>& x);

      /*!
      \brief Create Quad4 boundary cell
      */
      void NewQuad4Cell(Mesh& mesh, Facet* f, const std::vector<Point*>& x);

      /*!
      \brief Create boundary cell which contains more than 4 corner points
      */
      void NewArbitraryCell(Mesh& mesh, Facet* f, const std::vector<Point*>& x,
          const Core::FE::GaussIntegration& gp, const Core::LinAlg::Matrix<3, 1>& normal);

      /*!
      \brief Return the volume of the cell
       */
      double Volume() { return volume_; };

      void TestSurface();

      /*!
      \brief Return the number of Gauss points generated over the volumecell when tessellation is
      used
       */
      int NumGaussPoints(Core::FE::CellType shape);

      /*!
      \brief Write Geometry of volumecell together with gauss points produced from moment fitting
      method into GMSH output
       */
      void dump_gmsh_gauss_points_mom_fit(const std::vector<std::vector<double>>& gauspts);

      /*!
      \brief Write Geometry of volumecell together with gauss points produced from tessellation into
      GMSH output
       */
      void dump_gmsh_gauss_points_tessellation();

      /*!
      \brief Write Geometry of volumecell as lines into GMSH output
       */
      void DumpGmsh(std::ofstream& filename);

      /*!
      \brief Write geometry of volumecell as solid objects instead of lines
       */
      void DumpGmshSolid(std::ofstream& filename, Mesh& mesh);

      /*!
      \brief return the integration rule for this volumecell when MomentFitting or DirectDivergence
      is used
       */
      Teuchos::RCP<Core::FE::GaussPoints> get_gauss_rule() { return gp_; }

      /*!
      \brief replace the integration rule for this volumecell when MomentFitting or DirectDivergence
      is used
       */
      void SetGaussRule(Teuchos::RCP<Core::FE::GaussPoints> gps) { gp_ = gps; }

      /*!
      \brief Set the volume of this cell
       */
      void SetVolume(const double vol) { volume_ = vol; }

      /*!
      \brief Return the equation of the reference plane when DirectDivergence is used for volumecell
      treatment
       */
      const std::vector<double>& GetRefEqnPlane() { return ref_eqn_plane_; }

      /*!
      \brief Return Ids of all the points associated with this volumecell
       */
      const std::set<int>& VolumeCellPointIds();

      /*!
      \brief | Find Position of the Volumecell based on the orientation of the cut_sides
      */
      bool set_position_cut_side_based();

      /*!
      \brief Chcek whether this point is inside, outside or on boundary of this this volumecell.
      The output std::string will be either "outside", "inside" or "onBoundary"
       */
      std::string IsThisPointInside(Point* pt);

      /*!
      \brief Chcek whether this Global coordinate point is inside, outside or on boundary of this
      this volumecell. The output std::string will be either "outside", "inside" or "onBoundary"
       */
      std::string IsThisPointInside(Core::LinAlg::Matrix<3, 1>& x);

      /*!
      \brief Integrate specific predifined functions over the volumecell by integrating over each of
      the integration cells resulting from Tessellation. For MomentFitting and DIrectDivergence this
      can't be used
       */
      void integrate_specific_functions_tessellation();

      template <Core::FE::CellType distype>
      Teuchos::RCP<Core::FE::GaussPoints> create_projected(Core::Geo::Cut::IntegrationCell* ic);

      /*!
      \brief Returns whether this volumecell is negligibly small (used only in DirectDIvergence
      approach)
       */
      bool IsNegligiblySmall() { return this->is_negligible_small_; }

     private:
      void set_tet_points(const int* totet4, const std::vector<Point*>& points,
          std::vector<Point*>& tet4_points) const
      {
        for (int i = 0; i < 4; ++i)
        {
          tet4_points[i] = points[totet4[i]];
        }
      }

      /*!
      \brief Check whether the numbering of points need to be reversed. This is to ensure outward
      pointing normal for the boundarycells when triangulation is performed
       */
      bool to_reverse(const Core::Geo::Cut::Point::PointPosition posi,
          const std::vector<double>& parEqn, const std::vector<double>& facetEqn);

      /*!
      \brief return the Gauss points computed using moment fitting equations
      */
      Teuchos::RCP<Core::FE::GaussPoints> gauss_points_fitting();

      /*!
      \brief Generate internal gauss rule for every integration point on the facet when
      DirectDivergence method is used
       */
      Teuchos::RCP<Core::FE::GaussPoints> generate_internal_gauss_rule(
          Teuchos::RCP<Core::FE::GaussPoints>& gp);

      /// the element this is a part of
      Element* element_;

      /// position relative to cut surface
      Core::Geo::Cut::Point::PointPosition position_;

      /// defining facets
      plain_facet_set facets_;

      /// volume preserving simple shaped cells
      plain_integrationcell_set integrationcells_;

      /// boundary cells at cut surface with outward normals
      plain_boundarycell_set bcells_;

      /// dof set number of all element nodes
      std::vector<int> nodaldofset_;

      // point ids of a volume cell
      std::set<int> vcpoints_ids_;

      static int hex8totet4[5][4];
      static int wedge6totet4[3][4];
      static int pyramid5totet4[2][4];

      /// volume of this volumecell
      double volume_;

      /// Returns true if this volumecell is negligibly small (used only in DirectDivergence method)
      bool is_negligible_small_;

      /// store the Gauss point weights and location for volumecells when moment fitting equations
      /// are used to arrive at them
      Core::LinAlg::SerialDenseVector weights_;
      std::vector<std::vector<double>> gaus_pts_;

      /// store the Gauss point weights and location for boundarycells when moment fitting equations
      /// are used to arrive at them
      Core::LinAlg::SerialDenseVector bcellweights_;
      std::vector<std::vector<double>> bcellgaus_pts_;

      /// Equation of plane which contains the reference facet when using DirectDivergence
      std::vector<double> ref_eqn_plane_;

      /// Gauss rule for this volumecell. In case of DirectDivergence method, this just stores main
      /// Gauss points
      Teuchos::RCP<Core::FE::GaussPoints> gp_;
    };

  }  // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
