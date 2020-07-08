/*----------------------------------------------------------------------*/
/*! \file

\brief class that holds information about a mesh that is cut or about a cutmesh that cuts another
mesh

\level 3
 *------------------------------------------------------------------------------------------------*/

#ifndef CUT_MESH_H
#define CUT_MESH_H

#include <list>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyTraits.hpp>

#include "cut_boundingbox.H"
#include "cut_facet.H"

#include "../drt_inpar/inpar_cut.H"

#ifdef DEBUG
#include "cut_node.H"
#include "cut_edge.H"
#endif


namespace GEO
{
  // class BVTree;
  class SearchTree;

  namespace CUT
  {
    class Node;
    class Edge;
    class Side;
    class Element;

    class BoundingBox;

    class PointPool;
    class Options;

    class Point;
    class Line;
    class Facet;
    class VolumeCell;

    class BoundaryCell;
    class Point1BoundaryCell;
    class Line2BoundaryCell;
    class Tri3BoundaryCell;
    class Quad4BoundaryCell;
    class ArbitraryBoundaryCell;

    class IntegrationCell;
    class Line2IntegrationCell;
    class Tri3IntegrationCell;
    class Quad4IntegrationCell;
    class Hex8IntegrationCell;
    class Tet4IntegrationCell;
    class Wedge6IntegrationCell;
    class Pyramid5IntegrationCell;



    /*!
    \brief All the geometrical entities (mesh, volume, cut surface, nodes etc.) are contained in
    this.

      There is one background mesh and one mesh for the cut surface. These meshes
      have different objects but share the same points via the point pool.

      Mesh does the memory management for the whole thing. Therefore, all
      creation of cut library objects is done via the mesh.
     */
    class Mesh
    {
     public:
      /// constructor
      Mesh(Options& options, double norm = 1, Teuchos::RCP<PointPool> pp = Teuchos::null,
          bool cutmesh = false, int myrank = -1);

      /*========================================================================*/
      //! @name general create-routines for elements and sides
      /*========================================================================*/

      /// creates a new element, dependent on distype
      Element* CreateElement(
          int eid, const std::vector<int>& nids, DRT::Element::DiscretizationType distype);

      /// creates a new side, dependent on distype
      Side* CreateSide(
          int sid, const std::vector<int>& nids, DRT::Element::DiscretizationType distype);

      /*========================================================================*/
      //! @name Create-routines for elements
      /*========================================================================*/
      /// @{
      /// creates a new line2 element based on element id and node ids
      Element* CreateLine2(int eid, const std::vector<int>& nids);

      /// creates a new tri3 element based on element id and node ids
      Element* CreateTri3(int eid, const std::vector<int>& nids);

      /// creates a new quad4 element based on element id and node ids
      Element* CreateQuad4(int eid, const std::vector<int>& nids);

      /// creates a new tet4 element based on element id and node ids
      Element* CreateTet4(int eid, const std::vector<int>& nids);

      /// creates a new pyramid5 element based on element id and node ids
      Element* CreatePyramid5(int eid, const std::vector<int>& nids);

      /// creates a new wedge6 element based on element id and node ids
      Element* CreateWedge6(int eid, const std::vector<int>& nids);

      /// creates a new hex8 element based on element id and node ids
      Element* CreateHex8(int eid, const std::vector<int>& nids);
      /// @}

      /*========================================================================*/
      //! @name Create-routines for sides
      /*========================================================================*/
      /// @{
      /// creates a new tri3 side based on side id and node ids
      Side* CreateTri3Side(int sid, const std::vector<int>& nids);

      /// creates a new quad4 side based on side id and node ids
      Side* CreateQuad4Side(int sid, const std::vector<int>& nids);
      /// @}

      /*========================================================================*/
      //! @name Create-routines for points, lines, facets and volumecells
      /*========================================================================*/

      /// creates a new point, optional information about cut-edge and cut-side, and whether this is
      /// hapenning during
      // loading of the mesh
      Point* NewPoint(const double* x, Edge* cut_edge, Side* cut_side, double tolerance,
          double tol_scale = 1.0);

      /// creates a new line
      void NewLine(Point* p1, Point* p2, Side* cut_side1, Side* cut_side2, Element* cut_element,
          std::vector<Line*>* newlines = NULL);

      /// ?
      bool NewLinesBetween(const std::vector<Point*>& line, Side* cut_side1, Side* cut_side2,
          Element* cut_element, std::vector<Line*>* newlines = NULL);

      /// creates a new facet, consists of points, additional bool if it is a facet on a cutsurface
      Facet* NewFacet(const std::vector<Point*>& points, Side* side, bool cutsurface);

      /// creates a new volumecell, consists of facets
      VolumeCell* NewVolumeCell(const plain_facet_set& facets,
          const std::map<std::pair<Point*, Point*>, plain_facet_set>& volume_lines,
          Element* element);

      /*========================================================================*/
      //! @name Create-routines for 0D boundary cells
      /*========================================================================*/

      /// creates a new point1 boundary cell
      Point1BoundaryCell* NewPoint1Cell(
          VolumeCell* volume, Facet* facet, const std::vector<Point*>& points);

      /*========================================================================*/
      //! @name Create-routines for 1D boundary cells
      /*========================================================================*/

      Line2BoundaryCell* NewLine2Cell(
          VolumeCell* volume, Facet* facet, const std::vector<Point*>& points);

      /*========================================================================*/
      //! @name Create-routines for 2D boundary cells
      /*========================================================================*/

      /// creates a new tri3 boundary cell
      Tri3BoundaryCell* NewTri3Cell(
          VolumeCell* volume, Facet* facet, const std::vector<Point*>& points);

      /// creates a new quad4 boundary cell
      Quad4BoundaryCell* NewQuad4Cell(
          VolumeCell* volume, Facet* facet, const std::vector<Point*>& points);

      /// creates a new ??? boundary cell
      ArbitraryBoundaryCell* NewArbitraryCell(VolumeCell* volume, Facet* facet,
          const std::vector<Point*>& points, const DRT::UTILS::GaussIntegration& gaussRule,
          const LINALG::Matrix<3, 1>& normal);


      /*========================================================================*/
      //! @name Create-routines for 1D integration cells
      /*========================================================================*/

      /// creates a new line2 integration cell
      Line2IntegrationCell* NewLine2Cell(
          Point::PointPosition position, const std::vector<Point*>& points, VolumeCell* cell);

      /*========================================================================*/
      //! @name Create-routines for 2D integration cells
      /*========================================================================*/

      /// creates a new tri3 integration cell
      Tri3IntegrationCell* NewTri3Cell(
          Point::PointPosition position, const std::vector<Point*>& points, VolumeCell* cell);

      /// creates a new tri3 integration cell
      Quad4IntegrationCell* NewQuad4Cell(
          Point::PointPosition position, const std::vector<Point*>& points, VolumeCell* cell);


      /*========================================================================*/
      //! @name Create-routines for 3D integration cells
      /*========================================================================*/

      /// creates a new hex8 integration cell
      Hex8IntegrationCell* NewHex8Cell(
          Point::PointPosition position, const std::vector<Point*>& points, VolumeCell* cell);

      /// creates a new tet4 integration cell, based on points
      Tet4IntegrationCell* NewTet4Cell(
          Point::PointPosition position, const std::vector<Point*>& points, VolumeCell* cell);

      /// creates a new hex8 integration cell, based on xyz coordinates
      Tet4IntegrationCell* NewTet4Cell(
          Point::PointPosition position, const Epetra_SerialDenseMatrix& xyz, VolumeCell* cell);

      /// creates a new wedge6 integration cell
      Wedge6IntegrationCell* NewWedge6Cell(
          Point::PointPosition position, const std::vector<Point*>& points, VolumeCell* cell);

      /// creates a new pyramid5 integration cell
      Pyramid5IntegrationCell* NewPyramid5Cell(
          Point::PointPosition position, const std::vector<Point*>& points, VolumeCell* cell);


      /*========================================================================*/
      //! @name Basic routines to run the collision detection
      /*========================================================================*/

      /// build the static search tree for the collision detection in the self cut
      void BuildSelfCutTree();

      /// build the static search tree for the collision detection
      void BuildStaticSearchTree();

      /// detects if a side of the cut mesh possibly collides with an element of the background mesh
      void SearchCollisions(Mesh& cutmesh);


      /*========================================================================*/
      //! @name Basic routines to run the cut, creates a cloud of cut points
      /*========================================================================*/

      /// Cuts the background elements of the mesh with all the cut sides
      void Cut(Mesh& mesh, plain_element_set& elements_done);

      /// Cuts the background elements with this considered side
      void Cut(Side& side, const plain_element_set& done, plain_element_set& elements_done);

      /// Cuts the background elements with this levelset side
      void Cut(Side& side);

      /// Used in TetMeshIntersection, however, unclear for what use??!
      void RectifyCutNumerics();

      /// finds intersections between sides and edges
      void FindCutPoints();

      /*========================================================================*/
      //! @name Basic routines to create lines, facets and volumecells
      /*========================================================================*/

      /// create cut lines based on the point cloud
      void MakeCutLines();

      /// create facets based on the cut lines
      void MakeFacets();

      /// create volumecells based on created facets
      void MakeVolumeCells();


      /*========================================================================*/
      //! @name Basic routines to determine positions of nodes, facets, volumecells
      /*========================================================================*/

      /// find node positions and propagate the positions to facets, points and volumecells
      void FindNodePositions();

      /// ?
      void FindLSNodePositions();

      /// find facet positions for remaining facets, points, volumecells that have not been found
      /// using FindNodePositions()
      void FindFacetPositions();

      /// Check if there are nodes whose position is undecided (often the case in parallel), return
      /// whether undecided node positions available
      bool CheckForUndecidedNodePositions(
          std::map<int, int>& undecided_node, std::map<plain_int_set, int>& undecided_shadow_nodes);


      /*========================================================================*/
      //! @name Basic routines to determine nodal dofsets
      /*========================================================================*/

      /// still used???
      void FindNodalDOFSets(bool include_inner);


      /*========================================================================*/
      //! @name Basic routines to create integration cells and/or integration points
      /*========================================================================*/

      /// Execute Tessellation with QHULL for each element to generate integrationcells
      void CreateIntegrationCells(int count, bool tetcellsonly = false);

      /// Call the moment fitting method for each element to generate the Gaussian integration rule
      void MomentFitGaussWeights(bool include_inner, INPAR::CUT::BCellGaussPts Bcellgausstype);

      /// Call the DirectDivergence method for each element to generate the Gaussian integration
      /// rule
      void DirectDivergenceGaussRule(bool include_inner, INPAR::CUT::BCellGaussPts Bcellgausstype);


      /*========================================================================*/
      //! @name others ?
      /*========================================================================*/

      /// ?
      void RemoveEmptyVolumeCells();

      /// test if for all elements the element volume is equal to the volume of all integration
      /// cells
      void TestElementVolume(
          bool fatal, INPAR::CUT::VCellGaussPts VCellGP = INPAR::CUT::VCellGaussPts_Tessellation);

      /*!
      \brief Find the difference between the volume of background element and the sum of volume of
      all integration cells. There should be no difference between these two
       */
      void TestElementVolume(DRT::Element::DiscretizationType shape, Element& e, bool fatal,
          INPAR::CUT::VCellGaussPts VCellGP);


      /*========================================================================*/
      //! @name print statistics
      /*========================================================================*/

      /// ???
      void PrintCellStats();

      /// ???
      void Status();

      /// print all facets
      void PrintFacets();


      /*========================================================================*/
      //! @name GMSH output routines
      /*========================================================================*/

      /// Write full Gmsh Output
      void DumpGmsh(std::string name);

      /*!
      \brief Output information about the volume cell.
      If the cut has a level set side. Also write output for level set values and gradients.
       */
      void DumpGmshVolumeCells(std::string name, bool include_inner);

      /// ?
      void DumpGmshIntegrationCells(std::string name);

      /* writre boundary cells belonging to a volume cell with "pos" position
       * and their normals to the given file */
      void DumpGmshBoundaryCells(std::ofstream& file, Point::PointPosition pos);

      /// ?
      void DumpGmshVolumeCells(std::string name);

      /// DebugDump to call before runtime error!!!
      void DebugDump(GEO::CUT::Element* ele, std::string file = "", int line = -1);


      /*========================================================================*/
      //! @name Get routines for nodes, sides and elements
      /*========================================================================*/

      /// ? -> used in cut_tetmeshintersection
      void NewNodesFromPoints(std::map<Point*, Node*>& nodemap);

      /// get a map of node id and the pointer to the node
      void GetNodeMap(std::map<int, Node*>& nodemap);

      /// Returns the node with given id
      Node* GetNode(int nid) const;

      /*!
      \brief Returns the unique shadow node
      identified by given nids of a quad8 boundary side or all nodes of hex20 element
      for the inner shadow node
       */
      Node* GetNode(const plain_int_set& nids) const;

      /*!
      \brief If node with the given id exists return the node, else create a new node with
      given coordinates and levelset value
       */
      Node* GetNode(int nid, const double* xyz, double lsv = 0.0, double tolerance = 0.0);

#if 0
  Node* GetNode( int nid, Point * p, double lsv=0.0 );
#endif

      /// ?
      Node* GetNode(const plain_int_set& nids, const double* xyz, double lsv = 0.0);

      /// get the edge with begin node and end node
      Edge* GetEdge(Node* begin, Node* end);

      /// get the side that contains the nodes with the following node ids
      Side* GetSide(std::vector<int>& nids) const;

      /// ???
      Side* GetSide(int sid, const std::vector<int>& nids, const CellTopologyData* top_data);

#if 0
  Side* GetSide( int sid,
                 const std::vector<Point*> & points,
                 const CellTopologyData * top_data );
#endif

      /// ???
      Side* GetSide(int sid, const std::vector<Node*>& nodes, const CellTopologyData* top_data);

      /// Returns the element with given id
      Element* GetElement(int eid);

      /*!
      \brief  If element with the given id exists return the element, else create a new element
      with given node ids and details given in cell topology data
       */
      Element* GetElement(int eid, const std::vector<int>& nids, const CellTopologyData& top_data,
          bool active = true);
#if 0
  Element* GetSubElement( int parenteid, int subeid,
                       const std::vector<int> & nids,
                       const CellTopologyData & top_data,
                       bool active=true );

#endif

#if 0
  Element* GetElement( int eid,
                       const std::vector<Point*> & nids,
                       const CellTopologyData & top_data );
#endif

      /*! \brief Create a new element 1D/2D/3D element with given nodes.
       *
       *  All details of the element are in cell topology data. */
      GEO::CUT::Element* GetElement(
          int eid, const std::vector<Node*>& nodes, const CellTopologyData& top_data, bool active);

      /*! \brief Create a new element with desired element dimension
       *
       *  dim = 1: Create a new element 1D/line element with given nodes.
       *
       *  All details of the element are in cell topology data.
       *  This routine uses the elements (1-D) as edges (1-D) and sides (1-D).
       *
       *  dim = 1: Create a new element 2D/surface element with given nodes.
       *
       *  All details of the element are in cell topology data.
       *  This routine creates new sides (1-D) and uses them also as edges (1-D).
       *
       *  dim = 3: Create a new element 3D element with given nodes.
       *
       *  All details of the element (3-D) are in cell topology data.
       *  This routine creates also the corresponding edges (1-D) and sides (2-D) for
       *  3-dimensional elements. */
      template <unsigned dim>
      Element* GetElement(int eid, const std::vector<Node*>& nodes,
          const CellTopologyData& top_data, bool active = true);

      /*========================================================================*/
      //! @name Get routines for points, edges, volumecells
      /*========================================================================*/

      /// get the octTree based PointPool that contains all points of the current mesh
      Teuchos::RCP<PointPool> Points() { return pp_; }

      /// get a list of all volumecells
      const std::list<Teuchos::RCP<VolumeCell>>& VolumeCells() const { return cells_; }

      /// ???
      const std::map<plain_int_set, Teuchos::RCP<Edge>>& Edges() const { return edges_; }


      /*========================================================================*/
      //! @name others
      /*========================================================================*/

      /// check if xyz-coordinates lie within the mesh's bounding box
      bool WithinBB(const Epetra_SerialDenseMatrix& xyz);

      /// check if the element lies within the bounding box
      bool WithinBB(Element& element);

      //     Only used in cut_test_volume.cpp
      void CreateSideIds_CutTest(int lastid = 0);

      // only used for testing
      int CreateSideIdsAll_CutTest(int lastid = 0);

      //     Only used in cut_test_volume.cpp
      void AssignOtherVolumeCells_CutTest(const Mesh& other);

      /// return the options
      Options& CreateOptions() { return options_; }

      /// return the options
      Options& GetOptions() const { return options_; }

      /// ???
      void TestVolumeSurface();

      /// Test if the area of a cut facet is covered by the same area of boundary cells on both
      /// sides.
      void TestFacetArea(bool istetmeshintersection = false);


      /*========================================================================*/
      //! @name Routines which provides the interface to the selfcut
      /*========================================================================*/

      /// Returns all sides of the cutmesh
      const std::map<plain_int_set, Teuchos::RCP<Side>>& Sides() const { return sides_; }

      /// Returns of search tree all sides of the cutmesh
      const Teuchos::RCP<GEO::SearchTree>& SelfCutTree() const { return selfcuttree_; }

      /// Returns the bounding volumes of all sides of the cutmesh
      const std::map<int, LINALG::Matrix<3, 2>>& SelfCutBvs() const { return selfcutbvs_; }

      /// Returns the map of all sides of the cutmesh
      const std::map<int, Side*>& ShadowSides() const { return shadow_sides_; }

      /// Returns all nodes of the cutmesh
      const std::map<int, Teuchos::RCP<Node>>& Nodes() const { return nodes_; }

      /// Creates a new node in the cutmesh
      void GetNode(int nid, Node* node) { nodes_[nid] = Teuchos::rcp(node); }

      /// Creates a new edge in the cutmesh
      void GetEdge(plain_int_set eid, const Teuchos::RCP<Edge>& edge) { edges_[eid] = edge; }

      /// Erases a side of the cutmesh
      void EraseSide(plain_int_set sid) { sides_.erase(sid); }

      /// Erases a edge of the cutmesh
      void EraseEdge(plain_int_set eid) { edges_.erase(eid); }

      /// Erases a node of the cutmesh
      void EraseNode(int nid) { nodes_.erase(nid); }

      /// Move this Side from to the Storage container of the Mesh
      void MoveSidetoStorage(plain_int_set sid)
      {
        storagecontainer_sides_[sid] = sides_[sid];
        sides_.erase(sid);
      }

      /// Move this Node from to the Storage container of the Mesh
      void MoveNodetoStorage(int nid)
      {
        storagecontainer_nodes_[nid] = nodes_[nid];
        nodes_.erase(nid);
      }

     private:
      /*========================================================================*/
      //! @name private member functions
      /*========================================================================*/

      /// ?
      Edge* GetEdge(const plain_int_set& nids, const std::vector<Node*>& nodes,
          const CellTopologyData& edge_topology);

      /// ?
      Side* GetSide(int sid, const plain_int_set& nids, const std::vector<Node*>& nodes,
          const std::vector<Edge*>& edges, const CellTopologyData& side_topology);

      /// Create new line between the two given cut points that are in given two cut sides
      GEO::CUT::Line* NewLineInternal(
          Point* p1, Point* p2, Side* cut_side1, Side* cut_side2, Element* cut_element);


      /*========================================================================*/
      //! @name private member variables
      /*========================================================================*/

      /// options container
      Options& options_;

      /// mesh dependent point lookup norm
      double norm_;

      /// shared point storage with geometry based access (octtree)
      Teuchos::RCP<PointPool> pp_;

      /// bounding box of this mesh
      Teuchos::RCP<BoundingBox> bb_;

      /// (output) flag for cut mesh
      ///  TODO: Remove this!
      //     Only used in cut_test_volume.cpp
      bool cutmesh_;

      /// the spatial partitioning octree for a fast collision detection in the self cut
      Teuchos::RCP<GEO::SearchTree> selfcuttree_;

      /// the bounding volumes for a fast collision detection in the self cut
      std::map<int, LINALG::Matrix<3, 2>> selfcutbvs_;

      /// the spatial partitioning octree for a fast collision detection
      Teuchos::RCP<GEO::SearchTree> searchtree_;

      /// the bounding volumes for a fast collision detection
      std::map<int, LINALG::Matrix<3, 2>> boundingvolumes_;

      /*========================================================================*/
      //! @name Containers that hold all those mesh objects
      /*========================================================================*/

      /// Plain pointers are used within the library. Memory management is done here.
      std::list<Teuchos::RCP<Line>> lines_;
      std::list<Teuchos::RCP<Facet>> facets_;
      std::list<Teuchos::RCP<VolumeCell>> cells_;
      std::list<Teuchos::RCP<BoundaryCell>> boundarycells_;
      std::list<Teuchos::RCP<IntegrationCell>> integrationcells_;

      /// nodes by unique id, contains also shadow nodes with negative node-Ids
      /// Remark: the negative nids of shadow nodes are not unique over processors!
      std::map<int, Teuchos::RCP<Node>> nodes_;

      /// edges by set of nodal ids
      std::map<plain_int_set, Teuchos::RCP<Edge>> edges_;

      /// sides by set of nodal ids
      std::map<plain_int_set, Teuchos::RCP<Side>> sides_;

      /// elements by unique id
      std::map<int, Teuchos::RCP<Element>> elements_;

      /// internally generated nodes by nodal ids of element nodes
      /// there is at most one unique shadow node for each 2D side element (no for quad9 side, one
      /// for quad8 side, no for tri6 side, no for linear elements and so on) it is the center node
      /// of a quad8 side in case of a hex20 element, the eight nodes are used as key for the
      /// boundary shadow node the inner center node of hex20 element is a (kind of) unique inner!
      /// shadow node and is stored using all 20 nodes of the hex20 element as key for the map
      /// warning: in the context of the selfcut, nodes are erased for the cut mesh
      /// these nodes are not erased in shadow_nodes_
      /// so don't use the shadow_nodes_ of the cut mesh
      std::map<plain_int_set, Node*> shadow_nodes_;

      // unique map between int ( < 0 ) and element; int is not eid and not parentid
      // we need this for the search tree in the self cut
      // is not erased during self cut!!!
      std::map<int, Side*> shadow_sides_;

      /// new: unique map between int ( < 0 ) and element; int is not eid and not parentid
      std::map<int, Teuchos::RCP<Element>> shadow_elements_;

      /// A storage container for sides which should not interact with the CUT anymore (If we want
      /// to access these geometric objects still from outside)
      std::map<plain_int_set, Teuchos::RCP<Side>> storagecontainer_sides_;

      /// A storage container for nodes which should not interact with the CUT anymore (If we want
      /// to access these geometric objects still from outside)
      std::map<int, Teuchos::RCP<Node>> storagecontainer_nodes_;

      /// processor id --> required just for output!
      int myrank_;

      //@}
    };

    // instantiation of the function template specializations
    template <>
    Element* Mesh::GetElement<1>(
        int eid, const std::vector<Node*>& nodes, const CellTopologyData& top_data, bool active);
    template <>
    Element* Mesh::GetElement<2>(
        int eid, const std::vector<Node*>& nodes, const CellTopologyData& top_data, bool active);
    template <>
    Element* Mesh::GetElement<3>(
        int eid, const std::vector<Node*>& nodes, const CellTopologyData& top_data, bool active);


  }  // namespace CUT
}  // namespace GEO

#endif
