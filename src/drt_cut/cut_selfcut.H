/*----------------------------------------------------------------------*/
/*! \file

\brief class that provides all routines to handle cutsides which cut each other

\level 2

\maintainer Christoph Ager
 *------------------------------------------------------------------------------------------------*/

#ifndef CUT_SELFCUT_H
#define CUT_SELFCUT_H

#include "cut_point.H"


namespace GEO
{
  namespace CUT
  {
    class Point;
    class Node;
    class Edge;
    class Side;
    class Cycle;
    class Mesh;
    class MeshHandle;
    class BoundingBox;

    class SelfCut
    {
     public:
      /// constructor
      SelfCut(MeshHandle& cut_mesh_handle, int myrank);

      /*========================================================================*/
      //! @name Main routines to run the self cut
      /*========================================================================*/

      void PerformSelfCut();

      /// detects cutsides which cut each other by finding the respective cutpoints
      bool CollisionDetection();

      /// replaces cutted sides by creating new nodes, edges and sides
      void MeshIntersection();

      /// erases nodes, edges and sides which lies inside a structure body by locating there
      /// position
      void ElementSelection();

      /// repair the mesh from potential islands
      void MeshCorrection();

      /// represents the result by using text viewer or gmsh
      void StatusInspection();

      /*========================================================================*/
      //! @name Basic routines to run the collision detection
      /*========================================================================*/

      /// finds cutsides which possibly cuts the considered side
      void FindCuttingSides();

      /// if there are two nodes from cut sides are at the same position, we merge it into one node
      bool MergeCoincidingNodes(Side* side1, Side* side2);

      /// opearations to modify the nodal ids of sides and corresponding edges are performed
      void operationsForNodeMerging(
          std::vector<std::pair<const Node*, const Node*>> repl, bool initial);

      /// templated function to delete "nod" by "replwith" in both edge and side data-structures
      template <typename A, typename B>
      void ModifyEdgeOrSideMap(std::map<A, B>& data, int nod, int replwith);

      /// finds cutpoints by searching the intersection between the edges of one side with the other
      /// side and vice versa
      void FindSelfCutPoints();

      /// gets all cutted sides and there nodes and edges to store them as privat variables
      void GetSelfCutObjects();

      /*========================================================================*/
      //! @name Basic routines to run the mesh intersection
      /*========================================================================*/

      /// creates a node at every selfcutpoint with respect to the corresponding cutsides
      void CreateSelfCutNodes();

      /// creates a edge between two selfcutnodes with respect to the corresponding cutsides
      void CreateSelfCutEdges();

      /// finds triangles with respect to the self cut by using a pointgraph and a triangulation
      /// method
      void FindSelfCutTriangulation();

      /// creates triangular sides out of the self cut triangles
      void CreateSelfCutSides();

      /// erases all cutsides which are cut by another side
      void EraseCuttedSides();

      /// erases all edges which are cut by a cutside
      void EraseCuttedEdges();

      /// Is this edge connected to the backgound mesh
      bool ConnectedtoBackground(Edge* edge);

      /// Is this node connected to the backgound mesh
      bool ConnectedtoBackground(Node* node);

      /*========================================================================*/
      //! @name Basic routines to run the element selection
      /*========================================================================*/

      /// locates the position of nodes, edges and sides of a structure body with respect to the
      /// other bodys
      void DetermineSelfCutPosition();

      /// locates the position of nodes, edges and sides for special cases
      void PropagateSelfCutPosition();

      /// erases sides which lies inside a structure body by locating there position
      void EraseInsideSides();

      /// erases edges which lies inside a structure body by locating there position
      void EraseInsideEdges();

      /// erases nodes which lies inside a structure body by locating there position
      void EraseInsideNodes();

      /// construct the connectivity of the nodes toi find potential islands in the cut mesh
      void ConstructConnectivity();

      /// find the next node for the construction of the connectivity
      void NextNode(Node* node, plain_node_set& remainingnodes, int count);

      /// identify the islands
      void FindIslands();

      /// Get next Sides
      void NextSides(Side* cutside, Teuchos::RCP<GEO::CUT::BoundingBox>& tmpbb,
          // plain_side_set allselfcutsides,
          plain_side_set& selfcutsides, plain_side_set& islandsides, bool& IsIsland);

      /*========================================================================*/
      //! @name Basic routines to run the status inspection
      /*========================================================================*/

      /// Status of the cutted sides for text viewer
      void CuttedSideStatusText();

      /// Status of the cutmesh for text viewer
      void CutMeshStatusText();

      /// Status of one problematic side for text viewer
      void ErrorStatusText(Side& cutside);

      /// Status of the cutted sides for gmsh
      void CuttedSideStatusGmsh(const std::string& name);

      /// Status of the cutmesh in gmsh
      void WallGmsh(const std::string& name);

      /// Status of the selfcut objects in gmsh
      void SCObjectsGmsh(const std::string& name);

      /// Status of the cutmesh in gmsh for my CMGM (old)
      void SCmgmGmsh(const std::string& name);


      /// Status of all sides in separates files for gmsh
      void AllSingleGmsh(const std::string& location);

      /// Status of one problematic side for gmsh
      void ErrorGmsh(const std::string& name, Side& cutside);

      /*========================================================================*/
      //! @name Auxiliary routines to run the self cut
      /*========================================================================*/

      /// cuts the edges of the first cutside with the other cutside
      void PerformSelfCut(Side& cutside, Side& otherside, PointSet& selfcutpoints);

      /// gets the edges of the new created sides
      void GetSelfCutEdges(Side& cutside);

      /// checks if the direction of rotation of the cycle is correct
      bool CheckNormal(std::vector<double> cutsideplane, Cycle& maincycle);

      /// deletes all pointers which are pointing to a soon to be erased side
      void EraseSidePointer(Side& cutside, bool erase_nodepointers);

      /// erases a side
      void EraseSide(std::vector<plain_int_set>& cutsideids);

      /// erases a side which is in the unphysical part
      void EraseInsideSide(std::vector<plain_int_set>& cutsideids);

      /// deletes all pointers which are pointing to a soon to be erased edge
      void EraseEdgePointer(Edge& cutsideedge);

      /// erases an edge
      void EraseEdge(std::vector<plain_int_set>& cutsideedgeids);

      void PlotAcross();

      void PointPlotHead();

      void NodePlotHead();

      void EdgePlotHead();

      void SidePlotHead();

      void PointPlot(Point& point);

      void NodePlot(Node& node);

      void EdgePlot(Edge& edge);

      void SidePlot(Side& side);

      void SideNodePlot(Side& side);


     private:
      // My Processor Id --> Only for output purposes
      int myrank_;

      /// The Cut Mesh
      GEO::CUT::Mesh& mesh_;

      /// The Cut Meshhandle (Interface to the outer world)
      GEO::CUT::MeshHandle& meshhandle_;

      std::map<plain_int_set, Teuchos::RCP<Side>> selfcut_sides_;

      std::map<plain_int_set, Teuchos::RCP<Edge>> selfcut_edges_;

      std::map<int, Teuchos::RCP<Node>> selfcut_nodes_;

      std::map<int, plain_node_set> selfcut_connectivity_;

      /// characteristic interface element size
      double meshsizeparam_;
    };

  }  // namespace CUT
}  // namespace GEO

#endif
