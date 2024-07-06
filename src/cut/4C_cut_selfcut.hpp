/*----------------------------------------------------------------------*/
/*! \file

\brief class that provides all routines to handle cutsides which cut each other

\level 2

 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_SELFCUT_HPP
#define FOUR_C_CUT_SELFCUT_HPP

#include "4C_config.hpp"

#include "4C_cut_point.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Core::Geo
{
  namespace Cut
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

      void perform_self_cut();

      /// detects cutsides which cut each other by finding the respective cutpoints
      bool collision_detection();

      /// replaces cutted sides by creating new nodes, edges and sides
      void mesh_intersection();

      /// erases nodes, edges and sides which lies inside a structure body by locating there
      /// position
      void element_selection();

      /// repair the mesh from potential islands
      void mesh_correction();

      /// represents the result by using text viewer or gmsh
      void status_inspection();

      /*========================================================================*/
      //! @name Basic routines to run the collision detection
      /*========================================================================*/

      /// finds cutsides which possibly cuts the considered side
      void find_cutting_sides();

      /// if there are two nodes from cut sides are at the same position, we merge it into one node
      bool merge_coinciding_nodes(Side* keep, Side* replace);

      /// opearations to modify the nodal ids of sides and corresponding edges are performed
      void operations_for_node_merging(
          std::vector<std::pair<const Node*, const Node*>> repl, bool initial);

      /// templated function to delete "nod" by "replwith" in both edge and side data-structures
      template <typename A, typename B>
      void modify_edge_or_side_map(std::map<A, B>& data, int nod, int replwith);

      /// finds cutpoints by searching the intersection between the edges of one side with the other
      /// side and vice versa
      void find_self_cut_points();

      /// gets all cutted sides and there nodes and edges to store them as privat variables
      void get_self_cut_objects();

      /*========================================================================*/
      //! @name Basic routines to run the mesh intersection
      /*========================================================================*/

      /// creates a node at every selfcutpoint with respect to the corresponding cutsides
      void create_self_cut_nodes();

      /// creates a edge between two selfcutnodes with respect to the corresponding cutsides
      void create_self_cut_edges();

      /// finds triangles with respect to the self cut by using a pointgraph and a triangulation
      /// method
      void find_self_cut_triangulation();

      /// creates triangular sides out of the self cut triangles
      void create_self_cut_sides();

      /// erases all cutsides which are cut by another side
      void erase_cutted_sides();

      /// erases all edges which are cut by a cutside
      void erase_cutted_edges();

      /// Is this edge connected to the backgound mesh
      bool connectedto_background(Edge* edge);

      /// Is this node connected to the backgound mesh
      bool connectedto_background(Node* node);

      /*========================================================================*/
      //! @name Basic routines to run the element selection
      /*========================================================================*/

      /// locates the position of nodes, edges and sides of a structure body with respect to the
      /// other bodys
      void determine_self_cut_position();

      /// locates the position of nodes, edges and sides for special cases
      void propagate_self_cut_position();

      /// erases sides which lies inside a structure body by locating there position
      void erase_inside_sides();

      /// erases edges which lies inside a structure body by locating there position
      void erase_inside_edges();

      /// erases nodes which lies inside a structure body by locating there position
      void erase_inside_nodes();

      /// construct the connectivity of the nodes toi find potential islands in the cut mesh
      void construct_connectivity();

      /// find the next node for the construction of the connectivity
      void next_node(Node* node, plain_node_set& remainingnodes, int count);

      /// identify the islands
      void find_islands();

      /// Get next Sides
      void next_sides(Side* cutside, Teuchos::RCP<Core::Geo::Cut::BoundingBox>& tmpbb,
          // plain_side_set allselfcutsides,
          plain_side_set& selfcutsides, plain_side_set& islandsides, bool& IsIsland);

      /*========================================================================*/
      //! @name Basic routines to run the status inspection
      /*========================================================================*/

      /// Status of the cutted sides for text viewer
      void cutted_side_status_text();

      /// Status of the cutmesh for text viewer
      void cut_mesh_status_text();

      /// Status of one problematic side for text viewer
      void error_status_text(Side& cutside);

      /// Status of the cutted sides for gmsh
      void cutted_side_status_gmsh(const std::string& name);

      /// Status of the cutmesh in gmsh
      void wall_gmsh(const std::string& name);

      /// Status of the selfcut objects in gmsh
      void sc_objects_gmsh(const std::string& name);

      /// Status of the cutmesh in gmsh for my CMGM (old)
      void s_cmgm_gmsh(const std::string& name);


      /// Status of all sides in separates files for gmsh
      void all_single_gmsh(const std::string& location);

      /// Status of one problematic side for gmsh
      void error_gmsh(const std::string& name, Side& cutside);

      /*========================================================================*/
      //! @name Auxiliary routines to run the self cut
      /*========================================================================*/

      /// cuts the edges of the first cutside with the other cutside
      void perform_self_cut(Side& cutside, Side& otherside, PointSet& selfcutpoints);

      /// gets the edges of the new created sides
      void get_self_cut_edges(Side& cutside);

      /// checks if the direction of rotation of the cycle is correct
      bool check_normal(std::vector<double> cutsideplane, Cycle& maincycle);

      /// deletes all pointers which are pointing to a soon to be erased side
      void erase_side_pointer(Side& cutside, bool erase_nodepointers);

      /// erases a side
      void erase_side(std::vector<plain_int_set>& cutsideids);

      /// erases a side which is in the unphysical part
      void erase_inside_side(std::vector<plain_int_set>& cutsideids);

      /// deletes all pointers which are pointing to a soon to be erased edge
      void erase_edge_pointer(Edge& cutsideedge);

      /// erases an edge
      void erase_edge(std::vector<plain_int_set>& cutsideedgeids);

      void plot_across();

      void point_plot_head();

      void node_plot_head();

      void edge_plot_head();

      void side_plot_head();

      void point_plot(Point& point);

      void node_plot(Node& node);

      void edge_plot(Edge& edge);

      void side_plot(Side& side);

      void side_node_plot(Side& side);


     private:
      // My Processor Id --> Only for output purposes
      int myrank_;

      /// The Cut Mesh
      Core::Geo::Cut::Mesh& mesh_;

      /// The Cut Meshhandle (Interface to the outer world)
      Core::Geo::Cut::MeshHandle& meshhandle_;

      std::map<plain_int_set, Teuchos::RCP<Side>> selfcut_sides_;

      std::map<plain_int_set, Teuchos::RCP<Edge>> selfcut_edges_;

      std::map<int, Teuchos::RCP<Node>> selfcut_nodes_;

      std::map<int, plain_node_set> selfcut_connectivity_;

      /// characteristic interface element size
      double meshsizeparam_;
    };

  }  // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
