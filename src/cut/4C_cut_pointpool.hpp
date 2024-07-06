/*---------------------------------------------------------------------*/
/*! \file

\brief PointPool, stores a points in the cut and decides if points are merged or new points are
created

\level 3


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_POINTPOOL_HPP
#define FOUR_C_CUT_POINTPOOL_HPP

#include "4C_config.hpp"

#include "4C_cut_boundingbox.hpp"
#include "4C_cut_point.hpp"

#include <Teuchos_RCP.hpp>

#include <set>

FOUR_C_NAMESPACE_OPEN

namespace Core::Geo
{
  namespace Cut
  {
    class Point;
    class Edge;
    class Side;

    ///< Specifies the merge Strategy of the Pointpool
    enum PointpoolMergeStrategy
    {
      InitialLoad,
      SelfCutLoad,
      NormalCutLoad
    };

    /// OctTreeNode class that implements a single node of an OctTree
    /*!
     * Might need some speedup. (Uli)
     * this class is an OcTree-based pool that stores all points in tree-structure
     * Each node of this tree is represented by a bounding box and a RCPPointSet (a set of pointer
     * to nodes) that lie within this bounding box, and has 8 children nodes. Pointers of nodes are
     * contained on more levels of this tree, also the leafs of this tree is a node with the
     * additional info IsLeaf (no children) Points between 2 or more bounding boxes are strictly
     * related to one of these boxes
     */
    class OctTreeNode
    {
     public:
      /// empty constructor
      OctTreeNode(double norm) : bb_(Teuchos::rcp(BoundingBox::create())), norm_(norm) {}

      /*========================================================================*/
      //! @name create and get points
      /*========================================================================*/

      /// If a point with the coordinates "x" does not exists, it creates a new point
      /// correspondingly
      Point* new_point(const double* x, Edge* cut_edge, Side* cut_side, double tolerance,
          PointpoolMergeStrategy merge_strategy);

      /// Get the point with the specified coordinates "x" from the pointpool
      Point* get_point(const double* x, Edge* cut_edge, Side* cut_side, double tolerance,
          PointpoolMergeStrategy merge_strategy);

      /// Simply insert p into the pointpool and correspondingly modify the boundingbox size
      void add_point(const double* x, Teuchos::RCP<Point> p);


      /*========================================================================*/
      //! @name get info about the tree
      /*========================================================================*/

      /// how many points are contained in the current bounding box
      unsigned size() const { return points_.size(); }

      /// is the current tree-node a leaf of the tree?
      bool is_leaf() const { return nodes_[0] == Teuchos::null; }

      /*========================================================================*/
      //! @name collect adjacent objects
      /*========================================================================*/

      /// collect all edges
      void collect_edges(const BoundingBox& edgebox, plain_edge_set& edges);

      /// collect all sides
      void collect_sides(const BoundingBox& sidebox, plain_side_set& sides);

      /*!
      \brief collect all elements near the sidebox
             (unused, does not work properly when there is no point adjacent to elements in a tree's
      leaf, e.g. when side lies within an element)
      */
      void collect_elements(const BoundingBox& sidebox, plain_element_set& elements);


      /*========================================================================*/
      //! @name print routine and others
      /*========================================================================*/

      /// reset the Point::Position of outside points
      void reset_outside_points();

      /// print the tree at a given level
      void print(int level, std::ostream& stream);

      /// Get Points
      const RCPPointSet& get_points() { return points_; }

      // get the node where point with specified coordinate lies
      OctTreeNode* find_node(const double* coord, Point* p);

     private:
      /*========================================================================*/
      //! @name private routines
      /*========================================================================*/

      /// Create a point with the specified ID
      Teuchos::RCP<Point> create_point(
          unsigned newid, const double* x, Edge* cut_edge, Side* cut_side, double tolerance);

      /// get the leaf where the point with the given coordinates lies in
      OctTreeNode* leaf(const double* x);

      /// split the current boounding box (tree-node)
      void split(int level);


      /*========================================================================*/
      //! @name private class variables
      /*========================================================================*/

      /// children of the current OctTreeNode
      Teuchos::RCP<OctTreeNode> nodes_[8];

      /// set of points contained in the bounding box related to this tree-node
      RCPPointSet points_;

      /// bounding box related to this tree-node
      Teuchos::RCP<BoundingBox> bb_;

      /// point coordinates where the current bounding box is split into 8 children bboxes
      Core::LinAlg::Matrix<3, 1> splitpoint_;

      /// norm
      double norm_;
    };



    /*!
    \brief OctTree based storage container (the PointPool class)
     */
    class PointPool
    {
     public:
      /// constructor
      explicit PointPool(double norm = 1);

      /// create a new point in the point PointPool
      Point* new_point(const double* x, Edge* cut_edge, Side* cut_side, double tolerance)
      {
        switch (probdim_)
        {
          case 2:
          {
            // extend the point dimension
            double x_ext[3] = {0.0, 0.0, 0.0};
            std::copy(x, x + 2, x_ext);

            return tree_.new_point(
                x_ext, cut_edge, cut_side, get_tolerance(x_ext, tolerance), merge_strategy_);
          }
          case 3:
            return tree_.new_point(
                x, cut_edge, cut_side, get_tolerance(x, tolerance), merge_strategy_);
          default:
            FOUR_C_THROW("Unsupported problem dimension! (probdim = %d)", probdim_);
            exit(EXIT_FAILURE);
        }

        FOUR_C_THROW("Impossible to reach this line!");
        exit(EXIT_FAILURE);
      }

      /// get the pointer to a stored point if possible
      Point* get_point(const double* x, Edge* cut_edge, Side* cut_side, double tolerance)
      {
        return tree_.get_point(x, cut_edge, cut_side, get_tolerance(x, tolerance), merge_strategy_);
      }

      /// get default tolerance for the case that the tolerance of a point is set to 0!!!
      double get_tolerance(const double* x, double tolerance)
      {
        if (tolerance == 0.0)
        {
          Core::LinAlg::Matrix<3, 1> px(x);
          return (px.norm_inf() * POSITIONTOL);
        }
        else
          return tolerance;
      }

      /// get the size of the OctTree (total number of stored points)
      unsigned size() const { return tree_.size(); }

      /// Collects the edges that are fully or partially contained in the edgebox
      void collect_edges(const BoundingBox& edgebox, plain_edge_set& edges)
      {
        tree_.collect_edges(edgebox, edges);
      }

      /// Collects the sides that are fully or partially contained in the sidebox
      void collect_sides(const BoundingBox& sidebox, plain_side_set& sides)
      {
        tree_.collect_sides(sidebox, sides);
      }

      /// Collects the elements that are fully or partially contained in the sidebox
      void collect_elements(const BoundingBox& sidebox, plain_element_set& elements)
      {
        tree_.collect_elements(sidebox, elements);
      }

      /// reset Point::Position of all stored outside points
      void reset_outside_points() { tree_.reset_outside_points(); }

      /// print the OctTree
      void print(std::ostream& stream) { tree_.print(0, stream); }

      /// Get Points
      const RCPPointSet& get_points() { return tree_.get_points(); }

      void set_merge_strategy(PointpoolMergeStrategy strategy) { merge_strategy_ = strategy; }

     private:
      ///< Current Merge Strategy of the Pointpool
      PointpoolMergeStrategy merge_strategy_;

      ///< the PointPool class is represented by the first node in the OctTree
      OctTreeNode tree_;

      /// problem dimension
      const unsigned probdim_;
    };
  }  // namespace Cut
}  // namespace Core::Geo

FOUR_C_NAMESPACE_CLOSE

#endif
