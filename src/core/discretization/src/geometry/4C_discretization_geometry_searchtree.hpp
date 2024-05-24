/*----------------------------------------------------------------------*/
/*! \file

\brief provides a class with search tree with various search requests

\level 3


*----------------------------------------------------------------------*/
#ifndef FOUR_C_DISCRETIZATION_GEOMETRY_SEARCHTREE_HPP
#define FOUR_C_DISCRETIZATION_GEOMETRY_SEARCHTREE_HPP

#include "4C_config.hpp"

#include "4C_discretization_geometry_searchtree_nearestobject.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

#include <Teuchos_RCP.hpp>

#include <set>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace DRT
{
  class Discretization;
  class Element;
}  // namespace DRT

namespace CORE::LINALG
{
  class SerialDenseMatrix;
}

namespace CORE::GEO
{
  //! identifies tree type
  enum TreeType
  {
    OCTTREE,  ///< tree is a three-dimensional octtree
    QUADTREE
    ///< tree is a two-dimensional quadtree
  };

  //! identifies node types in a tree
  enum TreeNodeType
  {
    LEAF_NODE,  ///< indicates a leaf node (no further children)
    INNER_NODE
    ///< indicates an inner node (has children)
  };

  //! represents the whole data structure
  class SearchTree
  {
    //! data node for tree
   private:
    class TreeNode
    {
     private:
      //! no copy constructor and assignment operator wanted
      TreeNode(const TreeNode& old);
      TreeNode& operator=(const TreeNode& old);

      //! pointer to the parent node
      const TreeNode* const parent_;

      //! treedepth of this node
      const int treedepth_;

      //! is either STATE_LEAF_NODE or STATE_INNER_NODE
      TreeNodeType tree_node_type_;

      //! is either a OCTTREE or a QUADTREE
      const TreeType tree_type_;

      //! stores the label of the XFEM condition or 0 for fluid
      int label_;

      //! stores nearestObject
      CORE::GEO::NearestObject nearest_object_;

      //! axis aligned bounding box of this tree node
      const CORE::LINALG::Matrix<3, 2> node_box_;

      //! x-coord of the center of this treenode
      const double x_plane_coordinate_;

      //! y-coord of the center of this treenode
      const double y_plane_coordinate_;

      //! z-coord of the center of this treenode
      const double z_plane_coordinate_;

      //! treenode has 8 children (octtree) or 4 children (quadtree)
      std::vector<Teuchos::RCP<TreeNode>> children_;

      //! list of elements belonging to this treenode
      std::map<int, std::set<int>> element_list_;

      /*!
       \brief returns the node box of a child node
       \param childIndex           index of child node
       \return returns node box
       */
      CORE::LINALG::Matrix<3, 2> get_child_node_box(const int childIndex) const;

      /*!
       \brief returns the node box of a child node
       \param index           index of child node
       \param childNodeBox    child node box
       */
      void get_child_node_box(const int index, CORE::LINALG::Matrix<3, 2>& childNodeBox) const;

      /*!
       \brief returns the child node indices which overlaps with a given AABB
       \param AABB            AABB
       \param octants         vector of octantcs
       */
      void classify_xaabb(const CORE::LINALG::Matrix<3, 2>& AABB, std::vector<int>& octants) const;

      /*!
       \brief returns the child node indices which overlaps with a given AABB
       \param AABB            AABB
       \param octants         vector of octantcs
       */
      void classify_kdop(const CORE::LINALG::Matrix<9, 2>& KDOP, std::vector<int>& octants) const;

      /*!
       \brief returns the child node indices which overlaps with a given AABB
       \param AABB                AABB
       \return vector of childen indices
       */
      std::vector<int> classify_xaabb(const CORE::LINALG::Matrix<3, 2>& AABB) const;

      /*!
       \brief returns the index of the child node which overlaps with a given AABB
       \param index               index
       \param AABB                AABB
       \return vector of childen indices
       */
      bool classify_xaabb(int& index, const CORE::LINALG::Matrix<3, 2>& AABB) const;

      /*!
       \brief returns the index of the child node which overlaps with a given AABB
       \param index               index
       \param AABB                AABB
       \return vector of childen indices
       */
      bool classify_kdop(int& index, const CORE::LINALG::Matrix<9, 2>& KDOP) const;

      /*!
       \brief return child(ren) of this tree node in which the ele has to be inserted
       \param element              element
       \param currentpositions     current nodal positions in discretization
       \return vector of children ids
       */
      std::vector<int> classify_element(const DRT::Element* element,
          const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions) const;

      /*!
       \brief return child(ren) of this tree node in which the ele has to be inserted
       \param element              element
       \param xyze_element         coordinates of element
       \return vector of children ids
       */
      std::vector<int> classify_element(
          const DRT::Element* element, const CORE::LINALG::SerialDenseMatrix& xyze_element) const;

      /*!
       \brief return child(ren) of this tree node in which the ele has to be inserted
       \param element              element
       \param xyze_element         coordinates of element
       \return vector of children ids
       */
      std::vector<int> classify_element(const Teuchos::RCP<DRT::Element> element,
          const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions) const;

      /*!
       \brief return children, whose node box intersects with the circle with the given midpoint and
       radius \param radius              radius \param point               point \return child index
       */
      std::vector<int> classify_radius(
          const double radius, const CORE::LINALG::Matrix<3, 1>& point) const;

      /*!
       \brief create children of a treenode and insert elements
       \param dis                  discretization
       \param currentpositions     current nodal positions in discretization
       */
      void create_children(const DRT::Discretization& dis,
          const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions);

      /*!
       \brief create children of a treenode and insert elements
       a lot more efficient than the method above

       \param currentXAABBs        current elemental bounding boxes
       */
      void create_children(const std::map<int, CORE::LINALG::Matrix<3, 2>>& currentXAABBs);

      /*!
       \brief create children of a treenode and insert elements

       \param currentKDOPs        current elemental kdops
       */
      void create_children(const std::map<int, CORE::LINALG::Matrix<9, 2>>& currentKDOPs);


     public:
      /*!
       \brief constructor of tree node
       \param depth      depth of tree node
       \param nodeBox    node box
       \param parent     pointer to parent element
       */
      TreeNode(const TreeNode* const parent, const int depth,
          const CORE::LINALG::Matrix<3, 2>& nodeBox, const TreeType treeType);

      //! destructor
      virtual ~TreeNode() = default;

      /*!
       \brief checks, if this node has a parent node
       \return true if node has parent
       */
      bool hasParent() const
      {
        if (parent_ != nullptr)
          return true;
        else
          return false;
      };

      /*!
       \brief sets element list of a treenode
       \param elementsByLabel                elements sorted according XFEM label
       */
      void setElementList(const std::map<int, std::set<int>>& elementsByLabel);

      /*!
       \brief sets nearestObject  of a treenode
       \param nearestObject                 nearestObject
       */
      void setNearestObject(const CORE::GEO::NearestObject& nearestObject);

      /*!
       \brief returns tree node type INNER_NODE or LEAF_NODE
       \return returns tree node type INNER_NODE or LEAF_NODE
       */
      TreeNodeType getTreeNodeType() const { return tree_node_type_; };

      /*!
       \brief returns tree type OCTTREE or QUADTREE
       \return returns tree type OCTTREE or QUADTREE
       */
      TreeType getTreeType() const { return tree_type_; };

      /*!
       \brief return number of children for tree type
       \return return number of children for tree type
       */
      int getNumChildren() const;

      /*!
       \brief return pointer to the child node determined by the child index
       \param index   child node index
       \return retruns pointer to child node
       */
      Teuchos::RCP<CORE::GEO::SearchTree::TreeNode> getChild(const int index) const;

      /*!
       \brief return pointer to the parent node
       \return pointer to parent tree node
       */
      const TreeNode* getParent() const
      {
        if (this->hasParent()) return parent_;
        return nullptr;
      };

      /*!
       \brief returns elementList
       \return element list
       */
      const std::map<int, std::set<int>>& getElementList() const { return element_list_; };

      /*!
       \brief insert an element into the tree
       \param labelId              label id
       \param eleId                global ele id
       */
      void insertElement(const int labelId, const int eleId);

      /*!
       \brief returns a set of gids of nodes lying in a radius around a given point
       \param dis                  discretization
       \param currentpositions     current nodal positions in discretization
       \param point                point to be examined
       \param radius               radius
       \param label                label of structure the query point belongs to
       \return set of node gids
       */
      std::map<int, std::set<int>> search_elements_in_radius(const DRT::Discretization& dis,
          const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions,
          const CORE::LINALG::Matrix<3, 1>& point, const double radius, const int label);

      /*!
       \brief    build the static search tree for the collision detection
       \param currentBVs        map of all current AABBs for all elements
       */
      void build_static_search_tree(const std::map<int, CORE::LINALG::Matrix<3, 2>>& currentBVs);

      /*!
       \brief    build the static search tree for the collision detection
       \param currentBVs        map of all current 18-DOPs for all elements
       */
      void build_static_search_tree(const std::map<int, CORE::LINALG::Matrix<9, 2>>& currentBVs);

      /*!
       \brief returns a set of gids of elements whose AABB (bounding volume) intersect with
       the AABB of a given element (This is a query method for CONTACT-related search!)
       \param currentBVs          map of all current AABBs for all elements
       \param queryBV             18-DOP for considered element
       \param label               label ???
       \param collisions          ids of elements of overlapping current AABBs
       \return set of master contact element gids
       */
      void searchCollisions(const std::map<int, CORE::LINALG::Matrix<3, 2>>& currentBVs,
          const CORE::LINALG::Matrix<3, 2>& queryBV, const int label, std::set<int>& collisions);

      /*!
       \brief returns a set of gids of elements whose 18-DOP (bounding volume) intersect with
       the 18-DOP of a given element (This is a query method for CONTACT-related search!)
       \param currentBVs          map of all current 18-DOPs for all elements
       \param queryBV             18-DOP for considered element
       \param label               label ???
       \param collisions          ids of elements of overlapping current 18-DOPs
       \return set of master contact element gids
       */
      void searchCollisions(const std::map<int, CORE::LINALG::Matrix<9, 2>>& currentBVs,
          const CORE::LINALG::Matrix<9, 2>& queryBV, const int label, std::set<int>& collisions);
    };  // class TreeNode

    // class Tree
   public:
    /*!
     \brief tree constructor
     \param max_depth                  max_depth
     */
    SearchTree(const int max_depth  // 4 is reasonable, 5 is possible
    );

    //! destructor
    virtual ~SearchTree() = default;

    /*!
     \brief destroys the old tree if its exists and builds the root
     node of a new tree with a possibly different discretization

     \param nodebox              nodeBox
     \param elementsByLabel      elementsByLabel
     \param treetype             octtree or quadtree
     */
    void initializeTree(const CORE::LINALG::Matrix<3, 2>& nodeBox,
        const std::map<int, std::set<int>>& elementsByLabel, const TreeType treetype);

    /*!
     \brief destroys the old tree if its exists and builds the root node of a
     new tree with a possibly different discretization
     no label sorting
     \param nodebox              nodeBox
     \param dis                  discretization
     \param treetype             octtree or quadtree
     */
    void initializeTree(const CORE::LINALG::Matrix<3, 2>& nodeBox, const DRT::Discretization& dis,
        const TreeType treetype);

    void initializeTree(const CORE::LINALG::Matrix<3, 2>& nodeBox, const TreeType treetype);

    /*!
     \brief destroys the old tree if its exists and builds the root node of a
     new tree with elements; implemented for SlideALE problems where
     boundary elements are known
     \param nodebox              nodeBox
     \param elements             elements
     \param treetype             quadtree
     */
    void initialize_tree_slide_ale(const CORE::LINALG::Matrix<3, 2>& nodeBox,
        std::map<int, Teuchos::RCP<DRT::Element>>& elements, const TreeType treetype);

    void insertElement(const int eid);

    /*!
     \brief returns a set of gids of nodes lying in a radius around a given point for each object
     \param dis                  discretization
     \param currentpositions     current nodal positions in discretization
     \param point                point to be examined
     \param radius               radius
     \param label                label
     \return set of node gids
     */
    std::map<int, std::set<int>> search_elements_in_radius(const DRT::Discretization& dis,
        const std::map<int, CORE::LINALG::Matrix<3, 1>>& currentpositions,
        const CORE::LINALG::Matrix<3, 1>& point, const double radius, const int label);

    /*!
     \brief    build the static search tree for the collision detection
     \param currentBVs        map of all current AABBs for all elements
     */
    void build_static_search_tree(const std::map<int, CORE::LINALG::Matrix<3, 2>>& currentBVs);

    /*!
     \brief    build the static search tree for the collision detection
     \param currentBVs        map of all current 18-DOPs for all elements
     */
    void build_static_search_tree(const std::map<int, CORE::LINALG::Matrix<9, 2>>& currentBVs);

    /*!
     \brief returns a set of gids of elements whose AABB (bounding volume) intersect with
     the AABB of a given element (This is a query method for CONTACT-related search!)
     \param currentBVs          map of all current AABB for all elements
     \param queryBV           AABB for considered element
     \param label               label ???
     \return set of master contact element gids
     */
    void searchCollisions(const std::map<int, CORE::LINALG::Matrix<3, 2>>& currentBVs,
        const CORE::LINALG::Matrix<3, 2>& queryBV, const int label, std::set<int>& collisions);

    /*!
     \brief returns a set of gids of elements whose 18-DOP (bounding volume) intersect with
     the 18-DOP of a given element (This is a query method for CONTACT-related search!)
     \param currentKDOPs        map of all current 18-DOPs for all elements
     \param queryKDOP           18-DOP for considered element
     \param label               label ???
     \return set of master contact element gids
     */
    void searchCollisions(const std::map<int, CORE::LINALG::Matrix<9, 2>>& currentKDOPs,
        const CORE::LINALG::Matrix<9, 2>& queryKDOP, const int label, std::set<int>& contactEleIds);

   private:
    //! no copy constructor and assignment operator wanted
    SearchTree(const SearchTree& old);

    SearchTree& operator=(const SearchTree& old);

    //! maximum search depth
    const int max_depth_;

    //! pointer to the root of the tree
    Teuchos::RCP<CORE::GEO::SearchTree::TreeNode> tree_root_;
  };
  // class tree
}  // namespace CORE::GEO

FOUR_C_NAMESPACE_CLOSE

#endif
