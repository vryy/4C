/*----------------------------------------------------------------------*/
/*! \file

\brief provides a class with search tree with various search requests

\level 3


*----------------------------------------------------------------------*/
#ifndef SEARCHTREE_H
#define SEARCHTREE_H

#include <set>

#include <Teuchos_RCP.hpp>
#include "searchtree_nearestobject.H"
#include "../linalg/linalg_fixedsizematrix.H"

// forward declarations
namespace DRT
{
  class Discretization;
  class Element;
}  // namespace DRT

namespace LINALG
{
  class SerialDenseMatrix;
}

namespace GEO
{
  class InterfacePoint;

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
      TreeNodeType treeNodeType_;

      //! is either a OCTTREE or a QUADTREE
      const TreeType treeType_;

      //! stores the label of the XFEM condition or 0 for fluid
      int label_;

      //! stores nearestObject
      GEO::NearestObject nearestObject_;

      //! axis aligned bounding box of this tree node
      const LINALG::Matrix<3, 2> nodeBox_;

      //! x-coord of the center of this treenode
      const double xPlaneCoordinate_;

      //! y-coord of the center of this treenode
      const double yPlaneCoordinate_;

      //! z-coord of the center of this treenode
      const double zPlaneCoordinate_;

      //! treenode has 8 children (octtree) or 4 children (quadtree)
      std::vector<Teuchos::RCP<TreeNode>> children_;

      //! list of elements belonging to this treenode
      std::map<int, std::set<int>> elementList_;

      /*!
       \brief returns the node box of a child node
       \param childIndex           index of child node
       \return returns node box
       */
      LINALG::Matrix<3, 2> getChildNodeBox(const int childIndex) const;

      /*!
       \brief returns the node box of a child node
       \param index           index of child node
       \param childNodeBox    child node box
       */
      void getChildNodeBox(const int index, LINALG::Matrix<3, 2>& childNodeBox) const;

      /*!
       \brief returns the child node index, in which a point lies
       \param point                point
       \return index
       */
      int classifyPoint(const LINALG::Matrix<3, 1>& point) const;

      /*!
       \brief returns the child node indices which overlaps with a given AABB
       \param AABB            AABB
       \param octants         vector of octantcs
       */
      void classifyXAABB(const LINALG::Matrix<3, 2>& AABB, std::vector<int>& octants) const;

      /*!
       \brief returns the child node indices which overlaps with a given AABB
       \param AABB            AABB
       \param octants         vector of octantcs
       */
      void classifyKDOP(const LINALG::Matrix<9, 2>& KDOP, std::vector<int>& octants) const;

      /*!
       \brief returns the child node indices which overlaps with a given AABB
       \param AABB                AABB
       \return vector of childen indices
       */
      std::vector<int> classifyXAABB(const LINALG::Matrix<3, 2>& AABB) const;

      /*!
       \brief returns the index of the child node which overlaps with a given AABB
       \param index               index
       \param AABB                AABB
       \return vector of childen indices
       */
      bool classifyXAABB(int& index, const LINALG::Matrix<3, 2>& AABB) const;

      /*!
       \brief returns the index of the child node which overlaps with a given AABB
       \param index               index
       \param AABB                AABB
       \return vector of childen indices
       */
      bool classifyKDOP(int& index, const LINALG::Matrix<9, 2>& KDOP) const;

      /*!
       \brief return child(ren) of this tree node in which the ele has to be inserted
       \param element              element
       \param currentpositions     current nodal positions in discretization
       \return vector of children ids
       */
      std::vector<int> classifyElement(const DRT::Element* element,
          const std::map<int, LINALG::Matrix<3, 1>>& currentpositions) const;

      /*!
       \brief return child(ren) of this tree node in which the ele has to be inserted
       \param element              element
       \param xyze_element         coordinates of element
       \return vector of children ids
       */
      std::vector<int> classifyElement(
          const DRT::Element* element, const LINALG::SerialDenseMatrix& xyze_element) const;

      /*!
       \brief return child(ren) of this tree node in which the ele has to be inserted
       \param element              element
       \param xyze_element         coordinates of element
       \return vector of children ids
       */
      std::vector<int> classifyElement(const Teuchos::RCP<DRT::Element> element,
          const std::map<int, LINALG::Matrix<3, 1>>& currentpositions) const;

      /*!
       \brief return children, whose node box intersects with the circle with the given midpoint and
       radius \param radius              radius \param point               point \return child index
       */
      std::vector<int> classifyRadius(const double radius, const LINALG::Matrix<3, 1>& point) const;

      /*!
       \brief create children of a treenode and insert elements
       \param dis                  discretization
       \param currentpositions     current nodal positions in discretization
       */
      void createChildren(const DRT::Discretization& dis,
          const std::map<int, LINALG::Matrix<3, 1>>& currentpositions);

      /*!
       \brief create children of a treenode and insert elements
       \param masterelements      all given masterelements
       \param currentpositions    current nodal positions of masterelements
       */
      void createChildren(std::map<int, Teuchos::RCP<DRT::Element>>& masterelements,
          const std::map<int, LINALG::Matrix<3, 1>>& currentpositions);

      /*!
       \brief create children of a treenode and insert tree elements
       \param currentpositions     current nodal positions in discretization
       */
      void createChildren(const std::map<int, LINALG::Matrix<3, 1>>& currentpositions);

      /*!
       \brief create children of a treenode and insert elements
       a lot more efficient than the method above

       \param currentXAABBs        current elemental bounding boxes
       */
      void createChildren(const std::map<int, LINALG::Matrix<3, 2>>& currentXAABBs);

      /*!
       \brief create children of a treenode and insert elements

       \param currentKDOPs        current elemental kdops
       */
      void createChildren(const std::map<int, LINALG::Matrix<9, 2>>& currentKDOPs);

      /*!
       \brief set XFEM label of empty children
       \param dis                  discretization
       \param currentpositions     current nodal positions in discretization
       */
      void setXFEMLabelOfEmptyChildren(const DRT::Discretization& dis,
          const std::map<int, LINALG::Matrix<3, 1>>& currentpositions);

      /*!
       \brief set XFEM label of empty children
       \param triangleList         triangleList
       \param currentpositions     current nodal positions in discretization
       */
      void setXFEMLabelOfEmptyChildren(const std::vector<std::vector<int>>& triangleList,
          const std::vector<GEO::InterfacePoint>& pointList);

      /*!
       \brief set XFEM label and nearest object of empty children
       \param dis                  discretization
       \param currentpositions     current nodal positions in discretization
       */
      void setXFEMLabelAndNearestObjectOfEmptyChildren(const DRT::Discretization& dis,
          const std::map<int, LINALG::Matrix<3, 1>>& currentpositions);

     public:
      /*!
       \brief constructor of tree node
       \param depth      depth of tree node
       \param nodeBox    node box
       \param parent     pointer to parent element
       */
      TreeNode(const TreeNode* const parent, const int depth, const LINALG::Matrix<3, 2>& nodeBox,
          const TreeType treeType);

      //! destructor
      virtual ~TreeNode();

      /*!
       \brief clears node stet label to -1, treenodetype to LEAF_NODE and deletes children
       \return true if node has parent
       */
      void clear();

      /*!
       \brief checks, if this node has a parent node
       \return true if node has parent
       */
      bool hasParent() const
      {
        if (parent_ != NULL)
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
       \brief sets label of a treenode
       \param label                label
       */
      void setLabel(const int label);

      /*!
       \brief sets nearestObject  of a treenode
       \param nearestObject                 nearestObject
       */
      void setNearestObject(const GEO::NearestObject& nearestObject);

      /*!
       \brief returns tree node type INNER_NODE or LEAF_NODE
       \return returns tree node type INNER_NODE or LEAF_NODE
       */
      TreeNodeType getTreeNodeType() const { return treeNodeType_; };

      /*!
       \brief returns tree type OCTTREE or QUADTREE
       \return returns tree type OCTTREE or QUADTREE
       */
      TreeType getTreeType() const { return treeType_; };

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
      const Teuchos::RCP<GEO::SearchTree::TreeNode> getChild(const int index) const;

      /*!
       \brief return pointer to the parent node
       \return pointer to parent tree node
       */
      const TreeNode* getParent() const
      {
        if (this->hasParent()) return parent_;
        return NULL;
      };

      /*!
       \brief returns the center of this tree node
       \return returns the center of this tree node
       */
      LINALG::Matrix<3, 1> getCenterCoord() const;

      /*!
       \brief returns the node box of the tree node
       \return returns the node box of the tree node
       */
      const LINALG::Matrix<3, 2>& getNodeBox() const { return nodeBox_; };

      /*!
       \brief returns elementList
       \return element list
       */
      const std::map<int, std::set<int>>& getElementList() const { return elementList_; };

      /*!
       \brief insert an element into the tree
       \param labelId              label id
       \param eleId                global ele id
       */
      void insertElement(const int labelId, const int eleId);

      /*!
       \brief updates tree node
       \param AABB_old             AABB_old
       \param AABB_new             AABB_new
       */
      void updateTreeNode(
          const LINALG::Matrix<3, 2>& AABB_old, const LINALG::Matrix<3, 2>& AABB_new);

      /*!
       \brief returns the label 0 -> fluid 1 -> n solid
       \param dis                  discretization
       \param currentpositions     current nodal positions in discretization
       \param point                point to be examined
       \param nearestObject        nearestobject
       \return 0->fluid 1..n->Solid with XFEM-Coupling condition
       */
      int queryFSINearestObject(const DRT::Discretization& dis,
          const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
          const std::map<int, LINALG::Matrix<3, 2>>& currentXAABBs,
          const LINALG::Matrix<3, 1>& point, GEO::NearestObject& nearestobject);

      /*!
       \brief returns the label 0 -> fluid 1 -> n solid
       \param dis                  discretization
       \param currentpositions     current nodal positions in discretization
       \param point                point to be examined
       \return 0->fluid 1..n->Solid with XFEM-Coupling condition
       */
      int queryXFEMFSIPointType(const DRT::Discretization& dis,
          const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
          const std::map<int, LINALG::Matrix<3, 2>>& currentXAABBs,
          const LINALG::Matrix<3, 1>& point);

      /*!
       \brief fi intersection for contact with XFEM FSI
       \param triangleList                   discretization
       \param interfacePointList             list of interface interfacePointList
       \param currentpositions               current nodal positions in discretization
       \param currentXAABBs                  current XAABBs
       \param querypoint                     point to be examined
       \param querypointId                   point id in triangle list
       \return 0->fluid 1..n->Solid with XFEM-Coupling condition
       */
      void moveContactNodes(const std::vector<std::vector<int>>& triangleList,
          std::vector<GEO::InterfacePoint>& pointList,
          const std::map<int, LINALG::Matrix<3, 2>>& triangleXAABBs,
          const LINALG::Matrix<3, 1>& querypoint, const int querypointId,
          const int querypointLabel);

      /*!
       \brief returns a set of gids of nodes lying in a radius around a given point
       \param dis                  discretization
       \param currentpositions     current nodal positions in discretization
       \param point                point to be examined
       \param radius               radius
       \param label                label of structure the query point belongs to
       \return set of node gids
       */
      std::map<int, std::set<int>> searchElementsInRadius(const DRT::Discretization& dis,
          const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
          const LINALG::Matrix<3, 1>& point, const double radius, const int label);

      /*!
       \brief returns a vector of Ids of points closer than radius to the querypoint
       \param currentpositions     current positions of a number of points assigned to a unique Id
       \param querypoint           point around which should be searched
       \param radius               radius
       \return set of node gids
       */
      std::vector<int> searchPointsInRadius(
          const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
          const LINALG::Matrix<3, 1>& querypoint,  // point around which should be searched
          const double radius);

      /*!
       \brief   returns a set of elements whose XAABB s are
       intersecting with the XAABB of a given volume element

       \param element              given slave contact element
       \param xyze_element         coordinates of element
       \return set of node gids
       */
      void queryIntersectionCandidates(const std::map<int, LINALG::Matrix<3, 2>>& currentXAABBs,
          const LINALG::Matrix<3, 2>& eleXAABB, std::set<int>& elementset);

      /*!
       \brief   returns a set of elements whose XAABB s are
       intersecting with the XAABB of a given volume element and whose labels are different

       \param element              given slave contact element
       \param xyze_element         coordinates of element
       \param label                label of element to be examined
       \return set of node gids
       */
      void queryPotentialElements(const std::map<int, LINALG::Matrix<3, 2>>& currentXAABBs,
          const LINALG::Matrix<3, 2>& eleXAABB, std::map<int, std::set<int>>& elementset,
          const int label);

      /*!
       \brief   returns a vector of elements whose XAABB s are
       intersecting with the XAABB of a given volume element
       \param dis                  discretization
       \param currentpositions     current nodal positions in discretization
       \return set of node gids
       */
      void queryPotentialElements_Approx2(const DRT::Discretization& dis,
          const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
          const LINALG::Matrix<3, 2>& eleXAABB,
          const std::vector<LINALG::Matrix<3, 1>>& gaussPoints,
          std::map<int, std::map<int, GEO::NearestObject>>& potentialObjectsAtGP,
          const double cutoff, const int label, const int projectiontype);

      /*!
       \brief    build the static search tree for the collision detection
       \param currentBVs        map of all current AABBs for all elements
       */
      void buildStaticSearchTree(const std::map<int, LINALG::Matrix<3, 2>>& currentBVs);

      /*!
       \brief    build the static search tree for the collision detection
       \param currentBVs        map of all current 18-DOPs for all elements
       */
      void buildStaticSearchTree(const std::map<int, LINALG::Matrix<9, 2>>& currentBVs);

      /*!
       \brief returns a set of gids of elements whose AABB (bounding volume) intersect with
       the AABB of a given element (This is a query method for CONTACT-related search!)
       \param currentBVs          map of all current AABBs for all elements
       \param queryBV             18-DOP for considered element
       \param label               label ???
       \param collisions          ids of elements of overlapping current AABBs
       \return set of master contact element gids
       */
      void searchCollisions(const std::map<int, LINALG::Matrix<3, 2>>& currentBVs,
          const LINALG::Matrix<3, 2>& queryBV, const int label, std::set<int>& collisions);

      /*!
       \brief returns a set of gids of elements whose 18-DOP (bounding volume) intersect with
       the 18-DOP of a given element (This is a query method for CONTACT-related search!)
       \param currentBVs          map of all current 18-DOPs for all elements
       \param queryBV             18-DOP for considered element
       \param label               label ???
       \param collisions          ids of elements of overlapping current 18-DOPs
       \return set of master contact element gids
       */
      void searchCollisions(const std::map<int, LINALG::Matrix<9, 2>>& currentBVs,
          const LINALG::Matrix<9, 2>& queryBV, const int label, std::set<int>& collisions);

      /*!
       \brief print tree node
       \param max_depth           maximum depth of tree
       \param fc                  prefix node stream
       */
      void printTreeNode(const int max_depth, std::stringstream& fc) const;

      /*!
       \brief returns depth of subtree from the  tree node it is called (METRICS)
       \return tree depth
       */
      int getDepth() const;

      /*!
       \brief returns number of tree nodes in tree (METRICS)
       \return tree nodes
       */
      int getNumNodesInTree() const;

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
    virtual ~SearchTree();

    //! returns tree root node
    Teuchos::RCP<GEO::SearchTree::TreeNode> getRoot() const { return treeRoot_; };

    /*!
     \brief destroys the old tree if its exists and builds the root
     node of a new tree with a possibly different discretization

     \param nodebox              nodeBox
     \param elementsByLabel      elementsByLabel
     \param treetype             octtree or quadtree
     */
    void initializeTree(const LINALG::Matrix<3, 2>& nodeBox,
        const std::map<int, std::set<int>>& elementsByLabel, const TreeType treetype);

    /*!
     \brief destroys the old tree if its exists and builds the root node of a
     new tree with a possibly different discretization
     no label sorting
     \param nodebox              nodeBox
     \param dis                  discretization
     \param treetype             octtree or quadtree
     */
    void initializeTree(const LINALG::Matrix<3, 2>& nodeBox, const DRT::Discretization& dis,
        const TreeType treetype);

    void initializeTree(const LINALG::Matrix<3, 2>& nodeBox, const TreeType treetype);

    /*!
     \brief destroys the old tree if its exists and builds the root node of a
     new tree with a possibly different discretization;no label sorting;
     tree elements are not finite elements in discretization dis, but
     nodes in discretization dis
     \param nodebox              nodeBox
     \param dis                  discretization
     \param treetype             octtree or quadtree
     */
    void initializePointTree(const LINALG::Matrix<3, 2>& nodeBox,
        const std::map<int, LINALG::Matrix<3, 1>>& currentpositions, const TreeType treetype);

    /*!
     \brief destroys the old tree if its exists and builds the root node of a
     new tree with elements; implemented for SlideALE problems where
     boundary elements are known
     \param nodebox              nodeBox
     \param elements             elements
     \param treetype             quadtree
     */
    void initializeTreeSlideALE(const LINALG::Matrix<3, 2>& nodeBox,
        std::map<int, Teuchos::RCP<DRT::Element>>& elements, const TreeType treetype);

    void insertElement(const int eid);

    /*!
     \brief updates tree for labeled structures due to movements of these structures.
     For networks use initializeTree()
     \param dis                      discretization
     \param currentpositions_old     old current nodal positions in discretization
     \param currentpositions_new     new current nodal positions in discretization
     */
    void updateTree(const DRT::Discretization& dis,
        const std::map<int, LINALG::Matrix<3, 1>>& currentpositions_old,
        const std::map<int, LINALG::Matrix<3, 1>>& currentpositions_new);

    /*!
     \brief returns the label 0 -> fluid 1 -> n solid
     \param dis                  discretization
     \param currentpositions     current nodal positions in discretization
     \param point                point to be examined
     \param nearestObject        nearestobject
     \return 0->fluid 1..n->Solid with XFEM-Coupling condition
     */
    int queryFSINearestObject(const DRT::Discretization& dis,
        const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
        const std::map<int, LINALG::Matrix<3, 2>>& currentXAABBs, const LINALG::Matrix<3, 1>& point,
        GEO::NearestObject& nearestobject);

    /*!
     \brief returns the label 0 -> fluid 1 -> n solid
     \param dis                  discretization
     \param currentpositions     current nodal positions in discretization
     \param point                point to be examined
     \return 0->fluid 1..n->Solid with XFEM-Coupling condition
     */
    int queryXFEMFSIPointType(const DRT::Discretization& dis,
        const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
        const std::map<int, LINALG::Matrix<3, 2>>& currentXAABBs,
        const LINALG::Matrix<3, 1>& point);

    /*!
     \brief fix intersection for contact with XFEM FSI
     \param triangleList                   discretization
     \param interfacePointList             list of interface interfacePointList
     \param currentpositions               current nodal positions in discretization
     \param currentXAABBs                  current XAABBs
     \param querypoint                     point to be examined
     \param querypointId                   point id in triangle list
     */
    void moveContactNodes(const std::vector<std::vector<int>>& triangleList,
        std::vector<GEO::InterfacePoint>& pointList,
        const std::map<int, LINALG::Matrix<3, 2>>& triangleXAABBs,
        const LINALG::Matrix<3, 1>& querypoint, const int querypointId, const int querypointLabel);

    /*!
     \brief returns a set of gids of nodes lying in a radius around a given point for each object
     \param dis                  discretization
     \param currentpositions     current nodal positions in discretization
     \param point                point to be examined
     \param radius               radius
     \param label                label
     \return set of node gids
     */
    std::map<int, std::set<int>> searchElementsInRadius(const DRT::Discretization& dis,
        const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
        const LINALG::Matrix<3, 1>& point, const double radius, const int label);

    /*!
     \brief returns a vector of Ids of points in a radius around a given querypoint
     \param currentpositions     current positions of a number of points assigned to a unique Id
     \param querypoint           point around which should be searched
     \param radius               radius
     \return set of point Ids
     */
    std::vector<int> searchPointsInRadius(
        const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
        const LINALG::Matrix<3, 1>& querypoint,  // point around which should be searched
        const double radius);

    /*!
     \brief   returns a vector of elements whose XAABB s are
     intersecting with the XAABB of a given volume element
     \param dis                  discretization
     \param currentpositions     current nodal positions in discretization
     \return set of node gids
     */
    void queryIntersectionCandidates(const std::map<int, LINALG::Matrix<3, 2>>& currentXAABBs,
        const std::vector<LINALG::Matrix<3, 2>>& structure_AABBs,
        const LINALG::Matrix<3, 2>& eleXAABB, std::set<int>& elementset);

    /*!
     \brief   returns a vector of elements whose XAABB s are
     intersecting with the XAABB of a given volume element
     \param dis                  discretization
     \param currentpositions     current nodal positions in discretization
     \return set of node gids
     */
    void queryPotentialElements(const std::map<int, LINALG::Matrix<3, 2>>& currentXAABBs,
        const LINALG::Matrix<3, 2>& eleXAABB, std::map<int, std::set<int>>& elementset,
        const int label);

    /*!
     \brief   returns a vector of elements whose XAABB s are
     intersecting with the XAABB of a given volume element
     \param dis                  discretization
     \param currentpositions     current nodal positions in discretization
     \return set of node gids
     */
    void queryPotentialElements_Approx2(const DRT::Discretization& dis,
        const std::map<int, LINALG::Matrix<3, 1>>& currentpositions,
        const LINALG::Matrix<3, 2>& eleXAABB, const std::vector<LINALG::Matrix<3, 1>>& gaussPoints,
        std::map<int, std::map<int, GEO::NearestObject>>& potentialObjectsAtGP, const double cutoff,
        const int label, const int projectiontype);

    /*!
     \brief    build the static search tree for the collision detection
     \param currentBVs        map of all current AABBs for all elements
     */
    void buildStaticSearchTree(const std::map<int, LINALG::Matrix<3, 2>>& currentBVs);

    /*!
     \brief    build the static search tree for the collision detection
     \param currentBVs        map of all current 18-DOPs for all elements
     */
    void buildStaticSearchTree(const std::map<int, LINALG::Matrix<9, 2>>& currentBVs);

    /*!
     \brief returns a set of gids of elements whose AABB (bounding volume) intersect with
     the AABB of a given element (This is a query method for CONTACT-related search!)
     \param currentBVs          map of all current AABB for all elements
     \param queryBV           AABB for considered element
     \param label               label ???
     \return set of master contact element gids
     */
    void searchCollisions(const std::map<int, LINALG::Matrix<3, 2>>& currentBVs,
        const LINALG::Matrix<3, 2>& queryBV, const int label, std::set<int>& collisions);

    /*!
     \brief returns a set of gids of elements whose 18-DOP (bounding volume) intersect with
     the 18-DOP of a given element (This is a query method for CONTACT-related search!)
     \param currentKDOPs        map of all current 18-DOPs for all elements
     \param queryKDOP           18-DOP for considered element
     \param label               label ???
     \return set of master contact element gids
     */
    void searchCollisions(const std::map<int, LINALG::Matrix<9, 2>>& currentKDOPs,
        const LINALG::Matrix<9, 2>& queryKDOP, const int label, std::set<int>& contactEleIds);

    /*!
     \brief prints entire tree to gmsh file
     \param prefix               prefix
     \param step                 time step
     */
    void printTree(const std::string prefix, const int step) const;

    /*!
     \brief evaluates and prints tree metrics
     \param step                 time step
     */
    void evaluateTreeMetrics(const int step) const;

   private:
    //! no copy constructor and assignment operator wanted
    SearchTree(const SearchTree& old);

    SearchTree& operator=(const SearchTree& old);

    //! maximum search depth
    const int max_depth_;

    //! pointer to the root of the tree
    Teuchos::RCP<GEO::SearchTree::TreeNode> treeRoot_;
  };
  // class tree
}  // namespace GEO

#endif
