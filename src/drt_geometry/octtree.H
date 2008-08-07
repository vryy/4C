/*!
\file octtree.H

\brief provides a class with oct tree with various search requests

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
 */
#ifdef CCADISCRET

#ifndef XFSI_OCTTREE_H
#define XFSI_OCTTREE_H

#include "../drt_geometry/octtree_geometry_service.H"

namespace GEO
{

//! identifies node types in a tree
enum TreeNodeType {
    LEAF_NODE,    ///< indicates a leaf node (no further children)
    INNER_NODE    ///< indicates an inner node (has children)
  };


//! represents the whole data structure
class OctTree{
  
  
  //! data node for tree
private: class TreeNode{
     
private:

  //! pointer to the parent node
  const TreeNode* const parent_; 
    
  //! treedepth of this node
  const int treedepth_;
  
  //! is either STATE_LEAF_NODE or STATE_INNER_NODE
  TreeNodeType treeNodeType_;

  //! stores the label of the XFEM condition or 0 for fluid
  int label_;
  
  //! axis aligned bounding box of this tree node 
  const BlitzMat3x2 nodeBox_;
    
  //! x-coord of the center of this treenode
  const double xPlaneCoordinate_;

  //! y-coord of the center of this treenode
  const double yPlaneCoordinate_;

  //! z-coord of the center of this treenode
  const double zPlaneCoordinate_;

  //! treenode has 8 children (octree)
  std::vector<Teuchos::RCP<TreeNode> > children_;  
    
  //! list of elements belonging to this treenode
  std::map<int, std::set<int> > elementList_;


  /*!
   \brief returns the node box of a child node
   \param childIndex           index of child node
   \return returns node box
   */
  BlitzMat3x2 getChildNodeBox(
      const int childIndex                      
      ) const;

  /*!
   \brief returns the child node index, in which a point lies
   \param point                point
   \return index
   */ 
  const int classifyPoint(
      const BlitzVec3&    point
      ) const;


  /*!
   \brief returns the child node indices which overlaps with a given AABB
   \param AABB                AABB
   \return vector of childen indices
   */
  std::vector<int> classifyXAABB(
      const BlitzMat3x2&    AABB
      ) const;


  /*!
   \brief return child(ren) of this tree node in which the ele has to be inserted
   \param element              element
   \param currentpositions     current nodal positions in discretization
   \return vector of children ids 
   */
  std::vector<int> classifyElement(
      const DRT::Element*               element,
      const std::map< int,BlitzVec3 >&  currentpositions
      ) const; 


  /*!
   \brief create children of a treenode and insert elements
   \param dis                  discretization
   \param currentpositions     current nodal positions in discretization
   */
  void createChildren(
      const DRT::Discretization&          dis,    
      const std::map<int, BlitzVec3 >&    currentpositions
      );
  
  
public:
  /*!
   \brief constructor of tree node
   \param depth      depth of tree node
   \param nodeBox    node box
   \param parent     pointer to parent element
   */
  TreeNode(
      const TreeNode* const parent,
      const int depth,              
      const BlitzMat3x2& nodeBox
      );

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
  bool hasParent() const;

  
  /*!
   \brief sets element list of a treenode 
   \param elementsByLabel                elements sorted according XFEM label
   */
  void setElementList(
      const std::map<int, std::set<int> >& elementsByLabel
      );


  /*!
   \brief sets label of a treenode 
   \param label                label
   */
  void setLabel(
      const int label     
      );


  /*!
   \brief returns tree node type INNER_NODE or LEAF_NODE
   \return returns tree node type INNER_NODE or LEAF_NODE
   */
  const TreeNodeType getTreeNodeType() const;


  /*!
   \brief return pointer to the child node determined by the child index 
   \param index   child node index
   \return retruns pointer to child node
   */
  const Teuchos::RCP<TreeNode> getChild(
      const int index           
      ) const;


  /*!
   \brief return pointer to the parent node 
   \return pointer to parent tree node
   */
  const TreeNode* const getParent() const;


  /*!
   \brief returns the center of this tree node
   \return returns the center of this tree node
   */
  const BlitzVec3 getCenterCoord() const;


  /*!
   \brief returns the node box of the tree node 
   \return returns the node box of the tree node 
   */
  const BlitzMat3x2& getNodeBox() const;


  /*!
   \brief returns elementList
   \return element list
   */
  const std::map<int,std::set<int> > getElementList() const;


  /*!
   \brief insert an element into the tree
   \param labelId              label id
   \param eleId                global ele id
   */
  void insertElement(
      const int   labelId,
      const int   eleId                           
      );  

  
  /*!
   \brief updates tree node
   \param AABB_old             AABB_old
   \param AABB_new             AABB_new
   */
  void updateTreeNode(
      const BlitzMat3x2& AABB_old, 
      const BlitzMat3x2& AABB_new
      );

  
  /*!
  \brief returns the label 0 -> fluid 1 -> n solid
  \param dis                  discretization
  \param currentpositions     current nodal positions in discretization
  \param point                point to be examined
  \return 0->fluid 1..n->Solid with XFEM-Coupling condition
   */
  int queryXFEMFSIPointType(
      const DRT::Discretization&        dis,                   
      const std::map<int,BlitzVec3>&    currentpositions,  
      const BlitzVec3&                  point                               
      );
  
  
  /*!
  \brief print tree node
  \param max_depth           maximum depth of tree
  \param fc                  prefix node stream
   */
  void printTreeNode(
      const int       max_depth,
      stringstream&   fc
      ) const;
  

}; // class TreeNode
  
  

  // class Tree
public:

  //! tree constructor 
  OctTree(
      const int           max_depth,
      const BlitzMat3x2&  nodeBox
      );

  //! destructor
  virtual ~OctTree();


  /*!
  \brief destroys the old tree if its exists and builds the root 
         node of a new tree with a possibly different discretization

  \param elementsByLabel      elementsByLabel
   */
  void initializeTree(
      const std::map<int,std::set<int> >& elementsByLabel  	
      );


  /*!
  \brief destroys the old tree if its exists and builds the root node of a 
         new tree with a possibly different discretization
         no label sorting
  \param dis                  discretization
   */
  void initializeTree(
      const DRT::Discretization&        dis
      );


  /*!
  \brief updates tree for labeled structures due to movements of these structures.
         For networks use initializeTree()
  \param dis                      discretization
  \param currentpositions_old     old current nodal positions in discretization
  \param currentpositions_new     new current nodal positions in discretization
   */
  void updateTree(
      const DRT::Discretization& 	      dis,
      const std::map<int,BlitzVec3>& 	  currentpositions_old, 
      const std::map<int,BlitzVec3>& 	  currentpositions_new
  ); 


  /*!
  \brief returns the label 0 -> fluid 1 -> n solid
  \param dis                  discretization
  \param currentpositions     current nodal positions in discretization
  \param point                point to be examined
  \return 0->fluid 1..n->Solid with XFEM-Coupling condition
   */
  int queryXFEMFSIPointType(
      const DRT::Discretization&      dis,  			
      const std::map<int,BlitzVec3>&  currentpositions, 
      const BlitzVec3&                point
      );

  
  /*!
  \brief prints entire tree to gmsh file
  \param prefix               prefix
  \param step                 time step
   */
  void printTree(
      const string  prefix, 
      const int     step
      ) const;
  
  

private:
  //! maximum search depth
  const int max_depth_;

  //! bounding box of the whole tree
  const BlitzMat3x2 rootNodeBox_;      

  //! pointer to the root of the tree
  Teuchos::RCP< GEO::OctTree::TreeNode > treeRoot_;

}; 
} // namespace GEO


#endif
#endif
