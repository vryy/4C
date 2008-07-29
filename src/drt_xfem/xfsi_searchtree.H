/*!
\file xfsi_searchtree.H

\brief provides a class with search tree

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
 */
#ifdef CCADISCRET

#ifndef XFSI_SEARCHTREE_H
#define XFSI_SEARCHTREE_H

#include <list>
#include "geometry_service.H"
#include "xfem_condition.H"

namespace XFEM
{
  //! represents the whole data structure
  class XSearchTree {
    
    //! data node for tree
    private: class TreeNode {
      public: 
      //! identifies node types in a tree
      enum TreeNodeType {
        LEAF_NODE,    ///< indicates a leaf node (no further children)
        INNER_NODE,   ///< indicates an inner node (has children)
        ANY_NODE      ///< this is only used for recursive node counting
        };
      
      private:
      //! is either STATE_LEAF_NODE or STATE_INNER_NODE
      TreeNodeType TreeNodeType_;
      
      //! will be set if this node should not be refined any more.
      bool dont_refine_;
      
      //! if treenode is completly FLUID or SOLID it stores the label of the XFEM condition or 0 for fluid
      int labelID_;
      
      //! treedepth of this node
      const int actTreedepth_;
      
      //! treenode has 8 children (octree)
      TreeNode* children_[8];	
      
      //! pointer to the tree
      XSearchTree* tree_; 
      
      //! pointer to the tree
      const TreeNode* const parent_; 
      
      //! axis aligned bounding box of this treeNode
      const BlitzMat3x2 AABB_;
      
      //! x-coord of the center of this treenode
      const double XPlaneCoordinate_;
      
      //! y-coord of the center of this treenode
      const double YPlaneCoordinate_;
      
      //! z-coord of the center of this treenode
      const double ZPlaneCoordinate_;
      
      //! List of Elements belonging to this treenode (also for inner nodes!!!)
      list< const DRT::Element* > ElementList_;
      
      public:     
        //! standard constructor
        TreeNode();
        
        /*!
        \brief constructor of TreeNode-Object
        */
        TreeNode(
            const int Depth,          ///<depth of this node
            const BlitzMat3x2& aabb,  ///<AxisAlignedBoundingBox of this TreeNode (=Octant)
            XSearchTree* tree,        ///<Pointer to the tree this node belongs to
            const TreeNode* const parent          ///<pointer to the parent node
            );
        
        //! destructor
        virtual ~TreeNode();
        
        //! gives INNER_NODE or LEAF_NODE
        TreeNodeType getTreeNodeType() const;
        
        //! gives actTreedepth_
        int getDepth() const;
        
        //! gives number of nodes
        int getNodesInTree(
            const TreeNodeType type   ///< specifies which nodes are counted
            ) const;

        //! get child with index idx
        TreeNode* getChild(
            const int idx           ///< index of the child (1...8)
            ) const;

        //! get parent node
        const TreeNode* const getParent() const;

        //! gives the center coords of the AABB of this treenode
        BlitzVec3 getCenterCoord() const;

        //! gives a gmsh string to visualize the tree
        void printTree(
            stringstream& fc      ///< output is written to this stream
            ) const;
        
        //! get memory usage in byte
        int getMemoryUsage() const;

        //! gives the AxisAlignedBoundingBox of the treenode (which of course is the same as the "strict" BB)
        const BlitzMat3x2& getAABB() const;
        
        //! test if this node has a parent node (==treeRoot)
        bool hasParent() const;
        
        //! insert an element into the tree
        void insertElement(
            const DRT::Element* elem,                       ///< point element to insert
            const map<int,BlitzVec3>& currentpositions      ///< current position of element nodes
            );	
        
        /*!
         \brief Gives back the XFEM-label of the point X
         \return                             list of elements which are candidates, gives back an 
                                             empty list if point type (fluid/solid) is definitly known
         */
        int queryPointType(
            const DRT::Discretization& dis,                   ///<discretization containing elements 
            const std::map<int,BlitzVec3>& currentpositions,  ///< current position of elements in dis
            const BlitzVec3& X,                               ///<coords of point which will be examined
            int& closestElementID                             ///>XFEM-Label if point type is definitly known, has no defined value otherwise
            );
        
        //! gives back ElementList_ 
        list< const DRT::Element* > getElementList() const;
        
        //! recusivly searches for the nearest neighbour to pointcoords (NOT YET IMPLEMENTED)
        void queryNearestNeighbour(
            const BlitzVec3& pointcoords      ///< coord of point you want to know the nearest neighbours
            );       
        
        //! sets a treenode to be completly of one label (e.g. 0 for fluid in xfem-fsi-problem)
        void setLabel(
            const int label     ///< labelId you want to set
            );
               
  
      private:
        
        //! gives the octant(=child) of this treenode to which elem belongs (uses the XAABB of elem)
        list<int> classifyElement(
            const DRT::Element*       elem,
            const map<int,BlitzVec3>& currentpositions
            ) const; 
        
        //! gives the octant(=child-Idx) of this treenode to which AABB belongs
        list<int> classifyAABB(
             const BlitzMat3x2 AABB
            ) const;
        
        //! gives the octand(=child) of this treenode to which the point belongs 
        int classifyPoint(
            const BlitzVec3& pointcoords
            ) const;

        //! determines the XFEM-Label of an Point. The point should be located in the TreeNode calling the function
        int getXFEMLabelOfPointInTreeNode(
            const DRT::Discretization&          dis,   ///< discretization containing elements    
            const std::map<int,BlitzVec3>&      currentpositions, ///< current position of elements in dis
            const BlitzVec3&                    X     ///< coordinates of point to check      
            );
        
         //! gives the AABB of a child of this treenode
        BlitzMat3x2 getChildOctAABB(
            const int chldIdx                       ///< child index
            ) const;
        
        //! check XFEM-Label for TreeNodes which have no Elements
        void setLabelsForCandidateFreeChildren(
            const DRT::Discretization&          dis, ///< discretization containing elements 
            const std::map<int,BlitzVec3>&      currentpositions ///< current position of elements in dis
            );
        
        //! create children of a treenode and classify elements to these children 
        void createChildren(
            const std::map<int,BlitzVec3>& currentpositions///< current position of elements in dis
            );

        //! check percentage of overlap to decide if tree will be refined or not
        void doOverlapCheck(
            const std::map<int,BlitzVec3>& currentpositions///< current position of elements in dis
            );

      
    }; // class TreeNode
    
    private:
      //! maximum search depth
      const int MAX_TREEDEPTH_;
      
      //! if overlap of XXAABs of 2 elements is more than this percentage of the treenodes volume, tree wont refine
      const double MAX_OVERLAP_ALLOWED_;

      //! if only this number of elements is in a node they will be checked for overlaps with the treenode
      const unsigned int ELEMENTS_USED_FOR_OVERLAP_CHECK_;
      
      //! map with all label gids indexed by their element ids
      std::map<int, int> labelByElement_;
      
      //! bounding box of the whole tree
      const BlitzMat3x2 AABB_;      
      
      //! true bounding box of tree was given by an constructor argument
      const bool hasExternAABB_;
      
      //! true if tree has to be (re)initialized 
      bool TreeInit_;
      
      //! count nbr. of requests to the tree
      int searchRequests_;
      
      //! count nbr. of candidates
      int Candidates_;
      
      //! how deep did we have to search to answer a request
      int SearchLength_;
      
      //! whats the minimal searchpath we had to search to answer a request
      int MinSearchLength_;      
      
      //! pointer to the root of the tree
      XFEM::XSearchTree::TreeNode* treeRoot_;
      
      //! nodes which had to search in more than their own TreeNode
      int rootLevelSearches_;
      
      //! inserts an element into the tree
      void insertElement(
          const DRT::Element* elem,
          const map<int,BlitzVec3>& currentpositions
          );
      
    public:
      //! standard tree constructor
      XSearchTree();
      
      //! tree constructor defining an external bounding box, shoudl be used
      XSearchTree(
          const BlitzMat3x2 AABB
          );
      
      //! destructor
      virtual ~XSearchTree();
      
      //! gives actual tree depth
      int getDepth() const;
      
      //! gives number of total requests set on the tree
      int getTotalRequests() const;
      
      //! gives the mean path-length for a request
      double getMeanSearchlength() const;
      
      //! searches for the nearest neigbour to pointcoords (NOT YET IMPLEMENTED)
      void queryNearestNeighbour(
          const BlitzVec& pointCoords  ///< coords of point which will be examined
          );
      
      /*!
      \brief gives back the type of a point (fluid or solid): 
      \return                             0->fluid 1..n->Solid with XFEM-Coupling condition
      */
      int queryPointType(
          const DRT::Discretization& dis,  ///< discretization containing elements 
          const std::map<int,BlitzVec3>& currentpositions, ///< current position of elements in dis
          const BlitzVec3& pointCoords, ///< coords of point which will be examined
          int& closestElementId, ///> id of closest element
          double& distance ///> distance of closest element (<0 solid, >=0 fluid) 
          );
      
      //! destroys and rebuilds the whole tree
      void rebuild(
          const DRT::Discretization& dis, ///< discretization containing elements 
          const std::map<int,BlitzVec3>& currentpositions ///< current position of elements in dis
          );
      
      //! set maximum of overlap allowed for elements in a tree node, node wont build child nodes (=refine) this value is exeeded
      void setMaxOverlap(
          double maxO ///< value
          );

      //! get actual value of maximum allowed overlap
      double getMaxOverlap() const;
      
      //! get the maximum number of elements for which a overlap check will be performed      
      int getElementsForOverlapCheck() const{
        return ELEMENTS_USED_FOR_OVERLAP_CHECK_;
      }

      //! print tree metrics like requests, time, mean path length, memory usage, ...
      void printTreeMetrics(
          const int step ///< file postfix
          ) const;
      
      //! if the flag is set by this method, rebuild will be called within the next search request
      void setRebuildFlag();
      
      //! dumps a visualization of the tree in a gmsh file (tree0000<step>.pos)
      void printTree(
          const int step ///< file postfix 
          ) const;

      //! dumps a visualization of the tree in a gmsh file (<prefix>_tree0000<step>.pos)
      void printTree(
          const string prefix, ///< file prefix
          const int step ///< file postfix
          ) const;
     
      //! get memory usage
      int getMemoryUsage() const;
     
      //! the XFEM-Coupling-Condition-Label of a element by its gid
      int getLabelByElementID(
          const int gid ///< global elements id
          ) const;
      
  }; 
} // namespace XFEM


#endif
#endif
