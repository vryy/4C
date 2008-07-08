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
      
      //! nodes in tree can be leaf nodes or inner nodes
      static const int STATE_LEAF_NODE = 1;
      //! nodes in tree can be leaf nodes or inner nodes
      static const int STATE_INNER_NODE = 0;

      //! is either STATE_LEAF_NODE or STATE_INNER_NODE
      int State_;
      
      //! will be set if this node should not be refined any more.
      bool dont_refine_;
      
      //! if treenode is completly FLUID or SOLID it stores the label of the XFEM condition or 0 for fluid
      int labelID_;
      
      //! treedepth of this node
      int actTreedepth_;
      
      //! treenode has 8 children (octree)
      TreeNode* children_[8];	
      
      //! pointer to the tree
      XSearchTree* tree_; 
      
      //! axis aligned bounding box of this treeNode
      const BlitzMat3x2 AABB_;
      
      //! x-coord of the center of this treenode
      double XPlaneCoordinate_;
      
      //! y-coord of the center of this treenode
      double YPlaneCoordinate_;
      
      //! z-coord of the center of this treenode
      double ZPlaneCoordinate_;
      
      //! List of Elements belonging to this treenode (also for inner nodes!!!)
      list< const DRT::Element* > ElementList_;
      
      public:
        
        /*!
        \brief constructor of TreeNode-Object

        \param  aabb              : AxisAlignedBoundingBox of this TreeNode (=Octant)
        \param  tree              : Pointer to tree to which this Node belongs to
        */
        TreeNode(const int Depth, const BlitzMat3x2& aabb, XSearchTree* tree);
        virtual ~TreeNode();
        
        //! gives STATE_INNER_NODE or STATE_LEAF_NODE
        int getState();
        
        //! gives actTreedepth_
        int getDepth();
        
        //! gives number of nodes
        int getNodesInTree(int type);

        //! get child with index idx
        TreeNode* getChild(const int idx);
        
        //! insert Element
        void insertElement(
            const DRT::Element* elem,
            const map<int,BlitzVec3>& currentpositions
            );	
        
        /*!
         \brief Gibes candidates of surface elements in order to determine the point type of x_in, 
                gives back an empty list and labelId in lID if point type is definitly known

         \param  dis               (in)    : discretization containing elements 
         \param  currentpositions  (in)    : current position of elements in dis
         \param  pointCoords       (in)    : coords of point which will be examined
         \param  lID               (out)   : XFEM-Label if point type is definitly known, has no defined value otherwise
         \return                             list of elements which are candidates, gives back an 
                                             empty list if point type (fluid/solid) is definitly known
         */
        list< const DRT::Element* > queryPointType(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions, const BlitzVec3& pointcoords, int& lID);
        
        //! gives back ElementList_ 
        list< const DRT::Element* > getElementList();
        
        //! recusivly searches for the nearest neigbour to pointcoords (NOT YET IMPLEMENTED)
        void queryNearestNeighbour(const BlitzVec3& pointcoords);
        
        //! gives the center coords of the AABB of this treenode
        BlitzVec3 getCenterCoord();
        
        //! sets an treenode to be completly fluid (label should be set to 0)
        void setFluid(const int label);
        
        //! sets an treenode to be completly surrounded by a solid with given XFEM-coupling condidition
        void setSolid(const int label);
        
        //! gives a gmsh string to visualize the tree
        void printTree(stringstream& fc) const;
        
        //! get memory usage
        int getMemoryUsage() const;

        //! gives the AxisAlignedBoundingBox of the treenode (which of course is the same as the "strict" BB)
         const BlitzMat3x2& getAABB();
       
      private:
        
        //! gives the octand(=child) of this treenode to which elem belongs (uses the XAABB of elem)
        list<int> classifyElement(
            const DRT::Element*       elem,
            const map<int,BlitzVec3>& currentpositions
            ); 
        
        //! gives the octand(=child) of this treenode to which the point belongs 
        int classifyPoint(const BlitzVec3& pointcoords);
        
         //! gives the AABB of a child of this treenode
        BlitzMat3x2 getChildOctAABB(const int chldIdx);
      
    }; // class TreeNode
    
    private:
      
      //! maximum search depth
      static const int MAX_TREEDEPTH = 8;
      
      //! if overlap of XXAABs of 2 elements is more than this percentage of the treenodes volume, tree wont refine
      double MAX_OVERLAP_ALLOWED_;

      //! if only this number of elements is in a node they will be checked for overlaps with the treenode
      unsigned int ELEMENTS_USED_FOR_OVERLAP_CHECK_;
      
      //! true if tree has to be (re)initialized 
      bool TreeInit_;
      
      bool hasExternAABB_;
      const BlitzMat3x2 AABB_;
      
      //! count nbr. of requests to the tree
      int searchRequests_;
      
      //! count nbr. of candidates
      int Candidates_;
      
      //! how deep did we have to search to answer a request
      int SearchLength_;
      
      //! whats the minimal searchpath we had to search to answer a request
      int MinSearchLength_;
      
      //! map with all label gids indexed by their element ids
      std::map<int, int> labelByElement_;
      
      //! pointer to the root of the tree
      XFEM::XSearchTree::TreeNode* treeRoot_;
      
      //! inserts an element into the tree
      void insertElement(
          const DRT::Element* elem,
          const map<int,BlitzVec3>& currentpositions
          );
      
    public:
      XSearchTree();
      XSearchTree(const BlitzMat3x2 AABB);
      virtual ~XSearchTree();
      
      int getDepth();
      int getTotalRequests();
      double getMeanSearchlength();
      
      //! searches for the nearest neigbour to pointcoords (NOT YET IMPLEMENTED)
      void queryNearestNeighbour(const BlitzVec& pointCoords);
      
      /*!
      \brief gives back the type of a point (fluid or solid): 

      \param  dis               (in)    : discretization containing elements 
      \param  currentpositions  (in)    : current position of elements in dis
      \param  pointCoords       (in)    : coords of point which will be examined
      \return                             0->fluid 1..n->Solid with XFEM-Coupling condition
      */
      int queryPointType(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions, const BlitzVec3& pointCoords, int& closestElementId, double& distance);
      
      //! destroys and rebuilds the whole tree
      void rebuild(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions);
      
      void setMaxOverlap(double maxO);
      double getMaxOverlap();
      
      void setElementsForOverlapCheck(int maxO){
        ELEMENTS_USED_FOR_OVERLAP_CHECK_ = maxO;
      }
      int getElementsForOverlapCheck(){
        return ELEMENTS_USED_FOR_OVERLAP_CHECK_;
      }

      void printTreeMetrics(const int step) const;
      
      
      //! if the flag is set by this method, rebuild will be called within the next search request
      void setRebuildFlag();
      
      //! dumps a visualization of the tree in a gmsh file (tree0000<step>.pos)
      void printTree(const int step) const;

      //! dumps a visualization of the tree in a gmsh file (<prefiy>_tree0000<step>.pos)
      void printTree(const string prefix, const int step) const;

      
      //! get memory usage
      int getMemoryUsage() const;
     
      //! the XFEM-Coupling-Condition-Label of a element by its gid
      int getLabelByElementID(const int gid);
  }; 
  
} // namespace XFEM

#endif
#endif
