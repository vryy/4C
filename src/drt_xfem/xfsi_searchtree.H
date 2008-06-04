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
#include <bitset>
#include "../drt_lib/drt_element.H"
#include "intersection_service.H"
#include "geometry_service.H"
#include "enrichment.H"
#include "interface.H"
#include "physics.H"
#include "bele3.H"
#include "xfem.H"


namespace XFEM
{
  
  class XSearchTree {
    
    private: class TreeNode {
      
      static const int STATE_LEAF_NODE = 1;
      static const int STATE_INNER_NODE = 0;
      
      int State_;
      int labelID_;
      int actTreedepth_;
      TreeNode* children_[8];	
//      XSearchTree* tree_;
      
      BlitzMat3x2 AABB_; //axis aligned bounding box of this treeNode
      double XPlaneCoordinate_;
      double YPlaneCoordinate_;
      double ZPlaneCoordinate_;
      
      list< const DRT::Element* > ElementList_;
      
      public:
        TreeNode(int Depth, BlitzMat3x2 aabb);
        virtual ~TreeNode();
        
        int getState();
        int getDepth();
        TreeNode* getChild(const int idx);
        void insertElement(const DRT::Element* elem);	
        list< const DRT::Element* > queryPointType(const XFEM::InterfaceHandle& ih, const BlitzVec3& pointcoords, int& lID);
        list< const DRT::Element* > getElementList();
        void queryNearestNeighbour(const BlitzVec3& pointcoords);
        BlitzVec3& getCenterCoord();
        const BlitzMat3x2& getAABB();
        void setFluid(const int label);
        void setSolid(const int label);
        void printTree(stringstream& fc) const;
        
      private:
        list<int> classifyElement(const DRT::Element* elem); //give sector(s) to which elem belongs
        int classifyPoint(const BlitzVec3& pointcoords);
        BlitzMat3x2 getChildOctAABB(int chldIdx);
      
    }; // class TreeNode
    
    private:
      int tmp;
      static const int MAX_TREEDEPTH = 8;
      //	static const int AABB_Factor = 0;
      bool TreeInit_;
      int searchRequests_;
      int MeanSearchLength_;
      std::map<int,set<int> > labelsPerElementId_;
      XFEM::XSearchTree::TreeNode* treeRoot_;
      void insertElement(const DRT::Element* elem);
      void rebuild(const XFEM::InterfaceHandle& ih);
      void queryNearestNeighbour(const BlitzVec& pointCoords);
      
    public:
      XSearchTree();
      virtual ~XSearchTree();
      
      int getDepth();
      int getTotalRequests();
      int getMeanSearchlength();
      int queryPointType(const XFEM::InterfaceHandle& ih, const BlitzVec3& pointCoords);
      void setRebuildFlag();
      void printTree(const int step) const;
  } 
  static STree; // class xfsi_searchtree
  
} // namespace XFEM

#endif
#endif

//    private: class decoratedElement: public DRT::Element {
//      const DRT::Element* Element_;
//      
//      public:
//        decoratedElement(const DRT::Element* element);
//        virtual ~decoratedElement();
//        
//        int getID();
//        int getLabel();
//        list<BlitzVec3> getNodeCoords();    
//        
//      private:
//        
//    }; // class decoratedElement
//    

