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


#ifdef DEBUG 
#define DEBUGp(x) printf(x);
#else
#define DEBUGp(x)
#endif

#include <list>
#include <bitset>
#include <string>
#include "xfem.H"
#include "physics.H"
#include "enrichment.H"
#include "intersection_service.H"
#include "../drt_lib/drt_element.H"
#include "bele3.H"
#include <Teuchos_RefCountPtr.hpp>
#include "interface.H"

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
      XSearchTree* tree_;
      
      BlitzMat3x2 AABB_; //axis aligned bounding box of this treeNode
      double XPlaneCoordinate_;
      double YPlaneCoordinate_;
      double ZPlaneCoordinate_;
      
      list<RCP<DRT::Element> > ElementList_;
      
      public:
        TreeNode(int Depth, BlitzMat3x2 aabb, XSearchTree* tree);
        virtual ~TreeNode();
        
        int getState();
        int getDepth();
        TreeNode* getChild(const int idx);
        void insertElement(const RCP<DRT::Element> elem);	
        list<RCP<DRT::Element> > queryPointType(const BlitzVec3& pointcoords, int& lID, int& tmpL);
        list<RCP<DRT::Element> > getElementList();
        void queryNearestNeighbour(const BlitzVec3& pointcoords);
        BlitzVec3& getCenterCoord();
        const BlitzMat3x2& getAABB();
        void setFluid(const int label);
        void setSolid(const int label);
        void printTree(stringstream& fc) const;
        
      private:
        list<int> classifyElement(const RCP<DRT::Element> elem); //give sector(s) to which elem belongs
        int classifyPoint(const BlitzVec3& pointcoords);
        bool isContainedIn(const RCP<DRT::Element> elem);
        list<RCP<DRT::Element> > listFactory(const RCP<DRT::Element> elem);
        BlitzMat3x2 getChildOctAABB(int chldIdx);
        
        inline void classifyAndInsert(const RCP<DRT::Element> elem){		
          list<int> childIdx = classifyElement(elem);
          for (list<int>::const_iterator myIt = childIdx.begin(); myIt != childIdx.end(); myIt++){
            this->children_[*myIt-1]->insertElement(elem);
            BlitzMat3x2 ab = this->children_[*myIt-1]->getAABB();
            //					printf("inserted elem to AABB(%f,%f,%f,%f,%f,%f)\n", ab(0,0),ab(0,1),ab(1,0),ab(1,1),ab(2,0),ab(2,1));
          }
        }
        
    }; // class TreeNode
    
    private: class decoratedElement: public DRT::Element {
      const DRT::Element* Element_;
      
      public:
        decoratedElement(const DRT::Element* element);
        virtual ~decoratedElement();
        
        int getID();
        int getLabel();
        list<BlitzVec3> getNodeCoords(); 	
        
      private:
        
    }; // class decoratedElement
    
    private:
      int tmp;
      static const int MAX_TREEDEPTH = 8;
      //	static const int AABB_Factor = 0;
      bool TreeInit_;
      int searchRequests_;
      int MeanSearchLength_;
      XFEM::XSearchTree::TreeNode* treeRoot_;
      const XFEM::InterfaceHandle* ih_;
      std::map<int,set<int> > labelsPerElementId_;
      void insertElement(const RCP<DRT::Element> elem);
      void rebuild();
      BlitzMat3x2 getXAABBofDis();
      void queryNearestNeighbour(const BlitzVec& pointCoords);
      const DRT::Element* nearestNeighbourInList(const list<RCP<DRT::Element> > ElementList, const BlitzVec3& x, double& dist, int& tmpL);
      
    public:
      XSearchTree();
      virtual ~XSearchTree();
      
      int getDepth();
      int getTotalRequests();
      int getMeanSearchlength();
      //void elementDelete();
      int queryPointType(const XFEM::InterfaceHandle& ih, const BlitzVec3& pointCoords);
      void setRebuildFlag();
      void printTree(const int step) const;
  } 
  static STree; // class xfsi_searchtree
  
} // namespace XFEM

#endif
#endif
