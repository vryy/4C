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

class xfsi_searchtree {
	
	private: class TreeNode {
	        //! I don't have children
		static const int STATE_LEAF_NODE = 1;
		//! I have children
		static const int STATE_INNER_NODE = 0;

		int State;
		int labelID;
		int actTreedepth;
		TreeNode* children[8];	
		xfsi_searchtree* tree;
		
		BlitzMat3x2 AABB; //axis aligned bounding box of this treeNode
		double XPlaneCoordinate;
		double YPlaneCoordinate;
		double ZPlaneCoordinate;
		
		list<RCP<DRT::Element> > ElementList;

		public:
			TreeNode(int Depth, BlitzMat3x2 aabb, xfsi_searchtree* tree);
			virtual ~TreeNode();
			
			int getState();
			int getDepth();
			TreeNode* getChild(const int idx);
			void insertElement(const RCP<DRT::Element> elem);	
			list<RCP<DRT::Element> > queryPointType(const BlitzVec3& pointcoords, int& lID);
			list<RCP<DRT::Element> > getElementList();
			void queryNearestNeighbour(const BlitzVec3& pointcoords);
			BlitzVec3& getCenterCoord();
			const BlitzMat3x2& getAABB();
			void setFluid();
			void setSolid(int label);
			void printTree();
						
		private:
		    list<int> classifyElement(const RCP<DRT::Element> elem); //give sector(s) to which elem belongs
		    int classifyPoint(const BlitzVec3& pointcoords);
			bool isContainedIn(const RCP<DRT::Element> elem);
			list<RCP<DRT::Element> > listFactory(const RCP<DRT::Element> elem);
			BlitzMat3x2 getChildOctAABB(int chldIdx);
			
			inline void classifyAndInsert(const RCP<DRT::Element> elem){		
				list<int> childIdx = classifyElement(elem);
				for (list<int>::const_iterator myIt = childIdx.begin(); myIt != childIdx.end(); myIt++){
					this->children[*myIt-1]->insertElement(elem);
					BlitzMat3x2 ab = this->children[*myIt-1]->getAABB();
					printf("inserted elem to AABB(%f,%f,%f,%f,%f,%f)\n", ab(0,0),ab(0,1),ab(1,0),ab(1,1),ab(2,0),ab(2,1));
				}
			}
			
	}; // class TreeNode

private:
	static const int MAX_TREEDEPTH = 4;
	static const int AABB_Factor = 5;
	bool TreeInit;
	int searchRequests;
	int MeanSearchLength;
 	XFEM::xfsi_searchtree::TreeNode* treeRoot;
 	const XFEM::InterfaceHandle* ih;
 	std::map<int,set<int> > labelsPerElementId;
	void insertElement(const RCP<DRT::Element> elem);
	void rebuild();
 	BlitzMat3x2 getXAABBofDis();
	void queryNearestNeighbour(const BlitzVec& pointCoords);
	const DRT::Element* nearestNeighbourInList(const list<RCP<DRT::Element> > ElementList, const BlitzVec3& x, double& dist);
	
public:
	xfsi_searchtree();
	virtual ~xfsi_searchtree();

 	int getDepth();
 	int getTotalRequests();
 	int getMeanSearchlength();
	//void elementDelete();
	int queryPointType(const XFEM::InterfaceHandle& ih, const BlitzVec3& pointCoords);
	void setRebuildFlag();
	void printTree();
} 
static STree; // class xfsi_searchtree

} // namespace XFEM

#endif
#endif
