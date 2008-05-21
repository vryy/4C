/*!
\file xfsi_searchtree.cpp

\brief provides a class with search tree

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
#ifdef CCADISCRET
#include "xfsi_searchtree.H"
#include <typeinfo>
#include <string>

using namespace std;
using namespace XFEM;


XFEM::xfsi_searchtree::xfsi_searchtree() {
	this->TreeInit = false;
	ih = NULL;
	treeRoot=NULL;
}

XFEM::xfsi_searchtree::~xfsi_searchtree() {
	delete treeRoot;
}

void XFEM::xfsi_searchtree::setRebuildFlag() {
	TreeInit = false;
}

int XFEM::xfsi_searchtree::getTotalRequests(){
	return searchRequests;
}
int XFEM::xfsi_searchtree::getMeanSearchlength(){
	return MeanSearchLength;
}

int XFEM::xfsi_searchtree::queryPointType(const XFEM::InterfaceHandle& ih, const BlitzVec3& pointCoords) {
	searchRequests++;
	this->ih = &ih;
	if (!TreeInit)
		rebuild();
	// search candidates in tree
	list<RCP<DRT::Element> > candidates;
	int labID;
	candidates = treeRoot->queryPointType(pointCoords, labID);
	if (candidates.empty()){
		return labID;
	}
	// do excat routines
	int within = 0;
	list<RCP<DRT::Element> >::iterator myIt = candidates.begin();
	double dist = 0;
	const DRT::Element* closestEle =this->nearestNeighbourInList(candidates, pointCoords,dist);
	if (dist<0)
		within = 1;
	return within;
}

int XFEM::xfsi_searchtree::getDepth(){
	return treeRoot->getDepth();
}

void XFEM::xfsi_searchtree::insertElement(const RCP<DRT::Element> elem) {
	treeRoot->insertElement(elem);
}

void XFEM::xfsi_searchtree::rebuild() {
	// loop all elements on this processor
	if (treeRoot != NULL)
		delete treeRoot;
	
//	this->ih = ih;
	labelsPerElementId.clear();
	  for(std::map<int,set<int> >::const_iterator conditer = (*ih->elementsByLabel()).begin(); conditer!=(*ih->elementsByLabel()).end(); ++conditer)
	  {
	    for(std::set<int>::const_iterator eleid = conditer->second.begin(); eleid!=conditer->second.end(); ++eleid)
	    {
	      labelsPerElementId[*eleid].insert(conditer->first);
	    }
	  }
	
	
	//TODO: remove 2nd loop over data: getXAABBofDis also loops over elements of dis
	BlitzMat3x2 aabb =getXAABBofDis();
	
	// we want to cover the whole fluid domain, but we dont know its dimensions, 
	// therefor its estimated by expanding the bounding box of all solids by a
	// constant factor;
	aabb(0,0) = aabb(0,0)-AABB_Factor*abs(aabb(0,0));
	aabb(0,1) = aabb(0,1)+AABB_Factor*abs(aabb(0,1));
	aabb(1,0) = aabb(1,0)-AABB_Factor*abs(aabb(1,0));
	aabb(1,1) = aabb(1,1)+AABB_Factor*abs(aabb(1,1));
	aabb(2,0) = aabb(2,0)-AABB_Factor*abs(aabb(2,0));
	aabb(2,1) = aabb(2,1)+AABB_Factor*abs(aabb(2,1));

	treeRoot = new TreeNode(0,aabb,this);
	for (int i=0; i< ih->cutterdis()->NumMyRowElements(); ++i) {
		DEBUGp("insert dis-elem\n")
		RCP<DRT::Element> ele =  Teuchos::rcp(ih->cutterdis()->lRowElement(i));
		this->insertElement(ele);
	}
    DEBUGp("DEBUG_xfsi_searchtree: rebuilt treeRoot\n");
	TreeInit = true;

}

XFEM::xfsi_searchtree::TreeNode::TreeNode(int Depth,BlitzMat3x2 aabb, xfsi_searchtree* tree)
 {
	actTreedepth = Depth;
	this->AABB=aabb;
	XPlaneCoordinate = aabb(0,0)+(aabb(0,1)-aabb(0,0))/2.0;
	YPlaneCoordinate = aabb(1,0)+(aabb(1,1)-aabb(1,0))/2.0;
	ZPlaneCoordinate = aabb(2,0)+(aabb(2,1)-aabb(2,0))/2.0;
	this->State = STATE_LEAF_NODE;
	for (int i=0; i< 8; i++){
		children[i]=NULL;
	}
	this->tree = tree;
}

XFEM::xfsi_searchtree::TreeNode::~TreeNode() {
	for (int i=0; i< 8; i++){
		if (children[i]!=NULL) 
			delete children[i];
	}
}

int XFEM::xfsi_searchtree::TreeNode::getDepth(){
	if (State==STATE_LEAF_NODE){
		return actTreedepth;
	}
	else {
		int depth=actTreedepth;
		int tmp=0;
		for (int i=0; i< 8; i++){
			if (children[i] != NULL)
				tmp = children[i]->getDepth();
			else
				tmp =0;
			if (tmp>depth)
				depth=tmp;
		}
		return depth;
	}
	return actTreedepth;
}

void XFEM::xfsi_searchtree::TreeNode::setFluid(){
	this->labelID = 0;
}
void XFEM::xfsi_searchtree::TreeNode::setSolid(int label){
	this->labelID=label;
}

BlitzVec3& XFEM::xfsi_searchtree::TreeNode::getCenterCoord(){
	BlitzVec3* t = new BlitzVec3(this->XPlaneCoordinate, this->YPlaneCoordinate, this->ZPlaneCoordinate);
	return *t;
}

list<RCP<DRT::Element> > XFEM::xfsi_searchtree::TreeNode::queryPointType(const BlitzVec3& pointCoords, int& lID) {
//	printf("AABB(%f\t%f\t%f\t%f\t%f\t%f)\t", AABB(0,0),AABB(0,1),AABB(1,0),AABB(1,1),AABB(2,0),AABB(2,1));
//	printf("x_in(%f\t%f\t%f)\n",pointCoords(0), pointCoords(1),pointCoords(2));
	switch (this->State) {
	case STATE_INNER_NODE:				 // if inner node, do recursion 
//		printf("*");
//		printf("classified searchpoint to oct %d, so i will search there\n", classifyPoint(pointCoords)-1);
		return children[classifyPoint(pointCoords)-1]->queryPointType(pointCoords ,lID);
		break;
	case STATE_LEAF_NODE:	 // returns list of candidates
//		printf("/");
		if (actTreedepth >= MAX_TREEDEPTH){
//			printf("actTreedepth >= MAX_TREEDEPTH\n");
			lID = this->labelID;
			return ElementList;  // just return (maybe empty) list if leaf is at max depth, 
								 // or if no candidates (==fluid)
			}
		if (ElementList.empty()){
//			printf("empty ");
			lID = this->labelID;
			return ElementList;  // just return (maybe empty) list if leaf is at max depth, 
								 // or if no candidates (==fluid)
			}
		if (ElementList.size()>1) // dynamically grow tree
			{
			printf("more than 1 element in leaf -> growing tree\n");
			// create Octants
			for (int i=0; i<8; i++){
//				DEBUGp("creating child\n")
				BlitzMat3x2 chldAABB = getChildOctAABB(i+1);
				children[i] = new TreeNode(actTreedepth + 1, chldAABB, tree);
			}
			printf("created child nodes\n");
			// actual node becomes an inner tree node,
			// so we have to introduce one more tree-level
			for (list<RCP<DRT::Element> >::const_iterator myIt = this->ElementList.begin(); myIt != this->ElementList.end(); myIt++){
				classifyAndInsert(*myIt);
			}

			// if one of the created childs is empty, check if it is fluid or solid
			for (int i=0; i< 8; i++){
				if ((children[i]->getElementList()).empty()){
					double dista=0;
					BlitzVec3 c= children[i]->getCenterCoord();
					printf("centerX (%f,%f,%f)\n", c(0),c(1),c(2)); 
					const DRT::Element* closestEle = tree->nearestNeighbourInList(ElementList, children[i]->getCenterCoord(),dista);
//					printf("distance of empty octant center to next element is %f\n", dista);
					if (dista<0) {
						printf("solid octant with label %d\n", 1);
						int eleId = closestEle->Id();
						int l = tree->labelsPerElementId.find(eleId)->first;
						children[i]->setSolid(l);
					}
					else {
						printf("fluid octant\n");
						children[i]->setFluid();
					}
					// search for nearest Neigbour in ElementList and
					// check if child is fluid or solid
					// Element NN = NearestNeighbour(ElementList, child[i].center);
					// child[i] = fluid/solid
				}
				else {
					printf("non empty octant(#%d) with %d elements\n", i+1, children[i]->getElementList().size());
				}
			}
			// this node becomes an inner tree node
			this->State = STATE_INNER_NODE;
			this->ElementList.clear();
			// do recursion
			int l=0;
//			printf("classified searchpoint to oct %d, so i will search there", classifyPoint(pointCoords)-1);
			list<RCP<DRT::Element> > myL = children[classifyPoint(pointCoords)-1]->queryPointType(pointCoords, l);
			lID = l;
			return myL;
		}
		else  // if there is only one Element, just return it
			printf("returning candidate");
			return ElementList;
		break;
	}
	dserror("should not get here\n");
	return ElementList;
}



const DRT::Element* XFEM::xfsi_searchtree::nearestNeighbourInList(const list<RCP<DRT::Element> > ElementList, const BlitzVec3& x_in, double& dist)
{
    bool in_element = false;
    double min_ele_distance = 1.0e12;
    double distance = 1.0e12;
    const DRT::Element* closest_element;
    
    for (list<RCP<DRT::Element> >::const_iterator myIt = ElementList.begin(); myIt != ElementList.end(); myIt++)
    {
      const DRT::Element* cutterele = &**myIt;
      const BlitzMat xyze_cutter(getCurrentNodalPositions(cutterele, *ih->currentcutterpositions()));
      static BlitzVec2 eleCoord;
      static BlitzVec3 normal;
      in_element = XFEM::searchForNearestPointOnSurface(cutterele,xyze_cutter,x_in,eleCoord,normal,distance);
      if (abs(distance) < abs(min_ele_distance))
       {
         closest_element = cutterele;
         min_ele_distance = distance;
       }
    }
    dist = min_ele_distance;
    return closest_element;
}


void XFEM::xfsi_searchtree::TreeNode::insertElement(const RCP<DRT::Element> elem) {
	if ((actTreedepth >= XFEM::xfsi_searchtree::MAX_TREEDEPTH) || (this->State == STATE_LEAF_NODE) ) {
		this->ElementList.push_back(elem);
		printf("inserting element\n");
		this->State = STATE_LEAF_NODE;
	} else if(this->State == STATE_INNER_NODE) {
		classifyAndInsert(elem);
	} 
}

list<RCP<DRT::Element> > XFEM::xfsi_searchtree::TreeNode::getElementList(){
	return ElementList;
}

BlitzMat3x2 XFEM::xfsi_searchtree::TreeNode::getChildOctAABB(int i){
	BlitzMat3x2 chldAABB;
	if (i>4){
		chldAABB(0,0) = XPlaneCoordinate;
		chldAABB(0,1) = AABB(0,1);
	}
	else {
		chldAABB(0,0) = this->AABB(0,0);
		chldAABB(0,1) = XPlaneCoordinate;
	}
	if ((i==3) || (i==4) || (i==7) || (i==8)){
		chldAABB(1,0) = YPlaneCoordinate;
		chldAABB(1,1) = AABB(1,1);
	}
	else {
		chldAABB(1,0) = AABB(1,0);
		chldAABB(1,1) = YPlaneCoordinate;
	}
	if ((i%2)==0){
		chldAABB(2,0) = ZPlaneCoordinate;
		chldAABB(2,1) = AABB(2,1);
	}
	else {
		chldAABB(2,0) = AABB(2,0);
		chldAABB(2,1) = ZPlaneCoordinate;
	}		
	printf("created chldAABB(%f\t%f\t%f\t%f\t%f\t%f)\n", chldAABB(0,0),chldAABB(0,1),chldAABB(1,0),chldAABB(1,1),chldAABB(2,0),chldAABB(2,1));
	return chldAABB;
	
}


int XFEM::xfsi_searchtree::TreeNode::classifyPoint(const BlitzVec3& pointcoords) {
	int octIdx = 1;
	if (pointcoords(0) > XPlaneCoordinate)
		octIdx = octIdx + 4;
	if (pointcoords(1) > YPlaneCoordinate)
		octIdx = octIdx + 2;
	if (pointcoords(2) > ZPlaneCoordinate)
		octIdx = octIdx + 1;
	return octIdx;
}

list<int> XFEM::xfsi_searchtree::TreeNode::classifyElement(const RCP<DRT::Element> elem) {
	list<int> octants;
    const BlitzMat xyze(DRT::UTILS::PositionArrayBlitz(&*elem));
	const BlitzMat3x2 elemXAABB= XFEM::computeFastXAABB(&*elem, xyze);
	//bitset<8> PlaneBitSet; // every Bit stands for a plane
	//PlaneBitSet.reset();


	//TODO: eleganter lösen
/*	BlitzVec ful = BlitzVec(elemXAABB(0, 0), elemXAABB(1, 0), elemXAABB(2, 1));
	BlitzVec fur = BlitzVec(elemXAABB(0, 0), elemXAABB(1, 1), elemXAABB(1, 1));
	BlitzVec fll = BlitzVec(elemXAABB(0, 0), elemXAABB(1, 0), elemXAABB(1, 0));
	BlitzVec flr = BlitzVec(elemXAABB(0, 0), elemXAABB(1, 1), elemXAABB(1, 0));
	BlitzVec rul = BlitzVec(elemXAABB(0, 1), elemXAABB(1, 0), elemXAABB(2, 1));
	BlitzVec rur = BlitzVec(elemXAABB(0, 1), elemXAABB(1, 1), elemXAABB(1, 1));
	BlitzVec rll = BlitzVec(elemXAABB(0, 1), elemXAABB(1, 0), elemXAABB(1, 0));
	BlitzVec rlr = BlitzVec(elemXAABB(0, 1), elemXAABB(1, 1), elemXAABB(1, 0));*/

	if (elemXAABB(0, 1) > XPlaneCoordinate) {
		if (elemXAABB(1, 1) > YPlaneCoordinate) {
			if (elemXAABB(2, 1) > ZPlaneCoordinate)
			octants.push_back(8);
			if (elemXAABB(2, 0) <= ZPlaneCoordinate)
			octants.push_back(7);
		}
		if (elemXAABB(1, 0) <= YPlaneCoordinate) {
			if (elemXAABB(2, 1) > ZPlaneCoordinate)
			octants.push_back(6);
			if (elemXAABB(2, 0) <= ZPlaneCoordinate)
			octants.push_back(5);
		}
	}
	if (elemXAABB(0, 0) <= XPlaneCoordinate) {
		if (elemXAABB(1, 1) > YPlaneCoordinate) {
			if (elemXAABB(2, 1) > ZPlaneCoordinate)
				octants.push_back(4);
			if (elemXAABB(2, 0) <= ZPlaneCoordinate)
				octants.push_back(3);
		}
		if (elemXAABB(1, 0) <= YPlaneCoordinate) {
			if (elemXAABB(2, 1) > ZPlaneCoordinate)
				octants.push_back(2);
			if (elemXAABB(2, 0) <= ZPlaneCoordinate)
				octants.push_back(1);
		}

	}
	printf("classifiing elem(%f,%f,%f,%f,%f,%f)", elemXAABB(0,0),elemXAABB(0,1),elemXAABB(1,0),elemXAABB(1,1),elemXAABB(2,0),elemXAABB(2,1));
	printf("classified element to octants: ");
    for (list<int>::const_iterator myIt = octants.begin(); myIt != octants.end(); myIt++)
    {
	printf("%d ", *myIt);
    }
	printf("\n");
	return octants;

}

int XFEM::xfsi_searchtree::TreeNode::getState() {
	return State;
}

void XFEM::xfsi_searchtree::printTree() {
	printf("printing tree\n");
	this->treeRoot->printTree();
}

void XFEM::xfsi_searchtree::TreeNode::printTree() {
	string prefix("");
//	for (int i=0; i<ActTreedepth;i++) prefix=prefix+"\t";
//	printf(prefix + "act_depth: %d\n",actTreedepth);
//	if(this->State==STATE_INNER_NODE) printf((prefix+"*\n").c_str());
//	if(this->State==STATE_LEAF_NODE){
//		printf((prefix+"/").c_str());
//		printf("(%d)\n", this->ElementList.size());
//	}
	for (int j=0; j<8; j++){
		if (children[j]!=NULL)
			children[j]->printTree();
	}
}

const BlitzMat3x2& XFEM::xfsi_searchtree::TreeNode::getAABB(){
	return this->AABB;
}

XFEM::xfsi_searchtree::TreeNode* XFEM::xfsi_searchtree::TreeNode::getChild(int idx) {
	return children[idx-1];
}

bool XFEM::xfsi_searchtree::TreeNode::isContainedIn(const RCP<DRT::Element> element){
	// not yet implemented, even not necessary, because expensive 
	// routines are done anyway on candidates found by the octree
	return false;
}

BlitzMat3x2 XFEM::xfsi_searchtree::getXAABBofDis(){
	const int nsd = 3;
    BlitzMat3x2 XAABB;
    // first node

    const double* pos = ih->cutterdis()->lRowElement(0)->Nodes()[0]->X();
    for(int dim=0; dim<nsd; ++dim)
    {
        XAABB(dim, 0) = pos[dim] - TOL7;
        XAABB(dim, 1) = pos[dim] + TOL7;
    }

	for (int j=0; j< ih->cutterdis()->NumMyRowElements(); ++j) {
	    // remaining node
	    for(int i=1; i< ih->cutterdis()->lRowElement(j)->NumNode(); ++i)
	    {
	        const double* posEle = ih->cutterdis()->lRowElement(j)->Nodes()[i]->X();
	        for(int dim=0; dim<nsd; dim++)
	        {
	            XAABB(dim, 0) = std::min( XAABB(dim, 0), posEle[dim] - TOL7);
	            XAABB(dim, 1) = std::max( XAABB(dim, 1), posEle[dim] + TOL7);
	        }
	    }
	    
	    double maxDistance = fabs(XAABB(0,1) - XAABB(0,0));
	    for(int dim=1; dim<nsd; ++dim)
	       maxDistance = std::max(maxDistance, fabs(XAABB(dim,1)-XAABB(dim,0)) );
	    
	    // subtracts half of the maximal distance to minX, minY, minZ
	    // adds half of the maximal distance to maxX, maxY, maxZ 
	    const double halfMaxDistance = 0.5*maxDistance;
	    for(int dim=0; dim<nsd; ++dim)
	    {
	        XAABB(dim, 0) -= halfMaxDistance;
	        XAABB(dim, 1) += halfMaxDistance;
	    }   
	    
	    /*
	    printf("\n");
	    printf("axis-aligned bounding box:\n minX = %f\n minY = %f\n minZ = %f\n maxX = %f\n maxY = %f\n maxZ = %f\n", 
	              XAABB(0,0), XAABB(1,0), XAABB(2,0), XAABB(0,1), XAABB(1,1), XAABB(2,1));
	    printf("\n");
	    */

	
	}
	    
	    return XAABB;
}

///*----------------------------------------------------------------------*
// |  RQI:    searches the nearest point on a surface          u.may 02/08|
// |          element for a given point in physical coordinates           |
// *----------------------------------------------------------------------*/
//bool XFEM::searchForNearestPointOnSurfaceBound(
//    const DRT::Element*                     surfaceElement,
//    const BlitzMat&                         xyze_surfaceElement,
//    const BlitzVec3&                        physCoord,
//    BlitzVec2&                              xsi,
//    BlitzVec3&                              vector,
//    double&									distance)
//{
//	distance = -1.0;
//	vector = 0;
//	
//	CurrentToSurfaceElementCoordinates(surfaceElement, xyze_surfaceElement, physCoord, xsi);
//	
//	if (xsi(0)> 1) xsi(0)= 1;
//	if (xsi(0)<-1) xsi(0)=-1;
//	if (xsi(1)> 1) xsi(1)= 1;
//	if (xsi(1)<-1) xsi(1)=-1;
//	
//	const bool pointWithinElement = checkPositionWithinElementParameterSpace(xsi, surfaceElement->Shape());
//	
//	// normal vector at position xsi
//	static BlitzVec3 eleNormalAtXsi;
//	computeNormalToBoundaryElement(surfaceElement, xyze_surfaceElement, xsi, eleNormalAtXsi);
//	
//	BlitzVec3 x_surface_phys;
//	elementToCurrentCoordinates(surfaceElement, xyze_surfaceElement, xsi, x_surface_phys);
//	// vector pointing away from the surface towards physCoord
//	vector(0) = physCoord(0) - x_surface_phys(0);
//	vector(1) = physCoord(1) - x_surface_phys(1);
//	vector(2) = physCoord(2) - x_surface_phys(2);
//	// absolute distance between point and surface
//	distance = sqrt(vector(0)*vector(0) + vector(1)*vector(1) + vector(2)*vector(2));
//	// compute distance with sign
//	const double scalarproduct = eleNormalAtXsi(0)*vector(0) + eleNormalAtXsi(1)*vector(1) + eleNormalAtXsi(2)*vector(2);
//	const double teiler = Norm2(eleNormalAtXsi) * Norm2(vector);
//	const double cosphi = scalarproduct / teiler;
//	const double vorzeichen = cosphi/abs(cosphi);
//	distance *= vorzeichen;
//  
//	return pointWithinElement;
//}

#endif  // #ifdef CCADISCRET
