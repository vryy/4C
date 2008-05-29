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


XFEM::XSearchTree::XSearchTree() {
  TreeInit_  = false;
  ih_        = NULL;
  treeRoot_  = NULL;
}

XFEM::XSearchTree::~XSearchTree() {
  delete treeRoot_;
}

void XFEM::XSearchTree::setRebuildFlag() {
  TreeInit_ = false;
}

int XFEM::XSearchTree::getTotalRequests(){
  return searchRequests_;
}
int XFEM::XSearchTree::getMeanSearchlength(){
  return MeanSearchLength_;
}

int XFEM::XSearchTree::queryPointType(const XFEM::InterfaceHandle& ih, const BlitzVec3& pointCoords) {
  searchRequests_++;
  ih_ = &ih;
  if (!TreeInit_)
    rebuild();
  // search candidates in tree
  list<RCP<DRT::Element> > candidates;
  int labID;
  int tmpL =0;
  int tmpL1 =0;
  candidates = treeRoot_->queryPointType(pointCoords, labID, tmpL);
  if (candidates.empty()){
    //	return tmpL;
    return labID;
  }
  // do excat routines
  int within = 0;
  list<RCP<DRT::Element> >::iterator myIt = candidates.begin();
  double dist = 0;
  const DRT::Element* closestEle =this->nearestNeighbourInList(candidates, pointCoords,dist,tmpL1);
  if (dist<0){
    //	return tmpL;
    within = *(labelsPerElementId_.find(closestEle->Id())->second.begin());
  }
  else
    within = 0;
  return within;
}

int XFEM::XSearchTree::getDepth(){
  return treeRoot_->getDepth();
}

void XFEM::XSearchTree::insertElement(const RCP<DRT::Element> elem) {
  treeRoot_->insertElement(elem);
}

void XFEM::XSearchTree::rebuild() {
  // loop all elements on this processor
  if (treeRoot_ != NULL)
    delete treeRoot_;
  
  //  this->ih = ih;
  labelsPerElementId_.clear();
  for(std::map<int,set<int> >::const_iterator conditer = (*ih_->elementsByLabel()).begin(); conditer!=(*ih_->elementsByLabel()).end(); ++conditer)
  {
    for(std::set<int>::const_iterator eleid = conditer->second.begin(); eleid!=conditer->second.end(); ++eleid)
    {
      labelsPerElementId_[*eleid].insert(conditer->first);
    }
  }
  
  
  //TODO: remove 2nd loop over data: getXAABBofDis also loops over elements of dis
  BlitzMat3x2 aabb =getXAABBofDis();
  
  // we want to cover the whole fluid domain, but we dont know its dimensions, 
  // therefor its estimated by expanding the bounding box of all solids by a
  // constant factor;
  //  aabb(0,0) = aabb(0,0)-AABB_Factor*abs(aabb(0,0));
  //  aabb(0,1) = aabb(0,1)+AABB_Factor*abs(aabb(0,1));
  //  aabb(1,0) = aabb(1,0)-AABB_Factor*abs(aabb(1,0));
  //  aabb(1,1) = aabb(1,1)+AABB_Factor*abs(aabb(1,1));
  //  aabb(2,0) = aabb(2,0)-AABB_Factor*abs(aabb(2,0));
  //  aabb(2,1) = aabb(2,1)+AABB_Factor*abs(aabb(2,1));
  
  //  cout << "aabb=" << aabb(0,0) << ","<< aabb(0,1) << "," << aabb(1,0) << "," << aabb(1,1) << "," << aabb(2,0) << "," << aabb(2,1) << endl;
  treeRoot_ = new TreeNode(0,aabb,this);
  for (int i=0; i< ih_->cutterdis()->NumMyRowElements(); ++i) {
    //    DEBUGp("insert dis-elem\n")
    RCP<DRT::Element> ele =  Teuchos::rcp(ih_->cutterdis()->lRowElement(i));
    insertElement(ele);
  }
  //    DEBUGp("DEBUG_XFSI_searchtree: rebuilt treeRoot\n");
  TreeInit_ = true;
  
}

XFEM::XSearchTree::TreeNode::TreeNode(int Depth,BlitzMat3x2 aabb, XSearchTree* tree)
{
  actTreedepth_ = Depth;
  AABB_=aabb;
  XPlaneCoordinate_ = aabb(0,0)+(aabb(0,1)-aabb(0,0))/2.0;
  YPlaneCoordinate_ = aabb(1,0)+(aabb(1,1)-aabb(1,0))/2.0;
  ZPlaneCoordinate_ = aabb(2,0)+(aabb(2,1)-aabb(2,0))/2.0;
  State_ = STATE_LEAF_NODE;
  for (int i=0; i< 8; i++){
    children_[i]=NULL;
  }
  tree_ = tree;
}

XFEM::XSearchTree::TreeNode::~TreeNode() {
  for (int i=0; i< 8; i++){
    if (children_[i]!=NULL) 
      delete children_[i];
  }
}

int XFEM::XSearchTree::TreeNode::getDepth(){
  if (State_==STATE_LEAF_NODE){
    return actTreedepth_;
  }
  else {
    int depth=actTreedepth_;
    int tmp=0;
    for (int i=0; i< 8; i++){
      if (children_[i] != NULL)
        tmp = children_[i]->getDepth();
      else
        tmp =0;
      if (tmp>depth)
        depth=tmp;
    }
    return depth;
  }
  return actTreedepth_;
}

void XFEM::XSearchTree::TreeNode::setFluid(const int label){
  labelID_ = label;
}
void XFEM::XSearchTree::TreeNode::setSolid(const int label){
  labelID_ = label;
}

BlitzVec3& XFEM::XSearchTree::TreeNode::getCenterCoord(){
  BlitzVec3* t = new BlitzVec3(this->XPlaneCoordinate_, this->YPlaneCoordinate_, this->ZPlaneCoordinate_);
  return *t;
}

list<RCP<DRT::Element> > XFEM::XSearchTree::TreeNode::queryPointType(const BlitzVec3& pointCoords, int& lID, int& tmpL) {
  //  printf("AABB(%f\t%f\t%f\t%f\t%f\t%f)\t", AABB(0,0),AABB(0,1),AABB(1,0),AABB(1,1),AABB(2,0),AABB(2,1));
  //  printf("x_in(%f\t%f\t%f)\n",pointCoords(0), pointCoords(1),pointCoords(2));
  switch (State_) {
    case STATE_INNER_NODE:         // if inner node, do recursion 
      //    printf("*");
      //    printf("classified searchpoint to oct %d, so i will search there\n", classifyPoint(pointCoords)-1);
      return children_[classifyPoint(pointCoords)-1]->queryPointType(pointCoords ,lID, tmpL);
      break;
    case STATE_LEAF_NODE:   // returns list of candidates
      //    printf("/");
      if (actTreedepth_ >= MAX_TREEDEPTH){
        //      printf("actTreedepth >= MAX_TREEDEPTH\n");
        lID = labelID_;
        tmpL = actTreedepth_;
        return ElementList_;  // just return (maybe empty) list if leaf is at max depth, 
        // or if no candidates (==fluid)
      }
      if (ElementList_.empty()){
        //      printf("empty ");
        lID = labelID_;
        tmpL = -10-actTreedepth_;
        return ElementList_;  // just return (maybe empty) list if leaf is at max depth, 
        // or if no candidates (==fluid)
      }
      if (ElementList_.size()>1) // dynamically grow tree
      {
        //      printf("more than 1 element in leaf -> growing tree\n");
        // create Octants
        for (int i=0; i<8; i++){
          //        DEBUGp("creating child\n")
          BlitzMat3x2 chldAABB = getChildOctAABB(i+1);
          children_[i] = new TreeNode(actTreedepth_ + 1, chldAABB, tree_);
        }
        //      printf("created child nodes\n");
        // actual node becomes an inner tree node,
        // so we have to introduce one more tree-level
        for (list<RCP<DRT::Element> >::const_iterator myIt = ElementList_.begin(); myIt != ElementList_.end(); myIt++){
          classifyAndInsert(*myIt);
        }
        
        // if one of the created childs is empty, check if it is fluid or solid
        for (int i=0; i< 8; i++){
          if ((children_[i]->getElementList()).empty()){
            double dista=0;
            int tmpL1 = 0;
            BlitzVec3 c= children_[i]->getCenterCoord();
            //          printf("centerX (%f,%f,%f)\n", c(0),c(1),c(2)); 
            const DRT::Element* closestEle = tree_->nearestNeighbourInList(ElementList_, children_[i]->getCenterCoord(),dista, tmpL1);
            //          printf("distance of empty octant center to next element is %f\n", dista);
            if (dista<0) {
              int eleId = closestEle->Id();
              std::set<int>::const_iterator it =tree_->labelsPerElementId_.find(eleId)->second.begin();
              //            printf("solid octant with label %d\n", *it);
              children_[i]->setSolid(*it);
            }
            else {
              //            printf("fluid octant\n");
              children_[i]->setFluid(0);
            }
            // search for nearest Neigbour in ElementList and
            // check if child is fluid or solid
            // Element NN = NearestNeighbour(ElementList, child[i].center);
            // child[i] = fluid/solid
          }
          else {
            //          printf("non empty octant(#%d) with %d elements\n", i+1, children_[i]->getElementList().size());
          }
        }
        // this node becomes an inner tree node
        State_ = STATE_INNER_NODE;
        ElementList_.clear();
        // do recursion
        int l=0;
        //      printf("classified searchpoint to oct %d, so i will search there", classifyPoint(pointCoords)-1);
        list<RCP<DRT::Element> > myL = children_[classifyPoint(pointCoords)-1]->queryPointType(pointCoords, l, tmpL);
        lID = l;
        return myL;
      }
      else  // if there is only one Element, just return it
        tmpL = actTreedepth_;
      //      printf("returning candidate");
      return ElementList_;
      break;
  }
  dserror("should not get here\n");
  return ElementList_;
}



const DRT::Element* XFEM::XSearchTree::nearestNeighbourInList(const list<RCP<DRT::Element> > ElementList, const BlitzVec3& x_in, double& dist, int& tmpL)
{
  bool in_element = false;
  double min_ele_distance = 1.0e12;
  BlitzVec3 vectorX2minNode;
  const DRT::Element* closest_element;
  const DRT::Node* closest_node;
  bool foundNearSurfaceElement =false;
  
  for (list<RCP<DRT::Element> >::const_iterator myIt = ElementList.begin(); myIt != ElementList.end(); myIt++)
  {
    double distance = 1.0e12;
    const DRT::Element* cutterele = &**myIt;
    const BlitzMat xyze_cutter(getCurrentNodalPositions(cutterele, *ih_->currentcutterpositions()));
    static BlitzVec2 eleCoord;
    static BlitzVec3 normal;
    in_element = XFEM::searchForNearestPointOnSurface(cutterele,xyze_cutter,x_in,eleCoord,normal,distance);
    if (in_element && (abs(distance) < abs(min_ele_distance)))
    {
      closest_element = cutterele;
      min_ele_distance = distance;
      foundNearSurfaceElement =true;
    }
  }
  if (foundNearSurfaceElement)  
  {
    dist = min_ele_distance;
  }
  else 
  { 
    tmpL = -1;
    
    min_ele_distance = 1.0e12;
    //      printf("BRANCHING to node-distance ************************************************** \n");
    for (list<RCP<DRT::Element> >::const_iterator myIt2 = ElementList.begin(); myIt2 != ElementList.end(); myIt2++) {
      double distance = 1.0e12;
      const DRT::Element* cutterele = &**myIt2;       
      const int numnode = cutterele->NumNode();
      const DRT::Node*const* nodes = cutterele->Nodes();          
      for (int inode = 0; inode < numnode; ++inode)
      {
        BlitzVec3 vector;
        const DRT::Node* cutternode = nodes[inode]; 
        // node position in physical coordinates
        const BlitzVec3 x_node = ih_->currentcutterpositions()->find(cutternode->Id())->second;
        
        // vector pointing away from the node towards physCoord
        vector(0) = x_in(0) - x_node(0);
        vector(1) = x_in(1) - x_node(1);
        vector(2) = x_in(2) - x_node(2);
        // absolute distance between point and node
        distance = sqrt(vector(0)*vector(0) + vector(1)*vector(1) + vector(2)*vector(2));
        
        if (distance < min_ele_distance) {
          closest_node = cutternode;
          min_ele_distance = distance;
          vectorX2minNode=vector;
        }
      }          
    }
    //      	if (tmp!=closest_node->Id())
    //          printf("found nearest node with id %d (%f,%f,%f) for x_in (%f,%f,%f)", closest_node->Id(), closest_node->X()[0], closest_node->X()[1], closest_node->X()[2], x_in(0),x_in(1),x_in(2));
    // getElements corseponding to node
    // calculate normal at node by (vectorial) addition of element normals
    BlitzVec3 normal;
    normal(0)=0;normal(1)=0;normal(2)=0;
    DRT::Node*  node = ih_->cutterdis()->gNode(closest_node->Id());
    for(int j=0; j<node->NumElement();j++)
    {
      DRT::Element* surfaceElement = node->Elements()[j];
      BlitzMat xyze_surfaceElement(getCurrentNodalPositions(surfaceElement, *ih_->currentcutterpositions()));
      BlitzVec3 eleNormalAtXsi;
      BlitzVec2 xsi;
      CurrentToSurfaceElementCoordinates(surfaceElement, xyze_surfaceElement, node->X(), xsi);
      //        	printf("xsi (%f,%f)\n", xsi(0), xsi(1));
      // normal vector at position xsi
      computeNormalToBoundaryElement(surfaceElement, xyze_surfaceElement, xsi, eleNormalAtXsi);
      normal(0) = normal(0) +  eleNormalAtXsi(0);
      normal(1) = normal(1) +  eleNormalAtXsi(1);
      normal(2) = normal(2) +  eleNormalAtXsi(2);
    }
    //    	double norm = sqrt(normal(0)*normal(0)+normal(1)*normal(1)+normal(2)*normal(2));
    //      	if (tmp!=closest_node->Id())
    //       	  printf(", normal(%f,%f,%f)", normal(0),normal(1),normal(2));
    
    //      	normal(0) = normal(0)/norm;
    //        normal(1) = normal(1)/norm;
    //        normal(2) = normal(2)/norm;
    //      	printf("normed normal at node (%f,%f,%f)\n", normal(0),normal(1),normal(2));
    closest_element=node->Elements()[0];
    // compute distance with sign
    //      	if (tmp!=closest_node->Id())
    //           	  printf(", x2N: %f,%f,%f ", vectorX2minNode(0),vectorX2minNode(1),vectorX2minNode(2));
    const double scalarproduct = vectorX2minNode(0)*normal(0) + vectorX2minNode(1)*normal(1) + vectorX2minNode(2)*normal(2);
    const double vorzeichen = scalarproduct/abs(scalarproduct);
    //      	if (tmp!=closest_node->Id())
    //           	  printf(", vz: %f\n", vorzeichen);
    tmp=closest_node->Id();
    min_ele_distance *= vorzeichen;
    dist = min_ele_distance;
    //      	printf("end BRANCHING to node-distance ********************************************** \n");
  }
  return closest_element;
}

void XFEM::XSearchTree::TreeNode::insertElement(const RCP<DRT::Element> elem) {
  if ((actTreedepth_ >= XFEM::XSearchTree::MAX_TREEDEPTH) || (State_ == STATE_LEAF_NODE) ) {
    ElementList_.push_back(elem);
    //    printf("inserting element\n");
    State_ = STATE_LEAF_NODE;
  } else if(State_ == STATE_INNER_NODE) {
    classifyAndInsert(elem);
  } 
}

list<RCP<DRT::Element> > XFEM::XSearchTree::TreeNode::getElementList(){
  return ElementList_;
}

BlitzMat3x2 XFEM::XSearchTree::TreeNode::getChildOctAABB(int i){
  BlitzMat3x2 chldAABB;
  if (i>4){
    chldAABB(0,0) = XPlaneCoordinate_;
    chldAABB(0,1) = AABB_(0,1);
  }
  else {
    chldAABB(0,0) = this->AABB_(0,0);
    chldAABB(0,1) = XPlaneCoordinate_;
  }
  if ((i==3) || (i==4) || (i==7) || (i==8)){
    chldAABB(1,0) = YPlaneCoordinate_;
    chldAABB(1,1) = AABB_(1,1);
  }
  else {
    chldAABB(1,0) = AABB_(1,0);
    chldAABB(1,1) = YPlaneCoordinate_;
  }
  if ((i%2)==0){
    chldAABB(2,0) = ZPlaneCoordinate_;
    chldAABB(2,1) = AABB_(2,1);
  }
  else {
    chldAABB(2,0) = AABB_(2,0);
    chldAABB(2,1) = ZPlaneCoordinate_;
  }    
  //  printf("created chldAABB(%f\t%f\t%f\t%f\t%f\t%f)\n", chldAABB(0,0),chldAABB(0,1),chldAABB(1,0),chldAABB(1,1),chldAABB(2,0),chldAABB(2,1));
  return chldAABB;
  
}


int XFEM::XSearchTree::TreeNode::classifyPoint(const BlitzVec3& pointcoords) {
  int octIdx = 1;
  if (pointcoords(0) > XPlaneCoordinate_)
    octIdx = octIdx + 4;
  if (pointcoords(1) > YPlaneCoordinate_)
    octIdx = octIdx + 2;
  if (pointcoords(2) > ZPlaneCoordinate_)
    octIdx = octIdx + 1;
  return octIdx;
}

list<int> XFEM::XSearchTree::TreeNode::classifyElement(const RCP<DRT::Element> elem) {
  list<int> octants;
  const BlitzMat xyze(DRT::UTILS::PositionArrayBlitz(&*elem));
  const BlitzMat3x2 elemXAABB= XFEM::computeFastXAABB(&*elem, xyze);
  //bitset<8> PlaneBitSet; // every Bit stands for a plane
  //PlaneBitSet.reset();
  
  
  //TODO: eleganter lösen
  /*  BlitzVec ful = BlitzVec(elemXAABB(0, 0), elemXAABB(1, 0), elemXAABB(2, 1));
  BlitzVec fur = BlitzVec(elemXAABB(0, 0), elemXAABB(1, 1), elemXAABB(1, 1));
  BlitzVec fll = BlitzVec(elemXAABB(0, 0), elemXAABB(1, 0), elemXAABB(1, 0));
  BlitzVec flr = BlitzVec(elemXAABB(0, 0), elemXAABB(1, 1), elemXAABB(1, 0));
  BlitzVec rul = BlitzVec(elemXAABB(0, 1), elemXAABB(1, 0), elemXAABB(2, 1));
  BlitzVec rur = BlitzVec(elemXAABB(0, 1), elemXAABB(1, 1), elemXAABB(1, 1));
  BlitzVec rll = BlitzVec(elemXAABB(0, 1), elemXAABB(1, 0), elemXAABB(1, 0));
  BlitzVec rlr = BlitzVec(elemXAABB(0, 1), elemXAABB(1, 1), elemXAABB(1, 0));*/
  
  if (elemXAABB(0, 1) > XPlaneCoordinate_) {
    if (elemXAABB(1, 1) > YPlaneCoordinate_) {
      if (elemXAABB(2, 1) > ZPlaneCoordinate_)
        octants.push_back(8);
      if (elemXAABB(2, 0) <= ZPlaneCoordinate_)
        octants.push_back(7);
    }
    if (elemXAABB(1, 0) <= YPlaneCoordinate_) {
      if (elemXAABB(2, 1) > ZPlaneCoordinate_)
        octants.push_back(6);
      if (elemXAABB(2, 0) <= ZPlaneCoordinate_)
        octants.push_back(5);
    }
  }
  if (elemXAABB(0, 0) <= XPlaneCoordinate_) {
    if (elemXAABB(1, 1) > YPlaneCoordinate_) {
      if (elemXAABB(2, 1) > ZPlaneCoordinate_)
        octants.push_back(4);
      if (elemXAABB(2, 0) <= ZPlaneCoordinate_)
        octants.push_back(3);
    }
    if (elemXAABB(1, 0) <= YPlaneCoordinate_) {
      if (elemXAABB(2, 1) > ZPlaneCoordinate_)
        octants.push_back(2);
      if (elemXAABB(2, 0) <= ZPlaneCoordinate_)
        octants.push_back(1);
    }
    
  }
  //  printf("classifiing elem(%f,%f,%f,%f,%f,%f)", elemXAABB(0,0),elemXAABB(0,1),elemXAABB(1,0),elemXAABB(1,1),elemXAABB(2,0),elemXAABB(2,1));
  //  printf("classified element to octants: ");
  for (list<int>::const_iterator myIt = octants.begin(); myIt != octants.end(); myIt++)
  {
    //  printf("%d ", *myIt);
  }
  //  printf("\n");
  return octants;
  
}

int XFEM::XSearchTree::TreeNode::getState() {
  return State_;
}

void XFEM::XSearchTree::printTree(const int step) const{
  std::stringstream filename;
  std::stringstream fc;
  filename << "tree" << std::setw(5) << setfill('0') << step << ".pos";
  cout << "writing "<<filename.str()<<" ...";
  fc << "View \" " << "fsiOctree \" {" << endl;
  treeRoot_->printTree(fc);
  fc << "};" << endl;
  std::ofstream f_system(filename.str().c_str());
  f_system << fc.str();
  f_system.close();
  cout << " done" << endl;
}

void XFEM::XSearchTree::TreeNode::printTree(stringstream& fc) const{
  if (State_==STATE_INNER_NODE){
    for (int j=0; j<8; j++){
      if (children_[j]!=NULL)
        children_[j]->printTree(fc);
    }	  
  }
  else if (State_==STATE_LEAF_NODE)
  {
    BlitzMat XAABB(3,8);
    XAABB(0,0) = AABB_(0,0); XAABB(1,0) = AABB_(1,0);XAABB(2,0) = AABB_(2,0);
    XAABB(0,1) = AABB_(0,0); XAABB(1,1) = AABB_(1,1);XAABB(2,1) = AABB_(2,0);
    XAABB(0,2) = AABB_(0,0); XAABB(1,2) = AABB_(1,1);XAABB(2,2) = AABB_(2,1);
    XAABB(0,3) = AABB_(0,0); XAABB(1,3) = AABB_(1,0);XAABB(2,3) = AABB_(2,1);
    XAABB(0,4) = AABB_(0,1); XAABB(1,4) = AABB_(1,0);XAABB(2,4) = AABB_(2,0);
    XAABB(0,5) = AABB_(0,1); XAABB(1,5) = AABB_(1,1);XAABB(2,5) = AABB_(2,0);
    XAABB(0,6) = AABB_(0,1); XAABB(1,6) = AABB_(1,1);XAABB(2,6) = AABB_(2,1);
    XAABB(0,7) = AABB_(0,1); XAABB(1,7) = AABB_(1,0);XAABB(2,7) = AABB_(2,1);
    fc << IO::GMSH::cellWithScalarToString(DRT::Element::hex8, -actTreedepth_-10, XAABB)<< endl;
  }
  
}

const BlitzMat3x2& XFEM::XSearchTree::TreeNode::getAABB(){
  return AABB_;
}

XFEM::XSearchTree::TreeNode* XFEM::XSearchTree::TreeNode::getChild(int idx) {
  return children_[idx-1];
}

BlitzMat3x2 XFEM::XSearchTree::getXAABBofDis(){
  const int nsd = 3;
  BlitzMat3x2 XAABB;
  
  // XAABB of cutterElements
  // first node
  const double* pos = ih_->cutterdis()->lRowElement(0)->Nodes()[0]->X();
  for(int dim=0; dim<nsd; ++dim)
  {
    XAABB(dim, 0) = pos[dim] - TOL7;
    XAABB(dim, 1) = pos[dim] + TOL7;
  }
  for (int j=0; j< ih_->cutterdis()->NumMyRowElements(); ++j) {
    // remaining node
    for(int i=0; i< ih_->cutterdis()->lRowElement(j)->NumNode(); ++i)
    {
      const double* posEle = ih_->cutterdis()->lRowElement(j)->Nodes()[i]->X();
      for(int dim=0; dim<nsd; dim++)
      {
        XAABB(dim, 0) = std::min( XAABB(dim, 0), posEle[dim] - TOL7);
        XAABB(dim, 1) = std::max( XAABB(dim, 1), posEle[dim] + TOL7);
      }
      //      cout << "XAABB=" << XAABB(0,0) << ","<< XAABB(0,1) << "," << XAABB(1,0) << "," << XAABB(1,1) << "," << XAABB(2,0) << "," << XAABB(2,1) << endl;
    }  
  }
  
  // extend XAABB to xfem elements
  for (int j=0; j< ih_->xfemdis()->NumMyRowElements(); ++j) {
    // remaining node
    for(int i=0; i< ih_->xfemdis()->lRowElement(j)->NumNode(); ++i)
    {
      const double* posEle = ih_->xfemdis()->lRowElement(j)->Nodes()[i]->X();
      for(int dim=0; dim<nsd; dim++)
      {
        XAABB(dim, 0) = std::min( XAABB(dim, 0), posEle[dim] - TOL7);
        XAABB(dim, 1) = std::max( XAABB(dim, 1), posEle[dim] + TOL7);
      }
    }  
  }
  XAABB(0,0) = XAABB(0,0) - TOL7;
  XAABB(0,1) = XAABB(0,1) + TOL7;
  XAABB(1,0) = XAABB(1,0) - TOL7;
  XAABB(1,1) = XAABB(1,1) + TOL7;
  XAABB(2,0) = XAABB(2,0) - TOL7;
  XAABB(2,1) = XAABB(2,1) + TOL7;
  
  //  cout << "_XAABB=" << XAABB(0,0) << ","<< XAABB(0,1) << "," << XAABB(1,0) << "," << XAABB(1,1) << "," << XAABB(2,0) << "," << XAABB(2,1) << endl;
  return XAABB;
}

#endif  // #ifdef CCADISCRET
