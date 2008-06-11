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

using namespace std;
using namespace XFEM;


XFEM::XSearchTree::XSearchTree() {
  TreeInit_  = false;
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

int XFEM::XSearchTree::queryPointType(const RCP<DRT::Discretization> dis,const std::map<int,BlitzVec3>& currentpositions, const BlitzVec3& pointCoords) {
  searchRequests_++;
  if (!TreeInit_)
    rebuild(dis, currentpositions);
  if (dis->NumGlobalElements() == 0){ 
    return 0;
    }
  // search for candidates in tree
  list< const DRT::Element* > candidates;
  int labID;
  candidates = treeRoot_->queryPointType(dis, currentpositions, pointCoords, labID);
  if (candidates.empty()){
    return labID;
  }
  // do excat routines
  int within = 0;
  double dist = 0;
  list< const DRT::Element* >::iterator myIt = candidates.begin();
  const DRT::Element* closestEle = XFEM::nearestNeighbourInList(dis, currentpositions, candidates, pointCoords,dist);
  if (dist<0){
    within = labelByElement_[closestEle->Id()];
  }
  else
    within = 0;
  return within;
}

int XFEM::XSearchTree::getDepth(){
  return treeRoot_->getDepth();
}

void XFEM::XSearchTree::insertElement(const DRT::Element* elem) {
  treeRoot_->insertElement(elem);
}

void XFEM::XSearchTree::rebuild(const RCP<DRT::Discretization> dis,const std::map<int,BlitzVec3>& currentpositions) {
  if (treeRoot_ != NULL){
    delete treeRoot_;
  }
  BlitzMat3x2 aabb =getXAABBofDis(dis, currentpositions);
  treeRoot_ = new TreeNode(0,aabb, this);
  cout << "inserting new elements (" << dis->NumMyRowElements() << ")"<< endl;
  for (int i=0; i<dis->NumMyRowElements(); ++i) {
    const DRT::Element* ele =  dis->lRowElement(i);
    insertElement(ele);
  }
  TreeInit_ = true;
  
  std::map<int,set<int> >   elementsByLabel;
  cout<<"collectElementsByXFEMCouplingLabels" << endl;
  XFEM::CollectElementsByXFEMCouplingLabel(*dis, elementsByLabel);
  labelByElement_.clear();
  for(std::map<int,set<int> >::const_iterator conditer = elementsByLabel.begin(); conditer!=elementsByLabel.end(); ++conditer)
  {
    for(std::set<int>::const_iterator eleid = conditer->second.begin(); eleid!=conditer->second.end(); ++eleid)
    {
      labelByElement_[*eleid] = conditer->first;
    }
  }
  
}

XFEM::XSearchTree::TreeNode::TreeNode(int Depth,BlitzMat3x2 aabb, XSearchTree* tree)
{
  tree_ = tree;
  actTreedepth_ = Depth;
  AABB_=aabb;
  XPlaneCoordinate_ = aabb(0,0)+(aabb(0,1)-aabb(0,0))/2.0;
  YPlaneCoordinate_ = aabb(1,0)+(aabb(1,1)-aabb(1,0))/2.0;
  ZPlaneCoordinate_ = aabb(2,0)+(aabb(2,1)-aabb(2,0))/2.0;
  State_ = STATE_LEAF_NODE;
  for (int i=0; i< 8; i++){
    children_[i]=NULL;
  }
}

XFEM::XSearchTree::TreeNode::~TreeNode() {
//  cout << "destructor of treenode" << endl;
  for (int i=0; i< 8; i++){
//    cout << "deleting child  "  << i << endl;
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

list< const DRT::Element* > XFEM::XSearchTree::TreeNode::queryPointType(const RCP<DRT::Discretization> dis,const std::map<int,BlitzVec3>& currentpositions, const BlitzVec3& pointCoords, int& lID) {
  //  printf("AABB(%f\t%f\t%f\t%f\t%f\t%f)\t", AABB(0,0),AABB(0,1),AABB(1,0),AABB(1,1),AABB(2,0),AABB(2,1));
  //  printf("x_in(%f\t%f\t%f)\n",pointCoords(0), pointCoords(1),pointCoords(2));
  switch (State_) {
    case STATE_INNER_NODE:         // if inner node, do recursion 
      //    printf("*");
      //    printf("classified searchpoint to oct %d, so i will search there\n", classifyPoint(pointCoords)-1);
      return children_[classifyPoint(pointCoords)-1]->queryPointType(dis, currentpositions, pointCoords ,lID);
      break;
    case STATE_LEAF_NODE:   // returns list of candidates
      //    printf("/");
      if (actTreedepth_ >= MAX_TREEDEPTH){
        //      printf("actTreedepth >= MAX_TREEDEPTH\n");
        lID = labelID_;
        return ElementList_;  // just return (maybe empty) list if leaf is at max depth, 
        // or if no candidates (==fluid)
      }
      if (ElementList_.empty()){
        //      printf("empty ");
        lID = labelID_;
        return ElementList_;  // just return (maybe empty) list if leaf is at max depth, 
        // or if no candidates (==fluid)
      }
      if (ElementList_.size()>1) // dynamically grow tree
      {
        // create Octants
        for (int i=0; i<8; i++){
          BlitzMat3x2 chldAABB = getChildOctAABB(i+1);
          children_[i] = new TreeNode(actTreedepth_ + 1, chldAABB, tree_);
        }
        // actual node becomes an inner tree node,
        // so we have to introduce one more tree-level
        for (list< const DRT::Element* >::const_iterator myIt = ElementList_.begin(); myIt != ElementList_.end(); myIt++){
          list<int> childIdx = classifyElement(*myIt);
          for (list<int>::const_iterator myIt2 = childIdx.begin(); myIt2 != childIdx.end(); myIt2++){
            this->children_[*myIt2-1]->insertElement(*myIt);
            BlitzMat3x2 ab = this->children_[*myIt2-1]->getAABB();
            //  printf("inserted elem to AABB(%f,%f,%f,%f,%f,%f)\n", ab(0,0),ab(0,1),ab(1,0),ab(1,1),ab(2,0),ab(2,1));
           }
        }
        
        // if one of the created childs is empty, check if it is fluid or solid
        for (int i=0; i< 8; i++){
          if ((children_[i]->getElementList()).empty()){
            double dista=0;
            BlitzVec3 c= children_[i]->getCenterCoord();
            //          printf("centerX (%f,%f,%f)\n", c(0),c(1),c(2)); 
            const DRT::Element* closestEle = XFEM::nearestNeighbourInList(dis, currentpositions, ElementList_, children_[i]->getCenterCoord(),dista);
            //          printf("distance of empty octant center to next element is %f\n", dista);
            if (dista<0) {
              int eleId = closestEle->Id();
              children_[i]->setSolid(tree_->getLabelByElementID(eleId)); 
            }
            else {
              //            printf("fluid octant\n");
              children_[i]->setFluid(0);
            }
          }
          else {
            //          printf("non empty octant(#%d) with %d elements\n", i+1, children_[i]->getElementList().size());
          }
        }
        // this node becomes an inner tree node
        State_ = STATE_INNER_NODE;
        ElementList_.clear();
        // do recursion
        // cout << "classified searchpoint to oct "<< classifyPoint(pointCoords)-1 << " , so i will search there" << endl;
        list< const DRT::Element* > myL = children_[classifyPoint(pointCoords)-1]->queryPointType(dis, currentpositions, pointCoords, lID);
        return myL;
      }
      else  // if there is only one Element, just return it
        return ElementList_;
      break;
  }
  dserror("should not get here\n");
  return ElementList_;
}

void XFEM::XSearchTree::TreeNode::insertElement(const DRT::Element* elem) {
  if ((actTreedepth_ >= XFEM::XSearchTree::MAX_TREEDEPTH) || (State_ == STATE_LEAF_NODE) ) {
    ElementList_.push_back(elem);
//    cout << "inserted element" <<endl;
    State_ = STATE_LEAF_NODE;
  } else if(State_ == STATE_INNER_NODE) {
    list<int> childIdx = classifyElement(elem);
    for (list<int>::const_iterator myIt = childIdx.begin(); myIt != childIdx.end(); myIt++){
      this->children_[*myIt-1]->insertElement(elem);
//      BlitzMat3x2 ab = this->children_[*myIt-1]->getAABB();
//      printf("inserted elem to AABB(%f,%f,%f,%f,%f,%f)\n", ab(0,0),ab(0,1),ab(1,0),ab(1,1),ab(2,0),ab(2,1));
     }

  } 
}

list< const DRT::Element* > XFEM::XSearchTree::TreeNode::getElementList(){
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

list<int> XFEM::XSearchTree::TreeNode::classifyElement(const DRT::Element* elem) {
  list<int> octants;
  const BlitzMat xyze(DRT::UTILS::PositionArrayBlitz(&*elem));
  const BlitzMat3x2 elemXAABB= XFEM::computeFastXAABB(&*elem, xyze);
  
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
//    printf("classifiing elem(%f,%f,%f,%f,%f,%f)", elemXAABB(0,0),elemXAABB(0,1),elemXAABB(1,0),elemXAABB(1,1),elemXAABB(2,0),elemXAABB(2,1));
//    printf("classified element to octants: ");
//  for (list<int>::const_iterator myIt = octants.begin(); myIt != octants.end(); myIt++)
//  {
//      printf("%d ", *myIt);
//  }
//    printf("\n");
  return octants;
  
}

int XFEM::XSearchTree::TreeNode::getState() {
  return State_;
}

int XFEM::XSearchTree::getLabelByElementID(int gid){
  return labelByElement_[gid];
}

const BlitzMat3x2& XFEM::XSearchTree::TreeNode::getAABB(){
  return AABB_;
}

XFEM::XSearchTree::TreeNode* XFEM::XSearchTree::TreeNode::getChild(int idx) {
  return children_[idx-1];
}

void XFEM::XSearchTree::printTree(const int step) const{
  if (treeRoot_->getElementList().empty()){ 
    cout << "nothing to write, tree empty" << endl;
    return;
    }
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

#endif  // #ifdef CCADISCRET
