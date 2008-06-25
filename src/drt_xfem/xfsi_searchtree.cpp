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
#include "intersection_service.H"
#include "../drt_io/io_gmsh.H"
#include <Teuchos_TimeMonitor.hpp>

using namespace std;

XFEM::XSearchTree::XSearchTree() :
  MAX_OVERLAP_ALLOWED_(0.875),
  ELEMENTS_USED_FOR_OVERLAP_CHECK_(1),
  TreeInit_(false),
  hasExternAABB_(false),
  searchRequests_(0),
  Candidates_(0),
  SearchLength_(0),
  MinSearchLength_(INT_MAX),
  treeRoot_(NULL)
{
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
double XFEM::XSearchTree::getMeanSearchlength(){
  return SearchLength_/searchRequests_;
}

int XFEM::XSearchTree::queryPointType(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions, const BlitzVec3& pointCoords) {
  TEUCHOS_FUNC_TIME_MONITOR("XSearchTree - queryTime"); 
  searchRequests_++;
  if (dis.NumGlobalElements() == 0){ 
    return 0;
    }
  if (!TreeInit_)
    rebuild(dis, currentpositions);
  // search for candidates in tree
  int labID;
  const list< const DRT::Element* > candidates = treeRoot_->queryPointType(dis, currentpositions, pointCoords, labID);
  if (candidates.empty()){
    return labID;
  }
  // do excat routines
  #ifdef DEBUG
  Candidates_ = Candidates_ + candidates.size(); 
  #endif
  int within = 0;
  double dist = 0;
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

void XFEM::XSearchTree::insertElement(
    const DRT::Element* elem,
    const map<int,BlitzVec3>& currentpositions
    )
{
  treeRoot_->insertElement(elem,currentpositions);
}

void XFEM::XSearchTree::rebuild(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions) {
  if (treeRoot_ != NULL){
    delete treeRoot_;
  }
  if (hasExternAABB_) {
    treeRoot_ = new TreeNode(0,AABB_, this);
  }   
  else {
    const BlitzMat3x2 aabb(getXAABBofDis(dis, currentpositions));
    treeRoot_ = new TreeNode(0,aabb, this); 
  }  
  #ifdef DEBUG
  cout << endl <<"inserting new elements (" << dis.NumMyColElements() << ")"<< endl;
  
  std::stringstream filename;
  std::stringstream fc;
  int step = 0;
  filename << "xaabbs" << std::setw(5) << setfill('0') << step << ".pos";
  cout << endl << "writing... "<<filename.str()<<" ...";
  fc << "View \" " << "XAABB of Elements \" {" << endl;
  flush(cout);
  #endif

  
  for (int i=0; i<dis.NumMyColElements(); ++i) {
    insertElement(dis.lRowElement(i),currentpositions);
    
    #ifdef DEBUG
    const BlitzMat xyze(XFEM::getCurrentNodalPositions(dis.lRowElement(i), currentpositions));
    const BlitzMat3x2 elemXAABB = XFEM::computeFastXAABB(dis.lRowElement(i), xyze);
    BlitzMat XAABB(3,8);
    XAABB(0,0) = elemXAABB(0,0); XAABB(1,0) = elemXAABB(1,0);XAABB(2,0) = elemXAABB(2,0);
    XAABB(0,1) = elemXAABB(0,0); XAABB(1,1) = elemXAABB(1,1);XAABB(2,1) = elemXAABB(2,0);
    XAABB(0,2) = elemXAABB(0,0); XAABB(1,2) = elemXAABB(1,1);XAABB(2,2) = elemXAABB(2,1);
    XAABB(0,3) = elemXAABB(0,0); XAABB(1,3) = elemXAABB(1,0);XAABB(2,3) = elemXAABB(2,1);
    XAABB(0,4) = elemXAABB(0,1); XAABB(1,4) = elemXAABB(1,0);XAABB(2,4) = elemXAABB(2,0);
    XAABB(0,5) = elemXAABB(0,1); XAABB(1,5) = elemXAABB(1,1);XAABB(2,5) = elemXAABB(2,0);
    XAABB(0,6) = elemXAABB(0,1); XAABB(1,6) = elemXAABB(1,1);XAABB(2,6) = elemXAABB(2,1);
    XAABB(0,7) = elemXAABB(0,1); XAABB(1,7) = elemXAABB(1,0);XAABB(2,7) = elemXAABB(2,1);
    fc << IO::GMSH::cellWithScalarToString(DRT::Element::hex8, i+10, XAABB)<< endl;
    #endif

  }
  
  #ifdef DEBUG
    fc << "};" << endl;
    std::ofstream f_system(filename.str().c_str());
    f_system << fc.str();
    f_system.close();
    cout << " done" << endl;
  #endif
   
  TreeInit_ = true;
  
  std::map<int,set<int> >   elementsByLabel;
  XFEM::CollectElementsByXFEMCouplingLabel(dis, elementsByLabel);
  labelByElement_.clear();
  for(std::map<int,set<int> >::const_iterator conditer = elementsByLabel.begin(); conditer!=elementsByLabel.end(); ++conditer)
  {
    for(std::set<int>::const_iterator eleid = conditer->second.begin(); eleid!=conditer->second.end(); ++eleid)
    {
      labelByElement_[*eleid] = conditer->first;
    }
  }
  
}

XFEM::XSearchTree::TreeNode::TreeNode(const int Depth, const BlitzMat3x2& aabb, XSearchTree* tree) :
  State_(STATE_LEAF_NODE),
  dont_refine_(false),
  actTreedepth_(Depth),
  tree_(tree),
  AABB_(aabb)
{
  XPlaneCoordinate_ = aabb(0,0)+(aabb(0,1)-aabb(0,0))/2.0;
  YPlaneCoordinate_ = aabb(1,0)+(aabb(1,1)-aabb(1,0))/2.0;
  ZPlaneCoordinate_ = aabb(2,0)+(aabb(2,1)-aabb(2,0))/2.0;
  for (int i=0; i< 8; i++){
    children_[i]=NULL;
  }
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

int XFEM::XSearchTree::TreeNode::getNodesInTree(int type){
  int nodes=0;
  int tmp=0;
  if (State_==STATE_LEAF_NODE){
    switch (type) {
      case STATE_INNER_NODE:
        return 0;
      case STATE_LEAF_NODE:
        return 1;
      case -1:
        return 1;
    }
  }
  else {
    for (int i=0; i< 8; i++){
      if (children_[i] != NULL)
        tmp = children_[i]->getNodesInTree(type);
      else
        tmp =0;
      nodes = nodes + tmp;
    }
    switch (type) {
      case STATE_INNER_NODE:  nodes++;break;
      case STATE_LEAF_NODE:           break;
      case -1:                nodes++;break;
    }
  }
  return nodes;
}

void XFEM::XSearchTree::TreeNode::setFluid(const int label){
  labelID_ = label;
}
void XFEM::XSearchTree::TreeNode::setSolid(const int label){
  labelID_ = label;
}

BlitzVec3 XFEM::XSearchTree::TreeNode::getCenterCoord(){
  return BlitzVec3(this->XPlaneCoordinate_, this->YPlaneCoordinate_, this->ZPlaneCoordinate_);
}

list< const DRT::Element* > XFEM::XSearchTree::TreeNode::queryPointType(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions, const BlitzVec3& pointCoords, int& lID) {
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
        #ifdef DEBUG 
        tree_->MinSearchLength_ = std::min(tree_->MinSearchLength_, actTreedepth_);
        tree_->SearchLength_ = tree_->SearchLength_ + actTreedepth_+1;
        #endif
        return ElementList_;  // just return (maybe empty) list if leaf is at max depth, 
      }
      
      if (ElementList_.empty()){
        //      printf("empty ");
        lID = labelID_;
        #ifdef DEBUG 
        tree_->MinSearchLength_ = std::min(tree_->MinSearchLength_, actTreedepth_);
        tree_->SearchLength_ = tree_->SearchLength_ + actTreedepth_+1;
        #endif
        return ElementList_;  // just return empty list if there are no candidates (==fluid) 
      }
      
      if (dont_refine_){       //
        #ifdef DEBUG 
        tree_->MinSearchLength_ = std::min(tree_->MinSearchLength_, actTreedepth_);
        tree_->SearchLength_ = tree_->SearchLength_ + actTreedepth_+1;
        #endif
        return ElementList_;  // just return empty list if node was previously set not to refine any further  
      }

      if (ElementList_.size()>1) { // dynamically grow tree
        map<const DRT::Element*, list<int> > ElementClassification;
        if (ElementList_.size()<=tree_->ELEMENTS_USED_FOR_OVERLAP_CHECK_){
          list<BlitzMat3x2> AABBs;
          AABBs.clear();
          AABBs.push_back(AABB_);
          for (list< const DRT::Element* >::const_iterator myIt = ElementList_.begin(); myIt != ElementList_.end(); myIt++){
            const BlitzMat xyze(XFEM::getCurrentNodalPositions(*myIt,currentpositions));
            const BlitzMat3x2 boundingbox = XFEM::computeFastXAABB(*myIt, xyze, HIGHERORDER);
            AABBs.push_back(boundingbox);
            list<int> childIdx = classifyElement(*myIt,currentpositions);
            ElementClassification[*myIt]=childIdx;
          }
          const double overlap = XFEM::getOverlapArea(AABBs);
          const double AABBcoverage = overlap/XFEM::getArea(AABB_);
           if (AABBcoverage > tree_->MAX_OVERLAP_ALLOWED_){
            dont_refine_ = true;
            #ifdef DEBUG 
            tree_->MinSearchLength_ = std::min(tree_->MinSearchLength_, actTreedepth_);
            tree_->SearchLength_ = tree_->SearchLength_ + actTreedepth_+1;
            #endif
            return ElementList_;
          }
        }
        else {
          for (list< const DRT::Element* >::const_iterator myIt = ElementList_.begin(); myIt != ElementList_.end(); myIt++){
            list<int> childIdx = classifyElement(*myIt,currentpositions);
            ElementClassification[*myIt]=childIdx;
          }
        }
        // create Octants
        for (int i=0; i<8; i++){
          const BlitzMat3x2 chldAABB(getChildOctAABB(i+1));
          children_[i] = new TreeNode(actTreedepth_ + 1, chldAABB, tree_);
        }
        // actual node becomes an inner tree node,
        // so we have to introduce one more tree-level
        for (list< const DRT::Element* >::const_iterator myIt = ElementList_.begin(); myIt != ElementList_.end(); myIt++){
          const list<int> childIdx = ElementClassification[*myIt];
          for (list<int>::const_iterator myIt2 = childIdx.begin(); myIt2 != childIdx.end(); myIt2++){
            this->children_[*myIt2-1]->insertElement(*myIt,currentpositions);
            //const BlitzMat3x2 ab(this->children_[*myIt2-1]->getAABB());
            //  printf("inserted elem to AABB(%f,%f,%f,%f,%f,%f)\n", ab(0,0),ab(0,1),ab(1,0),ab(1,1),ab(2,0),ab(2,1));
          }
        }

        // if one of the created childs is empty, check if it is fluid or solid
        for (int i=0; i< 8; i++){
          if ((children_[i]->getElementList()).empty()){
            double dista=0;
            const BlitzVec3 c(children_[i]->getCenterCoord());
            //          printf("centerX (%f,%f,%f)\n", c(0),c(1),c(2)); 
            const DRT::Element* closestEle = XFEM::nearestNeighbourInList(dis, currentpositions, ElementList_, children_[i]->getCenterCoord(),dista);
            //          printf("distance of empty octant center to next element is %f\n", dista);
            if (dista<0) {
              const int eleId = closestEle->Id();
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
        // ElementList_.clear();
        // do recursion
        // cout << "classified searchpoint to oct "<< classifyPoint(pointCoords)-1 << " , so i will search there" << endl;
        return children_[classifyPoint(pointCoords)-1]->queryPointType(dis, currentpositions, pointCoords, lID);
      }
      else  // if there is only one Element, just return it
        #ifdef DEBUG 
        tree_->MinSearchLength_ = std::min(tree_->MinSearchLength_, actTreedepth_);
        tree_->SearchLength_ = tree_->SearchLength_ + actTreedepth_+1;
        #endif
        return ElementList_;
      break;
  }
  dserror("should not get here\n");
  return ElementList_;
}

void XFEM::XSearchTree::TreeNode::insertElement(
    const DRT::Element* elem,
    const map<int,BlitzVec3>& currentpositions
    ) {
  if ((actTreedepth_ >= XFEM::XSearchTree::MAX_TREEDEPTH) || (State_ == STATE_LEAF_NODE) ) {
    ElementList_.push_back(elem);
//    cout << "inserted element at depth " << actTreedepth_ <<endl;
    State_ = STATE_LEAF_NODE;
  } else if(State_ == STATE_INNER_NODE) {
    const list<int> childIdx(classifyElement(elem,currentpositions));
    for (list<int>::const_iterator myIt = childIdx.begin(); myIt != childIdx.end(); myIt++){
      this->children_[*myIt-1]->insertElement(elem,currentpositions);
//      BlitzMat3x2 ab = this->children_[*myIt-1]->getAABB();
//      printf("inserted elem to AABB(%f,%f,%f,%f,%f,%f)\n", ab(0,0),ab(0,1),ab(1,0),ab(1,1),ab(2,0),ab(2,1));
     }

  } 
}

list< const DRT::Element* > XFEM::XSearchTree::TreeNode::getElementList(){
  return ElementList_;
}

void XFEM::XSearchTree::setMaxOverlap(double maxO){
  MAX_OVERLAP_ALLOWED_ = maxO;
}

double XFEM::XSearchTree::getMaxOverlap(){
  return MAX_OVERLAP_ALLOWED_;
}


BlitzMat3x2 XFEM::XSearchTree::TreeNode::getChildOctAABB(const int i){
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
    octIdx += 4;
  if (pointcoords(1) > YPlaneCoordinate_)
    octIdx += 2;
  if (pointcoords(2) > ZPlaneCoordinate_)
    octIdx += 1;
  return octIdx;
}

list<int> XFEM::XSearchTree::TreeNode::classifyElement(
    const DRT::Element*       elem,
    const map<int,BlitzVec3>& currentpositions
    ) {
//  bool isAABBbiggerThanElemXAABB;
  list<int> octants;
  const BlitzMat xyze(XFEM::getCurrentNodalPositions(elem,currentpositions));
  const BlitzMat3x2 elemXAABB(XFEM::computeFastXAABB(elem, xyze, HIGHERORDER));
  
  if (elemXAABB(0, 1) > XPlaneCoordinate_) {
    if (elemXAABB(1, 1) > YPlaneCoordinate_) {
      if (elemXAABB(2, 1) > ZPlaneCoordinate_){
        octants.push_back(8);
//        do_refine[8] = isAABBbiggerThanElemXAABB;
      }
      if (elemXAABB(2, 0) <= ZPlaneCoordinate_){
        octants.push_back(7);
//        do_refine[7] = isAABBbiggerThanElemXAABB;
      }
    }
    if (elemXAABB(1, 0) <= YPlaneCoordinate_) {
      if (elemXAABB(2, 1) > ZPlaneCoordinate_){
        octants.push_back(6);
//        do_refine[6] = isAABBbiggerThanElemXAABB;
      }
      if (elemXAABB(2, 0) <= ZPlaneCoordinate_){
        octants.push_back(5);
//        do_refine[5] = isAABBbiggerThanElemXAABB;
      }
    }
  }
  if (elemXAABB(0, 0) <= XPlaneCoordinate_) {
    if (elemXAABB(1, 1) > YPlaneCoordinate_) {
      if (elemXAABB(2, 1) > ZPlaneCoordinate_){
        octants.push_back(4);
//        do_refine[4] = isAABBbiggerThanElemXAABB;
      }
      if (elemXAABB(2, 0) <= ZPlaneCoordinate_){
        octants.push_back(3);
//        do_refine[3] = isAABBbiggerThanElemXAABB;
      }
    }
    if (elemXAABB(1, 0) <= YPlaneCoordinate_) {
      if (elemXAABB(2, 1) > ZPlaneCoordinate_){
        octants.push_back(2);
//        do_refine[2] = isAABBbiggerThanElemXAABB;
      }
      if (elemXAABB(2, 0) <= ZPlaneCoordinate_){
        octants.push_back(1);
//        do_refine[1] = isAABBbiggerThanElemXAABB;
      }
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

int XFEM::XSearchTree::getLabelByElementID(const int gid){
  return labelByElement_[gid];
}

const BlitzMat3x2& XFEM::XSearchTree::TreeNode::getAABB(){
  return AABB_;
}

XFEM::XSearchTree::TreeNode* XFEM::XSearchTree::TreeNode::getChild(const int idx) {
  return children_[idx-1];
}

void XFEM::XSearchTree::printTree(const int step) const{
  cout << endl << "writing... ";
  if (!TreeInit_) {
    cout << "nothing to write, tree not initialized yet -> done" << endl;
    return;
  }
  if (treeRoot_->getElementList().empty()){ 
    cout << "nothing to write, tree empty -> done" << endl;
    return;
  }
  std::stringstream filename;
  std::stringstream fc;
  filename << "tree" << std::setw(5) << setfill('0') << step << ".pos";
  cout << " "<<filename.str()<<" ...";
  fc << "View \" " << "fsiOctree \" {" << endl;
  treeRoot_->printTree(fc);
  fc << "};" << endl;
  std::ofstream f_system(filename.str().c_str());
  f_system << fc.str();
  f_system.close();
  cout << " done" << endl;
}

int XFEM::XSearchTree::getMemoryUsage() const{
  int mem = treeRoot_->getMemoryUsage();
  mem = mem + sizeof this;
  return mem;
}

int XFEM::XSearchTree::TreeNode::getMemoryUsage() const{
  int mem=0;
  if (State_==STATE_LEAF_NODE){
    return sizeof this;
  }
  else {
    mem=sizeof this;
    int tmp=0;
    for (int i=0; i< 8; i++){
      if (children_[i] != NULL)
        tmp = children_[i]->getMemoryUsage();
      else
        tmp =0;
    mem = mem + tmp;
    }
  }
  return mem;
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

void XFEM::XSearchTree::printTreeMetrics(const int step) const{
#ifdef DEBUG 
  cout << "\t********************* TREE STATS ******************" << endl;
  if (!TreeInit_ ) {
    cout << "\tno stats, tree not initialized yet -> done" << endl;
    return;
  }
  if (treeRoot_->getElementList().empty()){ 
    cout << "\tno stats, tree is empty -> done" << endl;
    return;
  }
  cout.precision(3);
  cout << "\tnumber of requests      : " << searchRequests_ << endl;
  cout << "\tmean search path length : " << fixed << SearchLength_ / (double)searchRequests_ << endl;
  cout << "\tmin search path length  : " << MinSearchLength_+1 << endl;
  cout << "\tcandidates inspected    : " << Candidates_ << "\t(expensive routines)" << endl;
  cout << "\tcandidates per request  : " << fixed << Candidates_/(double)searchRequests_ << endl;
  cout << "\tfactor to linear        : " << fixed << treeRoot_->getElementList().size()/(Candidates_/(double)searchRequests_) << " (with respect to candidates)"<< endl;
  cout << "\telements in tree        : " << treeRoot_->getElementList().size() << endl;
  cout << "\tnodes in tree           : " << treeRoot_->getNodesInTree(-1) << endl;
  cout << "\tleaf nodes in tree      : " << treeRoot_->getNodesInTree(1) << endl;
  cout << "\ttree depth              : " << treeRoot_->getDepth() << " (max: "<< MAX_TREEDEPTH<<")"<< endl;
  cout << "\tmem size                : " << fixed << getMemoryUsage()/1024.0 << " kb"<<endl;  
  stringstream ts;
  Teuchos::TimeMonitor::summarize(ts);
  string s = ts.str().substr(ts.str().find("XSearchTree")+43,10);
  s = s.substr(0,s.find(" "));
  double time =  atof(s.data()); 
  cout << "\toverall answer time     : " <<fixed << time << " secs" << endl;
  cout << "\tmean answer time        : " << time/(double)searchRequests_ << " secs" << endl;
  cout << "\t***************************************************" << endl;
#else
  cout << "\tno statistics in fast mode, please use DEBUG flag." <<endl;
#endif
}



#endif  // #ifdef CCADISCRET
