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

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}
extern struct _FILES  allfiles;

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
  treeRoot_(NULL),
  nodesWithReTry_(0)
{
//  cout << " CREATING A TREE" << endl;
}

XFEM::XSearchTree::XSearchTree(const BlitzMat3x2 AABB) :
  MAX_OVERLAP_ALLOWED_(0.875),
  ELEMENTS_USED_FOR_OVERLAP_CHECK_(1),
  TreeInit_(false),
  hasExternAABB_(true),
  AABB_(AABB),
  searchRequests_(0),
  Candidates_(0),
  SearchLength_(0),
  MinSearchLength_(INT_MAX),
  treeRoot_(NULL),
  nodesWithReTry_(0)
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

int XFEM::XSearchTree::queryPointType(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions, const BlitzVec3& pointCoords, int& closestElementId, double& distance) {
  //  cout << " ASKING THE TREE" << endl;
  closestElementId = 0;
  TEUCHOS_FUNC_TIME_MONITOR("XSearchTree - queryTime"); 
  searchRequests_++;
  if (dis.NumGlobalElements() == 0){ 
    return 0;
    }
  if (!TreeInit_)
    rebuild(dis, currentpositions);
  // search for candidates in tree
  
  BlitzMat3x2 testAABB;
  if (hasExternAABB_)
    testAABB = AABB_;
  else
    testAABB = treeRoot_->getAABB();
  for (int i = 0; i<3;i++){
      if ( (pointCoords(i)< testAABB(i,0)) || (pointCoords(i)>testAABB(i,1)) ){
        dserror("Please specify a bounding box that covers all locations doing requests on the tree!");
      }
    }
  
  int closestElementID;
  return treeRoot_->queryPointType(dis, currentpositions, pointCoords, closestElementID);
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
    treeRoot_ = new TreeNode(0,AABB_, this, NULL);
  }   
  else {
    const BlitzMat3x2 aabb(getXAABBofDis(dis, currentpositions));
    treeRoot_ = new TreeNode(0,aabb, this, NULL); 
  }  
  #ifdef DEBUG
  cout << endl <<"inserting new elements (" << dis.NumMyColElements() << ")"<< endl;
  
  std::stringstream filename;
  std::stringstream fc;
  int step = 0;
  
  
  filename << allfiles.outputfile_kenner << "_xaabbs" << std::setw(5) << setfill('0') << step << ".pos";
  cout << endl << "writing... "<<filename.str()<<" ...";
  fc << "View \" " << "XAABB of Elements \" {" << endl;
  flush(cout);
  #endif

  
  for (int i=0; i<dis.NumMyColElements(); ++i) {
    insertElement(dis.lRowElement(i),currentpositions);
    
    #ifdef DEBUG
    const BlitzMat xyze(XFEM::getCurrentNodalPositions(dis.lRowElement(i), currentpositions));
    const BlitzMat3x2 elemXAABB = XFEM::computeFastXAABB(dis.lRowElement(i), xyze, HIGHERORDER);
    BlitzMat XAABB(3,8);
    XAABB(0,0) = elemXAABB(0,0); XAABB(1,0) = elemXAABB(1,0);XAABB(2,0) = elemXAABB(2,0);
    XAABB(0,1) = elemXAABB(0,0); XAABB(1,1) = elemXAABB(1,1);XAABB(2,1) = elemXAABB(2,0);
    XAABB(0,2) = elemXAABB(0,0); XAABB(1,2) = elemXAABB(1,1);XAABB(2,2) = elemXAABB(2,1);
    XAABB(0,3) = elemXAABB(0,0); XAABB(1,3) = elemXAABB(1,0);XAABB(2,3) = elemXAABB(2,1);
    XAABB(0,4) = elemXAABB(0,1); XAABB(1,4) = elemXAABB(1,0);XAABB(2,4) = elemXAABB(2,0);
    XAABB(0,5) = elemXAABB(0,1); XAABB(1,5) = elemXAABB(1,1);XAABB(2,5) = elemXAABB(2,0);
    XAABB(0,6) = elemXAABB(0,1); XAABB(1,6) = elemXAABB(1,1);XAABB(2,6) = elemXAABB(2,1);
    XAABB(0,7) = elemXAABB(0,1); XAABB(1,7) = elemXAABB(1,0);XAABB(2,7) = elemXAABB(2,1);
    fc << IO::GMSH::cellWithScalarToString(DRT::Element::hex8, dis.lRowElement(i)->Id(), XAABB)<< endl;
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

XFEM::XSearchTree::TreeNode::TreeNode(const int Depth, const BlitzMat3x2& aabb, XSearchTree* tree, TreeNode* parent) :
  State_(STATE_LEAF_NODE),
  dont_refine_(false),
  labelID_(-1),  
  actTreedepth_(Depth),
  tree_(tree),
  parent_(parent),
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

int XFEM::XSearchTree::TreeNode::queryPointType(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions, const BlitzVec3& pointCoords, int& closestElementID) {
  //  printf("AABB(%f\t%f\t%f\t%f\t%f\t%f)\t", AABB(0,0),AABB(0,1),AABB(1,0),AABB(1,1),AABB(2,0),AABB(2,1));
  //  printf("x_in(%f\t%f\t%f)\n",pointCoords(0), pointCoords(1),pointCoords(2));
  switch (State_) {
    case STATE_INNER_NODE:         // if inner node, do recursion 
      //    printf("*");
      //    printf("classified searchpoint to oct %d, so i will search there\n", classifyPoint(pointCoords)-1);
      return children_[classifyPoint(pointCoords)-1]->queryPointType(dis, currentpositions, pointCoords ,closestElementID);
      break;
    case STATE_LEAF_NODE:   // returns list of candidates
      //    printf("/");
      if (ElementList_.empty()){
        //      printf("empty ");
        closestElementID= -1;
        tree_->MinSearchLength_ = std::min(tree_->MinSearchLength_, actTreedepth_);
        tree_->SearchLength_ = tree_->SearchLength_ + actTreedepth_+1;
        return labelID_;  // just return empty list if there are no candidates (==fluid) 
      }
      
      if (actTreedepth_ >= MAX_TREEDEPTH){
        //      printf("actTreedepth >= MAX_TREEDEPTH\n");
        tree_->MinSearchLength_ = std::min(tree_->MinSearchLength_, actTreedepth_);
        tree_->SearchLength_ = tree_->SearchLength_ + actTreedepth_+1;
        return getXFEMLabelOfPointInTreeNode(dis, currentpositions, pointCoords); 
      }
           
      if (dont_refine_){       //
        tree_->MinSearchLength_ = std::min(tree_->MinSearchLength_, actTreedepth_);
        tree_->SearchLength_ = tree_->SearchLength_ + actTreedepth_+1;
        return getXFEMLabelOfPointInTreeNode(dis, currentpositions, pointCoords);
        // just return XFEM-Label if node was previously set not to refine any further  
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
          const double AABBcoverage = overlap/XFEM::getVolume(AABB_);
           if (AABBcoverage > tree_->MAX_OVERLAP_ALLOWED_){
            dont_refine_ = true;
            tree_->MinSearchLength_ = std::min(tree_->MinSearchLength_, actTreedepth_);
            tree_->SearchLength_ = tree_->SearchLength_ + actTreedepth_+1;
            return getXFEMLabelOfPointInTreeNode(dis, currentpositions, pointCoords);;
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
          children_[i] = new TreeNode(actTreedepth_ + 1, chldAABB, tree_, this);
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
            int XFEMLabel=-1;
            const BlitzVec3 childNodeCenter(children_[i]->getCenterCoord());
            XFEMLabel = getXFEMLabelOfPointInTreeNode(dis, currentpositions, childNodeCenter);
            if (XFEMLabel>0) {
              children_[i]->setSolid(XFEMLabel);
            }
            else {
//              cout << "fluid octant" << endl;
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
        return children_[classifyPoint(pointCoords)-1]->queryPointType(dis, currentpositions, pointCoords, closestElementID);
      }
      else  // if there is only one Element, just return it
        tree_->MinSearchLength_ = std::min(tree_->MinSearchLength_, actTreedepth_);
        tree_->SearchLength_ = tree_->SearchLength_ + actTreedepth_+1;
        return getXFEMLabelOfPointInTreeNode(dis, currentpositions, pointCoords);
      break;
  }
  dserror("should not get here\n");
  return -1;
}

int XFEM::XSearchTree::TreeNode::getXFEMLabelOfPointInTreeNode(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions, const BlitzVec3& pointCoords){
  if (ElementList_.size()==0) 
    dserror("should not have an empty list here!");
  double dista=0;
  const DRT::Element* closestEle;
  closestEle = XFEM::nearestNeighbourInList(dis, currentpositions, ElementList_, pointCoords,dista);
  double SearchRadius = fabs(dista);
  TreeNode* workingNode = this;
  double AABBradius;
   //cout << "haha" << endl;
  do {
    const BlitzMat3x2 AABB =workingNode->getAABB();
    for (int i=0;i<3;i++)
      AABBradius =max(AABBradius, fabs(AABB(0,1)-AABB(0,0)));
    tree_->Candidates_+=workingNode->getElementList().size();
    closestEle = XFEM::nearestNeighbourInList(dis, currentpositions, workingNode->getElementList(), pointCoords,dista);
    workingNode=workingNode->getParent();
  } while (workingNode!=NULL && workingNode->hasParent() && SearchRadius>AABBradius);
  
  if (dista<=0) {
    const int eleId = closestEle->Id();
//    cout << "classified as solid: "<<tree_->getLabelByElementID(eleId)<< endl;
    return (tree_->getLabelByElementID(eleId)); //solid
  }
  else {
//    cout << "classified as fluid"<< endl;
    return 0;  //fluid
  }
  dserror("should not get here");
  return -1;
}

list< const DRT::Element* > XFEM::XSearchTree::getCandidateList(const BlitzVec3& X,const double radius){
  BlitzMat3x2 AASearchRegion = XFEM::getAABBofSphere(X, radius);
  return treeRoot_->getCandidateList(AASearchRegion);
}

list< const DRT::Element* > XFEM::XSearchTree::TreeNode::getCandidateList(const BlitzMat3x2& AASearchRegion){
  if (State_ == STATE_LEAF_NODE || isContainedAinB(AASearchRegion, AABB_)){
    return ElementList_;
  }
  else // means INNER_NODE
  {  
    list< const DRT::Element* > ElemList;
    const list<int> childIdx(classifyAABB(AASearchRegion)); //which child nodes will be examined   
    for (list<int>::const_iterator myIt = childIdx.begin(); myIt != childIdx.end(); myIt++){
      list< const DRT::Element* > ElementList2Add =  children_[*myIt-1]->getCandidateList(AASearchRegion);
      ElemList.merge(ElementList2Add);
      ElemList.unique(compare_gId);
    }
    return ElemList;
  }
  dserror("should not get here");
}

void XFEM::XSearchTree::TreeNode::insertElement(
    const DRT::Element* elem,
    const map<int,BlitzVec3>& currentpositions
    ) {
  if ((actTreedepth_ >= XFEM::XSearchTree::MAX_TREEDEPTH) || (State_ == STATE_LEAF_NODE) ) {
    ElementList_.push_back(elem);
    ElementList_.sort(compare_gId);
//    cout << "inserted element at depth " << actTreedepth_ <<endl;
    State_ = STATE_LEAF_NODE;
  } else if(State_ == STATE_INNER_NODE) {
    ElementList_.push_back(elem);
    const list<int> childIdx(classifyElement(elem,currentpositions));
    for (list<int>::const_iterator myIt = childIdx.begin(); myIt != childIdx.end(); myIt++){
      this->children_[*myIt-1]->insertElement(elem,currentpositions);
//      BlitzMat3x2 ab = this->children_[*myIt-1]->getAABB();
//      printf("inserted elem to AABB(%f,%f,%f,%f,%f,%f)\n", ab(0,0),ab(0,1),ab(1,0),ab(1,1),ab(2,0),ab(2,1));
     }
    ElementList_.sort(compare_gId);
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
    ) 
{
  
  const BlitzMat xyze(XFEM::getCurrentNodalPositions(elem,currentpositions));
  const BlitzMat3x2 elemXAABB(XFEM::computeFastXAABB(elem, xyze, HIGHERORDER));
  
  return classifyAABB(elemXAABB);
  
}

list<int> XFEM::XSearchTree::TreeNode::classifyAABB(
    const BlitzMat3x2 AABB
    ) {
  list<int> octants;
 
  if (AABB(0, 1) > XPlaneCoordinate_) {
    if (AABB(1, 1) > YPlaneCoordinate_) {
      if (AABB(2, 1) > ZPlaneCoordinate_){
        octants.push_back(8);
//        do_refine[8] = isAABBbiggerThanElemXAABB;
      }
      if (AABB(2, 0) <= ZPlaneCoordinate_){
        octants.push_back(7);
//        do_refine[7] = isAABBbiggerThanElemXAABB;
      }
    }
    if (AABB(1, 0) <= YPlaneCoordinate_) {
      if (AABB(2, 1) > ZPlaneCoordinate_){
        octants.push_back(6);
//        do_refine[6] = isAABBbiggerThanElemXAABB;
      }
      if (AABB(2, 0) <= ZPlaneCoordinate_){
        octants.push_back(5);
//        do_refine[5] = isAABBbiggerThanElemXAABB;
      }
    }
  }
  if (AABB(0, 0) <= XPlaneCoordinate_) {
    if (AABB(1, 1) > YPlaneCoordinate_) {
      if (AABB(2, 1) > ZPlaneCoordinate_){
        octants.push_back(4);
//        do_refine[4] = isAABBbiggerThanElemXAABB;
      }
      if (AABB(2, 0) <= ZPlaneCoordinate_){
        octants.push_back(3);
//        do_refine[3] = isAABBbiggerThanElemXAABB;
      }
    }
    if (AABB(1, 0) <= YPlaneCoordinate_) {
      if (AABB(2, 1) > ZPlaneCoordinate_){
        octants.push_back(2);
//        do_refine[2] = isAABBbiggerThanElemXAABB;
      }
      if (AABB(2, 0) <= ZPlaneCoordinate_){
        octants.push_back(1);
//        do_refine[1] = isAABBbiggerThanElemXAABB;
      }
    }
    
  }
//    printf("classifiing AABB: ", elemXAABB(0,0),elemXAABB(0,1),elemXAABB(1,0),elemXAABB(1,1),elemXAABB(2,0),elemXAABB(2,1));
//    printf("classified AABB to octants: ");
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

XFEM::XSearchTree::TreeNode* XFEM::XSearchTree::TreeNode::getParent(){
  if (this->hasParent())
    return parent_;
  return NULL;
}

bool XFEM::XSearchTree::TreeNode::hasParent() const {
  if (parent_!=NULL)
    return true;
  return false;
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

void XFEM::XSearchTree::printTree(const string prefix, const int step) const{
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
  filename << prefix << "_tree" << std::setw(5) << setfill('0') << step << ".pos";
  cout << " "<<filename.str()<<" ...";
  fc << "View \" " << "fsiOctree \" {" << endl;  
  treeRoot_->printTree(fc);
  fc << "};" << endl;
  std::ofstream f_system(filename.str().c_str());
  f_system << fc.str();
  f_system.close();
  cout << " done" << endl;
}

void XFEM::XSearchTree::TreeNode::printTree(stringstream& fc) const{
  int factor = 1;
  if (actTreedepth_==-1)
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
      fc << IO::GMSH::cellWithScalarToString(DRT::Element::hex8, 0, XAABB)<< endl;
    }

  if (State_==STATE_INNER_NODE){
    for (int j=0; j<8; j++){
      if (children_[j]!=NULL)
        children_[j]->printTree(fc);
    }	  
  }
  else if (State_==STATE_LEAF_NODE)
  {
    if (labelID_<0)
      factor = 0; // more than one candidate in this octant
    else if (labelID_==0)
      factor = 1; // fluid octant
    else
      factor = 2; // solid octant
    BlitzMat XAABB(3,8);
    XAABB(0,0) = AABB_(0,0); XAABB(1,0) = AABB_(1,0);XAABB(2,0) = AABB_(2,0);
    XAABB(0,1) = AABB_(0,0); XAABB(1,1) = AABB_(1,1);XAABB(2,1) = AABB_(2,0);
    XAABB(0,2) = AABB_(0,0); XAABB(1,2) = AABB_(1,1);XAABB(2,2) = AABB_(2,1);
    XAABB(0,3) = AABB_(0,0); XAABB(1,3) = AABB_(1,0);XAABB(2,3) = AABB_(2,1);
    XAABB(0,4) = AABB_(0,1); XAABB(1,4) = AABB_(1,0);XAABB(2,4) = AABB_(2,0);
    XAABB(0,5) = AABB_(0,1); XAABB(1,5) = AABB_(1,1);XAABB(2,5) = AABB_(2,0);
    XAABB(0,6) = AABB_(0,1); XAABB(1,6) = AABB_(1,1);XAABB(2,6) = AABB_(2,1);
    XAABB(0,7) = AABB_(0,1); XAABB(1,7) = AABB_(1,0);XAABB(2,7) = AABB_(2,1);
    fc << IO::GMSH::cellWithScalarToString(DRT::Element::hex8, factor*tree_->MAX_TREEDEPTH+actTreedepth_, XAABB)<< endl;
  }
  
}

void XFEM::XSearchTree::printTreeMetrics(const int step) const{
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
  cout << "\tretry nodes             : " << nodesWithReTry_ << endl;
  cout << "\tmem size                : " << fixed << getMemoryUsage()/1024.0 << " kb"<<endl;  
  stringstream ts;
  Teuchos::TimeMonitor::summarize(ts);
  string s = ts.str().substr(ts.str().find("XSearchTree")+43,10);
  s = s.substr(0,s.find(" "));
  double time =  atof(s.data()); 
  cout << "\toverall answer time     : " <<fixed << time << " secs" << endl;
  cout << "\tmean answer time        : " << time/(double)searchRequests_ << " secs" << endl;
  cout << "\t***************************************************" << endl;
}

bool XFEM::compare_gId(const DRT::Element* A, const DRT::Element* B)
{
  if ( (A->Id())<(B->Id()) ) 
    return true;
  else 
    return false;
}


#endif  // #ifdef CCADISCRET
