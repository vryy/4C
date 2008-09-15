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
#include "../drt_geometry/intersection_service.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/standardtypes_cpp.H"
#include <Teuchos_TimeMonitor.hpp>

extern struct _FILES  allfiles;

using namespace std;

XFEM::XSearchTree::XSearchTree():
  MAX_TREEDEPTH_(8),
  MAX_OVERLAP_ALLOWED_(0.875),
  ELEMENTS_USED_FOR_OVERLAP_CHECK_(1),
  hasExternAABB_(true),
  TreeInit_(false),
  searchRequests_(0),
  Candidates_(0),
  SearchLength_(0),
  MinSearchLength_(INT_MAX),
  treeRoot_(NULL),
  rootLevelSearches_(0){
}

XFEM::XSearchTree::XSearchTree(
    const BlitzMat3x2 AABB,
    const int MAX_TREEDEPTH,
    const int ELEMENTS_USED_FOR_OVERLAP_CHECK,
    const double MAX_OVERLAP_ALLOWED):
      MAX_TREEDEPTH_(MAX_TREEDEPTH),
      MAX_OVERLAP_ALLOWED_(MAX_OVERLAP_ALLOWED),
      ELEMENTS_USED_FOR_OVERLAP_CHECK_(ELEMENTS_USED_FOR_OVERLAP_CHECK),
      AABB_(AABB),
      hasExternAABB_(true),
      TreeInit_(false),
      searchRequests_(0),
      Candidates_(0),
      SearchLength_(0),
      MinSearchLength_(INT_MAX),
      treeRoot_(NULL),
      rootLevelSearches_(0){
}

XFEM::XSearchTree::~XSearchTree() {
  delete treeRoot_;
}

void XFEM::XSearchTree::setRebuildFlag() {
  TreeInit_ = false;
}

int XFEM::XSearchTree::getTotalRequests() const{
  return searchRequests_;
}
double XFEM::XSearchTree::getMeanSearchlength() const{
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

int XFEM::XSearchTree::getDepth() const{
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
  for (int i=0; i<dis.NumMyColElements(); ++i) {
    insertElement(dis.lRowElement(i),currentpositions);
  }

    
  #ifdef DEBUG5
  cout << endl <<"inserting new elements (" << dis.NumMyColElements() << ")"<< endl;
  
  std::stringstream filename;
  std::stringstream fc;
  int step = 0;
  
  
  filename << allfiles.outputfile_kenner << "_xaabbs" << std::setw(5) << setfill('0') << step << ".pos";
  cout << endl << "writing... "<<filename.str()<<" ...";
  fc << "View \" " << "XAABB of Elements \" {" << endl;
  flush(cout);
  for (int i=0; i<dis.NumMyColElements(); ++i) {
    const BlitzMat xyze(DRT::UTILS::getCurrentNodalPositions(dis.lRowElement(i), currentpositions));
    const XFEM::EleGeoType eleGeoType(HIGHERORDER);
    checkRoughGeoType(*myIt, xyze, eleGeoType);
    const BlitzMat3x2 elemXAABB = XFEM::computeFastXAABB(dis.lRowElement(i), xyze, eleGeoType);
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
  }
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

XFEM::XSearchTree::TreeNode::TreeNode(
    const int Depth, 
    const BlitzMat3x2& aabb, 
    XSearchTree* tree, 
    const TreeNode* const parent) :
  TreeNodeType_(LEAF_NODE),
  dont_refine_(false),
  labelID_(-1),  
  actTreedepth_(Depth),
  tree_(tree),
  parent_(parent),
  AABB_(aabb),
  XPlaneCoordinate_(aabb(0,0)+(aabb(0,1)-aabb(0,0))/2.0),
  YPlaneCoordinate_(aabb(1,0)+(aabb(1,1)-aabb(1,0))/2.0),
  ZPlaneCoordinate_(aabb(2,0)+(aabb(2,1)-aabb(2,0))/2.0)
  {
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

int XFEM::XSearchTree::TreeNode::getDepth() const{
  if (TreeNodeType_==LEAF_NODE){
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

int XFEM::XSearchTree::TreeNode::getNodesInTree(TreeNodeType type) const{
  int nodes=0;
  int tmp=0;
  if (TreeNodeType_==LEAF_NODE){
    switch (type) {
      case INNER_NODE:
        return 0;
      case LEAF_NODE:
        return 1;
      case ANY_NODE:
        return 1;
    }
    return 0;
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
      case INNER_NODE:  nodes++;break;
      case LEAF_NODE:           break;
      case ANY_NODE:    nodes++;break;
    }
  }
  return nodes;
}

void XFEM::XSearchTree::TreeNode::setLabel(const int label){
  labelID_ = label;
}

BlitzVec3 XFEM::XSearchTree::TreeNode::getCenterCoord() const{
  return BlitzVec3(this->XPlaneCoordinate_, this->YPlaneCoordinate_, this->ZPlaneCoordinate_);
}

int XFEM::XSearchTree::TreeNode::queryPointType(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions, const BlitzVec3& pointCoords, int& closestElementID) {
  switch (TreeNodeType_) {
  case ANY_NODE: break;    // just to have this in enum
  case INNER_NODE:         // if inner node, do recursion 
    return children_[classifyPoint(pointCoords)-1]->queryPointType(dis, currentpositions, pointCoords ,closestElementID);
    break;
  case LEAF_NODE:   // returns list of candidates
    //    printf("/");
    if (ElementList_.empty()){
      //      printf("empty ");
      closestElementID= -1;
      tree_->MinSearchLength_ = std::min(tree_->MinSearchLength_, actTreedepth_);
      tree_->SearchLength_ = tree_->SearchLength_ + actTreedepth_+1;
      return labelID_;  // just return empty list if there are no candidates (==fluid) 
    }

    if (actTreedepth_ >= tree_->MAX_TREEDEPTH_){
      //      printf("actTreedepth >= tree_->MAX_TREEDEPTH_\n");
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

    if (ElementList_.size()==1) // if there is only one Element, just return it
    {
      tree_->MinSearchLength_ = std::min(tree_->MinSearchLength_, actTreedepth_);
      tree_->SearchLength_ = tree_->SearchLength_ + actTreedepth_+1;
      return getXFEMLabelOfPointInTreeNode(dis, currentpositions, pointCoords);
    }
    else // dynamically grow tree        
    { 
      if (ElementList_.size()<=tree_->ELEMENTS_USED_FOR_OVERLAP_CHECK_){
        doOverlapCheck(currentpositions);
        if (dont_refine_){
          tree_->MinSearchLength_ = std::min(tree_->MinSearchLength_, actTreedepth_);
          tree_->SearchLength_ = tree_->SearchLength_ + actTreedepth_+1;
          return getXFEMLabelOfPointInTreeNode(dis, currentpositions, pointCoords);;
        }
      }
      // create Octants and check if any of the created childs contains no elements 
      // and for empty children determine their XFEM-label
      createChildren(currentpositions);
      setLabelsForCandidateFreeChildren(dis, currentpositions);

      // search in apropriate child node
      return children_[classifyPoint(pointCoords)-1]->queryPointType(dis, currentpositions, pointCoords, closestElementID);
    }
    break;
  }
  dserror("should not get here\n");
  return -1;
}

void XFEM::XSearchTree::TreeNode::doOverlapCheck(const std::map<int,BlitzVec3>& currentpositions){
  list<BlitzMat3x2> AABBs; 
  AABBs.clear(); 
  AABBs.push_back(AABB_);
  for (list< const DRT::Element* >::const_iterator myIt = ElementList_.begin(); myIt != ElementList_.end(); myIt++){
    const BlitzMat xyze(GEO::getCurrentNodalPositions(*myIt,currentpositions));
    GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
    checkRoughGeoType(*myIt, xyze, eleGeoType);
    AABBs.push_back(GEO::computeFastXAABB(*myIt, xyze, eleGeoType));
  }
  const double AABBcoverage = XFEM::getOverlapArea(AABBs)/XFEM::getVolume(AABB_);
  if (AABBcoverage > tree_->MAX_OVERLAP_ALLOWED_){
    dont_refine_ = true;
  }

}

void XFEM::XSearchTree::TreeNode::createChildren(const std::map<int,BlitzVec3>& currentpositions){
  map<const DRT::Element*, list<int> > ElementClassification;
  
  for (list< const DRT::Element* >::const_iterator myIt = ElementList_.begin(); myIt != ElementList_.end(); myIt++){
    ElementClassification[*myIt]=classifyElement(*myIt,currentpositions);
  }

  for (int i=0; i<8; i++){
    children_[i] = new TreeNode(actTreedepth_ + 1, getChildOctAABB(i+1), tree_, this);
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
  // this node becomes an inner tree node
  TreeNodeType_ = INNER_NODE;
}


void XFEM::XSearchTree::TreeNode::setLabelsForCandidateFreeChildren(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions){
  // if one of the created childs is empty, check if it is fluid or solid
  for (int i=0; i< 8; i++){
    if ((children_[i]->getElementList()).empty()){
      int XFEMLabel=-1;
      const BlitzVec3 childNodeCenter(children_[i]->getCenterCoord());
      XFEMLabel = getXFEMLabelOfPointInTreeNode(dis, currentpositions, childNodeCenter);
      if (XFEMLabel>0) {
        children_[i]->setLabel(XFEMLabel);
      }
      else {
//              cout << "fluid octant" << endl;
        children_[i]->setLabel(0);
      }
    }
    else {
      //          printf("non empty octant(#%d) with %d elements\n", i+1, children_[i]->getElementList().size());
    }
  }

}

int XFEM::XSearchTree::TreeNode::getXFEMLabelOfPointInTreeNode(const DRT::Discretization& dis,const std::map<int,BlitzVec3>& currentpositions, const BlitzVec3& X){
  map<const DRT::Element*, double > squaredDistanceMap;
  squaredDistanceMap.clear();
  if (ElementList_.size()==0) 
    dserror("should not have an empty list here!");
  double distance=0.0;
  BlitzVec3 vectorX2minDistPoint(0.0,0.0,0.0);
  const DRT::Element* closestEle;
  const TreeNode* workingNode = this;
  XFEM::DistanceType distanceType;
  int tmpCandidateCounter=0;
  BlitzVec3 minDistPoint(0.0,0.0,0.0);
  BlitzMat3x2 workAABB;
  do {
    if (!workingNode->hasParent()){
      tree_->rootLevelSearches_++;  
    }
    tmpCandidateCounter = workingNode->getElementList().size();
    closestEle=XFEM::nearestNeighbourInList(dis, currentpositions, workingNode->getElementList(), X, squaredDistanceMap, minDistPoint, distance, distanceType);
    workAABB = workingNode->getAABB();
    if (XFEM::isContainedXinAABB(workAABB, minDistPoint)) {
      tree_->Candidates_+=workingNode->getElementList().size();
      if (distance<0){ 
        return (tree_->getLabelByElementID(closestEle->Id())); //solid
      }
      else {
        return 0;  //fluid
      }
    }
    workingNode=workingNode->getParent();
    
  } while (workingNode!=NULL);
  dserror("should not get here");
  return -1;
}

void XFEM::XSearchTree::TreeNode::insertElement(
    const DRT::Element* elem,
    const map<int,BlitzVec3>& currentpositions
    ) {
  if ((actTreedepth_ >= tree_->MAX_TREEDEPTH_) || (TreeNodeType_ == LEAF_NODE) ) {
    ElementList_.push_back(elem);
//    cout << "inserted element at depth " << actTreedepth_ <<endl;
    TreeNodeType_ = LEAF_NODE;
  } else if(TreeNodeType_ == INNER_NODE) {
    ElementList_.push_back(elem);
    const list<int> childIdx(classifyElement(elem,currentpositions));
    for (list<int>::const_iterator myIt = childIdx.begin(); myIt != childIdx.end(); myIt++){
      this->children_[*myIt-1]->insertElement(elem,currentpositions);
//      BlitzMat3x2 ab = this->children_[*myIt-1]->getAABB();
//      printf("inserted elem to AABB(%f,%f,%f,%f,%f,%f)\n", ab(0,0),ab(0,1),ab(1,0),ab(1,1),ab(2,0),ab(2,1));
     }
  } 
}

list< const DRT::Element* > XFEM::XSearchTree::TreeNode::getElementList() const{
    return ElementList_;
}

double XFEM::XSearchTree::getMaxOverlap() const{
  return MAX_OVERLAP_ALLOWED_;
}

BlitzMat3x2 XFEM::XSearchTree::TreeNode::getChildOctAABB(const int i) const{
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

int XFEM::XSearchTree::TreeNode::classifyPoint(const BlitzVec3& pointcoords) const{
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
    ) const
{
  
  const BlitzMat xyze(GEO::getCurrentNodalPositions(elem,currentpositions));
  GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
  checkRoughGeoType(elem, xyze, eleGeoType);
  const BlitzMat3x2 elemXAABB(GEO::computeFastXAABB(elem, xyze, eleGeoType));  
  return classifyAABB(elemXAABB);
  
}

list<int> XFEM::XSearchTree::TreeNode::classifyAABB(
    const BlitzMat3x2 AABB
    ) const 
{
  list<int> octants;
 
  if (AABB(0, 1) > XPlaneCoordinate_) {
    if (AABB(1, 1) > YPlaneCoordinate_) {
      if (AABB(2, 1) > ZPlaneCoordinate_){
        octants.push_back(8);
      }
      if (AABB(2, 0) <= ZPlaneCoordinate_){
        octants.push_back(7);
      }
    }
    if (AABB(1, 0) <= YPlaneCoordinate_) {
      if (AABB(2, 1) > ZPlaneCoordinate_){
        octants.push_back(6);
      }
      if (AABB(2, 0) <= ZPlaneCoordinate_){
        octants.push_back(5);
      }
    }
  }
  if (AABB(0, 0) <= XPlaneCoordinate_) {
    if (AABB(1, 1) > YPlaneCoordinate_) {
      if (AABB(2, 1) > ZPlaneCoordinate_){
        octants.push_back(4);
      }
      if (AABB(2, 0) <= ZPlaneCoordinate_){
        octants.push_back(3);
      }
    }
    if (AABB(1, 0) <= YPlaneCoordinate_) {
      if (AABB(2, 1) > ZPlaneCoordinate_){
        octants.push_back(2);
      }
      if (AABB(2, 0) <= ZPlaneCoordinate_){
        octants.push_back(1);
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

XFEM::XSearchTree::TreeNode::TreeNodeType XFEM::XSearchTree::TreeNode::getTreeNodeType() const{
  return TreeNodeType_;
}

int XFEM::XSearchTree::getLabelByElementID(
    const int gid) const
{
  return labelByElement_.find(gid)->second;
}

const BlitzMat3x2& XFEM::XSearchTree::TreeNode::getAABB() const
{
  return AABB_;
}

XFEM::XSearchTree::TreeNode* XFEM::XSearchTree::TreeNode::getChild(
    const int idx) const
{
  return children_[idx-1];
}

const XFEM::XSearchTree::TreeNode* const XFEM::XSearchTree::TreeNode::getParent() const
{
  if (this->hasParent())
    return parent_;
  return NULL;
}

bool XFEM::XSearchTree::TreeNode::hasParent() const 
{
  if (parent_!=NULL)
    return true;
  return false;
}

int XFEM::XSearchTree::getMemoryUsage() const
{
  int mem = treeRoot_->getMemoryUsage();
  mem = mem + sizeof this;
  return mem;
}

int XFEM::XSearchTree::TreeNode::getMemoryUsage() const
{
  int mem=0;
  if (TreeNodeType_==LEAF_NODE){
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

void XFEM::XSearchTree::printTree(const string prefix, const int step) const
{
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

void XFEM::XSearchTree::TreeNode::printTree(stringstream& fc) const
{
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

  if (TreeNodeType_==INNER_NODE){
    for (int j=0; j<8; j++){
      if (children_[j]!=NULL)
        children_[j]->printTree(fc);
    }
  }
  else if (TreeNodeType_==LEAF_NODE)
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
    fc << IO::GMSH::cellWithScalarToString(DRT::Element::hex8, factor*tree_->MAX_TREEDEPTH_+actTreedepth_, XAABB)<< endl;
  }
  
}

void XFEM::XSearchTree::printTreeMetricsFile(const int step) const
{
  if (!TreeInit_) {
     cout << "nothing to write, tree not initialized yet -> done" << endl;
     return;
   }
  std::stringstream filename;
  filename << allfiles.outputfile_kenner << "_treeMetrics_" << std::setw(5) << setfill('0') << step << ".txt";
  bool fileExists = false;
  std::fstream fin(filename.str().c_str(),ios::in);
  if( fin.is_open() )
  {
  fileExists=true;
  }
  fin.close();
  
  std::fstream f_system(filename.str().c_str(), ios::out | ios::app);
  stringstream filecontentToAdd;
  cout << "writing " << left << std::setw(50) <<filename.str()<<"...";
  filecontentToAdd.precision(10);
  if (!fileExists){
    filecontentToAdd << "searchRequests" <<  "\t"
                     << "searchLength" << "\t"
                     << "minSearchLength"  << "\t"
                     << "rootLevelSearches_" << "\t"
                     << "candidates" << "\t"
                     << "elements" << "\t"
                     << "nodes" << "\t"
                     << "leafnodes" << "\t"
                     << "depth" << "\t"
                     << "max_depth" << "\t"
                     << "mem" << "\t"
                     << "time"<< endl;
  }
  stringstream ts;
  Teuchos::TimeMonitor::summarize(ts);
  string s = ts.str().substr(ts.str().find("XSearchTree")+43,10);
  s = s.substr(0,s.find(" "));
  double time =  atof(s.data()); 
  filecontentToAdd << searchRequests_ <<  "\t"
                   << SearchLength_ << "\t"
                   << MinSearchLength_+1  << "\t"
                   << rootLevelSearches_ << "\t"
                   << Candidates_ << "\t"
                   << treeRoot_->getElementList().size()<<"\t"
                   << treeRoot_->getNodesInTree(XFEM::XSearchTree::TreeNode::ANY_NODE) << "\t"
                   << treeRoot_->getNodesInTree(XFEM::XSearchTree::TreeNode::LEAF_NODE) << "\t"
                   << treeRoot_->getDepth() << "\t"
                   << MAX_TREEDEPTH_ << "\t"
                   << getMemoryUsage() << "\t"
                   << fixed << time << endl;
  
  f_system << filecontentToAdd.str();
  f_system.close();
  
  cout << " done" << endl;

}

void XFEM::XSearchTree::printTreeMetrics(const int step) const
{
  cout << "\t********************* TREE STATS ******************" << endl;
  if (!TreeInit_ ) {
    cout << "\tno stats, tree not initialized yet -> done" << endl;
    return;
  }
  if (treeRoot_->getElementList().empty()){ 
    cout << "\tno stats, tree is empty -> done" << endl;
    return;
  }
  cout.precision(5);
  cout << "\tnumber of requests      : " << searchRequests_ << endl;
  cout << "\tmean search path length : " << fixed << SearchLength_ / (double)searchRequests_ << endl;
  cout << "\tmin search path length  : " << MinSearchLength_+1 << endl;
  cout << "\tsearches to root level  : " << rootLevelSearches_ << endl;
  cout << "\tcandidates inspected    : " << Candidates_ << "\t(expensive routines)" << endl;
  cout << "\tcandidates per request  : " << fixed << Candidates_/(double)searchRequests_ << endl;
  cout << "\tfactor to linear        : " << fixed << treeRoot_->getElementList().size()/(Candidates_/(double)searchRequests_) << " (with respect to candidates)"<< endl;
  cout << "\telements in tree        : " << treeRoot_->getElementList().size() << endl;
  cout << "\tnodes in tree           : " << treeRoot_->getNodesInTree(XFEM::XSearchTree::TreeNode::ANY_NODE) << endl;
  cout << "\tleaf nodes in tree      : " << treeRoot_->getNodesInTree(XFEM::XSearchTree::TreeNode::LEAF_NODE) << endl;
  cout << "\ttree depth              : " << treeRoot_->getDepth() << " (max: "<< MAX_TREEDEPTH_<<")"<< endl;
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
#endif  // #ifdef CCADISCRET
