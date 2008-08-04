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
#include "octtree.H"
#include "../drt_xfem/intersection_service.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/drt_utils.H"
#include <Teuchos_TimeMonitor.hpp>

extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
}
extern struct _FILES  allfiles;



/*----------------------------------------------------------------------*
 | constructor OctTree                                       u.may 07/08|
 *----------------------------------------------------------------------*/
GEO::OctTree::OctTree(
    const int           max_depth,
    const BlitzMat3x2&  nodeBox): 
      max_depth_(max_depth),
      rootNodeBox_(nodeBox),
      treeRoot_(null)
{
}
    


/*----------------------------------------------------------------------*
| destructor OctTree                                        u.may 07/08|
*----------------------------------------------------------------------*/
GEO::OctTree::~OctTree(){}



/*----------------------------------------------------------------------*
 | initialize or rebuild tree with possibly new              u.may 07/08| 
 | discretization, elements are sorted according to the given map       |
 *----------------------------------------------------------------------*/
void GEO::OctTree::initializeTree(
    const std::map<int,std::set<int> >&  elementsByLabel) 
{
  if (treeRoot_ != Teuchos::null)
    treeRoot_ = Teuchos::null;

  treeRoot_ = rcp(new TreeNode(Teuchos::null, max_depth_, rootNodeBox_)); 

  if(elementsByLabel.size()==0)
    dserror("please provide a filled elemens by label list or use overloaded method without elementsByLabel");

  // insert element map into tree root node
  treeRoot_->setElementList(elementsByLabel);
}



/*----------------------------------------------------------------------*
 | initialize or rebuild tree with possibly new              u.may 07/08| 
 | discretization, elements are taken unsorted from discretization      |
 *----------------------------------------------------------------------*/
void GEO::OctTree::initializeTree(
    const DRT::Discretization&      dis) 
{
  if (treeRoot_ != Teuchos::null)
    treeRoot_ = Teuchos::null;

  treeRoot_ = rcp( new TreeNode(Teuchos::null, max_depth_, rootNodeBox_)); 

  for (int i=0; i<dis.NumMyColElements(); ++i) 
    treeRoot_->insertElement(0,dis.lRowElement(i)->Id());
}


/*----------------------------------------------------------------------*
 | update tree                                               u.may 08/08|
 | only usefull for labeled structures,for networks use initializeTree  |
 *----------------------------------------------------------------------*/
void GEO::OctTree::updateTree(
    const DRT::Discretization& 	dis,
    const std::map<int,BlitzVec3>& 	currentpositions_old, 
    const std::map<int,BlitzVec3>& 	currentpositions_new) 
{

  std::vector< BlitzMat3x2 > AABBs_old = GEO::computeXAABBForLabeledStructures(dis, currentpositions_old, treeRoot_->getElementList());
  std::vector< BlitzMat3x2 > AABBs_new = GEO::computeXAABBForLabeledStructures(dis, currentpositions_new, treeRoot_->getElementList());

  for(unsigned int i = 0; i < AABBs_new.size(); i++)
    treeRoot_->updateTreeNode(AABBs_old[i], AABBs_new[i]);
}



/*----------------------------------------------------------------------*
 | returns xfem label of point                               u.may 07/08|
 *----------------------------------------------------------------------*/
int GEO::OctTree::queryPointType(
    const DRT::Discretization& 	     dis,
    const std::map<int,BlitzVec3>& 	 currentpositions, 
    const BlitzVec3& 		             point) 
{
  //  cout << " ASKING THE TREE" << endl;
  TEUCHOS_FUNC_TIME_MONITOR("OctTree - queryTime"); 

  return treeRoot_->queryPointType(dis, currentpositions, point);
}



/*----------------------------------------------------------------------*
 | c-tor TreeNode                                            u.may 07/08|
 *----------------------------------------------------------------------*/
GEO::OctTree::TreeNode::TreeNode(
    const Teuchos::RCP<GEO::OctTree::TreeNode>  parent,
    const int                                   depth, 
    const BlitzMat3x2&                          nodeBox):
    parent_(parent),
    treedepth_(depth),
    treeNodeType_(LEAF_NODE),
    label_(-1),
    nodeBox_(nodeBox),
    xPlaneCoordinate_( (nodeBox_(0,0)+0.5*(nodeBox_(0,1)-nodeBox_(0,0))) ),
    yPlaneCoordinate_( (nodeBox_(1,0)+0.5*(nodeBox_(1,1)-nodeBox_(1,0))) ),
    zPlaneCoordinate_( (nodeBox_(2,0)+0.5*(nodeBox_(2,1)-nodeBox_(2,0))) )
{    
  std::vector< Teuchos::RCP<TreeNode> > children_(8, Teuchos::null);
  
}
                  



/*====================================================================*/
  /* bounding box topology*/
  /*--------------------------------------------------------------------*/
  /* parameter coordinates (x,y,z) of nodes
   * child 0: (minX, minY, minZ)
   * child 1: (maxX, minY, minZ)
   * child 2: (maxX, maxY, minZ)
   * child 3: (minX, maxY, minZ)
   * child 4: (minX, minY, maxZ)
   * child 5: (maxX, minY, maxZ)
   * child 6: (maxX, maxY, maxZ)
   * child 7: (minX, maxY, maxZ)
   * 
   *                      z
   *                      |           
   *             4========|================7
   *           //|        |               /||
   *          // |        |              //||
   *         //  |        |             // ||
   *        //   |        |            //  ||
   *       //    |        |           //   ||
   *      //     |        |          //    ||
   *     //      |        |         //     ||
   *     5=========================6       ||
   *    ||       |        |        ||      ||
   *    ||       |        o--------||---------y
   *    ||       |       /         ||      ||
   *    ||       0------/----------||------3
   *    ||      /      /           ||     //
   *    ||     /      /            ||    //
   *    ||    /      /             ||   //
   *    ||   /      /              ||  //
   *    ||  /      /               || //
   *    || /      x                ||//
   *    ||/                        ||/
   *     1=========================2
   *
   */
  /*====================================================================*/





/*----------------------------------------------------------------------*
 | d-tor TreeNode                                            u.may 07/08|
 *----------------------------------------------------------------------*/
GEO::OctTree::TreeNode::~TreeNode() {}



/*----------------------------------------------------------------------*
 | clears node and deletes chlidren                          u.may 08/08|
 *----------------------------------------------------------------------*/
void GEO::OctTree::TreeNode::clear()
{
  treeNodeType_ = LEAF_NODE;
  label_ = -1; 

  children_.clear();
  children_.assign(8,Teuchos::null);
}


/*----------------------------------------------------------------------*
 | has parent element                                        peder 07/08|
 *----------------------------------------------------------------------*/
bool GEO::OctTree::TreeNode::hasParent() const 
{
  if (parent_!=Teuchos::null)
    return true;
  return false;
}



/*----------------------------------------------------------------------*
 | sets element list                                         u.may 07/08|
 *----------------------------------------------------------------------*/
void GEO::OctTree::TreeNode::setElementList(
   const std::map<int, std::set<int> >& elementsByLabel)
{
   elementList_ = elementsByLabel;
}



/*----------------------------------------------------------------------*
 | set label                                                 peder 07/08|
 *----------------------------------------------------------------------*/
void GEO::OctTree::TreeNode::setLabel(
	const int label)
{
  label_ = label;
}



/*----------------------------------------------------------------------*
 | returns pointer to parent element                         peder 07/08|
 *----------------------------------------------------------------------*/
Teuchos::RCP<GEO::OctTree::TreeNode> GEO::OctTree::TreeNode::getParent() const
{
  if (this->hasParent())
    return parent_;
  return null;
}



/*----------------------------------------------------------------------*
 | get center of treenode                                  peder   07/08|
 *----------------------------------------------------------------------*/
BlitzVec3 GEO::OctTree::TreeNode::getCenterCoord() const
{
  return BlitzVec3(this->xPlaneCoordinate_, this->yPlaneCoordinate_, this->zPlaneCoordinate_);
}



/*----------------------------------------------------------------------*
 | get map of elements                                     peder   07/08|
 *----------------------------------------------------------------------*/
std::map<int,std::set<int> > GEO::OctTree::TreeNode::getElementList() const
{
  return elementList_;
}



/*----------------------------------------------------------------------*
 | get type of tree node                                   peder   07/08|
 *----------------------------------------------------------------------*/
GEO::TreeNodeType GEO::OctTree::TreeNode::getTreeNodeType() const
{
  return treeNodeType_;
}



/*----------------------------------------------------------------------*
 | get node box                                            peder   07/08|
 *----------------------------------------------------------------------*/
const BlitzMat3x2& GEO::OctTree::TreeNode::getNodeBox() const
{
  return nodeBox_;
}



/*----------------------------------------------------------------------*
 | get child                                               peder   07/08|
 *----------------------------------------------------------------------*/
Teuchos::RCP<GEO::OctTree::TreeNode> GEO::OctTree::TreeNode::getChild(
    const int index) const
{
  //TODO change to index
  // was index-1 before
  return children_[index];
}



/*----------------------------------------------------------------------*
 | get node box of a child specified by index              peder   07/08|
 *----------------------------------------------------------------------*/
BlitzMat3x2 GEO::OctTree::TreeNode::getChildNodeBox(
	const int index) const
{
  BlitzMat3x2 childNodeBox;
  if (index>3){
    childNodeBox(0,0) = xPlaneCoordinate_;
    childNodeBox(0,1) = nodeBox_(0,1);
  }
  else 
  {
    childNodeBox(0,0) = nodeBox_(0,0);
    childNodeBox(0,1) = xPlaneCoordinate_;
  }

  if ((index==2) || (index==3) || (index==6) || (index==7))
  {
    childNodeBox(1,0) = yPlaneCoordinate_;
    childNodeBox(1,1) = nodeBox_(1,1);
  }
  else 
  {
    childNodeBox(1,0) = nodeBox_(1,0);
    childNodeBox(1,1) = yPlaneCoordinate_;
  }

  if ((index%2)==1){
    childNodeBox(2,0) = zPlaneCoordinate_;
    childNodeBox(2,1) = nodeBox_(2,1);
  }
  else 
  {
    childNodeBox(2,0) = nodeBox_(2,0);
    childNodeBox(2,1) = zPlaneCoordinate_;
  }    
  //  printf("created chldAABB(%f\t%f\t%f\t%f\t%f\t%f)\n",
  //  chldAABB(0,0),chldAABB(0,1),chldAABB(1,0),chldAABB(1,1),chldAABB(2,0),chldAABB(2,1));
  return childNodeBox;
  
}



/*----------------------------------------------------------------------*
 | insert element in tree node                              u.may  07/08|
 *----------------------------------------------------------------------*/
void GEO::OctTree::TreeNode::insertElement(
    const int label,
    const int eleId) 
{
    elementList_[label].insert(eleId);
}



/*----------------------------------------------------------------------*
 | create children                                         u.may   08/08|
 *----------------------------------------------------------------------*/
void GEO::OctTree::TreeNode::createChildren(
    const DRT::Discretization&      dis,
    const std::map<int,BlitzVec3>& 	currentpositions)
{
  for(int index = 0; index < 8; index++)
  {
    
    children_[index] = rcp(new TreeNode(rcp(this), (treedepth_-1), getChildNodeBox(index)));

    // insert elements into child node
    for (std::map<int, std::set<int> >::const_iterator labelIter = elementList_.begin(); labelIter != elementList_.end(); labelIter++)
      for (std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
        std::vector<int> elementClassification = classifyElement(dis.gElement(*eleIter),currentpositions);
        for(unsigned int count = 0; count < elementClassification.size(); count++)
          if(elementClassification[count] == index)
            children_[index]->insertElement(labelIter->first,*eleIter);
      }
  
    // if one of the created children is empty, set label
    if ((children_[index]->getElementList()).empty())
    {
      const BlitzVec3 childNodeCenter(children_[index]->getCenterCoord());
      children_[index]->setLabel(getXFEMLabel(dis, currentpositions, childNodeCenter, elementList_));
    }
  }
  // this node becomes an inner tree node
  treeNodeType_ = INNER_NODE;
}


// TODO child IDS fixed here
/*----------------------------------------------------------------------*
 | classifiy point in node                                  peder   07/08|
 *----------------------------------------------------------------------*/
int GEO::OctTree::TreeNode::classifyPoint(
    const BlitzVec3&   point) const
{
  int childIndex = 0;
  if (point(0) > xPlaneCoordinate_)
    childIndex += 4;
  if (point(1) > yPlaneCoordinate_)
    childIndex += 2;
  if (point(2) > zPlaneCoordinate_)
    childIndex += 1;
  return childIndex;
}


// TODO child IDS fixed here
/*----------------------------------------------------------------------*
 | classifiy AABB in node                                  peder   07/08|
 *----------------------------------------------------------------------*/
std::vector<int> GEO::OctTree::TreeNode::classifyXAABB(
    const BlitzMat3x2&   AABB
    ) const 
{
  std::vector<int> octants;

  if (AABB(0, 1) > xPlaneCoordinate_) 
  {
    if (AABB(1, 1) > yPlaneCoordinate_) 
    {
      if (AABB(2, 1) > zPlaneCoordinate_)
        octants.push_back(7);

      if (AABB(2, 0) <= zPlaneCoordinate_)
        octants.push_back(6);

    }
    if (AABB(1, 0) <= yPlaneCoordinate_) 
    {
      if (AABB(2, 1) > zPlaneCoordinate_)
        octants.push_back(5);

      if (AABB(2, 0) <= zPlaneCoordinate_)
        octants.push_back(4);

    }
  }

  if (AABB(0, 0) <= xPlaneCoordinate_) 
  {
    if (AABB(1, 1) > yPlaneCoordinate_) 
    {
      if (AABB(2, 1) > zPlaneCoordinate_)
        octants.push_back(3);

      if (AABB(2, 0) <= zPlaneCoordinate_)
        octants.push_back(2);

    }
    if (AABB(1, 0) <= yPlaneCoordinate_) 
    {
      if (AABB(2, 1) > zPlaneCoordinate_)
        octants.push_back(1);

      if (AABB(2, 0) <= zPlaneCoordinate_)
        octants.push_back(0);

    }
  }
  return octants;
}



/*----------------------------------------------------------------------*
 | classifiy element in node                               peder   07/08|
 *----------------------------------------------------------------------*/
std::vector<int> GEO::OctTree::TreeNode::classifyElement(
    const DRT::Element*       		  element,
    const std::map<int,BlitzVec3>& 	currentpositions
    ) const
{
  const BlitzMat xyze(DRT::UTILS::getCurrentNodalPositions(element,currentpositions));
  XFEM::EleGeoType eleGeoType(XFEM::HIGHERORDER);
  GEO::checkRoughGeoType(element, xyze, eleGeoType);
  const BlitzMat3x2 elemXAABB(XFEM::computeFastXAABB(element, xyze, eleGeoType));  
  return classifyXAABB(elemXAABB);
}



/*----------------------------------------------------------------------*
 | update tree node                                        u.may   08/08|
 *----------------------------------------------------------------------*/
void GEO::OctTree::TreeNode::updateTreeNode(
   const BlitzMat3x2& AABB_old, 
   const BlitzMat3x2& AABB_new)
{
  for(int i = 0; i < 8; i++)
  {
    if(children_[i] != Teuchos::null)
    {
      if( GEO::inSameNodeBox(AABB_old, AABB_new, children_[i]->getNodeBox()) )
        children_[i]->updateTreeNode(AABB_old, AABB_new);
      else
        children_[i]->clear();
    }
  }
} 



/*----------------------------------------------------------------------*
 | return xfem label for point (interface method)          peder   07/08|
 *----------------------------------------------------------------------*/
int GEO::OctTree::TreeNode::queryPointType(
    const DRT::Discretization& 	     dis,
    const std::map<int,BlitzVec3>& 	 currentpositions, 
    const BlitzVec3& 		             point) 
{
  switch (treeNodeType_) 
  {
    case INNER_NODE:
    {       
      return children_[classifyPoint(point)]->queryPointType(dis, currentpositions, point);
      break;
    }
    case LEAF_NODE:   
    {
      if (elementList_.empty())
        return label_;

      // max depth reached, counts reverse
      if (treedepth_ <= 0 || elementList_.size()==1)
        return GEO::getXFEMLabel(dis, currentpositions, point, elementList_); 

      // dynamically grow tree otherwise, create children and set label for empty children
      createChildren(dis, currentpositions);
      // search in apropriate child node
      return children_[classifyPoint(point)]->queryPointType(dis, currentpositions, point);
      break;
    }
    default:
      dserror("should not get here\n");
  }
  return -1;
}


#endif  // #ifdef CCADISCRET

