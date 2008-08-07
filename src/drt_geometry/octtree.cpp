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
#include "intersection_service.H"
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
    const int           max_depth): 
      max_depth_(max_depth),
      treeRoot_(Teuchos::null)      
{}
    


/*----------------------------------------------------------------------*
| destructor OctTree                                        u.may 07/08|
*----------------------------------------------------------------------*/
GEO::OctTree::~OctTree(){}



/*----------------------------------------------------------------------*
 | initialize or rebuild tree with possibly new              u.may 07/08| 
 | discretization, elements are sorted according to the given map       |
 *----------------------------------------------------------------------*/
void GEO::OctTree::initializeTree(
    const BlitzMat3x2&                    nodeBox,
    const std::map<int,std::set<int> >&   elementsByLabel) 
{
  
  
  treeRoot_ = Teuchos::null;
 
  // TODO initialize if map is empty
  treeRoot_ = rcp(new TreeNode(NULL, max_depth_, nodeBox)); 

  // insert element map into tree root node
  if(elementsByLabel.size()>0)
    treeRoot_->setElementList(elementsByLabel);
  
}



/*----------------------------------------------------------------------*
 | initialize or rebuild tree with possibly new              u.may 07/08| 
 | discretization, elements are taken unsorted from discretization      |
 *----------------------------------------------------------------------*/
void GEO::OctTree::initializeTree(
    const BlitzMat3x2&              nodeBox,
    const DRT::Discretization&      dis) 
{

  treeRoot_ = rcp( new TreeNode(NULL, max_depth_, nodeBox)); 

  // inserts all elements in a map with key -1
  for (int i=0; i<dis.NumMyColElements(); ++i) 
    treeRoot_->insertElement(-1, dis.lRowElement(i)->Id());
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
  if(treeRoot_ == Teuchos::null)
    dserror("tree is not yet initialized !!!");
  
  std::vector< BlitzMat3x2 > AABBs_old = GEO::computeXAABBForLabeledStructures(dis, currentpositions_old, treeRoot_->getElementList());
  std::vector< BlitzMat3x2 > AABBs_new = GEO::computeXAABBForLabeledStructures(dis, currentpositions_new, treeRoot_->getElementList());

  for(unsigned int i = 0; i < AABBs_new.size(); i++)
    treeRoot_->updateTreeNode(AABBs_old[i], AABBs_new[i]);
}



/*----------------------------------------------------------------------*
 | returns xfem label of point                               u.may 07/08|
 *----------------------------------------------------------------------*/
int GEO::OctTree::queryXFEMFSIPointType(
    const DRT::Discretization& 	     dis,
    const std::map<int,BlitzVec3>& 	 currentpositions, 
    const BlitzVec3& 		             point) 
{
  //  cout << " ASKING THE TREE" << endl;
  TEUCHOS_FUNC_TIME_MONITOR("OctTree - queryTime");
  
  if(treeRoot_ == Teuchos::null)
      dserror("tree is not yet initialized !!!");

  if(!treeRoot_->getElementList().empty())
    return treeRoot_->queryXFEMFSIPointType(dis, currentpositions, point);
  else 
    return 0;
}



/*----------------------------------------------------------------------*
 | returns nodes in the radius of a given point              u.may 07/08|
 *----------------------------------------------------------------------*/
std::set<int> GEO::OctTree::searchNodesInRadius(
    const DRT::Discretization& 	     dis,
    const std::map<int,BlitzVec3>&   currentpositions, 
    const BlitzVec3& 		             point,
    const double                     radius, 
    const int                        label) 
{
  //  cout << " ASKING THE TREE" << endl;
  TEUCHOS_FUNC_TIME_MONITOR("OctTree - queryTime");

  std::set<int> nodeset;

  if(treeRoot_ == Teuchos::null)
      dserror("tree is not yet initialized !!!");

  if(!treeRoot_->getElementList().empty())
    nodeset = treeRoot_->searchNodesInRadius(dis, currentpositions, point, radius, label);
  else
    dserror("element list is empty");

  return nodeset;
}


/*----------------------------------------------------------------------*
 | print tree to gmsh file                                 peder   07/08|
 *----------------------------------------------------------------------*/
void GEO::OctTree::printTree(
    const string  prefix, 
    const int     step) const
{
  cout << endl << "writing... ";
  if (treeRoot_ == Teuchos::null) {
    cout << "nothing to write, tree not initialized yet -> done" << endl;
    return;
  }
  if (treeRoot_->getElementList().empty()){ 
    cout << "nothing to write, tree empty -> done" << endl;
    return;
  }
  
  std::stringstream filename;
  std::stringstream node_string;
  filename << prefix << "_octtree" << std::setw(5) << setfill('0') << step << ".pos";
  cout << " " << filename.str() << " ...";
  node_string << "View \" " << "XFEM_FSI_Octree \" {" << endl;  
  treeRoot_->printTreeNode(max_depth_,node_string);
  node_string << "};" << endl;
  std::ofstream f_system(filename.str().c_str());
  f_system << node_string.str();
  f_system.close();
  cout << " done" << endl;
}



/*----------------------------------------------------------------------*
 | c-tor TreeNode                                            u.may 07/08|
 *----------------------------------------------------------------------*/
GEO::OctTree::TreeNode::TreeNode(
    const TreeNode* const                       parent,
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
  children_.assign(8, Teuchos::null);
  
}
                  

/*====================================================================*/
  /* octtree */
  /*--------------------------------------------------------------------*/
  /* numbering of children
   * child 0: 
   * child 1: 
   * child 2: 
   * child 3: 
   * child 4:
   * child 5: 
   * child 6: 
   * child 7: 
   * 
   *                      z
   *                      |           
   *             1========|================3
   *           //|        |               /||
   *          // |        |              //||
   *         //  |        |             // ||
   *        //   |        |            //  ||
   *       //    |        |           //   ||
   *      //     |        |          //    ||
   *     //      |        |         //     ||
   *     5=========================7       ||
   *    ||       |        |        ||      ||
   *    ||       |        o--------||---------y
   *    ||       |       /         ||      ||
   *    ||       0------/----------||------2
   *    ||      /      /           ||     //
   *    ||     /      /            ||    //
   *    ||    /      /             ||   //
   *    ||   /      /              ||  //
   *    ||  /      /               || //
   *    || /      x                ||//
   *    ||/                        ||/
   *     4=========================6
   *
   */
  /*====================================================================*/


/*  int index = 0;
    if (point(0) > xPlaneCoordinate_)
      index += 4;
    if (point(1) > yPlaneCoordinate_)
      octIdx += 2;
    if (point(2) > zPlaneCoordinate_)
      index += 1;
*/

    
    

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
  if (parent_!=NULL)
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
const GEO::OctTree::TreeNode* const GEO::OctTree::TreeNode::getParent() const
{
  if (this->hasParent())
    return parent_;
  
  return NULL;
}



/*----------------------------------------------------------------------*
 | get center of treenode                                  peder   07/08|
 *----------------------------------------------------------------------*/
const BlitzVec3 GEO::OctTree::TreeNode::getCenterCoord() const
{
  return BlitzVec3(this->xPlaneCoordinate_, this->yPlaneCoordinate_, this->zPlaneCoordinate_);
}



/*----------------------------------------------------------------------*
 | get map of elements                                     peder   07/08|
 *----------------------------------------------------------------------*/
const std::map<int,std::set<int> > GEO::OctTree::TreeNode::getElementList() const
{
  return elementList_;
}



/*----------------------------------------------------------------------*
 | get type of tree node                                   peder   07/08|
 *----------------------------------------------------------------------*/
const GEO::TreeNodeType GEO::OctTree::TreeNode::getTreeNodeType() const
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
const Teuchos::RCP<GEO::OctTree::TreeNode> GEO::OctTree::TreeNode::getChild(
    const int index) const
{
  return children_[index];
}



/*----------------------------------------------------------------------*
 | get node box of a child specified by index              peder   07/08|
 *----------------------------------------------------------------------*/
BlitzMat3x2 GEO::OctTree::TreeNode::getChildNodeBox(
	const int index) const
{
  BlitzMat3x2 childNodeBox;
  
  // determine z-coordinates
  if (index>3)
  {
    childNodeBox(0,0) = xPlaneCoordinate_;
    childNodeBox(0,1) = nodeBox_(0,1);
  }
  else 
  {
    childNodeBox(0,0) = nodeBox_(0,0);
    childNodeBox(0,1) = xPlaneCoordinate_;
  }

  // determine y-coordinates
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

  // determine z-coordinates
  if ((index==1) || (index==3) || (index==5) || (index==7))
  {
    childNodeBox(2,0) = zPlaneCoordinate_;
    childNodeBox(2,1) = nodeBox_(2,1);
  }
  else 
  {
    childNodeBox(2,0) = nodeBox_(2,0);
    childNodeBox(2,1) = zPlaneCoordinate_;
  }    
  //  printf("created chldAABB(%f\t%f\t%f\t%f\t%f\t%f)\n", childNodeBox(0,0),childNodeBox(0,1),childNodeBox(1,0),childNodeBox(1,1),childNodeBox(2,0),childNodeBox(2,1));
  return childNodeBox;
  
}



/*----------------------------------------------------------------------*
 | insert element in tree node                              u.may  07/08|
 *----------------------------------------------------------------------*/
void GEO::OctTree::TreeNode::insertElement(
    const int   label,
    const int   eleId) 
{
    elementList_[label].insert(eleId);
}



/*----------------------------------------------------------------------*
 | create children                                         u.may   08/08|
 *----------------------------------------------------------------------*/
void GEO::OctTree::TreeNode::createChildren(
    const DRT::Discretization&      	dis,
    const std::map<int,BlitzVec3>& 	  currentpositions)
{
  // create empty children
  for(int index = 0; index < 8; index++)
    children_[index] = rcp(new TreeNode(this, (treedepth_-1), getChildNodeBox(index) ));
  
  
  // insert elements into child node
  for (std::map<int, std::set<int> >::const_iterator labelIter = elementList_.begin(); labelIter != elementList_.end(); labelIter++)
    for (std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
    {
      std::vector<int> elementClassification = classifyElement(dis.gElement(*eleIter),currentpositions);
      for(unsigned int count = 0; count < elementClassification.size(); count++)
        children_[elementClassification[count]]->insertElement(labelIter->first,*eleIter);
    }
  
  for(int index = 0; index < 8; index++)
  {
    // if one of the created children is empty, set label immediately
    if ((children_[index]->getElementList()).empty())
    {
      const BlitzVec3 childNodeCenter(children_[index]->getCenterCoord());
      // xfem label has to be computed on this level because child is empty
      children_[index]->setLabel(getXFEMLabel(dis, currentpositions, childNodeCenter, elementList_));
    }
  }
  // this node becomes an inner tree node
  treeNodeType_ = INNER_NODE;
}



/*----------------------------------------------------------------------*
 | classifiy point in node                                  peder   07/08|
 *----------------------------------------------------------------------*/
const int GEO::OctTree::TreeNode::classifyPoint(
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



/*----------------------------------------------------------------------*
 | classifiy AABB in node                                  u.may   07/08|
 *----------------------------------------------------------------------*/
std::vector<int> GEO::OctTree::TreeNode::classifyXAABB(
    const BlitzMat3x2&   AABB
    ) const 
{
  
  // collect all children the XAABB is lying in
  // use tolerances such that it is ensured that no child is left behind :-) 
  // XAABB s which are lying between plane-tol and plane + tol are collected in each of the two chlidren
  
  std::vector<int> octants;

  // check max_x greater than x-plane
  if (AABB(0, 1) > (xPlaneCoordinate_ - GEO::TOL7) ) 
  {
    // check max_y greater than y-plane
    if (AABB(1, 1) > (yPlaneCoordinate_ - GEO::TOL7) ) 
    {
      // check max_z greater than z-plane
      if (AABB(2, 1) > (zPlaneCoordinate_ - GEO::TOL7) )
        octants.push_back(7);

      // check min_z less than z-plane
      if (AABB(2, 0) < (zPlaneCoordinate_ + GEO::TOL7) )
        octants.push_back(6);
    }
    
    // check min_y less than y-plane
    if (AABB(1, 0) < ( yPlaneCoordinate_ + GEO::TOL7) ) 
    {
      // check min_z less than z-plane
      if (AABB(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
        octants.push_back(4);
      
      // check max_z greater than z-plane
      if (AABB(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
        octants.push_back(5);

    }
  }

  // check min_x less than x-plane
  if (AABB(0, 0) < ( xPlaneCoordinate_ + GEO::TOL7) ) 
  {
    // check min_y less than y-plane
    if (AABB(1, 0) < ( yPlaneCoordinate_ + GEO::TOL7) ) 
    {
      // check min_z less than z-plane
      if (AABB(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
        octants.push_back(0);
      
      // check max_z greater than z-plane
      if (AABB(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
        octants.push_back(1);
    }
    
    // check max_y greater than y-plane
    if (AABB(1, 1) > ( yPlaneCoordinate_ - GEO::TOL7) ) 
    {
      // check max_z greater than z-plane
      if (AABB(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
        octants.push_back(3);
      
      // check min_z less than z-plane
      if (AABB(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
        octants.push_back(2);

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
  GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
  GEO::checkRoughGeoType(element, xyze, eleGeoType);
  const BlitzMat3x2 elemXAABB(GEO::computeFastXAABB(element, xyze, eleGeoType));  
  return classifyXAABB(elemXAABB);
}




/*----------------------------------------------------------------------*
 | classifiy element in node                               u.may   08/08|
 *----------------------------------------------------------------------*/
int GEO::OctTree::TreeNode::classifyRadius(
    const double      radius,
    const BlitzVec3&	point
    ) const
{
  BlitzMat3x2 radiusXAABB;  

  for(int dim = 0; dim < 3; dim++)
  {
    radiusXAABB(dim,0) = ( point(dim) - radius) - GEO::TOL7;
    radiusXAABB(dim,1) = ( point(dim) + radius) + GEO::TOL7;
  }
  return classifyAABBCompletelyInNodeBox(radiusXAABB);
}



/*----------------------------------------------------------------------*
 | checks if a AABB is completely in a child node          u.may   08/08|
 *----------------------------------------------------------------------*/
int GEO::OctTree::TreeNode::classifyAABBCompletelyInNodeBox(
    const BlitzMat3x2&  AABB
    ) const
{
  int childindex = -1;
  
  int indexMin = 0;
  if(AABB(0,0) > xPlaneCoordinate_)
    indexMin += 4;
  if(AABB(1,0) > yPlaneCoordinate_)
    indexMin += 2;
  if(AABB(2,0) > zPlaneCoordinate_)
    indexMin += 1;
  
  int indexMax = 0;
  if(AABB(0,1) > xPlaneCoordinate_)
    indexMax += 4;
  if(AABB(1,1) > yPlaneCoordinate_)
    indexMax += 2;
  if(AABB(2,1) > zPlaneCoordinate_)
    indexMax += 1;
    
  
  if(indexMin == indexMax)
    childindex = indexMin;
  
  return childindex;
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
 | return xfem label for point (interface method)          u.may   07/08|
 *----------------------------------------------------------------------*/
int GEO::OctTree::TreeNode::queryXFEMFSIPointType(
    const DRT::Discretization& 	     dis,
    const std::map<int,BlitzVec3>& 	 currentpositions, 
    const BlitzVec3& 		             point) 
{
  switch (treeNodeType_) 
  {
    case INNER_NODE:
    {       
      return children_[classifyPoint(point)]->queryXFEMFSIPointType(dis, currentpositions, point);
      break;
    }
    case LEAF_NODE:   
    {
      if (elementList_.empty())
        return label_;

      // max depth reached, counts reverse
      if (treedepth_ <= 0 || (elementList_.size()==1 && (elementList_.begin()->second).size() == 1) )
        return GEO::getXFEMLabel(dis, currentpositions, point, elementList_); 

      // dynamically grow tree otherwise, create children and set label for empty children
      createChildren(dis, currentpositions);
      // search in apropriate child node
      return children_[classifyPoint(point)]->queryXFEMFSIPointType(dis, currentpositions, point);
      break;
    }
    default:
      dserror("should not get here\n");
  }
  return -1;
}



/*----------------------------------------------------------------------*
 | returns nodes in the radius of a given point              u.may 08/08|
 *----------------------------------------------------------------------*/
std::set<int> GEO::OctTree::TreeNode::searchNodesInRadius(
    const DRT::Discretization& 	     dis,
    const std::map<int,BlitzVec3>&   currentpositions, 
    const BlitzVec3& 		             point,
    const double                     radius, 
    const int 			                 label) 
{
  
  std::set<int> nodeset;

  switch (treeNodeType_) 
  {
    case INNER_NODE:
    {       
      const int childindex = classifyRadius(radius, point);
      if(childindex != -1) // child node found which encloses AABB so step down
        return children_[childindex]->searchNodesInRadius(dis, currentpositions, point, radius, label);
      else
        return GEO::getNodeSetInRadius(dis, currentpositions, point, radius, label, elementList_); 
        
      break;
    }
    case LEAF_NODE:   
    {
      if(elementList_.empty())
        return nodeset;

      // max depth reached, counts reverse
      if (treedepth_ <= 0 || (elementList_.size()==1 && (elementList_.begin()->second).size() == 1) )
        return GEO::getNodeSetInRadius(dis, currentpositions, point, radius, label, elementList_); 

      // dynamically grow tree otherwise, create children and set label for empty children
      // search in apropriate child node     
      const int childindex = classifyRadius(radius, point);
      if(childindex != -1)  // child node found
      {
        createChildren(dis, currentpositions);
        return children_[childindex]->searchNodesInRadius(dis, currentpositions, point, radius, label);
      }
      else // AABB does not fit into a single child node box anymore so don t refine any further
        return GEO::getNodeSetInRadius(dis, currentpositions, point, radius, label, elementList_); 
      break;
    }
    default:
      dserror("should not get here\n");
  }
  return nodeset;

}



/*----------------------------------------------------------------------*
 | print tree node to gmsh file                            peder   07/08|
 *----------------------------------------------------------------------*/
void GEO::OctTree::TreeNode::printTreeNode(
    const int       max_depth,
    stringstream&   fc) const
{
  if(treedepth_== max_depth)
  {
    BlitzMat printBox(3,8);
    printBox(0,0) = nodeBox_(0,0); printBox(1,0) = nodeBox_(1,0); printBox(2,0) = nodeBox_(2,0);
    printBox(0,1) = nodeBox_(0,0); printBox(1,1) = nodeBox_(1,1); printBox(2,1) = nodeBox_(2,0);
    printBox(0,2) = nodeBox_(0,0); printBox(1,2) = nodeBox_(1,1); printBox(2,2) = nodeBox_(2,1);
    printBox(0,3) = nodeBox_(0,0); printBox(1,3) = nodeBox_(1,0); printBox(2,3) = nodeBox_(2,1);
    printBox(0,4) = nodeBox_(0,1); printBox(1,4) = nodeBox_(1,0); printBox(2,4) = nodeBox_(2,0);
    printBox(0,5) = nodeBox_(0,1); printBox(1,5) = nodeBox_(1,1); printBox(2,5) = nodeBox_(2,0);
    printBox(0,6) = nodeBox_(0,1); printBox(1,6) = nodeBox_(1,1); printBox(2,6) = nodeBox_(2,1);
    printBox(0,7) = nodeBox_(0,1); printBox(1,7) = nodeBox_(1,0); printBox(2,7) = nodeBox_(2,1);
    fc << IO::GMSH::cellWithScalarToString(DRT::Element::hex8, 0, printBox)<< endl;
  }

  if(treeNodeType_== GEO::INNER_NODE)
  {
    for (int i=0; i<8; i++)
      if (children_[i] != Teuchos::null)
        children_[i]->printTreeNode(max_depth, fc);
   
  }
  else if (treeNodeType_== GEO::LEAF_NODE)
  {
    int factor = -1;
    
    if (label_ < 0)
      factor = 0; // more than one candidate in this octant
    else if (label_ == 0)
      factor = 1; // fluid octant
    else
      factor = 2; // solid octant
    BlitzMat printBox(3,8);
    printBox(0,0) = nodeBox_(0,0); printBox(1,0) = nodeBox_(1,0); printBox(2,0) = nodeBox_(2,0);
    printBox(0,1) = nodeBox_(0,0); printBox(1,1) = nodeBox_(1,1); printBox(2,1) = nodeBox_(2,0);
    printBox(0,2) = nodeBox_(0,0); printBox(1,2) = nodeBox_(1,1); printBox(2,2) = nodeBox_(2,1);
    printBox(0,3) = nodeBox_(0,0); printBox(1,3) = nodeBox_(1,0); printBox(2,3) = nodeBox_(2,1);
    printBox(0,4) = nodeBox_(0,1); printBox(1,4) = nodeBox_(1,0); printBox(2,4) = nodeBox_(2,0);
    printBox(0,5) = nodeBox_(0,1); printBox(1,5) = nodeBox_(1,1); printBox(2,5) = nodeBox_(2,0);
    printBox(0,6) = nodeBox_(0,1); printBox(1,6) = nodeBox_(1,1); printBox(2,6) = nodeBox_(2,1);
    printBox(0,7) = nodeBox_(0,1); printBox(1,7) = nodeBox_(1,0); printBox(2,7) = nodeBox_(2,1);
    fc << IO::GMSH::cellWithScalarToString(DRT::Element::hex8, factor + treedepth_ + max_depth, printBox)<< endl;
  }
}


#endif  // #ifdef CCADISCRET

