 /*!
\file searchtree.cpp

\brief provides a class with search tree

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
 */
#ifdef CCADISCRET
#include "searchtree.H"
#include "intersection_service.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_lib/standardtypes_cpp.H"
#include <Teuchos_TimeMonitor.hpp>



/*----------------------------------------------------------------------*
 | constructor SearchTree                                    u.may 07/08|
 *----------------------------------------------------------------------*/
GEO::SearchTree::SearchTree(
    const int           max_depth): 
    max_depth_(max_depth),
    treeRoot_(Teuchos::null)      
{}
    


/*----------------------------------------------------------------------*
| destructor SearchTree                                     u.may 07/08|
*----------------------------------------------------------------------*/
GEO::SearchTree::~SearchTree(){}



/*----------------------------------------------------------------------*
 | initialize or rebuild tree with possibly new              u.may 07/08| 
 | discretization, elements are sorted according to the given map       |
 *----------------------------------------------------------------------*/
void GEO::SearchTree::initializeTree(
    const LINALG::Matrix<3,2>&            nodeBox,
    const std::map<int,std::set<int> >&   elementsByLabel,
    const TreeType                        treetype) 
{
  
  
  treeRoot_ = Teuchos::null;
 
  // TODO initialize if map is empty
  treeRoot_ = rcp(new TreeNode(NULL, max_depth_, nodeBox, treetype)); 

  // insert element map into tree root node
  if(elementsByLabel.size()>0)
    treeRoot_->setElementList(elementsByLabel);
  
}



/*----------------------------------------------------------------------*
 | initialize or rebuild tree with possibly new              u.may 07/08| 
 | discretization, elements are taken unsorted from discretization      |
 *----------------------------------------------------------------------*/
void GEO::SearchTree::initializeTree(
    const LINALG::Matrix<3,2>&  nodeBox,
    const DRT::Discretization&  dis,
    const TreeType              treetype) 
{

  treeRoot_ = Teuchos::null;
  treeRoot_ = rcp( new TreeNode(NULL, max_depth_, nodeBox, treetype)); 
  
  // inserts all elements in a map with key -1 and global id
  for (int i=0; i<dis.NumMyColElements(); i++)
    treeRoot_->insertElement(-1, dis.lColElement(i)->Id());
}



/*----------------------------------------------------------------------*
 | update tree                                               u.may 08/08|
 | only usefull for labeled structures,for networks use initializeTree  |
 *----------------------------------------------------------------------*/
void GEO::SearchTree::updateTree(
    const DRT::Discretization& 	dis,
    const std::map<int,LINALG::Matrix<3,1> >& 	currentpositions_old, 
    const std::map<int,LINALG::Matrix<3,1> >& 	currentpositions_new) 
{
  if(treeRoot_ == Teuchos::null)
    dserror("tree is not yet initialized !!!");
  
  std::vector< LINALG::Matrix<3,2> > AABBs_old = GEO::computeXAABBForLabeledStructures(dis, currentpositions_old, treeRoot_->getElementList());
  std::vector< LINALG::Matrix<3,2> > AABBs_new = GEO::computeXAABBForLabeledStructures(dis, currentpositions_new, treeRoot_->getElementList());

  for(unsigned int i = 0; i < AABBs_new.size(); i++)
    treeRoot_->updateTreeNode(AABBs_old[i], AABBs_new[i]);
}



/*----------------------------------------------------------------------*
 | returns xfem label and nearest object to point            u.may 07/08|
 *----------------------------------------------------------------------*/
int GEO::SearchTree::queryFSINearestObject(
    const DRT::Discretization&                  dis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const LINALG::Matrix<3,1>&                  point,
    GEO::NearestObject&                         nearestobject)
{
  
  TEUCHOS_FUNC_TIME_MONITOR("GEO::SearchTree - queryTime");
  
  if(treeRoot_ == Teuchos::null)
      dserror("tree is not yet initialized !!!");

  if(!treeRoot_->getElementList().empty())
    return treeRoot_->queryFSINearestObject(dis, currentpositions, point, nearestobject);
  else 
    return 0;
}



/*----------------------------------------------------------------------*
 | returns xfem label of point                               u.may 07/08|
 *----------------------------------------------------------------------*/
int GEO::SearchTree::queryXFEMFSIPointType(
    const DRT::Discretization&                  dis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions,
    const LINALG::Matrix<3,1>&                  point)
{
  
  TEUCHOS_FUNC_TIME_MONITOR("GEO::SearchTree - queryTime");
  
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
std::map<int,std::set<int> > GEO::SearchTree::searchElementsInRadius(
    const DRT::Discretization& 	                dis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions, 
    const LINALG::Matrix<3,1>&                  point,
    const double                                radius, 
    const int                                   label) 
{
  TEUCHOS_FUNC_TIME_MONITOR("SearchTree - queryTime");

  std::map<int,std::set<int> > nodeset;

  if(treeRoot_ == Teuchos::null)
      dserror("tree is not yet initialized !!!");

  if(!(treeRoot_->getElementList().empty()))
    nodeset = treeRoot_->searchElementsInRadius(dis, currentpositions, point, radius, label);
  else
    dserror("element list is empty");

  return nodeset;
}



/*----------------------------------------------------------------------*
 | returns a vector of elements whose XAABB s are            u.may 09/08|
 | intersecting with the XAABB of a given volume element                |
 *----------------------------------------------------------------------*/
std::vector<int> GEO::SearchTree::queryIntersectionCandidates(
    const DRT::Discretization&                  dis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions, 
    DRT::Element*                               element,
    const LINALG::SerialDenseMatrix&            xyze_element) 
{
  TEUCHOS_FUNC_TIME_MONITOR("SearchTree - queryTime");

  std::vector<int> elementset;

  if(treeRoot_ == Teuchos::null)
      dserror("tree is not yet initialized !!!");

  if(!(treeRoot_->getElementList().empty()))
    elementset = treeRoot_->queryIntersectionCandidates(dis, currentpositions, element, xyze_element);

  return elementset;
}



/*----------------------------------------------------------------------*
 | returns intersection elements (CONTACT)                    popp 07/08|
 *----------------------------------------------------------------------*/
std::vector<int> GEO::SearchTree::searchIntersectionElements(
    const DRT::Discretization&                  dis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions, 
    DRT::Element*                               element) 
{
  TEUCHOS_FUNC_TIME_MONITOR("SearchTree - queryTime");

  std::vector<int> elementset;

  if(treeRoot_ == Teuchos::null)
      dserror("tree is not yet initialized !!!");

  if(!treeRoot_->getElementList().empty())
    elementset = treeRoot_->searchIntersectionElements(dis, currentpositions, element);
  else
    dserror("element list is empty");

  return elementset;
}



/*----------------------------------------------------------------------*
 | print tree to gmsh file                                 peder   07/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::printTree(
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
 | evaluate tree metrics (METRICS)                         peder   07/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::evaluateTreeMetrics(
    const int     step) const
{
  cout << "\t********************* TREE METRICS ******************" << endl;

  if (treeRoot_ == Teuchos::null)
  { 
    cout << "tree root was not initialized yet, nothing to print" << endl;
    return;
  }
  
  cout.precision(5);
  cout << "\tnumber tree nodes        : " << treeRoot_->getNumNodesInTree() << endl;
  cout << "\ttree depth               : " << max_depth_-treeRoot_->getDepth() << " (max: "<< max_depth_<<")"<< endl;
  cout << "\t***************************************************" << endl;
  return;
}



/*----------------------------------------------------------------------*
 | c-tor TreeNode                                            u.may 07/08|
 *----------------------------------------------------------------------*/
GEO::SearchTree::TreeNode::TreeNode(
    const TreeNode* const                       parent,
    const int                                   depth, 
    const LINALG::Matrix<3,2> &                 nodeBox,
    const TreeType                              treeType):
    parent_(parent),
    treedepth_(depth),
    treeNodeType_(LEAF_NODE),
    treeType_(treeType),
    label_(-1),
    nodeBox_(nodeBox),
    xPlaneCoordinate_( (nodeBox_(0,0)+0.5*(nodeBox_(0,1)-nodeBox_(0,0))) ),
    yPlaneCoordinate_( (nodeBox_(1,0)+0.5*(nodeBox_(1,1)-nodeBox_(1,0))) ),
    zPlaneCoordinate_( (nodeBox_(2,0)+0.5*(nodeBox_(2,1)-nodeBox_(2,0))) )
{    
  children_.assign(getNumChildren(), Teuchos::null);
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
   *             4========|================6¨
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
   *     1=========================3
   *
   */
  /*====================================================================*/


/*  int index = 0;
    if (point(0) > xPlaneCoordinate_)
      index += 1;
    if (point(1) > yPlaneCoordinate_)
      octIdx += 2;
    if (point(2) > zPlaneCoordinate_)
      index += 4;
*/

    
    

/*----------------------------------------------------------------------*
 | d-tor TreeNode                                            u.may 07/08|
 *----------------------------------------------------------------------*/
GEO::SearchTree::TreeNode::~TreeNode() {}



/*----------------------------------------------------------------------*
 | clears node and deletes chlidren                          u.may 08/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::clear()
{
  treeNodeType_ = LEAF_NODE;
  label_ = -1; 

  children_.clear();
  children_.assign(getNumChildren(),Teuchos::null);
}



/*----------------------------------------------------------------------*
 | has parent element                                        peder 07/08|
 *----------------------------------------------------------------------*/
bool GEO::SearchTree::TreeNode::hasParent() const 
{
  if (parent_!=NULL)
    return true;
  return false;
}



/*----------------------------------------------------------------------*
 | sets element list                                         u.may 07/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::setElementList(
   const std::map<int, std::set<int> >& elementsByLabel)
{
   elementList_ = elementsByLabel;
}



/*----------------------------------------------------------------------*
 | set label                                                 peder 07/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::setLabel(
	const int label)
{
  label_ = label;
}



/*----------------------------------------------------------------------*
 | set nearest object                                        u.may 09/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::setNearestObject(
      const GEO::NearestObject&   nearestObject)
{
  nearestObject_ = nearestObject;
}



/*----------------------------------------------------------------------*
 | returns pointer to parent element                         peder 07/08|
 *----------------------------------------------------------------------*/
const GEO::SearchTree::TreeNode* const GEO::SearchTree::TreeNode::getParent() const
{
  if (this->hasParent())
    return parent_;
  
  return NULL;
}



/*----------------------------------------------------------------------*
 | get center of treenode                                  peder   07/08|
 *----------------------------------------------------------------------*/
const LINALG::Matrix<3,1> GEO::SearchTree::TreeNode::getCenterCoord() const
{
  LINALG::Matrix<3,1> centerCoord;

  centerCoord(0) = this->xPlaneCoordinate_;
  centerCoord(1) = this->yPlaneCoordinate_;
  centerCoord(2) = this->zPlaneCoordinate_;
  return centerCoord;
}



/*----------------------------------------------------------------------*
 | get map of elements                                     peder   07/08|
 *----------------------------------------------------------------------*/
const std::map<int,std::set<int> > GEO::SearchTree::TreeNode::getElementList() const
{
  return elementList_;
}



/*----------------------------------------------------------------------*
 | get type of tree node                                   peder   07/08|
 *----------------------------------------------------------------------*/
const GEO::TreeNodeType GEO::SearchTree::TreeNode::getTreeNodeType() const
{
  return treeNodeType_;
}



/*----------------------------------------------------------------------*
 | get type of tree                                        u.may   08/08|
 *----------------------------------------------------------------------*/
const GEO::TreeType GEO::SearchTree::TreeNode::getTreeType() const
{
  return treeType_;
}



/*----------------------------------------------------------------------*
 | get node box                                            peder   07/08|
 *----------------------------------------------------------------------*/
const LINALG::Matrix<3,2>& GEO::SearchTree::TreeNode::getNodeBox() const
{
  return nodeBox_;
}



/*----------------------------------------------------------------------*
 | get number of children                                  u.may   08/08|
 *----------------------------------------------------------------------*/
int GEO::SearchTree::TreeNode::getNumChildren(
    ) const
{
  if(treeType_ == OCTTREE)
    return 8;
  else if(treeType_ == QUADTREE)
    return 4;
  else
    dserror("treetype does not exist");
  
  return -1;
}



/*----------------------------------------------------------------------*
 | get child                                               peder   07/08|
 *----------------------------------------------------------------------*/
const Teuchos::RCP<GEO::SearchTree::TreeNode> GEO::SearchTree::TreeNode::getChild(
    const int index) const
{
  return children_[index];
}



/*----------------------------------------------------------------------*
 | get node box of a child specified by index              peder   07/08|
 *----------------------------------------------------------------------*/
LINALG::Matrix<3,2> GEO::SearchTree::TreeNode::getChildNodeBox(
	const int index) const
{
  LINALG::Matrix<3,2> childNodeBox;
  
  // determine x-coordinates
  if ((index==1) || (index==3) || (index==5) || (index==7))
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
  
  // determine z coordinates
  if (index>3)
  {
    childNodeBox(2,0) = zPlaneCoordinate_;
    childNodeBox(2,1) = nodeBox_(2,1);
  }
  else 
  {
    if(treeType_ == OCTTREE)
    {
      childNodeBox(2,0) = nodeBox_(2,0);
      childNodeBox(2,1) = zPlaneCoordinate_;
    }
    else
    {
      childNodeBox(2,0) = 0.0;
      childNodeBox(2,1) = 0.0;
    }
  }
  //  printf("created chldAABB(%f\t%f\t%f\t%f\t%f\t%f)\n", childNodeBox(0,0),childNodeBox(0,1),childNodeBox(1,0),childNodeBox(1,1),childNodeBox(2,0),childNodeBox(2,1));
  return childNodeBox;
}



/*----------------------------------------------------------------------*
 | insert element in tree node                              u.may  07/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::insertElement(
    const int   label,
    const int   eleId) 
{
    elementList_[label].insert(eleId);
}



/*----------------------------------------------------------------------*
 | create children                                         u.may   08/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::createChildren(
    const DRT::Discretization&      	             dis,
    const std::map<int,LINALG::Matrix<3,1> >& 	   currentpositions)
{
  // create empty children
  for(int index = 0; index < getNumChildren(); index++)
    children_[index] = rcp(new TreeNode(this, (treedepth_-1), getChildNodeBox(index), treeType_));
  
  
  // insert elements into child node
  for (std::map<int, std::set<int> >::const_iterator labelIter = elementList_.begin(); labelIter != elementList_.end(); labelIter++)
  {
    for (std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
    {
      std::vector<int> elementClassification = classifyElement(dis.gElement(*eleIter),currentpositions);
      for(unsigned int count = 0; count < elementClassification.size(); count++)
      {
        children_[elementClassification[count]]->insertElement(labelIter->first,*eleIter);
      }
    }
  }
  // this node becomes an inner tree node
  treeNodeType_ = INNER_NODE;
}
  
 
/*----------------------------------------------------------------------*
 | set XFEM label of empty children                        u.may   08/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::setXFEMLabelOfEmptyChildren(
    const DRT::Discretization&                    dis,
    const std::map<int,LINALG::Matrix<3,1> >&     currentpositions)
{    
  for(int index = 0; index < getNumChildren(); index++)
  {
    // if one of the created children is empty, set label immediately
    if ((children_[index]->getElementList()).empty())
    {
      const LINALG::Matrix<3,1>  childNodeCenter(children_[index]->getCenterCoord());
      // xfem label has to be computed on this level because child is empty
      children_[index]->setLabel(getXFEMLabel(dis, currentpositions, childNodeCenter, elementList_));
    }
  }
}



/*----------------------------------------------------------------------*
 | set XFEM label and nearest object of empty children     u.may   08/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::setXFEMLabelAndNearestObjectOfEmptyChildren(
    const DRT::Discretization&                    dis,
    const std::map<int,LINALG::Matrix<3,1> >&     currentpositions)
{    
  for(int index = 0; index < getNumChildren(); index++)
  {
    // if one of the created children is empty, set label immediately
    if ((children_[index]->getElementList()).empty())
    {
      const LINALG::Matrix<3,1> childNodeCenter(children_[index]->getCenterCoord());
      // xfem label has to be computed on this level because child is empty
      GEO::NearestObject nearestObject;
      int label = getXFEMLabelAndNearestObject(dis, currentpositions, childNodeCenter, elementList_, nearestObject);
      children_[index]->setLabel(label);
      children_[index]->setNearestObject(nearestObject);
    }
  }
}



/*----------------------------------------------------------------------*
 | classifiy point in node                                  peder   07/08|
 *----------------------------------------------------------------------*/
const int GEO::SearchTree::TreeNode::classifyPoint(
    const LINALG::Matrix<3,1>&   point) const
{
  
  int childIndex = 0;
  if (point(0) > xPlaneCoordinate_)
    childIndex += 1;
  if (point(1) > yPlaneCoordinate_)
    childIndex += 2;
  if(treeType_ == OCTTREE)
    if (point(2) > zPlaneCoordinate_)
      childIndex += 4;
  
  return childIndex;
}



/*----------------------------------------------------------------------*
 | classifiy AABB in node                                  u.may   07/08|
 *----------------------------------------------------------------------*/
std::vector<int> GEO::SearchTree::TreeNode::classifyXAABB(
    const LINALG::Matrix<3,2>&   AABB
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
      if(treeType_ == OCTTREE)
      {
        // check max_z greater than z-plane
        if (AABB(2, 1) > (zPlaneCoordinate_ - GEO::TOL7) )
          octants.push_back(7);

        // check min_z less than z-plane
        if (AABB(2, 0) < (zPlaneCoordinate_ + GEO::TOL7) )
          octants.push_back(3);
      }
      
      else if (treeType_ == QUADTREE)
        octants.push_back(3);
    }
    
    // check min_y less than y-plane
    if (AABB(1, 0) < ( yPlaneCoordinate_ + GEO::TOL7) ) 
    {
      if(treeType_ == OCTTREE)
      {
        // check max_z greater than z-plane
        if (AABB(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
          octants.push_back(5);
              
        // check min_z less than z-plane
        if (AABB(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
          octants.push_back(1);
      }
      
      else if(treeType_ == QUADTREE)
        octants.push_back(1);
    }
  }

  // check min_x less than x-plane
  if (AABB(0, 0) < ( xPlaneCoordinate_ + GEO::TOL7) ) 
  {
    // check min_y less than y-plane
    if (AABB(1, 0) < ( yPlaneCoordinate_ + GEO::TOL7) ) 
    {
      
      if(treeType_ == OCTTREE)
      {
        // check min_z less than z-plane
        if (AABB(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
          octants.push_back(0);
      
        // check max_z greater than z-plane
        if (AABB(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
          octants.push_back(4);
      }
      
      else if(treeType_ == QUADTREE)
        octants.push_back(0);
    }
    
    // check max_y greater than y-plane
    if (AABB(1, 1) > ( yPlaneCoordinate_ - GEO::TOL7) ) 
    {       
      if(treeType_ == OCTTREE)
      {
        // check max_z greater than z-plane
        if (AABB(2, 1) > ( zPlaneCoordinate_ - GEO::TOL7) )
          octants.push_back(6);
      
        // check min_z less than z-plane
        if (AABB(2, 0) < ( zPlaneCoordinate_ + GEO::TOL7) )
          octants.push_back(2);
      }
      
      else if(treeType_ == QUADTREE)
        octants.push_back(2);
    }
  }
  return octants;
}



/*----------------------------------------------------------------------*
 | classifiy element in node                               peder   07/08|
 *----------------------------------------------------------------------*/
std::vector<int> GEO::SearchTree::TreeNode::classifyElement(
    const DRT::Element*       		              element,
    const std::map<int,LINALG::Matrix<3,1> >& 	currentpositions
    ) const
{
  const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(element,currentpositions));
  GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
  GEO::checkRoughGeoType(element, xyze, eleGeoType);
  const LINALG::Matrix<3,2> elemXAABB(GEO::computeFastXAABB(element, xyze, eleGeoType));  
  return classifyXAABB(elemXAABB);
}



/*----------------------------------------------------------------------*
 | classifiy element in tree node                          u.may   07/08|
 *----------------------------------------------------------------------*/
std::vector<int> GEO::SearchTree::TreeNode::classifyElement(
        const DRT::Element*               element,
        const LINALG::SerialDenseMatrix&  xyze_element    
        ) const 
{
  GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
  GEO::checkRoughGeoType(element, xyze_element , eleGeoType);
  const LINALG::Matrix<3,2> elemXAABB(GEO::computeFastXAABB(element, xyze_element , eleGeoType));  
  return classifyXAABB(elemXAABB);
}



/*----------------------------------------------------------------------*
 | classifiy element in node                               u.may   08/08|
 *----------------------------------------------------------------------*/
std::vector<int> GEO::SearchTree::TreeNode::classifyRadius(
    const double                radius,
    const LINALG::Matrix<3,1>&	point
    ) const
{
  LINALG::Matrix<3,2> radiusXAABB;  

  for(int dim = 0; dim < 3; dim++)
  {
    radiusXAABB(dim,0) = ( point(dim) - radius) - GEO::TOL7;
    radiusXAABB(dim,1) = ( point(dim) + radius) + GEO::TOL7;
  }
  
  if(treeType_ == QUADTREE)
  {
    radiusXAABB(2,0) = 0.0;
    radiusXAABB(2,1) = 0.0;
  }
   
  return classifyXAABB(radiusXAABB);
}



/*----------------------------------------------------------------------*
 | update tree node                                        u.may   08/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::updateTreeNode(
   const LINALG::Matrix<3,2>& AABB_old, 
   const LINALG::Matrix<3,2>& AABB_new)
{
  for(int i = 0; i < getNumChildren(); i++)
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
 | and nearest object                                                   |
 *----------------------------------------------------------------------*/
int GEO::SearchTree::TreeNode::queryFSINearestObject(
    const DRT::Discretization&                      dis,
    const std::map<int,LINALG::Matrix<3,1> >&       currentpositions, 
    const LINALG::Matrix<3,1>&                      point,
    GEO::NearestObject&                             nearestObject
    ) 
{
  switch (treeNodeType_) 
  {
    case INNER_NODE:
    {       
      return children_[classifyPoint(point)]->queryFSINearestObject(dis, currentpositions, point, nearestObject);
      break;
    }
    case LEAF_NODE:   
    {
      if (elementList_.empty())
      {
        nearestObject = nearestObject_;
        return label_;
      }
      
      // max depth reached, counts reverse
      if (treedepth_ <= 0 || (elementList_.size()==1 && (elementList_.begin()->second).size() == 1) )
      {
        // nearest object refers only to the nearest object found in this particular tree node
        int xfemLabel = GEO::getXFEMLabelAndNearestObject(dis, currentpositions, point, elementList_, nearestObject);
        
        const TreeNode* workingNode = this;
        while(!GEO::pointInMinCircleInTreeNode(nearestObject.getPhysCoord(), point, workingNode->getNodeBox(), (!workingNode->hasParent()) ))
        {
          if(!workingNode->hasParent())
            dserror("this treenode has no parent");
          
          workingNode = workingNode->getParent();
          xfemLabel = GEO::getXFEMLabelAndNearestObject(dis, currentpositions, point, workingNode->getElementList(), nearestObject);
        }
        return xfemLabel;
      }
      //return GEO::getXFEMLabelAndNearestObject(dis, currentpositions, point, elementList_, nearestObject);

      // dynamically grow tree otherwise, create children and set label for empty children
      createChildren(dis, currentpositions);
      setXFEMLabelOfEmptyChildren(dis, currentpositions);
      // search in apropriate child node
      return children_[classifyPoint(point)]->queryFSINearestObject(dis, currentpositions, point, nearestObject);
      break;
    }
    default:
      dserror("should not get here\n");
  }
  return -1;
}



/*----------------------------------------------------------------------*
 | return xfem label for point (interface method)          u.may   07/08|
 *----------------------------------------------------------------------*/
int GEO::SearchTree::TreeNode::queryXFEMFSIPointType(
    const DRT::Discretization&                    dis,
    const std::map<int,LINALG::Matrix<3,1> >&     currentpositions, 
    const LINALG::Matrix<3,1>&                    point
    ) 
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
      {
        // nearest object refers only to the nearest object found in this particular tree node
        GEO::NearestObject nearestObjectInNode;
        int xfemLabel = GEO::getXFEMLabelAndNearestObject(dis, currentpositions, point, elementList_, nearestObjectInNode);
        
        const TreeNode* workingNode = this;
        while(!GEO::pointInTreeNode(nearestObjectInNode.getPhysCoord(), workingNode->nodeBox_))
        {
          if(!workingNode->hasParent())
            dserror("this treenode has no parent");
          
          workingNode = workingNode->getParent();
          xfemLabel = GEO::getXFEMLabelAndNearestObject(dis, currentpositions, point, workingNode->getElementList(), nearestObjectInNode);
        }
        return xfemLabel;
      }
      
      // dynamically grow tree otherwise, create children and set label for empty children
      createChildren(dis, currentpositions);
      setXFEMLabelOfEmptyChildren(dis, currentpositions);
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
std::map<int,std::set<int> > GEO::SearchTree::TreeNode::searchElementsInRadius(
    const DRT::Discretization& 	                  dis,
    const std::map<int,LINALG::Matrix<3,1> >&     currentpositions, 
    const LINALG::Matrix<3,1>&                    point,
    const double                                  radius, 
    const int                                     label) 
{  
  std::map<int,std::set<int> > eleMap;

  switch (treeNodeType_) 
  {
    case INNER_NODE:
    {       
      const vector<int> childindex = classifyRadius(radius, point);
      if(childindex.size() < 1)
        dserror("no child found\n");
      else if (childindex.size() ==1) // child node found which encloses AABB so step down
        return children_[childindex[0]]->searchElementsInRadius(dis, currentpositions, point, radius, label);
      else
        return GEO::getElementsInRadius(dis, currentpositions, point, radius, label, elementList_); 
      break;
    }
    case LEAF_NODE:   
    {
      if(elementList_.empty())
        return eleMap;

      // max depth reached, counts reverse
      if (treedepth_ <= 0 || (elementList_.size()==1 && (elementList_.begin()->second).size() == 1) )
        return GEO::getElementsInRadius(dis, currentpositions, point, radius, label, elementList_); 

      // dynamically grow tree otherwise, create children and set label for empty children
      // search in apropriate child node 
      const vector<int> childindex = classifyRadius(radius, point);
      if(childindex.size() < 1)
        dserror("no child found\n");
      else if (childindex.size() ==1) // child node found which encloses AABB so refine further
      { 
        createChildren(dis, currentpositions);
        return children_[childindex[0]]->searchElementsInRadius(dis, currentpositions, point, radius, label);
      }
      else // AABB does not fit into a single child node box anymore so don t refine any further
        return GEO::getElementsInRadius(dis, currentpositions, point, radius, label, elementList_); 
      break;
    }
    default:
      dserror("should not get here\n");
  }
  return eleMap;
}



/*----------------------------------------------------------------------*
 | returns a set of elements whose XAABB s are               u.may 09/08|
 | intersecting with the XAABB of a given volume element                |
 *----------------------------------------------------------------------*/
std::vector<int> GEO::SearchTree::TreeNode::queryIntersectionCandidates(
    const DRT::Discretization&                  dis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions, 
    DRT::Element*                               element,
    const LINALG::SerialDenseMatrix&            xyze_element)
{
  std::vector< int > elementset;
  
  switch(treeNodeType_) 
  {
    case INNER_NODE:
    {     
      const vector<int> childindex = classifyElement(element,xyze_element);
      if(childindex.size() < 1)
        dserror("no child found\n");
      else if (childindex.size() ==1)
        return children_[childindex[0]]->queryIntersectionCandidates(dis, currentpositions, element,xyze_element);
      else
        return GEO::getIntersectionCandidates(dis, currentpositions, element, elementList_); 
        
      break;
    }
    case LEAF_NODE:   
    {
      if(elementList_.empty())
        return elementset;

      // max depth reached, counts reverse
      if (treedepth_ <= 0 || (elementList_.begin()->second).size() == 1)
        return GEO::getIntersectionCandidates(dis, currentpositions, element, elementList_); 
  
      // dynamically grow tree otherwise, create children and set label for empty children
      // search in apropriate child node  
      const vector<int> childindex = classifyElement(element, xyze_element);
      if((int)childindex.size() < 1)
        dserror("no child found\n");
      else if ((int)childindex.size() == 1)
      {
        createChildren(dis, currentpositions);
        return children_[childindex[0]]->queryIntersectionCandidates(dis, currentpositions, element, xyze_element);
      }
      else
        return GEO::getIntersectionCandidates(dis, currentpositions, element, elementList_); 

      break;
    }
    default:
      dserror("should not get here\n");
  }
  return elementset;
}



/*----------------------------------------------------------------------*
 | returns intersection elements (CONTACT)                    popp 07/08|
 *----------------------------------------------------------------------*/
std::vector<int> GEO::SearchTree::TreeNode::searchIntersectionElements(
    const DRT::Discretization&                  dis,
    const std::map<int,LINALG::Matrix<3,1> >&   currentpositions, 
    DRT::Element*                               element) 
{
  std::vector<int> elementset;
  
  switch (treeNodeType_) 
  {
    case INNER_NODE:
    {     
      const vector<int> childindex = classifyElement(element,currentpositions);
      if(childindex.size() < 1)
        dserror("no child found\n");
      else if (childindex.size() ==1)
        return children_[childindex[0]]->searchIntersectionElements(dis, currentpositions, element);
      else
        return GEO::getIntersectionElements(dis, currentpositions, element, elementList_); 
        
      break;
    }
    case LEAF_NODE:   
    {
      if(elementList_.empty())
        return elementset;

      // max depth reached, counts reverse
      if (treedepth_ <= 0 || (elementList_.begin()->second).size() == 1)
        return GEO::getIntersectionElements(dis, currentpositions, element, elementList_);
  
      // dynamically grow tree otherwise, create children and set label for empty children
      // search in apropriate child node  
      const vector<int> childindex = classifyElement(element,currentpositions);
      if((int)childindex.size() < 1)
        dserror("no child found\n");
      else if ((int)childindex.size() == 1)
      {
        createChildren(dis, currentpositions);
        return children_[childindex[0]]->searchIntersectionElements(dis, currentpositions, element);
      }
      else
        return GEO::getIntersectionElements(dis, currentpositions, element, elementList_); 

      break;
    }
    default:
      dserror("should not get here\n");
  }
  return elementset;
}


/*----------------------------------------------------------------------*
 | print tree node to gmsh file                            peder   07/08|
 *----------------------------------------------------------------------*/
void GEO::SearchTree::TreeNode::printTreeNode(
    const int       max_depth,
    stringstream&   fc) const
{
  if(treedepth_== max_depth)
  {
    LINALG::Matrix<3,8> printBox(true);
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
    for (int i=0; i < getNumChildren(); i++)
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
    LINALG::Matrix<3,8> printBox(true);
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



/*----------------------------------------------------------------------*
 | get depth of tree node (METRICS)                        peder   07/08|
 *----------------------------------------------------------------------*/
int GEO::SearchTree::TreeNode::getDepth() const
{
  int depth = -1;
  
  if (treeNodeType_==LEAF_NODE)
    depth = treedepth_;
  else 
  {
    int tmp_depth = -1;
    depth = treedepth_;
    for (int i = 0; i < getNumChildren(); i++)
    {
      if (children_[i] != Teuchos::null)
        tmp_depth = children_[i]->getDepth();
      else 
        dserror("inner node has to have leave nodes");
      
      if(tmp_depth < depth)
        depth = tmp_depth;
    }
  }
  return depth;
}



/*----------------------------------------------------------------------*
 | get number of tree nodes in tree (METRICS)              u.may   08/08|
 *----------------------------------------------------------------------*/
int GEO::SearchTree::TreeNode::getNumNodesInTree() const
{
  int numTreeNodes = 0;
  
  switch(treeNodeType_)
  { 
    case LEAF_NODE:
    {
      numTreeNodes = 1;
      break;
    }
    case INNER_NODE:
    {
      for(int i = 0; i < getNumChildren(); i++)
        numTreeNodes += children_[i]->getNumNodesInTree();
  
      break;
    }
    default:
      dserror("wrong tree node type");
  }
  return numTreeNodes;
}


#endif  // #ifdef CCADISCRET

